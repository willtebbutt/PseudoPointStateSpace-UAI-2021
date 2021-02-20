println("housing")

using Revise

using Pkg
Pkg.activate(".");
Pkg.resolve();
Pkg.instantiate();

using ArgParse
using Bijectors
using Clustering
using CSV
using BenchmarkTools
using DrWatson
using FileIO
using Formatting
using GADM
using GeoDatasets
using Images
using LinearAlgebra
using LineSearches
using Optim
using ParameterHandling
using Plots
using Random
using Stheno
using TemporalGPs
using TexTables
using Zygote

using ParameterHandling: positive, deferred, flatten
using TemporalGPs: approx_posterior_marginals, RegularInTime, Separable
using Zygote: dropgrad

const exp_dir_name = joinpath("icml-2021", "house_data");
const data_dir = joinpath(datadir(), exp_dir_name);
const metadata_dir = joinpath(data_dir, "meta");
const training_output_dir = joinpath(data_dir, "training");
const prediction_output_dir = joinpath(data_dir, "predictions");
const results_dir = joinpath(projectdir(), "results", exp_dir_name);

include(joinpath(srcdir(), "housing_data.jl"))

# Pull out the period of time specified from the full data set, sort it, and return it.
function generate_raw_data_set(all_data, first_date::Date, last_date::Date)
    data = all_data[first_date .< all_data.date .< last_date, :]
    return sort(data, :date)
end

# Log and normalise the data set, and return the mean / std used to transform it.
function log_normalised_price(data, train_indices)
    log_price = log.(data.price)
    return (log_price .- mean(log_price[train_indices])) ./ std(log_price[train_indices])
end

# Take a vector of non-sequential `Int`s and construct a vector of sequential `Int`s with
# the same start and end numbers. Return two addition vectors, one containing the indices of
# the elements of the returned vector that are newly added, and the other the indices of
# the existing `Int`s.
function fill_in_empties(x::AbstractVector{Int}, y::AbstractVector, fill_val)
    return vcat(
        map((y, d) -> vcat([y], fill(fill_val, d - 1)), y, diff(x))...,
        [y[end]],
    )
end

smse(predicted, observed) = mean(abs2, filter(!ismissing, mean.(predicted) .- observed))

rsmse(predicted, observed) = sqrt(smse(predicted, observed))

function pplp(predicted, observed)
    present_idx = findall(!ismissing, observed)
    return mean(logpdf.(predicted[present_idx], observed[present_idx]))
end

# Encodes whether or not a property is a new-build as a binary variable.
encode_new_build(new_build_flags::Vector) = map(x -> x == "N" ? 0.0 : 1.0, new_build_flags)

# Converts property types to a 1-hot encoding.
function encode_property_type(property_types::Vector{String})
    ks = unique(property_types)
    vals = eachindex(ks)
    mapping = Dict([k => val for (k, val) in zip(ks, vals)])

    feature_mat = fill(0.0, length(ks), length(property_types))
    for n in eachindex(property_types)
        feature_mat[mapping[property_types[n]], n] = 1.0
    end
    return feature_mat
end

function train_indices_name(start_date::Date, finish_date::Date, N_train::Int)
    return joinpath(data_dir, "$start_date-$finish_date-$N_train.jld2")
end

function maybe_create_train_indices(
    rng::AbstractRNG, start_date::Date, finish_date::Date, N_train::Int, N_total::Int,
)
    fname = train_indices_name(start_date, finish_date, N_train)

    if !isfile(fname)
        @info "CREATING NEW TEST INDICES FOR DATES $start_date-$finish_date-$N_train."
        train_indices = randperm(rng, N_total)[1:N_train]
        FileIO.save(fname, Dict("train_indices" => train_indices))
    else
        @info "NOT CREATING NEW TEST INDICES AS THEY ALREADY EXIST."
    end
end

function get_train_indices(start_date::Date, finish_date::Date, N_train::Int)
    return FileIO.load(
        train_indices_name(start_date, finish_date, N_train),
    )["train_indices"]
end

function generate_data_set(
    data, first_date::Date, last_date::Date, train_indices::Vector{Int},
)
    # Mark train and test data.
    data.is_train = fill(false, size(data, 1))
    data.is_train[train_indices] .= true

    # Standardise the data using only statistics from the training data.
    data.normalised_price = log_normalised_price(data, train_indices)

    # Get new build encoding.
    new_build = encode_new_build(collect(data.new_build))

    # Get property type encoding.
    property_type = encode_property_type(collect(data.property_type))

    data.new_build_enc = new_build
    data.property_type_enc = ColVecs(property_type)

    day_1 = minimum(data.date)
    N_days = (maximum(data.date) - day_1).value + 1
    ts = RegularSpacing(1.0, 1.0, N_days)
    data.integer_day = map(d -> (d - day_1).value + 1, data.date)

    train_data = build_processed_dataset(filter(row -> row.is_train, data), ts)
    test_data = build_processed_dataset(filter(row -> !row.is_train, data), ts)
    return train_data, test_data
end

function build_processed_dataset(data::DataFrame, ts::RegularSpacing)

    # Group the training data by day / integer date.
    daily_data = map(collect(groupby(data, :date))) do d
        X = collect(hcat(d.lat, d.lon)')
        return (
            t = only(unique(d.integer_day)),
            x = ColVecs(X),
            y = collect(d.normalised_price),
        )
    end

    # Construct dictionary from the vector.
    daily_data_dict = Dict([d.t => d for d in daily_data])

    # Define fallbacks for what should be filled in if no data is provided.
    default = (
        x = ColVecs(collect(reshape(first(first(daily_data).x), :, 1))),
        y = [missing],
    )

    # Construct inputs.
    x_space = map(t -> get(daily_data_dict, t, default).x, ts)
    inputs = RegularInTime(ts, x_space)

    # Construct outputs.
    outputs = vcat(map(t -> get(daily_data_dict, t, default).y, ts)...)

    return (inputs=inputs, outputs=outputs)
end

# Utility functionality for working with `ParameterHandling` a bit more conveniently.
function wrap_objective(f, θ)
    θ_flat, unflatten = flatten(θ)
    f_flat(θ_flat) = f(ParameterHandling.value(unflatten(θ_flat)))
    return f_flat, θ_flat, unflatten
end

function ParameterHandling.flatten(x::ColVecs)
    X_flat, unflatten_X = flatten(x.X)
    unflatten_ColVecs(X_flat) = ColVecs(unflatten_X(X_flat))
    return X_flat, unflatten_ColVecs
end

# Make it possible to parse strings from the command line into Vector{Int}s.
function ArgParse.parse_item(::Type{Vector{Int}}, x::String)
    @assert x[1] == '[' && x[end] == ']'
    return CSV.read(IOBuffer(x[2:prevind(x, end)]); header=false, transpose=true).Column1
end

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--run-experiments"
            help = "true / false - run the experiments or not"
            arg_type = Bool
            default = true
        "--M"
            help = "number of pseudo-points"
            arg_type = Int
            default = 50
    end
    return parse_args(s)
end

const args = parse_commandline()

# Just load the big data set once and re-use it lots. It's quite time consuming to load it
# back up lots of times, particularly when iterating on this script.
const all_data = @time generate_housing_data(data_dir);





# Restrict data set to consider only London. These coordinates pretty much restrict us to
# the M25 corridor.
london_lons = (-0.55, 0.31);
london_lats = (51.27, 51.72);

london_lons = (-6, 2);
london_lats = (50, 56);

struct GridBox{Tlats<:Tuple{<:Real, <:Real}, Tlons<:Tuple{<:Real, <:Real}}
    lats::Tlats
    lons::Tlons
end

function in_box(box::GridBox, lat::Real, lon::Real)
    return box.lats[1] <= lat && box.lats[2] > lat &&
        box.lons[1] <= lon && box.lons[2] > lon
end

# Filter observational data and restrict to a box around Greater London.
london_data = filter(
    row -> in_box(GridBox(london_lats, london_lons), row.lat, row.lon) &&
        row.property_type == "F",
    all_data,
);

# Define locations at which to compute the posterior marginals.
N_base = 50; # related to how densely we plot the points.
N_lons = Int(floor(N_base * (london_lons[2] - london_lons[1])));
N_lats = Int(floor(N_base * (london_lats[2] - london_lats[1])));
const london_lats_range = range(london_lats[1], london_lats[2]; length=N_lats);
const london_lons_range = range(london_lons[1], london_lons[2]; length=N_lons);
const lat_lon_pairs = tuple.(london_lats_range, london_lons_range');

BLAS.set_num_threads(8)


# Set proportion of data to be used as training.
const train_prop = 0.5





#
# Specify experiments to run.
#

mkpath(metadata_dir)
const settings_fname = joinpath(metadata_dir, "settings.bson")

let
    wsave(
        settings_fname,
        :settings => dict_list(
            Dict(
                # The total number of pseudo-points to use per time step.
                :M => [75],

                # Kernels to try, specified using the `ParameterHandling.jl` interface.
                :kernel_info => [
                    (
                        name = "separable",
                        value = deferred(
                            (λ_time, λ_space, s) -> Separable(
                                s * stretch(EQ(), λ_space),
                                stretch(Matern32(), λ_time),
                            ),
                            positive(1e-2),
                            fill(positive(1.0), 2),
                            positive(1.0),
                        ),
                    ),
                    (
                        name = "sum-separable",
                        value = deferred(
                            (λs_time, λs_space, ss) -> begin
                                k1 = Separable(
                                    ss[1] * stretch(EQ(), λs_space[1]),
                                    stretch(Matern32(), λs_time[1]),
                                )
                                k2 = Separable(
                                    ss[2] * stretch(EQ(), λs_space[2]),
                                    stretch(Matern32(), λs_time[2]),
                                )
                                return k1 + k2
                            end,
                            (positive(1e-3), positive(1e-1)),
                            (
                                fill(positive(1.0), 2),
                                fill(positive(5.0), 2),
                            ),
                            (positive(0.7), positive(0.3)),
                        ),
                    ),
                ],

                # The dates between which the training data is selected.
                :date_range => [
                    (start=Date(2010), stop=Date(2020)),
                ]
            ),
        ),
    );
end









#
# Run experiments. This is split into two sections:
# 1. training: this takes by far the longest amount of time. The entire model is saved as an
#     output from this step.
# 2. prediction: loads the results from step 1 and uses them to make predictions. Not run in
#     the training block because you don't want to have to re-run the training procedure to
#     make slightly different kinds of predictions. Not in the results block because it can
#     take a while, which is less than helpful when tweaking plots.
#

if args["run-experiments"]

println("Running experiments")

mkpath(training_output_dir)

let
    println("Loading housing data")

    settings = load(settings_fname)[:settings];

    for (n, setting) in enumerate(settings)

        println("Experiment $n / $(length(settings))")
        display(setting)
        println()

        # Generate the data set.
        first_date = setting[:date_range].start
        last_date = setting[:date_range].stop

        london_data_local = generate_raw_data_set(london_data, first_date, last_date);
        N_total = size(london_data_local, 1);
        N_train = Int(floor(train_prop * N_total));
        maybe_create_train_indices(
            MersenneTwister(123456), first_date, last_date, N_train, N_total,
        )

        train_indices = get_train_indices(first_date, last_date, N_train);
        train_data, test_data = generate_data_set(
            london_data_local, first_date, last_date, train_indices,
        );

        inputs = train_data.inputs
        outputs = train_data.outputs
        println("Total number of training data is $(length(outputs)).")
        println("Total number of test data is $(length(test_data.outputs))")

        # Compute pseudo-point initialisation.
        x_r_full = ColVecs(hcat(getfield.(inputs.vs, :X)...))
        Random.seed!(123456)
        clustering_result = kmeans(x_r_full.X, setting[:M])
        Z_r = clustering_result.centers

        # Specify the initial parameters in structured form.
        initial_parameters = (
            k = setting[:kernel_info].value,
            σ² = bounded(0.5, 1e-2, 2.0),
            Z_r = fixed(Z_r),
        )

        # Create the objective function, wrapping it in `ParameterHandling.jl` stuff to
        # make it play nicely with Optim.
        objective, flat_initial_parameters, unflatten = wrap_objective(
            θ -> begin
                f = to_sde(GP(θ.k, GPC()))
                return -elbo(f(dropgrad(inputs), θ.σ²), dropgrad(outputs), ColVecs(θ.Z_r))
            end,
            initial_parameters,
        )

        # Construct the optimiser.
        training_results = Optim.optimize(
            objective,
            θ -> Zygote.gradient(objective, θ)[1],
            flat_initial_parameters,
            LBFGS(
                alphaguess = LineSearches.InitialStatic(scaled=true),
                linesearch = LineSearches.BackTracking(),
                m=50,
            ),
            # The numbers involved in my objective function are generally quite large,
            # so looser tolerances are appropriate.
            Optim.Options(
                x_tol = 1e-4,
                f_tol = 1e-8,
                g_tol = 1e-1,
                f_calls_limit = 1_000,
                g_calls_limit = 100,
                iterations = 1_000,
                show_trace = true,
                extended_trace=true,
                show_every=10,
            );
            inplace=false,
        )

        wsave(
            joinpath(
                training_output_dir,
                savename(
                    Dict(
                        :M => setting[:M],
                        :kernel => setting[:kernel_info].name,
                        :start_date => first_date,
                    ),
                    "bson",
                ),
            ),
            Dict(
                :setting => setting,
                :training_results => training_results,
                :unflatten => unflatten,
            ),
        )
    end
end

end










#
# Make predictions at locations that we wish to view them, and compute RMSE on data etc.
#

println("Generating predictions.")

mkpath(prediction_output_dir)

let
    trained_models = collect_results(training_output_dir);

    for (n, trained_model_) in enumerate(eachrow(trained_models))

        println("n, $(trained_model_.path)")

        # Defensive copy. Something in here is mutating for some unknown reason...
        trained_model = deepcopy(trained_model_);

        # Generate the data set.
        first_date = trained_model.setting[:date_range].start
        last_date = trained_model.setting[:date_range].stop

        # Generate the dataset.
        london_data_local = generate_raw_data_set(london_data, first_date, last_date);
        N_total = size(london_data_local, 1);
        N_train = Int(floor(train_prop * N_total));
        train_indices = get_train_indices(first_date, last_date, N_train);
        train_data, test_data = generate_data_set(
            london_data_local, first_date, last_date, train_indices,
        );
        println("Total number of training data is $(length(train_data.outputs)).")
        println("Total number of test data is $(length(test_data.outputs))")

        # Build the model.
        unflatten = trained_model.unflatten
        results = trained_model.training_results
        learned_parameters = ParameterHandling.value(unflatten(results.minimizer))

        f = to_sde(GP(learned_parameters.k, GPC()))
        fx = f(train_data.inputs, learned_parameters.σ²);

        # Make predictions at entire grid of points at which we might consider predicting.
        get_vector(x::Tuple) = [x[1], x[2]]
        x_r_pr = ColVecs(
            vcat(
                hcat(get_vector.(lat_lon_pairs)...),
                # repeat([0.0, 0.0, 1.0, 0.0, 0.0, 0.0], 1, length(lat_lon_pairs)),
            ),
        );

        # Compute posterior marginals for the sake of plotting.
        z_r = ColVecs(learned_parameters.Z_r)
        pos = length(TemporalGPs.get_time(train_data.inputs)) - 100
        @show pos, length(TemporalGPs.get_time(train_data.inputs))
        preds = @time approx_posterior_marginals(
            dtc, fx, train_data.outputs, z_r, x_r_pr, pos,
        )

        # Compute posterior marginals at the training data and the test data.
        preds_train = @time approx_posterior_marginals(
            dtc, fx, train_data.outputs, z_r, train_data.inputs,
        )
        preds_test = @time approx_posterior_marginals(
            dtc, fx, train_data.outputs, z_r, test_data.inputs,
        )

        # Write the results.
        wsave(
            joinpath(
                prediction_output_dir,
                savename(
                    Dict(
                        :M => trained_model.setting[:M],
                        :kernel => trained_model.setting[:kernel_info].name,
                        :start_date => first_date,
                    ),
                    "bson",
                ),
            ),
            Dict(
                "trained_model" => trained_model,
                "preds_train" => preds_train,
                "preds_test" => preds_test,
                "preds_extrapolate" => preds, 
                "x_r_pr" => x_r_pr,
                "lat_lon_pairs" => lat_lon_pairs,
                "smse_train" => rsmse(preds_train, train_data.outputs),
                "smse_test" => rsmse(preds_test, test_data.outputs),
                "pplp_train" => pplp(preds_train, train_data.outputs),
                "pplp_test" => pplp(preds_test, test_data.outputs),
            ),
        )
    end
end










#
# Visualise results.
#



# Load shoreline data.
# const shore_data = GeoDatasets.gshhg('f');

my_format(x) = format(round(x; sigdigits=3); stripzeros=true)

let
    gr()

    results = collect_results(prediction_output_dir)

    name_map = Dict(
        "separable" => "Separable",
        "sum-separable" => "Sum-Separable",
    );

    # Compute some additional fields necessary for table creation.
    results.M = map(x -> x.setting[:M], results.trained_model);
    results.start_date = map(x -> x.setting[:date_range].start, results.trained_model);
    results.rsmse_test_string = map(my_format, results.smse_test);
    results.npplp_test_string = map(my_format, -results.pplp_test);
    results.row_keys = map(
        x -> "($(x.setting[:M])) - " * name_map[x.setting[:kernel_info].name],
        results.trained_model,
    );

    # Construct tables for each training window.
    for group in groupby(results, :start_date)

        headers = [
            "RSMSE",
            "NPPLP",
        ]
        col_data = [
            group.rsmse_test_string,
            group.npplp_test_string,
        ]
        cols = [
            TableCol(header, collect(group.row_keys), collect(data)) for
            (header, data) in zip(headers, col_data)
        ]
        results_tex_table = hcat(cols...)

        # Generate paper-ready tables of results, comparing each model.
        start_date = only(unique(group.start_date))
        open(joinpath(results_dir, "performance_results_date=$start_date.txt"), "w") do io
            write(io, to_ascii(results_tex_table))
        end
        open(joinpath(results_dir, "performance_results_date=$start_date.tex"), "w") do io
            write(io, to_tex(results_tex_table))
        end
    end

    box = GridBox(london_lats, london_lons);

    # London boroughs info, including geometries.
    # boroughs = GADM.get("GBR", "England", "Greater London"; children=true)[2];
    boroughs = vcat(
        GADM.get("GBR", "England"; children=true)[2].geom,
        GADM.get("GBR", "Wales"; children=true)[2].geom,
    );

    for result in eachrow(results)

        trained_model = result.trained_model;
        preds_mat = reshape(result.preds_extrapolate, N_lats, :);

        # Generate the dataset.
        first_date = trained_model.setting[:date_range].start
        last_date = trained_model.setting[:date_range].stop
        london_data_local = generate_raw_data_set(london_data, first_date, last_date);
        N_total = size(london_data_local, 1);
        N_train = Int(floor(train_prop * N_total));
        train_indices = get_train_indices(first_date, last_date, N_train);
        train_data, test_data = generate_data_set(
            london_data_local, first_date, last_date, train_indices,
        );

        x_r = train_data.inputs.vs[end-20];

        unflatten = trained_model.unflatten
        results = trained_model.training_results
        learned_parameters = ParameterHandling.value(unflatten(results.minimizer))
        z_r = ColVecs(learned_parameters.Z_r)

        colour_scheme = :viridis;

        let
            plt = plot(
                xlims=london_lons,
                ylims=london_lats,
                xticks=nothing,
                yticks=nothing,
            )

            heatmap!(
                plt, london_lons_range, london_lats_range, mean.(preds_mat);
                label="", c=cgrad(colour_scheme), right_margin=6Plots.mm,
            );

            # scatter!(
            #     plt, getindex.(x_r, 2), first.(x_r);
            #     label="Sold Houses",
            #     markersize=5,
            #     marker=:square,
            #     color=colorant"#A95AA1", # Purple
            # )
            scatter!(
                plt, getindex.(z_r, 2), first.(z_r);
                label="z",
                markersize=4,
                color=colorant"#F5793A", # Orange
            )

            # Plot the Greater London boroughs.
            for b in boroughs
                plot!(plt, b; fillcolor=nothing)
            end

            mkpath(results_dir)
            savefig(
                plt,
                joinpath(
                    results_dir,
                    savename(
                        Dict(
                            :M => trained_model.setting[:M],
                            :kernel => trained_model.setting[:kernel_info].name,
                            :start_date => first_date,
                        ),
                        "posterior_mean.png",
                    ),
                )
            )
        end

        let
            plt = plot(
                xlims=london_lons,
                ylims=london_lats,
                xticks=nothing,
                yticks=nothing,
            )

            heatmap!(
                london_lons_range, london_lats_range, std.(preds_mat);
                label="", c=colour_scheme, right_margin=3Plots.mm,
            )

            # scatter!(plt, getindex.(x_r, 2), first.(x_r); label="Sold Houses", color=:blue)
            # scatter!(plt, getindex.(z_r, 2), first.(z_r); label="z", α=0.8, color=:red)

            # scatter!(
            #     plt, getindex.(x_r, 2), first.(x_r);
            #     label="Sold Houses",
            #     markersize=5,
            #     marker=:square,
            #     color=colorant"#A95AA1", # Purple
            # )
            scatter!(
                plt, getindex.(z_r, 2), first.(z_r);
                label="z",
                markersize=4,
                color=colorant"#F5793A", # Orange
            )

            # Plot the Greater London boroughs.
            for b in boroughs
                plot!(plt, b; fillcolor=nothing)
            end

            mkpath(results_dir)
            savefig(
                plt,
                joinpath(
                    results_dir,
                    savename(
                        Dict(
                            :M => trained_model.setting[:M],
                            :kernel => trained_model.setting[:kernel_info].name,
                            :start_date => first_date,
                        ),
                        "posterior_std.png",
                    ),
                )
            )
        end
    end
end












# Functionality to construct a land-sea mask that's actually at high resolution using
# ArchGDAL's functionality. It's rather slow, so is probably the kind of thing that you
# really only want to run once. Will do this later on when closer to the deadline if time
# permits.
# using GeometricalPredicates
# geom = vcat(
#     GADM.get("GBR", "England"; children=true)[2].geom,
#     GADM.get("GBR", "Wales"; children=true)[2].geom,
# )[1];
# geom = GADM.get("GBR", "England").geom;
# london_lons = (-6, 2);
# london_lats = (50, 56);

# pt = ArchGDAL.createpoint(-1.23, 52.24);

# ArchGDAL.within.(geom, Ref(pt))

# ArchGDAL.contains.(geom, Ref(pt))

# function point_map(lat_lon_pairs, geom)
#     return map(
#         ((lat, lon), ) -> ArchGDAL.contains(geom, ArchGDAL.createpoint(lon, lat)),
#         lat_lon_pairs,
#     )
# end

# pm = point_map(lat_lon_pairs, geom[1])


# boroughs = vcat(
#     GADM.get("GBR", "England"; children=true)[2].geom,
#     GADM.get("GBR", "Wales"; children=true)[2].geom,
# );
