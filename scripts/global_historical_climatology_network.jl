println("GHCN")

using Pkg
Pkg.activate(".");
Pkg.instantiate();

using ArgParse
using Bijectors
using Clustering
using CSV
using DataFrames
using Dates
using BenchmarkTools
using DrWatson
using FileIO
using Formatting
using GADM
using GeoDatasets
using GHCNData
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

using TemporalGPs:
    approx_posterior_marginals,
    get_time,
    get_space,
    RectilinearGrid,
    Separable

const exp_dir_name = joinpath("global_historical_climatology_network");
const data_dir = joinpath(datadir(), exp_dir_name);
const metadata_dir = joinpath(data_dir, "meta");
const training_output_dir = joinpath(data_dir, "training");
const prediction_output_dir = joinpath(data_dir, "predictions");
const results_dir = joinpath(projectdir(), "results", exp_dir_name);

# Washington + a bit of Canada.
const lon_interval = (-125, -122);
const lat_interval = (47, 49);

# Time intervals to consider.
const time_interval = (Date(1980), Date(2020));
const element = "TMAX"

function make_tmax_data_path(lat_interval, lon_interval, time_interval)
    return joinpath(
        data_dir,
        "tmax_data-$(lat_interval[1])-$(lat_interval[2])" *
            "-$(lon_interval[1])-$(lon_interval[2])" *
            "-$(time_interval[1])-$(time_interval[2])" *
            ".jld2",
    )
end

raw_tmax_data_path = make_tmax_data_path(lat_interval, lon_interval, time_interval);

raw_tmax_data, station_lat_lon_pairs, IDs = if isfile(raw_tmax_data_path)

    @info "DATA ALREADY SERIALISED"

    # If the data has been serialised already, just load it.
    data = FileIO.load(raw_tmax_data_path)

    data["raw_tmax_data"], data["station_lat_lon_pairs"], data["IDs"]
else

    @info "SERIALISING RAW DATA"

    # Otherwise load the raw data, serialise it, and load the serialised version.
    raw_tmax_data, station_lat_lon_pairs, IDs = GHCNData.select_data(
        time_interval, lat_interval, lon_interval, element,
    )

    # Serialise the data.
    FileIO.save(
        raw_tmax_data_path,
        Dict(
            "raw_tmax_data" => raw_tmax_data,
            "station_lat_lon_pairs" => station_lat_lon_pairs,
            "IDs" => IDs,
        ),
    )

    # Load the data and return it.
    data = FileIO.load(raw_tmax_data_path)
    data["raw_tmax_data"], data["station_lat_lon_pairs"], data["IDs"]
end

const station_lats = first.(station_lat_lon_pairs);
const station_lons = last.(station_lat_lon_pairs);

# Get the elevations associated with each station.
station_metadata = load_station_metadata();
elevation_dict = Dict([row.ID => row.ELEVATION for row in eachrow(station_metadata)]);
elevations = [elevation_dict[ID] for ID in IDs];

# TMAX data is in 10ths of a degree C. Remove obvious outliers.
const tmax_data = map(
    x -> (x !== missing && (-30 < x < 60)) ? x : missing,
    raw_tmax_data ./ 10,
);
const tmax_dates = time_interval[1]:Day(1):(time_interval[2]-Day(1))

# Determine locations at which to make predictions.
const N_base = 250;
const N_lons = N_base * (lon_interval[2] - lon_interval[1]);
const N_lats = N_base * (lat_interval[2] - lat_interval[1]);
const uk_lats_range = range(lat_interval[1], lat_interval[2]; length=N_lons);
const uk_lons_range = range(lon_interval[1], lon_interval[2]; length=N_lats);
const lat_lon_grid = tuple.(uk_lats_range, uk_lons_range');


# Look at some basic statistics of the data.
@show minimum(filter(!ismissing, tmax_data));
@show maximum(filter(!ismissing, tmax_data));
@show prop_missing = count(ismissing, tmax_data) / length(tmax_data);

# Construct the spatial locations of the stations in a GP-friendly format.
# Combine with the elevation data.
get_vector(x::Tuple{A, B} where {A, B}) = [x[1], x[2]]
get_vector(x::Tuple{A, B, C} where {A, B, C}) = [x[1], x[2], x[3]]

lat_lon_mat = hcat(get_vector.(station_lat_lon_pairs)...);
space_inputs_mat = Float64.(vcat(lat_lon_mat, elevations' ./ 100));
const space_inputs = ColVecs(space_inputs_mat);

# Utility functionality for working with `ParameterHandling` a bit more conveniently.
function wrap_objective(f, θ)
    θ_flat, unflatten = flatten(θ)
    function f_flat(θ_flat)
        return f(ParameterHandling.value(unflatten(θ_flat)))
    end
    return f_flat, θ_flat, unflatten
end

function construct_training_outputs(test_idx::Vector{Int}, outputs::AbstractVector)
    train_outputs = deepcopy(outputs)
    train_outputs[test_idx] .= missing
    return test_idx, train_outputs
end

function test_indices_name(start_date::Date, finish_date::Date, Nte::Int)
    return joinpath(data_dir, "$start_date-$finish_date-$Nte.jld2")
end

function maybe_create_test_indices(
    rng::AbstractRNG,
    start_date::Date,
    finish_date::Date,
    Nte::Int,
    outputs::AbstractVector,
)
    fname = test_indices_name(start_date, finish_date, Nte)

    if !isfile(fname)
        @warn "CREATING NEW TEST INDICES FOR DATES $start_date-$finish_date-$Nte."
        present_idx = findall(!ismissing, outputs)
        test_indices = present_idx[randperm(rng, length(present_idx))[1:Nte]]
        FileIO.save(fname, Dict("test_indices" => test_indices))
    else
        @warn "NOT CREATING NEW TEST INDICES AS THEY ALREADY EXIST."
    end
end

function get_test_indices(start_date::Date, finish_date::Date, Nte::Int)
    fname = test_indices_name(start_date, finish_date, Nte)
    return FileIO.load(fname)["test_indices"]
end

# Make it possible to parse strings from the command line into Vector{Int}s.
function ArgParse.parse_item(::Type{Vector{Int}}, x::String)
    @assert x[1] == '[' && x[end] == ']'
    return CSV.read(IOBuffer(x[2:prevind(x, end)]); header=false, transpose=true).Column1
end

function bad_experiment_type_error(experiment_type::String)
    throw(error("Unknown experiment type: $experiment_type."))
end

function sor_data(inputs::RectilinearGrid, outputs::AbstractVector, M::Int)

    # Compute the K-means clustering of the data.
    println(size(get_space(inputs).X))
    cluster_indices = if size(get_space(inputs).X, 2) > M
        println("Using subset of data.")
        Random.seed!(123456)
        clustering_result = kmeans(get_space(inputs).X, M)
        cluster_indices = fill(0, M)
        for (n, a) in enumerate(assignments(clustering_result))
            if cluster_indices[a] == 0
                cluster_indices[a] = n
            end
        end
        cluster_indices
    else
        println("Using all data!")
        collect(1:M)
    end

    # Subset the SoR data.
    inputs_sor = RectilinearGrid(get_space(inputs)[cluster_indices], get_time(inputs))
    outputs_sor = collect(vec(
        reshape(outputs, length(get_space(inputs)), :)[cluster_indices, :],
    ))
    return inputs_sor, outputs_sor
end

# Default behaviour is Matern32 -- Matern52 experiments were added afterwards.
function build_kernel(::Val{:separable}, params::NamedTuple)
    return Separable(
        stretch(EQ(), params.λ_space),
        params.σ^2 * stretch(Matern32(), params.λ_time),
    )
end

function build_kernel(::Val{:sum_separable}, params::Vector{<:NamedTuple})
    length(params) != 2 && throw(error("Expected two sets of params."))
    k1 = build_kernel(Val(:separable), params[1])
    k2 = build_kernel(Val(:separable), params[2])
    return k1 + k2
end

function build_kernel(::Val{:separable_matern_52}, params::NamedTuple)
    return Separable(
        stretch(EQ(), params.λ_space),
        params.σ^2 * stretch(Matern52(), params.λ_time),
    )
end

function build_kernel(::Val{:sum_separable_matern_52}, params::Vector{<:NamedTuple})
    length(params) != 2 && throw(error("Expected two sets of params."))
    k1 = build_kernel(Val(:separable_matern_52), params[1])
    k2 = build_kernel(Val(:separable_matern_52), params[2])
    return k1 + k2
end

function build_kernel(params::NamedTuple{(:kernel_type, :params)})
    return build_kernel(params.kernel_type, params.params)
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
            default = 60
    end
    return parse_args(s)
end

const args = parse_commandline()

BLAS.set_num_threads(8);

# Need to be able to handle `Val`s for dispatch reasons.
function ParameterHandling.flatten(::Type{T}, ::Val{V}) where {T<:Real, V}
    v = T[]
    unflatten_to_Val(::Vector{T}) = Val{V}()
    return v, unflatten_to_Val
end

ParameterHandling.value(v::Val) = v





#
# Specify experiments to run.
#

const settings_fname = joinpath(metadata_dir, "settings.bson")

let
    mkpath(metadata_dir)

    separable_kernel_info_matern_32 = (
        name = "separable-Matern-32",
        value = (
            kernel_type = Val(:separable),
            params = (
                σ = positive(1.0),
                λ_space = fill(positive(1.0), 3),
                λ_time = positive(1e-2),
            ),
        ),
    )

    sum_separable_kernel_info_matern_32 = (
        name = "sum-separable-Matern32",
        value = (
            kernel_type = Val(:sum_separable),
            params = [
                (
                    σ = positive(0.7),
                    λ_space = fill(positive(1.0), 3),
                    λ_time = positive(1e-2),
                ),
                (
                    σ = positive(0.3),
                    λ_space = fill(positive(5.0), 3),
                    λ_time = positive(1e-1),
                ),
            ],
        ),
    )

    separable_kernel_info_matern_52 = (
        name = "separable-Matern-52",
        value = (
            kernel_type = Val(:separable_matern_52),
            params = (
                σ = positive(1.0),
                λ_space = fill(positive(1.0), 3),
                λ_time = positive(1e-2),
            ),
        ),
    )

    sum_separable_kernel_info_matern_52 = (
        name = "sum-separable-Matern52",
        value = (
            kernel_type = Val(:sum_separable_matern_52),
            params = [
                (
                    σ = positive(0.7),
                    λ_space = fill(positive(1.0), 3),
                    λ_time = positive(1e-2),
                ),
                (
                    σ = positive(0.3),
                    λ_space = fill(positive(5.0), 3),
                    λ_time = positive(1e-1),
                ),
            ],
        ),
    )

    wsave(
        settings_fname,
        :settings => vcat(

            dict_list(Dict(
                # Use the subset of regressors approximation.
                :experiment_type => ["SoD"],

                # Try both kinds of kernel.
                :kernel_info => [sum_separable_kernel_info_matern_52],

                # The total number of spatial locations to use.
                # :M => [1, 2, 5],
                :M => [10, 20, 50],
                # :M => [99],

                # The first date for the training data.
                :start_date => [Date(2000)],
            )),

            dict_list(Dict(
                # Use the new pseudo-point approximation.
                :experiment_type => ["pseudo-point"],

                # The total number of pseudo-points to use per time step.
                # :M => [1, 2, 5],
                :M => [10, 20, 50],
                # :M => [50],

                # Kernels to try, specified using the `ParameterHandling.jl` interface.
                :kernel_info => [sum_separable_kernel_info_matern_52],

                # The first date for the training data.
                :start_date => [Date(2000)],
            ),),
        ),
    )

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

# if args["run-experiments"]

println("Running experiments")

mkpath(training_output_dir)

let
    settings = load(settings_fname)[:settings];

    for (n, setting) in enumerate(settings)

        println("Experiment $n / $(length(settings))")
        display(setting)
        println()

        experiment_type = setting[:experiment_type]

        # Compute the observation times.
        first_date = setting[:start_date]
        first_date_index = only(findall(tmax_dates .== first_date))

        time_inputs_as_dates = tmax_dates[first_date_index:end]
        time_inputs_as_range = RegularSpacing(0.0, 1.0, length(time_inputs_as_dates))
        inputs = RectilinearGrid(space_inputs, time_inputs_as_range)
        outputs_raw = collect(vec(tmax_data[:, first_date_index:end]))

        # If a train-test split hasn't already been created, create one and write it to
        # disk so that it is permanent. SHOULD BE COMMITTED IN THE GIT HISTORY!
        N_obs = length(filter(!ismissing, outputs_raw))
        Nte = Int(floor(0.1 * N_obs))
        rng = MersenneTwister(123456)
        finish_date = time_interval[2]
        maybe_create_test_indices(rng, first_date, finish_date, Nte, outputs_raw)

        # Construct the training data set.
        test_idx = get_test_indices(first_date, finish_date, Nte)
        test_idx, outputs_train_raw = construct_training_outputs(test_idx, outputs_raw)

        # Standardise the training data.
        m = mean(filter(!ismissing, outputs_train_raw))
        s = std(filter(!ismissing, outputs_train_raw))
        outputs_train = (outputs_train_raw .- m) ./ s
        @show typeof(outputs_train)
        println("Total number of training data is $(length(filter(!ismissing, outputs_train))).")
        println("Total number of data is $(length(filter(!ismissing, outputs_raw)))")
        @show mean(filter(!ismissing, outputs_train)), std(filter(!ismissing, outputs_train))

        # Create the objective function, wrapping it in `ParameterHandling.jl` stuff to
        # make it play nicely with Optim.
        objective, flat_initial_parameters, unflatten = if experiment_type == "SoD"

            inputs_sor, outputs_train_sor = sor_data(inputs, outputs_train, setting[:M])

            # Pick one point from each cluster.
            wrap_objective(
                θ -> begin
                    k = build_kernel(θ.k.kernel_type, θ.k.params)
                    f = to_sde(GP(k, GPC()))
                    fx = f(dropgrad(inputs_sor), θ.σ²)
                    return -logpdf(fx, dropgrad(outputs_train_sor))
                end,
                (
                    k = setting[:kernel_info].value,
                    σ² = bounded(0.5, 1e-2, 2.0),
                ),
            )
        elseif experiment_type == "pseudo-point"

            # Compute pseudo-point initialisation. Because Clustering doesn't let you
            # provide a seed to the kmeans++ initialisation, it's necessary to set the
            # global seed.
            Random.seed!(123456)
            clustering_result = kmeans(space_inputs.X, setting[:M])
            Z_r = clustering_result.centers

            wrap_objective(
                θ -> begin
                    k = build_kernel(θ.k.kernel_type, θ.k.params)
                    f = to_sde(GP(k, GPC()))
                    fx = f(dropgrad(inputs), θ.σ²)
                    z_r = ColVecs(θ.Z_r)
                    return -elbo(fx, dropgrad(outputs_train), z_r)
                end,
                (
                    k = setting[:kernel_info].value,
                    σ² = bounded(0.5, 1e-2, 2.0),
                    Z_r = fixed(Z_r),
                ),
            )
        else
            bad_experiment_type_error(experiment_type)
        end

        # Optimise the parameters.
        training_results = Optim.optimize(
            θ -> TemporalGPs.time_ad(Val(:disabled), "objective", objective, θ),
            θ -> Zygote.gradient(objective, θ)[1],
            flat_initial_parameters,
            LBFGS(
                alphaguess = LineSearches.InitialStatic(scaled=true, alpha=1),
                linesearch = LineSearches.BackTracking(),
                m=50,
            ),
            Optim.Options(
                x_tol = 1e-4,
                f_tol = 1e-8,
                g_tol = 1e-1,
                f_calls_limit = 1_000,
                g_calls_limit = 1_000,
                iterations = 1_000,
                show_trace = true,
                extended_trace = true,
                store_trace=true,
                show_every=10,
            );
            inplace=false,
        )

        display(training_results)
        println()

        wsave(
            joinpath(
                training_output_dir,
                savename(
                    if experiment_type == "exact"
                        Dict(
                            :experiment_type => experiment_type,
                            :kernel => setting[:kernel_info].name,
                            :start_date => first_date,
                        )
                    else
                        Dict(
                            :experiment_type => experiment_type,
                            :M => setting[:M],
                            :kernel => setting[:kernel_info].name,
                            :start_date => first_date,
                        )
                    end,
                    "jld2",
                ),
            ),
            Dict(
                "setting" => setting,
                "training_results" => training_results,
                "unflatten" => unflatten,
                "test_idx" => test_idx,
            ),
        )
    end
end

# end





#
# Make predictions at locations that we wish to view them, and compute RMSE on data etc.
#

println("Generating predictions.")

mkpath(prediction_output_dir)

let

    include(joinpath(srcdir(), "elevation_data.jl"))

    trained_models = collect_results(training_output_dir);

    for (n, trained_model_) in enumerate(eachrow(trained_models))

        println("n, $(trained_model_.path)")
        println(n)

        # Defensive copy. Something in here is mutating for some unknown reason.
        trained_model = deepcopy(trained_model_);

        # Construct the inputs.
        first_date = trained_model.setting[:start_date]
        first_date_index = only(findall(tmax_dates .== first_date))
        time_inputs_as_dates = tmax_dates[first_date_index:end]
        time_inputs_as_range = RegularSpacing(0.0, 1.0, length(time_inputs_as_dates))
        inputs = RectilinearGrid(space_inputs, time_inputs_as_range);

        # Reconstruct the training data.
        outputs_raw = collect(vec(tmax_data[:, first_date_index:end]))
        test_idx = trained_model.test_idx
        @show typeof(test_idx)
        _, outputs_train_raw = construct_training_outputs(test_idx, outputs_raw)

        # Standardise the training data.
        m = mean(filter(!ismissing, outputs_train_raw))
        s = std(filter(!ismissing, outputs_train_raw))
        outputs_train = (outputs_train_raw .- m) ./ s
        @show typeof(outputs_train)
        println("Total number of data is $(length(outputs_train)).")

        # Build the model.
        unflatten = trained_model.unflatten
        results = trained_model.training_results
        learned_parameters = ParameterHandling.value(unflatten(results.minimizer))

        k = build_kernel(learned_parameters.k)
        f = to_sde(GP(k, GPC()))

        # Get the elevation for each point in the lat-lon grid and return a new data set.
        lat_lon_el_grid = map(lat_lon_grid) do lat_lon
            lat, lon = lat_lon
            elevation = nearest_point(elevation_data_mat, lat, lon) / 100
            return (lat, lon, elevation)
        end

        # Make predictions at entire grid of points at which we might consider predicting.
        x_r_pr = ColVecs(hcat(get_vector.(lat_lon_el_grid)...))
        x_r = space_inputs
        pos = length(time_inputs_as_range) - 300

        # target_times = [1800, 3600];
        # trace_times = map(x -> x.metadata["time"], results.trace);
        # trace_indices_raw = map(x -> findfirst(trace_times .> x), target_times)
        # trace_indices = map(x -> x === nothing ? 1 : x, trace_indices_raw)

        experiment_type = trained_model.setting[:experiment_type]
        (preds_mat, station_preds) = if experiment_type ∈ ["SoD", "SoR"]

            inputs_sor, outputs_train_sor = sor_data(
                inputs, outputs_train, trained_model.setting[:M],
            )
            inputs_space_sor = get_space(inputs_sor)

            # Compute SoR posterior marginals at requested locations and time.
            # For some reason, my exact inference algorithm is broken at the time of
            # writing, but using the approximate inference algorithm with z_r = x_r yields
            # SoR inference.
            σ² = learned_parameters.σ²
            fx = f(inputs_sor, σ²)
            preds = @time approx_posterior_marginals(
                dtc, fx, outputs_train_sor, inputs_space_sor, x_r_pr, pos,
            )
            preds_mat = reshape(preds, N_lons, :)
            preds_mat = Normal.(mean.(preds_mat), std.(preds_mat) .+ sqrt(σ²))

            # Compute posterior marginals at the training data.
            fx = f(inputs_sor, learned_parameters.σ²)
            station_preds = @time approx_posterior_marginals(
                dtc, fx, outputs_train_sor, inputs_space_sor, x_r,
            )
            station_preds = Normal.(mean.(station_preds), std.(station_preds) .+ sqrt(σ²))

            (preds_mat, station_preds)
        elseif experiment_type == "pseudo-point"
            z_r = ColVecs(learned_parameters.Z_r)
            σ² = learned_parameters.σ²
            fx = f(inputs, σ²)

            # Compute approx posterior marginals at requested locations and time.
            preds = @time approx_posterior_marginals(
                dtc, fx, outputs_train, z_r, x_r_pr, pos,
            )
            preds_mat = reshape(preds, N_lons, :)
            preds_mat = Normal.(mean.(preds_mat), std.(preds_mat) .+ sqrt(σ²))

            # Compute approx. posterior marginals at the training data.
            station_preds = @time approx_posterior_marginals(
                dtc, fx, outputs_train, z_r, x_r,
            )
            station_preds = Normal.(mean.(station_preds), std.(station_preds) .+ sqrt(σ²))

            preds_mat, station_preds
        else
            bad_experiment_type_error(experiment_type)
        end

        # Compute the RSMSEs on the training data.
        train_idx = setdiff(eachindex(outputs_raw), test_idx)
        preds_train = station_preds # Test data is missing from training, so doesn't matter.
        preds_test = station_preds[test_idx]
        outputs_test = (outputs_raw[test_idx] .- m) ./ s
        @show size(preds_train), size(preds_test), size(outputs_test), size(outputs_train)
        smse_train = mean(abs2, filter(!ismissing, outputs_train .- mean.(preds_train)))
        smse_test = mean(abs2, filter(!ismissing, outputs_test .- mean.(preds_test)))

        # Compute the Posterior Predictive Log Probs (PPLP).
        present_train_idx = findall(!ismissing, outputs_train)
        present_test_idx = findall(!ismissing, outputs_test)
        @show size(present_train_idx), size(present_test_idx)
        pplp_train = mean(logpdf.(
            preds_train[present_train_idx],
            outputs_train[present_train_idx],
        ))
        pplp_test = mean(logpdf.(
            preds_test[present_test_idx],
            outputs_test[present_test_idx],
        ))

        # Write the results.
        @time wsave(
            joinpath(
                prediction_output_dir,
                savename(
                    if experiment_type == "exact"
                        Dict(
                            :experiment_type => experiment_type,
                            :kernel => trained_model.setting[:kernel_info].name,
                            :start_date => first_date,
                        )
                    else
                        Dict(
                            :experiment_type => experiment_type,
                            :M => trained_model.setting[:M],
                            :kernel => trained_model.setting[:kernel_info].name,
                            :start_date => first_date,
                        )
                    end,
                    "jld2",
                ),
            ),
            Dict(
                "trained_model" => trained_model,
                "preds_train" => preds_train,
                "preds_test" => preds_test,
                "preds_extrapolate" => preds_mat, 
                "x_r_pr" => x_r_pr,
                "lat_lon_pairs" => lat_lon_grid,
                "smse_train" => smse_train,
                "smse_test" => smse_test,
                "pplp_train" => pplp_train,
                "pplp_test" => pplp_test,
            ),
        )
    end

end











#
# Visualise results.
#

# Load shoreline data.
# const shore_data = GeoDatasets.gshhg('f');

struct GridBox{Tlats<:Tuple{<:Real, <:Real}, Tlons<:Tuple{<:Real, <:Real}}
    lats::Tlats
    lons::Tlons
end

const box = GridBox(lat_interval, lon_interval);

function in_box(box::GridBox, lat::Real, lon::Real)
    return box.lats[1] <= lat && box.lats[2] > lat &&
        box.lons[1] <= lon && box.lons[2] > lon
end

let

    results = @time collect_results(prediction_output_dir);

    _, us_states = GADM.get("USA"; children=true);
    state_map = Dict(
        map(
            ((state_name, geom), ) -> state_name => geom,
            zip(us_states.NAME_1, us_states.geom),
        ),
    );

    _, can_states = GADM.get("CAN"; children=true);
    can_state_map = Dict(
        map(
            ((state_name, geom), ) -> state_name => geom,
            zip(can_states.NAME_1, can_states.geom),
        ),
    );

    country_geometries = [
        # GADM.get("AUT"),
        # GADM.get("FIN"),
        # GADM.get("DNK"),
        # GADM.get("ITA"),
        # GADM.get("LIE"),
        # GADM.get("NLD"),
        # GADM.get("NOR"),
        # GADM.get("POL"),
        # GADM.get("SWE"),
        # GADM.get("CHE"),
        # GADM.get("DEU"),
        (geom=state_map["Washington"],),
        (geom=can_state_map["British Columbia"],),
    ];

    # Compute some additional columns that will be used later.
    let
        optim_results = map(x -> x.training_results, results.trained_model)
        results.iters = Optim.iterations.(optim_results)
        results.timings = Optim.time_run.(optim_results)
        results.time_per_iter = results.timings ./ results.iters
        results.trace = getproperty.(optim_results, :trace)
    end;

    gr()

    for result in eachrow(results)

        preds_mat = result.preds_extrapolate
        trained_model = deepcopy(result.trained_model)

        # Construct the inputs.
        first_date = trained_model.setting[:start_date]
        first_date_index = only(findall(tmax_dates .== first_date))
        time_inputs_as_dates = tmax_dates[first_date_index:end]
        time_inputs_as_range = RegularSpacing(0.0, 1.0, length(time_inputs_as_dates))
        inputs = RectilinearGrid(space_inputs, time_inputs_as_range);

        # Get the pseudo-point locations.
        unflatten = trained_model.unflatten
        training_results = trained_model.training_results
        learned_parameters = ParameterHandling.value(unflatten(training_results.minimizer))

        @show learned_parameters.k
        @show result.smse_train, result.smse_test, result.pplp_train, result.pplp_test
        @show training_results.minimum

        experiment_type = trained_model.setting[:experiment_type]
        colour_scheme = :viridis;

        # Plot the posterior mean.
        let
            plt = plot(
                xlims=lon_interval,
                ylims=lat_interval,
                xticks=nothing,
                yticks=nothing,
                legend=:bottomleft,
            )

            heatmap!(
                plt, uk_lons_range, uk_lats_range, mean.(preds_mat);
                label="", c=cgrad(colour_scheme), right_margin=6Plots.mm,
            )

            # Plot national borders in the relevant area of Europe.
            for country in country_geometries
                plot!(plt, country.geom; fillcolor=nothing, label="")
            end

            scatter!(
                plt, getindex.(space_inputs, 2), getindex.(space_inputs, 1);
                label="Stations Locations",
                # color=:blue,
                markersize=5,
                marker=:square,
                color=colorant"#A95AA1",
            )

            if experiment_type == "pseudo-point"
                z_r = ColVecs(learned_parameters.Z_r)
                scatter!(
                    plt, getindex.(z_r, 2), getindex.(z_r, 1);
                    label="z",
                    # color=:red,
                    markersize=4,
                    color=colorant"#F5793A", # Orange
                )
            end

            # for (lons, lats) in shore_data

            #     # Filter out stuff outside of the box.
            #     points_in_box = filter(
            #         ((lat, lon), ) -> in_box(box, lat, lon),
            #         collect(zip(lats, lons)),
            #     )

            #     if length(points_in_box) > 0
            #         lats_f = first.(points_in_box)
            #         lons_f = last.(points_in_box)

            #         # Plot stuff that's in the box.
            #         plot!(plt, lons_f, lats_f; label="", color=:black, alpha=1.0)
            #     end
            # end

            mkpath(results_dir)
            savefig(
                plt,
                joinpath(
                    results_dir,
                    savename(
                        if experiment_type == "SoD"
                            Dict(
                                :experiment_type => experiment_type,
                                :M => trained_model.setting[:M],
                                :kernel => trained_model.setting[:kernel_info].name,
                                :start_date => first_date,
                            )
                        elseif experiment_type == "pseudo-point"
                            Dict(
                                :experiment_type => experiment_type,
                                :M => trained_model.setting[:M],
                                :kernel => trained_model.setting[:kernel_info].name,
                                :start_date => first_date,
                            )
                        else
                            bad_experiment_type_error(experiment_type)
                        end,
                        "posterior_mean.png",
                    ),
                )
            )
        end

        # Plot the posterior standard deviation.
        let
            plt = plot(
                xlims=lon_interval,
                ylims=lat_interval,
                xticks=nothing,
                yticks=nothing,
            )

            heatmap!(
                uk_lons_range, uk_lats_range, std.(preds_mat);
                label="", c=colour_scheme, right_margin=3Plots.mm,
            )

            # Plot national borders in the relevant area of Europe.
            for country in country_geometries
                plot!(plt, country.geom; fillcolor=nothing, label="", linecolor=:black)
            end

            scatter!(
                plt, getindex.(space_inputs, 2), getindex.(space_inputs, 1);
                label="Stations Locations",
                markersize=5,
                marker=:square,
                color=colorant"#A95AA1", # Purple
                # color=:blue,
            )

            if experiment_type == "pseudo-point"
                z_r = ColVecs(learned_parameters.Z_r)
                scatter!(
                    plt, getindex.(z_r, 2), getindex.(z_r, 1);
                    label="z",
                    # color=:red,
                    markersize=4,
                    color=colorant"#F5793A", # Orange
                )
            end

            # for (lons, lats) in shore_data

            #     # Filter out stuff outside of the box.
            #     points_in_box = filter(
            #         ((lat, lon), ) -> in_box(box, lat, lon),
            #         collect(zip(lats, lons)),
            #     )

            #     if length(points_in_box) > 0
            #         lats_f = first.(points_in_box)
            #         lons_f = last.(points_in_box)

            #         # Plot stuff that's in the box.
            #         plot!(plt, lons_f, lats_f; label="", color=:black, alpha=1.0)
            #     end
            # end

            mkpath(results_dir)
            savefig(
                plt,
                joinpath(
                    results_dir,
                    savename(
                        if experiment_type == "SoD"
                            Dict(
                                :experiment_type => experiment_type,
                                :M => trained_model.setting[:M],
                                :kernel => trained_model.setting[:kernel_info].name,
                                :start_date => first_date,
                            )
                        elseif experiment_type == "pseudo-point"
                            Dict(
                                :experiment_type => experiment_type,
                                :M => trained_model.setting[:M],
                                :kernel => trained_model.setting[:kernel_info].name,
                                :start_date => first_date,
                            )
                        else
                            bad_experiment_type_error(experiment_type)
                        end,
                        "posterior_std.png",
                    ),
                )
            )
        end
    end

    #
    # Build TexTable to pretty print the results.
    #

    my_format(x) = format(round(x; sigdigits=3); stripzeros=true)

    name_map = Dict(
        "SoD" => "SoD",
        "exact" => "Exact",
        "pseudo-point" => "P-P",
        # "separable" => "Separable",
        # "sum-separable" => "Sum-Separable",
        "separable-Matern-52" => "Separable",
        "sum-separable-Matern52" => "Sum-Separable",
    );

    # Group everything by date.
    results_df_unsorted = DataFrame(
        start_date = map(x -> x.setting[:start_date], results.trained_model),
        lml_approx = map(
            x -> my_format(x.training_results.minimum),
            results.trained_model,
        ),
        rsmse_test = sqrt.(results.smse_test),
        npplp_test = -results.pplp_test,
        rsmse_test_string = map(my_format, sqrt.(results.smse_test)),
        npplp_test_string = map(my_format, -results.pplp_test),
        row_keys = map(
            x -> begin
                M = x.setting[:M]
                experiment_name = name_map[x.setting[:experiment_type]] * "($M)"
                kernel_name = name_map[x.setting[:kernel_info].name]
                return experiment_name * " - " * kernel_name
            end,
            results.trained_model,
        ),
        M = map(x -> x.setting[:M], results.trained_model),
        experiment_type = map(x -> x.setting[:experiment_type], results.trained_model),
        kernel_name = map(x -> x.setting[:kernel_info].name, results.trained_model),
        run_time = map(x -> Optim.time_run(x.training_results), results.trained_model),
        time_per_iter = results.time_per_iter,
    );

    results_df = sort(results_df_unsorted, :M);

    for group in groupby(results_df, :start_date)

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



    pgfplotsx();

    # Wessel's colour pallete that is friendly to those suffering colour-blindness.
    colours = [
        colorant"#4BA6FB", # Modified blue (original was #85C0F9)
        colorant"#A95AA1", # Pink
        colorant"#000000", # Black
        colorant"#F5793A", # Orange
    ]
    # colours = [:black, :red, :blue, :cyan]
    shapes = [:utriangle, :diamond, :square, :circle, :cross]

    # Generate plot. This should be roughly half a page width.
    width_inches = 2.8
    width_px = width_inches * Plots.DPI
    height_px = width_px / 2


    # Build speed-accuracy graphs for each start date in the results.
    for date_group in groupby(results_df, :start_date)
        rsmse_plot = plot(
            size=(width_px, height_px),
            xscale=:log10,
            legend=:topright,
            background_color_legend=RGBA(1, 1, 1, 0),
            foreground_color_legend=RGBA(1, 1, 1, 0),
        )
        npplp_plot = plot(
            size=(width_px, height_px),
            xscale=:log10,
            legend=nothing,
            background_color_legend=RGBA(1, 1, 1, 0),
            foreground_color_legend=RGBA(1, 1, 1, 0),
        )
        start_date = only(unique(date_group.start_date))
        for (n, group) in enumerate(groupby(date_group, [:experiment_type, :kernel_name]))
            experiment_type = only(unique(group.experiment_type))
            kernel_name = only(unique(group.kernel_name))
            label = name_map[experiment_type] * ", " * name_map[kernel_name]
            plot!(
                rsmse_plot, group.run_time, group.rsmse_test;
                xlabel="Run Time (s)",
                ylabel="RSMSE",
                label=label,
                color=colours[n],
                marker=shapes[n],
                markerstrokealpha=0,
                markersize=3,
                markercolor=colours[n],
                tickfontsize=nothing,
                legendfontsize=nothing,
            )
            plot!(
                npplp_plot, group.run_time, group.npplp_test,
                xlabel="Run Time (s)",
                ylabel="NPPLP",
                label=label,
                color=colours[n],
                marker=shapes[n],
                markerstrokealpha=0,
                markersize=3,
                markercolor=colours[n],
                tickfontsize=nothing,
                legendfontsize=nothing,
            )
        end

        savefig(rsmse_plot, joinpath(results_dir, "rsmse-tradeoff_date-$start_date.tikz"))
        savefig(npplp_plot, joinpath(results_dir, "pplp-tradeoff_date-$start_date.tikz"))
    end

    # Build time-per-iteration graphs.
    let
        for date_group in groupby(results_df, :start_date)
            start_date = only(unique(date_group.start_date))
            plt = plot(
                yscale=:log10,
                size=(width_px, height_px),
                legend=:bottomright,
                background_color_legend=RGBA(1, 1, 1, 0),
                foreground_color_legend=RGBA(1, 1, 1, 0),
            )
            for (n, group) in enumerate(groupby(date_group, [:experiment_type, :kernel_name]))
                experiment_type = only(unique(group.experiment_type))
                kernel_name = only(unique(group.kernel_name))
                label = 
                plot!(
                    plt, group.M, group.time_per_iter;
                    xlabel="M",
                    ylabel="Time Per Iteration",
                    label=name_map[experiment_type] * ", " * name_map[kernel_name],
                    color=colours[n],
                    marker=shapes[n],
                    markerstrokealpha=0,
                    markersize=2,
                    markercolor=colours[n],
                    tickfontsize=nothing,
                    legendfontsize=6,
                )
            end
            savefig(plt, joinpath(results_dir, "time-per-iteration-$start_date.tikz"))
        end
    end

end
