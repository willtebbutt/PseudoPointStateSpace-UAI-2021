#
# This is the script to generate the result for my AABI submission.
# It will evolve over the next few days as I populate the results.
#

println("pseudo_point_synthetic_experiments_regular_obs")

using Revise

using Pkg
Pkg.activate(".");
Pkg.resolve();
Pkg.instantiate();

using ArgParse
using CSV

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
            default = false
        "--N-space-per-time"
            help = "number of points to consider in space per time point"
            arg_type = Int
            default = 50
        "--N-missing-per-time"
            help = "number of points to consider in space per time point"
            arg_type = Int
            default = 5
        "--N-time"
            help = "number of time points to use -- expects multiple settings"
            arg_type = Vector{Int}
            default = [
                10, 25, 50,
                100, 250, 500,
                1_000, 2_500, 5_000,
                10_000, 25_000, 50_000,
                100_000,
            ]
        "--M-per-time"
            help = "number of pseudo-points per time point -- expects multiple settings. 0 means exact"
            arg_type = Vector{Int}
            default = [0, 5, 10, 20]
        "--BLAS-threads"
            help = "number of threads to let BLAS use"
            arg_type = Int
            default = 1
    end
    return parse_args(s)
end

            # default = [10, 25, 50, 100, 250, 500, 1_000, 2_500, 5_000, 10_000, 25_000, 50_000, 100_000]

args = parse_commandline()

using BenchmarkTools
using DrWatson
using LinearAlgebra
using Random
using Stheno
using TemporalGPs

BLAS.set_num_threads(args["BLAS-threads"])

const exp_dir_name = "pseudo_point_scaling_regular_obs"


#
# Specify some functionality required to run the experiments.
#

function separable_kernel(σ::Real, λ_space::Real, λ_time::Real)
    k_space = σ^2 * stretch(EQ(), λ_space)
    k_time = stretch(Matern52(), λ_time)
    return TemporalGPs.Separable(k_space, k_time)
end



#
# Specify experiments to run.
#

# Ensure that output directory for results exists.
const data_dir = joinpath(datadir(), exp_dir_name)
mkpath(data_dir)

# Ensure that directory for experiment specification exists.
mkpath(joinpath(data_dir, "meta"))

wsave(
    joinpath(data_dir, "meta", "settings.bson"),
    :settings => dict_list(
        Dict(
            # The number of kernels to include in the summation.
            :model => [
                (
                    name = "separable",
                    θ = (
                        σ = 0.92,
                        λ_space = 1.1,
                        λ_time = 0.83,
                        S_min = 0.1,
                    ),
                ),
            ],

            :N_missing_per_time => args["N-missing-per-time"],

            # The total number of points in space to consider.
            :N_space_per_time => args["N-space-per-time"],

            # The total number of points in time to consider.
            :Ns_time => args["N-time"],

            # The total number of pseudo-points per time point.
            :Ms_per_time => args["M-per-time"],
        ),
    ),
)



#
# Run experiments.
#

if args["run-experiments"]

println("Running experiments")

let
    settings = load(joinpath(data_dir, "meta", "settings.bson"))[:settings]

    for (n, setting) in enumerate(settings)
        println("Starting experiment $n / $(length(settings))")

        # Set seed for reproducibility.
        rng = MersenneTwister(123456)

        # Extract experiment parameters.
        model_name = setting[:model].name
        θ = setting[:model].θ
        T = setting[:Ns_time]
        N_space_per_time = setting[:N_space_per_time]
        N_missing_per_time = setting[:N_missing_per_time]
        M_per_time = setting[:Ms_per_time]

        # Specify input locations. We place these irregularly through time.
        t = sort(rand(rng, T)) .* 50
        t = RegularSpacing(0.0, 50.0, T) # THIS SHOULD CHANGE TO SOMETHING SENSIBLE!
        r = collect(range(0.0, 5.0; length=N_space_per_time))
        x = TemporalGPs.RectilinearGrid(r, t);
        x_naive_full = collect(x)

        # Specify inducing inputs. They are regularly-spaced between 0 and 10.
        z_r = collect(range(0.0, 10.0; length=M_per_time))
        z_naive = collect(TemporalGPs.RectilinearGrid(z_r, t))

        # Specify observation noise.
        S_full = [Diagonal(rand(rng, N_space_per_time) .+ θ.S_min) for _ in 1:T]

        # Generate entire data set exactly / naively.
        k = separable_kernel(θ.σ, θ.λ_space, θ.λ_time)
        f = GP(k, GPC())

        # Generate entire data set using state-space model.
        fx = to_sde(f)(x, Diagonal(vcat(diag.(S_full)...)))
        y = rand(rng, TemporalGPs.build_lgssm(fx))

        # Remove a couple of data points at each time point directly from
        # x_naive and y_naive.
        missing_idxs = map(t -> randperm(rng, N_space_per_time)[1:N_missing_per_time], 1:T)
        x_naive_vecs = map(1:T) do t
            x_naive_full[(t - 1) * N_space_per_time + 1:t * N_space_per_time]
        end
        x_naive_missing_vecs = map(eachindex(missing_idxs)) do n
            present_idx = deleteat!(collect(1:N_space_per_time), sort(missing_idxs[n]))
            return x_naive_vecs[n][present_idx]
        end
        y_naive_missing_vecs = map(eachindex(missing_idxs)) do n
            present_idx = deleteat!(collect(1:N_space_per_time), sort(missing_idxs[n]))
            return y[n][present_idx]
        end
        S_naive_missing_vecs = map(eachindex(missing_idxs)) do n
            present_idx = deleteat!(collect(1:N_space_per_time), sort(missing_idxs[n]))
            return diag(S_full[n])[present_idx]
        end

        x_naive = vcat(x_naive_missing_vecs...)
        y_naive = vcat(y_naive_missing_vecs...)
        S_naive = vcat(S_naive_missing_vecs...)

        # Modify S where there is missing data.
        S = map(eachindex(missing_idxs)) do n
            s = copy(diag(S_full[n]))
            s[missing_idxs[n]] .= 1e9
            return Diagonal(s)
        end
        S_full_diag = Diagonal(vcat(diag.(S)...))
        y_long = vcat(y...)

        display((length(x_naive), length(y_naive)))
        println()

        # Naive log marginal likelihood.
        if M_per_time == 0 && length(y_naive) <= 10_000
            fx_naive = f(x_naive, S_naive)
            println("Performing exact inference.")
            wsave(
                joinpath(
                    data_dir,
                    savename(
                        Dict(
                            :k => model_name,
                            :N_space_per_time => N_space_per_time,
                            :N_missing_per_time => N_missing_per_time,
                            :N_time => T,
                            :M_per_time => M_per_time,
                            :method => :exact,
                        ),
                        "bson",
                    ),
                ),
                Dict(
                    :setting => setting,
                    :lml_result => logpdf(fx_naive, y_naive),
                    :timing => (@benchmark logpdf($fx_naive, $y_naive)),
                    :method => :exact,
                ),
            )
        end

        # State-Space log marginal likelihood.
        if M_per_time == 0
            println("State-space log marginal likelihood.")

            # This accounts for some extra terms associated with missing data.
            fudge_factor = 0.5 * T * N_missing_per_time * log(2π * 1e9)

            wsave(
                joinpath(
                    data_dir,
                    savename(
                        Dict(
                            :k => model_name,
                            :N_space_per_time => N_space_per_time,
                            :N_missing_per_time => N_missing_per_time,
                            :N_time => T,
                            :M_per_time => M_per_time,
                            :method => :exact_state_space,
                        ),
                        "bson",
                    ),
                ),
                Dict(
                    :setting => setting,
                    :lml_result => logpdf(to_sde(f)(x, S_full_diag), y_long) + fudge_factor,
                    :timing => (@benchmark logpdf(to_sde($f)($x, $S_full_diag), $y_long)),
                    :method => :exact_state_space,
                ),
            )
        end

        # Naive ELBO computation. Really just here for verification purposes.
        if M_per_time > 0 && M_per_time * T < 2_000
            fx_naive = f(x_naive, S_naive)
            println("Using the naive ELBO")
            wsave(
                joinpath(
                    data_dir,
                    savename(
                        Dict(
                            :k => model_name,
                            :N_space_per_time => N_space_per_time,
                            :N_missing_per_time => N_missing_per_time,
                            :N_time => T,
                            :M_per_time => M_per_time,
                            :method => :naive_elbo,
                        ),
                        "bson",
                    ),
                ),
                Dict(
                    :setting => setting,
                    :lml_result => elbo(fx_naive, y_naive, f(z_naive, 1e-9)),
                    :timing => (@benchmark elbo($fx_naive, $y_naive, $f($z_naive, 1e-9))),
                    :method => :naive_elbo,
                ),
            )
        end

        # State-Space ELBO.
        if M_per_time > 0
            println("Using the state-space ELBO")

            # This accounts for some extra terms associated with missing data.
            fudge_factor = 0.5 * T * N_missing_per_time * log(2π * 1e9)

            wsave(
                joinpath(
                    data_dir,
                    savename(
                        Dict(
                            :k => model_name,
                            :N_space_per_time => N_space_per_time,
                            :N_missing_per_time => N_missing_per_time,
                            :N_time => T,
                            :M_per_time => M_per_time,
                            :method => :state_space_elbo,
                        ),
                        "bson",
                    ),
                ),
                Dict(
                    :setting => setting,
                    :lml_result => elbo(to_sde(f)(x, S_full_diag), y_long, z_r) + fudge_factor,
                    :timing => (@benchmark elbo(to_sde($f)($x, $S_full_diag), $y_long, $z_r)),
                    :method => :state_space_elbo,
                ),
                # Dict(
                #     :setting => setting,
                #     :lml_result => elbo(k, x, y, z_r, S) + fudge_factor,
                #     :timing => (@benchmark elbo($k, $x, $y, $z_r, $S)),
                #     :method => :state_space_elbo,
                # ),
            )
        end
    end
end

end # run-experiments



#
# Summarise results. Generate output (plots and tables).
#

using ColorTypes, DataFrames, LaTeXStrings, PGFPlotsX, Plots, TexTables
pgfplotsx();

Plots.default(
    tickfontsize=nothing,
    legendfontsize=nothing,
    guidefontsize=nothing,
);

# Use default document font-size if nothing provided.
function Plots.pgfx_font(fontsize::Nothing, thickness_scaling = 1, font = "\\selectfont")
    return string("{", font, "}")
end

# Ensure that results directory exists.
const results_dir = joinpath(projectdir(), "results", exp_dir_name)
mkpath(results_dir)

let
    # Load the results and augment columns with extra info from `settings` column.
    results_unsorted = DataFrame(
        map(eachrow(collect_results(data_dir))) do row
            merge(
                copy(row),
                (
                    N_time = row.setting[:Ns_time],
                    M_per_time = row.setting[:Ms_per_time],
                ),
            )
        end
    )

    results = sort(results_unsorted, [:method, :M_per_time]; rev=true)

    # Each inference method gets its own marker shape.
    marker_shapes = Dict(
        :exact => :square,
        :exact_state_space => :utriangle,
        :naive_elbo => :circle,
        :state_space_elbo => :circle,
    )

    # Each inference method gets its own colour.
    colours = Dict(
        :exact => :black,
        :exact_state_space => :green,
        :naive_elbo => :red,
        :state_space_elbo => :blue,
    )

    method_name_map = Dict(
        :exact => "exact",
        :exact_state_space => "exact (sde)",
        :naive_elbo => "elbo (naive)",
        :state_space_elbo => "elbo",
    )

    # Each M gets its own line style. This means that the maximum value of `M` is 5.
    # We really shouldn't be using that many anyway though...
    line_styles = [:solid, :solid, :dash, :dot, :dashdot, :dashdotdot]
    Ms_per_time = sort(unique(results.M_per_time))
    line_styles_map = Dict(
        M => style for (M, style) in zip(Ms_per_time, line_styles[eachindex(Ms_per_time)])
    )

    # Locations to save the results to.
    save_dirs = [
        results_dir,
        "/Users/willtebbutt/Dropbox/University/PhD/first_year_report/chapter_3/pseudo_points"
    ]

    # Ensure that save_dirs exist.
    foreach(mkpath, save_dirs)

    let
        width_inches = 2.8
        width_px = width_inches * Plots.DPI
        height_px = width_px

        plt_lml = plot(
            xlabel="T",
            ylabel="ELBO / T",
            legend=:bottomright,
            legendfontsize=9,
            background_color_legend=RGBA(1, 1, 1, 0.6),
            size=(width_px, height_px),
            grid = (1, 0.4),
        );
        plt_comp_time = plot(
            xlabel="T",
            ylabel="Compute Time (s)",
            legend=:bottomright,
            legendfontsize=9,
            background_color_legend=RGBA(1, 1, 1, 0.6),
            size=(width_px, height_px),
            grid = (1, 0.4),
        );
        for group in groupby(results, [:method, :M_per_time])
            method_name = only(unique(group.method))
            M_per_time = only(unique(group.M_per_time))
            if method_name == :exact && M_per_time != first(Ms_per_time)
                continue
            end
            if method_name == :naive_elbo
                continue
            end

            idx = sortperm(group.N_time)
            label = if method_name == :exact
                method_name_map[method_name]
            elseif method_name == :exact_state_space
                method_name_map[method_name]
            else
                L"M_\tau=%$M_per_time"
            end

            Ts = group.N_time[idx]

            plot!(plt_lml, Ts, group.lml_result[idx] ./Ts ;
                markershape=marker_shapes[method_name],
                color=colours[method_name],
                linestyle=line_styles_map[M_per_time],
                markerstrokewidth=0,
                markersize=3,
                label=label,
                xscale=:log10,
            )

            min_time = getfield.(minimum.(group.timing), :time)[idx] ./ 1e9
            plot!(plt_comp_time, Ts, min_time;
                markershape=marker_shapes[method_name],
                color=colours[method_name],
                linestyle=line_styles_map[M_per_time],
                markersize=3,
                markerstrokewidth=0,
                label=label,
                yscale=:log10,
                xscale=:log10,
            )
        end
        foreach(
            save_dir -> begin
                savefig(plt_lml, joinpath(save_dir, "regular_lml_plot.tikz"))
                savefig(plt_comp_time, joinpath(save_dir, "regular_timing_plot.tikz"))
                savefig(
                    plot(plt_lml, plt_comp_time; layout=(2, 1)),
                    joinpath(save_dir, "regular_lml_timing_plot.tikz"),
                )
            end,
            save_dirs,
        )
    end
end
