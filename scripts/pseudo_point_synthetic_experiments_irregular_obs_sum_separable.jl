#
# This script was used to generate results for the camera-ready AABI submission.
#

println("pseudo_point_synthetic_experiments_irregular_obs_sum_separable")

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
            default = 10
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


args = parse_commandline()

using BenchmarkTools
using DrWatson
using LinearAlgebra
using Random
using Stheno
using TemporalGPs

BLAS.set_num_threads(args["BLAS-threads"])

const exp_dir_name = "pseudo_point_scaling_irregular_obs_sum_separable"



#
# Specify some functionality required to run the experiments.
#

function separable_kernel(σ::Real, λ_space::Real, λ_time::Real)
    k_space = σ^2 * stretch(EQ(), λ_space)
    k_time = stretch(Matern52(), λ_time)
    return TemporalGPs.Separable(k_space, k_time)
end

function sum_separable_kernel(σ::Real, λ_space::Real, λ_time::Real)
    k1 = separable_kernel(σ, λ_space, λ_time)
    k2 = separable_kernel(1.1 * σ, 0.3 * λ_space, 1.05 * λ_time)
    return k1 + k2
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
                    name = "sum-separable",
                    θ = (
                        σ = 0.92,
                        λ_space = 1.1,
                        λ_time = 0.83,
                        S_min = 0.1,
                    ),
                    ctor = sum_separable_kernel,
                ),
            ],

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
        @show model_name
        θ = setting[:model].θ
        ctor = setting[:model].ctor
        T = setting[:Ns_time]
        N_space_per_time = setting[:N_space_per_time]
        M_per_time = setting[:Ms_per_time]

        # Specify input locations. We place these irregularly through time.
        t = sort(rand(rng, T)) .* 50
        r = [rand(rng, N_space_per_time) * 10 for _ in 1:T]
        x = TemporalGPs.RegularInTime(t, r)
        x_naive = collect(x)

        # Specify inducing inputs. They are regularly-spaced between 0 and 10.
        z_r = collect(range(0.0, 10.0; length=M_per_time))
        z_naive = collect(TemporalGPs.RectilinearGrid(z_r, t))

        # Specify observation noise.
        S = [Diagonal(rand(rng, N_space_per_time) .+ θ.S_min) for _ in 1:T]
        S_naive = vcat(diag.(S)...)

        # Generate data exactly / naively.
        k = ctor(θ.σ, θ.λ_space, θ.λ_time)
        f = GP(k, GPC())

        # Hack together generation.
        N_per_slice = 5_000
        ys_naive = map(1:div(length(x_naive), N_per_slice)) do j
            idx = N_per_slice * (j - 1) + 1:min(N_per_slice * j, length(x_naive))
            return rand(rng, f(x_naive[idx], S_naive[idx]))
        end
        ys_naive = if rem(length(x_naive), N_per_slice) == 0
            ys_naive
        else
            idx = length(x_naive) - rem(length(x_naive), N_per_slice) + 1:length(x_naive)
            vcat(
                ys_naive,
                rand(rng, f(x_naive[idx], S_naive[idx]))
            )
        end
        y_naive = vcat(ys_naive...)
        y = TemporalGPs.match_to(y_naive, x)

        display((length(x_naive), length(y_naive)))
        println()
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

        # State-Space ELBO
        println("Using the state-space ELBO")
        if M_per_time > 0
            wsave(
                joinpath(
                    data_dir,
                    savename(
                        Dict(
                            :k => model_name,
                            :N_space_per_time => N_space_per_time,
                            :N_time => T,
                            :M_per_time => M_per_time,
                            :method => :state_space_elbo,
                        ),
                        "bson",
                    ),
                ),
                Dict(
                    :setting => setting,
                    :lml_result => elbo(to_sde(f)(x, S_naive), y_naive, z_r),
                    :timing => (@benchmark elbo(to_sde($f)($x, $S_naive), $y_naive, $z_r)),
                    :method => :state_space_elbo,
                ),
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
        :naive_elbo => :circle,
        :state_space_elbo => :circle,
    )

    # Each inference method gets its own colour.
    colours = Dict(
        :exact => :black,
        :naive_elbo => :red,
        :state_space_elbo => :blue,
    )

    method_name_map = Dict(
        :exact => "exact",
        :naive_elbo => "elbo (naive)",
        :state_space_elbo => "elbo"
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

        # width_pts = 246.0  # or any other value
        # inches_per_points = 1.0 / 72.27
        # width_inches = width_pts *inches_per_points
        width_inches = 2.8
        width_px = width_inches * Plots.DPI  # or  width_inches*DPI
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
        )
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
                savefig(
                    plt_lml,
                    joinpath(save_dir, "irregular_lml_plot_sum_separable.tikz"),
                )
                savefig(
                    plt_comp_time,
                    joinpath(save_dir, "irregular_timing_plot_sum_separable.tikz"),
                )
                savefig(
                    plot(plt_lml, plt_comp_time; layout=(2, 1)),
                    joinpath(save_dir, "irregular_lml_timing_plot_sum_separable.tikz"),
                )
            end,
            save_dirs,
        )
    end
end
