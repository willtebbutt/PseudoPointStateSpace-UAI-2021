#
# This is a script to generate some results for the AABI-2020 submission. It will produce
# a couple of graphs that will both provide an explanation of the experimental set up, and
# show an approximate posterior mean (because why not :shrug:).
#
# This uses my reduced scripting format that only includes a single experiment block /
# results analysis / plotting block. There's no need for the usual first block because there
# are only two things generated. It's a bit of a mush, but I was pushed for time.
#

println("pseudo_point_synthetic_posterior_plots")

using Revise

using Pkg
Pkg.activate(".");
Pkg.resolve();
Pkg.instantiate();

using ColorTypes
using DataFrames
using BenchmarkTools
using DrWatson
using FileIO
using ImageMagick
using Images
using LaTeXStrings
using LinearAlgebra
using PGFPlotsX
using Plots
using Random
using Stheno
using TemporalGPs
using TexTables

# Configure Plots.
pgfplotsx();
Plots.default(
    tickfontsize=nothing,
    legendfontsize=nothing,
);

# Use default document font-size if nothing provided. This is a hack.
function Plots.pgfx_font(fontsize::Nothing, thickness_scaling = 1, font = "\\selectfont")
    return string("{", font, "}")
end

BLAS.set_num_threads(1)

const exp_dir_name = "pseudo_point_synthetic_posterior_plots"

# A kernel used throughout the experiments.
function separable_kernel(σ::Real, λ_space::Real, λ_time::Real)
    k_space = σ^2 * stretch(EQ(), λ_space)
    k_time = stretch(Matern52(), λ_time)
    return TemporalGPs.Separable(k_space, k_time)
end

# Ensure that output directory for results exists.
const data_dir = joinpath(datadir(), exp_dir_name)
mkpath(data_dir)
mkpath(joinpath(data_dir, "meta"))

# Ensure that results directory exists.
const results_dir = joinpath(projectdir(), "results", exp_dir_name)
mkpath(results_dir)

# Locations to save the plots to.
save_dirs = [
    results_dir,
    "/Users/willtebbutt/Dropbox/University/PhD/first_year_report/chapter_3/pseudo_points"
]
foreach(mkpath, save_dirs)


# Run experiments.
let

    let
        # Set seed for reproducibility.
        rng = MersenneTwister(123456)

        # Extract experiment parameters.
        θ = (σ = 0.92, λ_space = 1.1, λ_time = 0.83, S_min = 0.1)
        T = 100
        N_space_per_time = 50
        N_missing_per_time = 5
        M_per_time = 10

        # Specify input locations. We place these irregularly through time.
        t = sort(rand(rng, T)) .* 50
        t = RegularSpacing(0.0, 0.2, T)
        r = collect(range(0.0, 5.0; length=N_space_per_time))
        r_pred = collect(range(minimum(r), maximum(r); length=100))
        x = TemporalGPs.RectilinearGrid(r, t);
        x_naive_full = collect(x)

        # Specify inducing inputs. They are regularly-spaced between 0 and 10.
        z_r = collect(range(0.0, maximum(r); length=M_per_time))
        z_naive = collect(TemporalGPs.RectilinearGrid(z_r, t))

        # Specify observation noise.
        S_full = [Diagonal(rand(rng, N_space_per_time) .+ θ.S_min) for _ in 1:T]
        S_full_diag = Diagonal(vcat(diag.(S_full)...))

        # Generate entire data set exactly / naively.
        k = separable_kernel(θ.σ, θ.λ_space, θ.λ_time)
        f = GP(k, GPC())
        fx = to_sde(f)(x, S_full_diag)

        # Generate entire data set using state-space model.
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
        S_diag = vcat(diag.(S)...)

        # Perform approximate inference using the state-space approach.
        post_marginals_sde = TemporalGPs.approx_posterior_marginals(
            dtc, to_sde(f)(x, S_diag), vcat(y...), z_r, r_pred,
        )

        # Perform approximate inference naively.
        fx = f(x_naive, S_naive)
        f_post_approx = f | Stheno.PseudoObs(fx ← y_naive, f(z_naive, 1e-9))
        x_pred_naive = collect(TemporalGPs.RectilinearGrid(r_pred, t))
        post_marginals_naive = marginals(f_post_approx(x_pred_naive))

        approx_post_mean_sde = reshape(mean.(post_marginals_sde), length(r_pred), :);
        approx_post_mean_naive = reshape(mean.(post_marginals_naive), length(r_pred), :);
        @assert maximum(abs.(approx_post_mean_sde .- approx_post_mean_naive)) < 1e-7

        approx_post_std_sde = reshape(std.(post_marginals_sde), length(r_pred), :);
        approx_post_std_naive = reshape(std.(post_marginals_naive), length(r_pred), :);
        @assert maximum(abs.(approx_post_std_sde .- approx_post_std_naive)) < 1e-7

        # Generate plot. This should be roughly half a page width.
        width_inches = 2.8
        width_px = width_inches * Plots.DPI
        height_px = width_px

        plt_post = plot(
            xlabel="time",
            ylabel="space",
            legend=:bottomright,
            legendfontsize=9,
            background_color_legend=RGBA(1, 1, 1, 1),
            size=(width_px, height_px),
            grid = (0, 0),
            colorbar=nothing,
            xlims=(minimum(t), maximum(t)),
            ylim=(minimum(r_pred), maximum(r_pred)),
        );

        scatter!(plt_post, fill(0.0, length(z_r)), z_r;
            color=:black,
            label="Pseudo-Inputs",
        );
        scatter!(plt_post, last.(x_naive), first.(x_naive);
            markersize=0.1,
            color=:black,
            label="Observation Inputs",
        )
        savefig(plt_post, joinpath(first(save_dirs), "post_rectlinear.tikz"));

        # This image will have to be manually added onto the plot after it's been made.
        mean_min = minimum(approx_post_mean_sde);
        mean_max = maximum(approx_post_mean_sde);
        img = colorview(Gray, (approx_post_mean_sde .- mean_min) ./ (mean_max - mean_min));
        FileIO.save(joinpath(first(save_dirs), "post_rectlinear.jpg"), img);
    end

    let
        # Set seed for reproducibility.
        rng = MersenneTwister(123456)

        # Extract experiment parameters.
        θ = (σ = 0.92, λ_space = 1.1, λ_time = 0.83, S_min = 0.1)
        T = 100
        N_space_per_time = 10
        M_per_time = 10

        # Specify input locations. We place these irregularly through time.
        t = sort(rand(rng, T)) .* 50;
        t = RegularSpacing(0.0, 0.3, T);
        r = [rand(rng, N_space_per_time) * 10 for _ in 1:T];
        r_pred = collect(range(0.0, 10.0; length=100))
        x = TemporalGPs.RegularInTime(t, r);
        x_naive = collect(x);

        # Specify inducing inputs. They are regularly-spaced between 0 and 10.
        z_r = collect(range(0.0, 10.0; length=M_per_time));
        z_naive = collect(TemporalGPs.RectilinearGrid(z_r, t));

        # Specify observation noise.
        S = [Diagonal(rand(rng, N_space_per_time) .+ θ.S_min) for _ in 1:T];
        S_naive = vcat(diag.(S)...);

        # Generate data exactly / naively.
        k = separable_kernel(θ.σ, θ.λ_space, θ.λ_time);
        f = GP(k, GPC());

        # Hack together generation.
        N_per_slice = 5_000;
        ys_naive = map(1:div(length(x_naive), N_per_slice)) do j
            idx = N_per_slice * (j - 1) + 1:min(N_per_slice * j, length(x_naive))
            return rand(rng, f(x_naive[idx], S_naive[idx]))
        end;
        ys_naive = if rem(length(x_naive), N_per_slice) == 0
            ys_naive
        else
            idx = length(x_naive) - rem(length(x_naive), N_per_slice) + 1:length(x_naive)
            vcat(
                ys_naive,
                rand(rng, f(x_naive[idx], S_naive[idx]))
            )
        end;
        y_naive = vcat(ys_naive...);
        y = TemporalGPs.match_to(y_naive, x);


        # Perform approximate inference using the state-space approach.
        post_marginals_sde = TemporalGPs.approx_posterior_marginals(
            dtc, to_sde(f)(x, S_naive), y_naive, z_r, r_pred,
        )

        # Perform approximate inference naively.
        fx = f(x_naive, S_naive);
        f_post_approx = f | Stheno.PseudoObs(fx ← y_naive, f(z_naive, 1e-9));
        x_pred_naive = collect(TemporalGPs.RectilinearGrid(r_pred, t));
        post_marginals_naive = marginals(f_post_approx(x_pred_naive));

        approx_post_mean_sde = reshape(mean.(post_marginals_sde), length(r_pred), :);
        approx_post_mean_naive = reshape(mean.(post_marginals_naive), length(r_pred), :);
        @assert maximum(abs.(approx_post_mean_sde .- approx_post_mean_naive)) < 1e-7

        approx_post_std_sde = reshape(std.(post_marginals_sde), length(r_pred), :);
        approx_post_std_naive = reshape(std.(post_marginals_naive), length(r_pred), :);
        @assert maximum(abs.(approx_post_std_sde .- approx_post_std_naive)) < 1e-7

        # Generate plot. This should be roughly half a page width.
        width_inches = 2.8
        width_px = width_inches * Plots.DPI
        height_px = width_px

        plt_post = plot(
            xlabel="time",
            ylabel="space",
            legend=:bottomright,
            legendfontsize=9,
            background_color_legend=RGBA(1, 1, 1, 1),
            size=(width_px, height_px),
            grid = (0, 0),
            colorbar=nothing,
            xlims=(minimum(t), maximum(t)),
            ylim=(minimum(r_pred), maximum(r_pred)),
        );

        scatter!(plt_post, fill(0.0, length(z_r)), z_r;
            color=:black,
            label="Pseudo-Inputs",
        );
        scatter!(plt_post, last.(x_naive), first.(x_naive);
            markersize=0.1,
            color=:black,
            label="Observation Inputs",
        );
        savefig(plt_post, joinpath(first(save_dirs), "post_irregular.tikz"));

        # This image will have to be manually added onto the plot after it's been made.
        mean_min = minimum(approx_post_mean_sde);
        mean_max = maximum(approx_post_mean_sde);
        img = colorview(Gray, (approx_post_mean_sde .- mean_min) ./ (mean_max - mean_min));
        FileIO.save(joinpath(first(save_dirs), "post_irregular.jpg"), img);
    end
end
