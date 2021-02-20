julia scripts/pseudo_point_synthetic_experiments_irregular_obs.jl --run-experiments=true
julia scripts/pseudo_point_synthetic_experiments_regular_obs.jl --run-experiments=true
julia scripts/pseudo_point_synthetic_posterior_plots.jl
julia scripts/pseudo_point_synthetic_experiments_irregular_obs_sum_separable.jl --run-experiments=true
julia scripts/pseudo_point_synthetic_experiments_regular_obs_sum_separable.jl --run-experiments=true
julia scripts/global_historical_climatology_network.jl --run-experiments=true
julia scripts/apartment_prices.jl --run-experiments=true
