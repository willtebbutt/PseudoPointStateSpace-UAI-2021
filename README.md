# Combining Pseudo-Point and State Space Approximations for Sum-Separable Gaussian Processes

Implements the experiments in [1].
This repo is helpful for the sake of reproducing the experiments, but a better starting
point for using the approximation developed in [1] is [TemporalGPs.jl](https://github.com/JuliaGaussianProcesses/TemporalGPs.jl).

All experiments should be run from the top level of this repo -- don't navigate into the scripts directory.

You must have a working version of Julia installed to run this code -- experiments were run using version 1.5.3, a copy of which can be obtained from https://julialang.org/downloads/ .
I recommend running these experiments using version 1.5.3.

There are no additional binary dependencies other than Julia itself.

This project comes equipped with a `Manifest.toml`, which is a specification of the exact
commit of each dependency which was used in this project.
Consequently, when running any of the scripts below, you should wind up with precisely the
same versions of all dependencies originally used to produce the results from the paper.

## Benchmarking
Run the `scripts/run_benchmarking_experiments.sh` script.
Results will be output to `results/benchmarking`.

## Global Historical Climatology Network

Firstly, you'll need to get hold of the NASADEM (NASA Digital Elevation Map) data set.
First sign up for an account here: https://urs.earthdata.nasa.gov/users/new to get credentials.
Then navigate to the `scripts/nasa_data_grabbing` directory and run the
`get_all_of_nasadem.sh` script.
Running this script will prompt you for the credentials you obtained from NASA.
It will download the data associated with the region studied in this work.
It will produce a number of zip files, which you should unzip into `scripts/nasa_data_grabbing/nasa_data`.

Once the above has been achieved, run
```
julia scripts/global_historical_climatology_network.jl
```
to run all experiments and reproduce results.
Running this script will prompt you to download a large amount of data from NOAA's
Global Historical Climatology Network daily dataset, which is necessary to run these
experiments.

Results will be output into `results/global_historical_climatology_network`.

## Apartment Prices

Run
```
julia scripts/apartment_prices.jl
```
Results will be output into `results/apartment_prices`.
Running this script will prompt you to download a dataset from HM Land Registry database,
and a database from Camden council containing a mapping between postcodes and lat-lon
coordinates.


## Something Doesn't Work?

Please open an issue.
