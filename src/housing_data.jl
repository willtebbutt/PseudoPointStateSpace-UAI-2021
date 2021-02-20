using CSV
using DataDeps
using DataFrames
using Dates
using JLD2

# Obtained link from this webpage:
# https://www.gov.uk/government/statistical-data-sets/price-paid-data-downloads
register(DataDep(
    "uk_house_price_data",
    "This is the UK house price transaction data. Roughly 3.6GB CSV file.",
    "http://prod.publicdata.landregistry.gov.uk.s3-website-eu-west-1.amazonaws.com" *
        "/pp-complete.csv",
))

register(DataDep(
    "uk_postode_lat_lon",
    "UK postcode vs lat-lon mappings.",
    "https://opendata.camden.gov.uk/api/views/77ra-mbbn/rows.csv?accessType=DOWNLOAD",
))

function generate_housing_data(data_dir)

    data_fname = joinpath(data_dir, "housing.jld2")
    if !isfile(data_fname)

        @warn "SERIALISING DATA"

        house_data_file_name = datadep"uk_house_price_data/pp-complete.csv"

        postcode_file_name =
            datadep"uk_postode_lat_lon/National_Statistics_Postcode_Lookup_UK_Coordinates.csv"

        # Load price data.
        raw_data = CSV.read(house_data_file_name, DataFrame; header=false);

        named_data = dropmissing(DataFrame(
            price=raw_data.Column2,
            date=map(x -> parse(Date, x[1:10]), raw_data.Column3),
            post_code=map(x -> x === missing ? x : filter(!isspace, x), raw_data.Column4),
            property_type=raw_data.Column5,
            new_build=raw_data.Column6,
            estate_type=raw_data.Column7,
        ));

        # Load lat lon data.
        raw_lat_lon = CSV.read(postcode_file_name, DataFrame);

        lat_lon = DataFrame(
            post_code=filter.(!isspace, raw_lat_lon[!, "Postcode 1"]),
            lat=raw_lat_lon.Latitude,
            lon=raw_lat_lon.Longitude,
        );

        # Merge the two data sets.
        all_data = innerjoin(named_data, lat_lon; on=:post_code);

        # Serialise the data.
        FileIO.save(data_fname, Dict("data" => all_data))
    end

    # Load the from the serialised format.
    return FileIO.load(data_fname)["data"]
end
