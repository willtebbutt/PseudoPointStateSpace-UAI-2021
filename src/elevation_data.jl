using ArchGDAL
using CSV
using DataFrames
using Geodesy
using Images
using ImageTransformations
using Plots

# File-naming-parsing functionality.
get_north(x) = x[1] == 'n' ? parse(Int, x[2:3]) : -parse(Int, x[2:3])
get_east(x) = x[4] == 'e' ? parse(Int, x[5:7]) : -parse(Int, x[5:7])

# Is x inside the bounding box defined by lats / lons?
is_inside(x, lats, lons) = lats[1] <= x.lat <= lats[2] && lons[1] <= x.lon <= lons[2]

# Set these to define the region that the grid covers. Must be set large enough to cover
# all of the elevation that you might be interested in.

# Old GB + Ireland + a-bit-of-France bbox.
# const lons = (-11, 2);
# const lats = (49, 60);

# # Central Europe bbox.
# const lons = (5, 15);
# const lats = (45, 55);

const lons = (-125, -122);
const lats = (47, 49);

# const lon_interval = (-125, -122);
# const lat_interval = (47, 49);



# The NASA data comprises square images of size img_side_length x img_side_length.
# ratio is used to choose lowering the resolution of these images so that the combined grid
# actually fits in memory. A Europe-wide full-resolution image uses about 17GB of memory.
img_side_length = 3601;
ratio = 1 / 8;

# Apply the imresize operation to a fake full-resolution imaget to determine the size of the
# lower-resolution images.
new_size = size(imresize(randn(img_side_length, img_side_length); ratio=ratio));
img_height = new_size[1]
img_width = new_size[2]

# Pre-allocate for the joined-up image data.
elevation_data_mat = fill(
    0.0,
    img_width * (lons[2] - lons[1] + 1),
    img_height * (lats[2] - lats[1] + 1),
);

# Load each of the NASA data files, and plonk them in the correct location in the matrix.
fnames = filter(x -> x[9:11] == "hgt", readdir("scripts/nasa_data_grabbing/nasa_data"));

for fname in fnames
    lat = get_north(fname)
    lon = get_east(fname)
    if is_inside(LLA(lat, lon, 0), lats, lons)

        start_row = img_height * (lat - lats[1]) + 1
        row_range = start_row:(start_row + img_height - 1)

        start_col = img_width * (lon - lons[1]) + 1
        col_range = start_col:(start_col + img_width - 1)

        data = ArchGDAL.read(joinpath("scripts", "nasa_data_grabbing", "nasa_data", fname));

        band = ArchGDAL.getband(data, 1);
        small_band = imresize(band[:, :]; ratio=ratio)

        elevation_data_mat[col_range, row_range] .= reverse(small_band; dims=2)
    end
end;

pixel_lons = range(
    lons[1];
    step=(lons[2] - lons[1] + 1) / size(elevation_data_mat, 1),
    length=size(elevation_data_mat, 1),
);
pixel_lats = range(
    lats[1];
    step=(lats[2] - lats[1] + 1) / size(elevation_data_mat, 2),
    length=size(elevation_data_mat, 2),
);


# Functionality to determine the value of the closest pixel to the requested lat-lon.
function nearest_point(elevation_data, lat::Real, lon::Real)
    return elevation_data[
        nearest_point_idx(pixel_lons, lon),
        nearest_point_idx(pixel_lats, lat),
    ]
end

function nearest_point_idx(xs::StepRangeLen, x::Real)
    return convert(Int, ceil((x - first(xs) + 1e-9) / (last(xs) - first(xs)) * length(xs)))
end



# # Plot the matrix as an image.
# min_val = minimum(elevation_data_mat);
# max_val = maximum(elevation_data_mat);
# img = Gray.((elevation_data_mat .- min_val) ./ (max_val - min_val));
# Plots.plot(reverse(img; dims=2)')

# log_img = log.(Float64.(img) .+ 1);
# min_log_val = minimum(log_img);
# max_log_val = maximum(log_img);
# new_img = Gray.((log_img .- min_log_val) ./ (max_log_val - min_log_val));
# Plots.plot(reverse(new_img; dims=2)')



# heatmap(log.(elevation_data_mat .- minimum(elevation_data_mat) .+ 1))
