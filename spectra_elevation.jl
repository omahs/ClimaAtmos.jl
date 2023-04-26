using NCDatasets
using DelimitedFiles
using Plots
using ImageFiltering
using Interpolations
using ClimaCore
using ClimaCore.InputOutput
using ClimaCore: Geometry, Meshes, Domains, Topologies, Spaces
using ClimaCoreTempestRemap
import ClimaCoreSpectra as CCS

# Load Topography Data
data = NCDataset("/Users/asridhar/Research/Codes/Topography/ETOPO1_Ice_g_gdal.grd")

# Unpack information
x_range = data["x_range"][:] 
y_range = data["y_range"][:]
z_range = data["z_range"][:]
spacing = data["spacing"][:]
dimension = data["dimension"][:]
elevation = data["z"][:]
lon = x_range[1]:spacing[1]:x_range[2]
lat = y_range[1]:spacing[2]:y_range[2]
nlon = Int64(dimension[1])
nlat = Int64(dimension[2])
zlevels = reshape(elevation, nlon, nlat)
zlevels[zlevels .<= 0.0] .= 0.0
zlevels = reverse(zlevels, dims=2)
zlevels_coarse = zlevels[1:50:end-1, 1:50:end-1]
lat_coarse = lat[1:100:end-1]
lon_coarse = lon[1:100:end-1]

# Apply Smoothing

plotlist = [];

target_lon = -180:1:180
target_lat = -90:1:90

add_dim(x::Array) = reshape(x, (size(x)...,1))
var_grid = zlevels_coarse;
add_dim(var_grid)

### Fourier Transforms
### Assume surface level

    num_lev = 1
    num_lat = length(lat_coarse)
    num_lon = length(lon_coarse)
    num_fourier = Int(num_lon)
    # get number of positive Fourier coefficients incl. 0
    if mod(num_lon, 2) == 0 # even
        num_pfourier = div(num_lon, 2)
    else # odd
        num_pfourier = div(num_lon, 2) + 1
    end
    zon_spectrum = zeros(FT, num_pfourier, num_lat, num_lev)
    freqs = zeros(FT, num_pfourier, num_lat, num_lev)
    for k in 1:num_lev
        for j in 1:num_lat
            # compute fft frequencies for each latitude
            x = lon ./ 180 .* π
            dx = (lon[2] - lon[1]) ./ 180 .* π

            freqs_ = FFTW.fftfreq(num_fourier, 1.0 / dx) # 0,+ve freq,-ve freqs (lowest to highest)
            freqs[:, j, k] = freqs_[1:num_pfourier] .* 2.0 .* π

            # compute the fourier coefficients for all latitudes
            fourier = FFTW.fft(FT.(var_grid[:, j, k])) # e.g. vcos_grid, ucos_grid
            fourier = (fourier / num_fourier)

            # convert to energy spectra
            zon_spectrum[1, j, k] =
                zon_spectrum[1, j, k] +
                weight[k] * fourier[1] .* conj(fourier[1])

            for m in 2:num_pfourier
                zon_spectrum[m, j, k] =
                    zon_spectrum[m, j, k] +
                    2 * weight[k] * fourier[m] * conj(fourier[m]) # factor 2 for neg freq contribution
            end
        end
    end


Plots.plot(
  new_wn[2:end-1,1,1],
  new_freq[2:end-1,1,1],
  xaxis = (:log, "wavenumber"),
  yaxis = (:log, "nm spectrum"),
  marker=:o,
  label = "raw_data"
)
Plots.plot!(
            new_wn[2:end-1,1,1],
            new_wn[2:end-1,1,1].^(-1) .* 10^8,
            ls=:dash, label = "nm⁻¹"
)

### Fourier Transforms

### REMAP Functions
#fspace = axes(Y.f)
#cspace = axes(Y.c)
#
#datafile_cc = "newtest.nc"
#ᶠz_surface = Fields.coordinate_field(Y.f.w).z
#
#NCDataset(datafile_cc, "c") do nc
#  def_space_coord(nc, cspace, type = "cgll")
#  def_space_coord(nc, fspace, type = "cgll")
#  nc_time = def_time_coord(nc)
#  nc_zsfc = defVar(nc, "zsurface", FT, fspace, ("time",))
#  nc_u = defVar(nc, "u", FT, cspace, ("time",))
#  nc_v = defVar(nc, "v", FT, cspace, ("time",))
#  nc_w = defVar(nc, "w", FT, cspace, ("time",))
#  nc_time[1] = FT(0)
#  nc_u[:, 1] = @. Geometry.project(Geometry.UAxis(), p.ᶜu).components.data.:1
#  nc_v[:, 1] = @. Geometry.project(Geometry.VAxis(), p.ᶜu).components.data.:1
#  nc_w[:, 1] = @. Geometry.project(Geometry.WAxis(), p.ᶜu).components.data.:1
#  nc_zsfc[:, 1] = Fields.coordinate_field(Y.f).z
#end
#
## write out our cubed sphere mesh
#meshfile_cc = "mesh_cubedsphere.g"
#write_exodus(meshfile_cc, fspace.horizontal_space.topology)
#
## write out RLL mesh
#nlat = 90
#nlon = 180
#meshfile_rll = "mesh_rll.g"
#rll_mesh(meshfile_rll; nlat = nlat, nlon = nlon)
#
## construct overlap mesh
#meshfile_overlap = "mesh_overlap.g"
#overlap_mesh(meshfile_overlap, meshfile_cc, meshfile_rll)
#
## construct remap weight file
#weightfile = "remap_weights.nc"
#remap_weights(
#    weightfile,
#    meshfile_cc,
#    meshfile_rll,
#    meshfile_overlap;
#    in_type = "cgll",
#    in_np = Spaces.Quadratures.degrees_of_freedom(Spaces.quadrature_style(fspace)),
#)
#
## apply remap
#datafile_rll = "data_rll.nc"
#apply_remap(datafile_rll, datafile_cc, weightfile, ["zsurface", "u", "v", "w"])
