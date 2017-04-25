"Returns empty ocean type for initialization purposes."
function createEmptyOcean()
    return Ocean(false,

                 zeros(1),

                 zeros(1,1),
                 zeros(1,1),
                 zeros(1,1),
                 zeros(1,1),

                 zeros(1),
                 zeros(1),

                 zeros(1,1,1,1),
                 zeros(1,1,1,1),
                 zeros(1,1,1,1),
                 zeros(1,1,1,1))
end

"""
Read ocean NetCDF files generated by MOM6 from disk and return as `Ocean` data 
structure.

# Arguments
* `velocity_file::String`: Path to NetCDF file containing ocean velocities, 
    etc., (e.g. `prog__####_###.nc`).
* `grid_file::String`: Path to NetCDF file containing ocean super-grid 
    information (typically `INPUT/ocean_hgrid.nc`).
"""
function readOceanNetCDF(velocity_file::String, grid_file::String)

    time, u, v, h, e, zl, zi = readOceanStateNetCDF(velocity_file)
    xh, yh, xq, yq = readOceanGridNetCDF(grid_file)

    if size(u[:,:,1,1]) != size(xq) || size(v[:,:,1,1]) != size(xq) ||
        size(xq) != size(yq)
        error("size mismatch between velocities and grid
              (u: $(size(u[:,:,1,1])), v: $(size(v[:,:,1,1])),
              xq: $(size(xq)), yq: $(size(yq)))")
    end

    ocean = Ocean([grid_file, velocity_file],

                  time,

                  xq,
                  yq,
                  xh,
                  yh,

                  zl,
                  zi,

                  u,
                  v,
                  h,
                  e)
    return ocean
end

"""
Read NetCDF file with ocean state generated by MOM6 (e.g.  `prog__####_###.nc` 
or `########.ocean_month.nc`) from disk and return time stamps, velocity fields, 
layer thicknesses, interface heights, and vertical coordinates.

# Returns
* `time::Array{Float, 1}`: Time [s]
* `u::Array{Float, 2}`: Cell corner zonal velocity [m/s],
    dimensions correspond to placement in `[xq, yq, zl, time]`
* `v::Array{Float, 2}`: Cell corner meridional velocity [m/s],
    dimensions correspond to placement in `[xq, yq, zl, time]`
* `h::Array{Float64, 2}`: layer thickness [m], dimensions correspond to 
    placement in `[xh, yh, zl, time]`
* `e::Array{Float64, 2}`: interface height relative to mean sea level [m],  
    dimensions correspond to placement in `[xh, yh, zi, time]`
* `zl::Array{Float64, 1}`: layer target potential density [kg m^-3]
* `zi::Array{Float64, 1}`: interface target potential density [kg m^-3]
"""
function readOceanStateNetCDF(filename::String)

    if !isfile(filename)
        error("$(filename) could not be opened")
    end

    u_staggered::Array{float, 4} = NetCDF.ncread(filename, "u")
    v_staggered::Array{float, 4} = NetCDF.ncread(filename, "v")
    u, v = interpolateOceanVelocitiesToCorners(u_staggered, v_staggered)

    time = NetCDF.ncread(filename, "time")*24.*60.*60.
    h = NetCDF.ncread(filename, "h")
    e = NetCDF.ncread(filename, "e")

    zl = NetCDF.ncread(filename, "zl")
    zi = NetCDF.ncread(filename, "zi")

    return time, u, v, h, e, zl, zi
end

"""
Read NetCDF file with ocean *supergrid* information generated by MOM6 (e.g.  
`ocean_hrid.nc`) from disk and return as `Ocean` data structure.  This file is 
located in the simulation `INPUT/` subdirectory.

# Returns
* `xh::Array{Float, 2}`: Longitude for cell centers [deg]
* `yh::Array{Float, 2}`: Latitude for cell centers [deg]
* `xq::Array{Float, 2}`: Longitude for cell corners [deg]
* `yq::Array{Float, 2}`: Latitude for cell corners [deg]
"""
function readOceanGridNetCDF(filename::String)

    if !isfile(filename)
        error("$(filename) could not be opened")
    end
    x::Array{float, 2} = NetCDF.ncread(filename, "x")
    y::Array{float, 2} = NetCDF.ncread(filename, "y")

    xh = x[2:2:end, 2:2:end]
    yh = y[2:2:end, 2:2:end]

    xq = x[1:2:end, 1:2:end]
    yq = y[1:2:end, 1:2:end]

    return xh, yh, xq, yq
end

"""
Convert gridded data from Arakawa-C type (decomposed velocities at faces) to 
Arakawa-B type (velocities at corners) through interpolation.
"""
function interpolateOceanVelocitiesToCorners(u_in::Array{float, 4},
                                             v_in::Array{float, 4})

    if size(u_in) != size(v_in)
        error("size of u_in ($(size(u_in))) must match v_in ($(size(v_in)))")
    end

    nx, ny, nz, nt = size(u_in)
    #u = Array{float}(nx+1, ny+1, nz, nt)
    #v = Array{float}(nx+1, ny+1, nz, nt)
    u = zeros(nx+1, ny+1, nz, nt)
    v = zeros(nx+1, ny+1, nz, nt)
    for i=1:nx
        for j=1:ny
            if j < ny - 1
                u[i, j, :, :] = (u_in[i, j, :, :] + u_in[i, j+1, :, :])/2.
            else
                u[i, j, :, :] = u_in[i, j, :, :]
            end
            if i < nx - 1
                v[i, j, :, :] = (v_in[i, j, :, :] + v_in[i+1, j, :, :])/2.
            else
                v[i, j, :, :] = v_in[i, j, :, :]
            end
        end
    end
    return u, v
end

"""
Ocean data is containted in `Ocean` type at discrete times (`Ocean.time`).  This 
function performs linear interpolation between time steps to get the approximate 
ocean state at any point in time.  If the `Ocean` data set only contains a 
single time step, values from that time are returned.
"""
function interpolateOceanState(ocean::Ocean, t::float)
    if length(ocean.time) == 1
        return ocean.u, ocean.v, ocean.h, ocean.e
    elseif t < ocean.time[1] || t > ocean.time[end]
        error("selected time (t = $(t)) is outside the range of time steps in 
              the ocean data")
    end

    i = 1
    rel_time = 0.
    while i < length(ocean.time)
        if ocean.time[i+1] < t
            i += 1
            continue
        end

        dt = ocean.time[i+1] - ocean.time[i]
        rel_time = (t - ocean.time[i])/dt
        if rel_time < 0. || rel_time > 1.
            error("time bounds error")
        end
        break
    end

    return ocean.u[:,:,:,i]*(1. - rel_time + ocean.u[:,:,:,i+1]*rel_time),
        ocean.v[:,:,:,i]*(1. - rel_time) + ocean.v[:,:,:,i+1]*rel_time,
        ocean.h[:,:,:,i]*(1. - rel_time) + ocean.h[:,:,:,i+1]*rel_time,
        ocean.e[:,:,:,i]*(1. - rel_time) + ocean.e[:,:,:,i+1]*rel_time
end

"""
Add Stokes-type drag from velocity difference between ocean and ice floe.
"""
function addOceanDrag!(simulation::Simulation)
    if !simulation.ocean.id
        error("no ocean data read")
    end

    u, v, h, e = interpolateOceanState(simulation.ocean, simulation.time)

end
