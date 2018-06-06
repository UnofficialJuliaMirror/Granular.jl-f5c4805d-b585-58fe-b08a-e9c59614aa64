import Compat
using Compat.Test
using Compat.LinearAlgebra

hasNetCDF = false
if VERSION < v"0.7.0-alpha"
    if typeof(Compat.Pkg.installed("NetCDF")) == VersionNumber
        import NetCDF
        hasNetCDF = true
    end
else
    import Pkg
    if haskey(Pkg.installed(), "NetCDF")
        import NetCDF
        hasNetCDF = true
    end
end
if !hasNetCDF
    Compat.@warn "Package NetCDF not found. " *
         "Ocean/atmosphere grid read not supported. " * 
         "If required, install NetCDF and its " *
         "requirements with `Pkg.add(\"NetCDF\")`."
end

export createEmptyOcean
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
                 zeros(1,1,1,1),
                 Array{Array{Int, 1}}(undef, 1, 1),
                 zeros(1,1),
                 1, 1, 1, 1,
                 false, [0.,0.,0.], [1.,1.,1.], [1,1,1], [1.,1.,1.])
end

export readOceanNetCDF
"""
Read ocean NetCDF files generated by MOM6 from disk and return as `Ocean` data 
structure.

# Arguments
* `velocity_file::String`: path to NetCDF file containing ocean velocities, 
    etc., (e.g. `prog__####_###.nc`).
* `grid_file::String`: path to NetCDF file containing ocean super-grid 
    information (typically `INPUT/ocean_hgrid.nc`).
* `regular_grid::Bool=false`: `true` if the grid is regular (all cells
    equal and grid is Cartesian) or `false` (default).
"""
function readOceanNetCDF(velocity_file::String, grid_file::String;
                         regular_grid::Bool=false)

    if !hasNetCDF
        Compat.@warn "Package NetCDF not found. " *
            "Ocean/atmosphere grid read not supported. " * 
             "Please install NetCDF and its " *
             "requirements with `Pkg.add(\"NetCDF\")`."
    else

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
                      e,
                      Array{Array{Int, 1}}(undef, size(xh, 1), size(xh, 2)),
                      zeros(size(xh)),
                      1, 1, 1, 1,

                      false, [0.,0.,0.], [1.,1.,1.], [1,1,1], [1.,1.,1.]
                     )
        return ocean
    end
end

export readOceanStateNetCDF
"""
Read NetCDF file with ocean state generated by MOM6 (e.g.  `prog__####_###.nc` 
or `########.ocean_month.nc`) from disk and return time stamps, velocity fields, 
layer thicknesses, interface heights, and vertical coordinates.

# Returns
* `time::Vector{Float64}`: Time [s]
* `u::Array{Float64, 2}`: Cell corner zonal velocity [m/s],
    dimensions correspond to placement in `[xq, yq, zl, time]`
* `v::Array{Float64, 2}`: Cell corner meridional velocity [m/s],
    dimensions correspond to placement in `[xq, yq, zl, time]`
* `h::Array{Float64, 2}`: layer thickness [m], dimensions correspond to 
    placement in `[xh, yh, zl, time]`
* `e::Array{Float64, 2}`: interface height relative to mean sea level [m],  
    dimensions correspond to placement in `[xh, yh, zi, time]`
* `zl::Vector{Float64}`: layer target potential density [kg m^-3]
* `zi::Vector{Float64}`: interface target potential density [kg m^-3]
"""
function readOceanStateNetCDF(filename::String)

    if !hasNetCDF
        Compat.@warn "Package NetCDF not found. " *
            "Ocean/atmosphere grid read not supported. " * 
             "Please install NetCDF and its " *
             "requirements with `Pkg.add(\"NetCDF\")`."
    else

        if !isfile(filename)
            error("$(filename) could not be opened")
        end

        u_staggered = convert(Array{Float64, 4}, NetCDF.ncread(filename, "u"))
        v_staggered = convert(Array{Float64, 4}, NetCDF.ncread(filename, "v"))
        u, v = interpolateOceanVelocitiesToCorners(u_staggered, v_staggered)

        time = convert(Vector{Float64},
                       NetCDF.ncread(filename, "time") .* 24. * 60. * 60.)
        h = convert(Array{Float64, 4}, NetCDF.ncread(filename, "h"))
        e = convert(Array{Float64, 4}, NetCDF.ncread(filename, "e"))

        zl = convert(Vector{Float64}, NetCDF.ncread(filename, "zl"))
        zi = convert(Vector{Float64}, NetCDF.ncread(filename, "zi"))

        return time, u, v, h, e, zl, zi
    end
end

export readOceanGridNetCDF
"""
Read NetCDF file with ocean *supergrid* information generated by MOM6 (e.g.  
`ocean_hrid.nc`) from disk and return as `Ocean` data structure.  This file is 
located in the simulation `INPUT/` subdirectory.

# Returns
* `xh::Array{Float64, 2}`: Longitude for cell centers [deg]
* `yh::Array{Float64, 2}`: Latitude for cell centers [deg]
* `xq::Array{Float64, 2}`: Longitude for cell corners [deg]
* `yq::Array{Float64, 2}`: Latitude for cell corners [deg]
"""
function readOceanGridNetCDF(filename::String)

    if !hasNetCDF
        Compat.@warn "Package NetCDF not found. " *
            "Ocean/atmosphere grid read not supported. " * 
             "Please install NetCDF and its " *
             "requirements with `Pkg.add(\"NetCDF\")`."
    else

        if !isfile(filename)
            error("$(filename) could not be opened")
        end
        x = convert(Array{Float64, 2}, NetCDF.ncread(filename, "x"))
        y = convert(Array{Float64, 2}, NetCDF.ncread(filename, "y"))

        xh = x[2:2:end, 2:2:end]
        yh = y[2:2:end, 2:2:end]

        xq = x[1:2:end, 1:2:end]
        yq = y[1:2:end, 1:2:end]

        return xh, yh, xq, yq
    end
end

export interpolateOceanVelocitiesToCorners
"""
Convert gridded data from Arakawa-C type (decomposed velocities at faces) to 
Arakawa-B type (velocities at corners) through interpolation.
"""
function interpolateOceanVelocitiesToCorners(u_in::Array{Float64, 4}, 
                                             v_in::Array{Float64, 4})

    if size(u_in) != size(v_in)
        error("size of u_in ($(size(u_in))) must match v_in ($(size(v_in)))")
    end

    nx, ny, nz, nt = size(u_in)
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

export interpolateOceanState
"""
Ocean data is containted in `Ocean` type at discrete times (`Ocean.time`).  This 
function performs linear interpolation between time steps to get the approximate 
ocean state at any point in time.  If the `Ocean` data set only contains a 
single time step, values from that time are returned.
"""
function interpolateOceanState(ocean::Ocean, t::Float64)
    if length(ocean.time) == 1
        return ocean.u, ocean.v, ocean.h, ocean.e
    elseif t < ocean.time[1] || t > ocean.time[end]
        error("selected time (t = $(t)) is outside the range of time steps " *
              "in the ocean data")
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

    return ocean.u[:,:,:,i]*(1. - rel_time) + ocean.u[:,:,:,i+1]*rel_time,
        ocean.v[:,:,:,i]*(1. - rel_time) + ocean.v[:,:,:,i+1]*rel_time,
        ocean.h[:,:,:,i]*(1. - rel_time) + ocean.h[:,:,:,i+1]*rel_time,
        ocean.e[:,:,:,i]*(1. - rel_time) + ocean.e[:,:,:,i+1]*rel_time
end

export createRegularOceanGrid
"""

    createRegularOceanGrid(n, L[, origo, time, name,
                           bc_west, bc_south, bc_east, bc_north])

Initialize and return a regular, Cartesian `Ocean` grid with `n[1]` by `n[2]` 
cells in the horizontal dimension, and `n[3]` vertical cells.  The cell corner 
and center coordinates will be set according to the grid spatial dimensions 
`L[1]`, `L[2]`, and `L[3]`.  The grid `u`, `v`, `h`, and `e` fields will contain 
one 4-th dimension matrix per `time` step.  Sea surface will be at `z=0.` with 
the ocean spanning `z<0.`.  Vertical indexing starts with `k=0` at the sea 
surface, and increases downwards.

# Arguments
* `n::Vector{Int}`: number of cells along each dimension [-].
* `L::Vector{Float64}`: domain length along each dimension [m].
* `origo::Vector{Float64}`: domain offset in each dimension [m] (default =
    `[0.0, 0.0]`).
* `time::Vector{Float64}`: vector of time stamps for the grid [s].
* `name::String`: grid name (default = `"unnamed"`).
* `bc_west::Integer`: grid boundary condition for the grains.
* `bc_south::Integer`: grid boundary condition for the grains.
* `bc_east::Integer`: grid boundary condition for the grains.
* `bc_north::Integer`: grid boundary condition for the grains.
"""
function createRegularOceanGrid(n::Vector{Int},
                                L::Vector{Float64};
                                origo::Vector{Float64} = zeros(2),
                                time::Vector{Float64} = zeros(1),
                                name::String = "unnamed",
                                bc_west::Integer = 1,
                                bc_south::Integer = 1,
                                bc_east::Integer = 1,
                                bc_north::Integer = 1)

    xq = repeat(Compat.range(origo[1], stop=origo[1] + L[1],
                             length=n[1] + 1),
                outer=[1, n[2] + 1])
    yq = repeat(Compat.range(origo[2], stop=origo[2] + L[2],
                             length=n[2] + 1)',
                outer=[n[1] + 1, 1])

    dx = L./n
    xh = repeat(Compat.range(origo[1] + .5*dx[1],
                             stop=origo[1] + L[1] - .5*dx[1],
                             length=n[1]),
                outer=[1, n[2]])
    yh = repeat(Compat.range(origo[2] + .5*dx[2],
                             stop=origo[2] + L[2] - .5*dx[2],
                             length=n[2])',
                outer=[n[1], 1])

    zl = -Compat.range(.5*dx[3], stop=L[3] - .5*dx[3], length=n[3])
    zi = -Compat.range(0., stop=L[3], length=n[3] + 1)

    u = zeros(n[1] + 1, n[2] + 1, n[3], length(time))
    v = zeros(n[1] + 1, n[2] + 1, n[3], length(time))
    h = zeros(n[1] + 1, n[2] + 1, n[3], length(time))
    e = zeros(n[1] + 1, n[2] + 1, n[3], length(time))

    return Ocean(name,
                 time,
                 xq, yq,
                 xh, yh,
                 zl, zi,
                 u, v, h, e,
                 Array{Array{Int, 1}}(undef, size(xh, 1), size(xh, 2)),
                 zeros(size(xh)),
                 bc_west, bc_south, bc_east, bc_north,
                 true, origo, L, n, dx)
end

export addOceanDrag!
"""
Add drag from linear and angular velocity difference between ocean and all ice 
floes.
"""
function addOceanDrag!(simulation::Simulation)
    if typeof(simulation.ocean.input_file) == Bool
        error("no ocean data read")
    end

    u, v, h, e = interpolateOceanState(simulation.ocean, simulation.time)
    uv_interp = Vector{Float64}(undef, 2)
    sw = Vector{Float64}(undef, 2)
    se = Vector{Float64}(undef, 2)
    ne = Vector{Float64}(undef, 2)
    nw = Vector{Float64}(undef, 2)

    for grain in simulation.grains

        if !grain.enabled
            continue
        end

        i, j = grain.ocean_grid_pos
        k = 1

        x_tilde, y_tilde = getNonDimensionalCellCoordinates(simulation.ocean,
                                                            i, j,
                                                            grain.lin_pos)
        x_tilde = clamp(x_tilde, 0., 1.)
        y_tilde = clamp(y_tilde, 0., 1.)

        bilinearInterpolation!(uv_interp, u, v, x_tilde, y_tilde, i, j, k, 1)
        applyOceanDragToGrain!(grain, uv_interp[1], uv_interp[2])
        applyOceanVorticityToGrain!(grain,
                                      curl(simulation.ocean, x_tilde, y_tilde,
                                           i, j, k, 1, sw, se, ne, nw))
    end
    nothing
end

export applyOceanDragToGrain!
"""
Add Stokes-type drag from velocity difference between ocean and a single ice 
floe.
"""
function applyOceanDragToGrain!(grain::GrainCylindrical,
                                  u::Float64, v::Float64)
    freeboard = .1*grain.thickness  # height above water
    ρ_o = 1000.   # ocean density
    draft = grain.thickness - freeboard  # height of submerged thickness

    drag_force = ρ_o * π *
    (2.0*grain.ocean_drag_coeff_vert*grain.areal_radius*draft + 
        grain.ocean_drag_coeff_horiz*grain.areal_radius^2.0) *
        ([u, v] - grain.lin_vel)*norm([u, v] - grain.lin_vel)

    grain.force += drag_force
    grain.ocean_stress = drag_force/grain.horizontal_surface_area
    nothing
end

export applyOceanVorticityToGrain!
"""
Add Stokes-type torque from angular velocity difference between ocean and a 
single grain.  See Eq. 9.28 in "Introduction to Fluid Mechanics" by Nakayama 
and Boucher, 1999.
"""
function applyOceanVorticityToGrain!(grain::GrainCylindrical, 
                                       ocean_curl::Float64)
    freeboard = .1*grain.thickness  # height above water
    ρ_o = 1000.   # ocean density
    draft = grain.thickness - freeboard  # height of submerged thickness

    grain.torque +=
        π * grain.areal_radius^4. * ρ_o * 
        (grain.areal_radius/5. * grain.ocean_drag_coeff_horiz + 
        draft * grain.ocean_drag_coeff_vert) * 
        abs(.5 * ocean_curl - grain.ang_vel) * (.5 * ocean_curl - grain.ang_vel)
    nothing
end

export compareOceans
"""
    compareOceans(ocean1::Ocean, ocean2::Ocean)

Compare values of two `Ocean` objects using the `Base.Test` framework.
"""
function compareOceans(ocean1::Ocean, ocean2::Ocean)

    @test ocean1.input_file == ocean2.input_file
    @test ocean1.time ≈ ocean2.time

    @test ocean1.xq ≈ ocean2.xq
    @test ocean1.yq ≈ ocean2.yq

    @test ocean1.xh ≈ ocean2.xh
    @test ocean1.yh ≈ ocean2.yh

    @test ocean1.zl ≈ ocean2.zl
    @test ocean1.zi ≈ ocean2.zi

    @test ocean1.u ≈ ocean2.u
    @test ocean1.v ≈ ocean2.v
    @test ocean1.h ≈ ocean2.h
    @test ocean1.e ≈ ocean2.e

    if isassigned(ocean1.grain_list, 1)
        @test ocean1.grain_list == ocean2.grain_list
    end
    nothing
end
