import WriteVTK
import NetCDF
hasJLD = false
if typeof(Pkg.installed("JLD")) == VersionNumber
    import JLD
    hasJLD = true
end
using Compat.DelimitedFiles

# load Homebrew/WinRPM for gnuplot and imagemagick
if is_apple() 
    using Homebrew
end

if is_windows() 
    using WinRPM
end

## IO functions

export writeSimulation
"""
    writeSimulation(simulation::Simulation;
                         filename::String="",
                         folder::String=".",
                         verbose::Bool=true)

Write all content from `Simulation` to disk in JDL format.  If the `filename` 
parameter is not specified, it will be saved to a subdirectory under the current 
directory named after the simulation identifier `simulation.id`.
"""
function writeSimulation(simulation::Simulation;
                         filename::String="",
                         folder::String=".",
                         verbose::Bool=true)
    if !hasJLD
        warn("Package JLD not found. Simulation save/read not supported. " * 
             "Please install JLD and its " *
             "requirements with `Pkg.add(\"JLD\")`.")
    else
        if filename == ""
            folder = folder * "/" * simulation.id
            mkpath(folder)
            filename = string(folder, "/", simulation.id, ".",
                              simulation.file_number, ".jld")
        end

        JLD.save(filename, "simulation", simulation)

        if verbose
            info("simulation written to $filename")
        end
    end
    nothing
end

export readSimulation
"""
    readSimulation(filename::String="";
                   verbose::Bool=true)

Return `Simulation` content read from disk using the JDL format.

# Arguments
* `filename::String`: path to file on disk containing the simulation
    information.
* `verbose::Bool=true`: confirm to console that the file has been read.
"""
function readSimulation(filename::String;
                         verbose::Bool=true)
    if !hasJLD
        warn("Package JLD not found. Simulation save/read not supported. " * 
             "Please install JLD and its " *
             "requirements with `Pkg.add(\"JLD\")`.")
        nothing
    else
        return JLD.load(filename, "simulation")
        if verbose
            info("Read simulation from $filename")
        end
    end
end
"""
    readSimulation(simulation::Simulation;
                   step::Integer = -1,
                   verbose::Bool = true)

Read the simulation state from disk and return as new simulation object.

# Arguments
* `simulation::Simulation`: use the `simulation.id` to determine the file name
    to read from, and read information from the file into this object.
* `step::Integer=-1`: attempt to read this output file step. At its default
    value (`-1`), the function will try to read the latest file, determined by
    calling [`readSimulationStatus`](@ref).
* `verbose::Bool=true`: confirm to console that the file has been read.
"""
function readSimulation(simulation::Simulation;
                         step::Integer = -1,
                         verbose::Bool = true)
    if !hasJLD
        warn("Package JLD not found. Simulation save/read not supported. " * 
             "Please install JLD and its " *
             "requirements with `Pkg.add(\"JLD\")`.")
        nothing
    else
        if step == -1
            step = readSimulationStatus(simulation)
        end
        filename = string(simulation.id, "/", simulation.id, ".$step.jld")
        if verbose
            info("Read simulation from $filename")
        end
        return JLD.load(filename, "simulation")
    end
end

export writeSimulationStatus
"""
    writeSimulationStatus(simulation::Simulation;
                          folder::String=".",
                          verbose::Bool=false)

Write current simulation status to disk in a minimal txt file.
"""
function writeSimulationStatus(simulation::Simulation;
                               folder::String=".",
                               verbose::Bool=false)
    folder = folder * "/" * simulation.id
    mkpath(folder)
    filename = string(folder, "/", simulation.id, ".status.txt")

    writedlm(filename, [simulation.time
                        simulation.time/simulation.time_total*100.
                        float(simulation.file_number)])
    if verbose
        info("Wrote status to $filename")
    end
    nothing
end

export readSimulationStatus
"""
    readSimulationStatus(simulation_id[, folder, verbose])

Read the current simulation status from disk (`<sim.id>/<sim.id>.status.txt`)
and return the last output file number.

# Arguments
* `simulation_id::String`: the simulation identifying string.
* `folder::String="."`: the folder in which to search for the status file.
* `verbose::Bool=true`: show simulation status in console.
"""
function readSimulationStatus(simulation_id::String;
                              folder::String=".",
                              verbose::Bool=true)

    folder = folder * "/" * simulation_id
    filename = string(folder, "/", simulation_id, ".status.txt")

    data = readdlm(filename)
    if verbose
        info("$simulation_id:\n" *
             "  time:             $(data[1]) s\n" *
             "  complete:         $(data[2])%\n" *
             "  last output file: $(Int(round(data[3])))\n")
    end
    return Int(round(data[3]))
"""
    readSimulationStatus(simulation[, folder, verbose])

Read the current simulation status from disk (`<sim.id>/<sim.id>.status.txt`)
and return the last output file number.

# Arguments
* `simulation::Simulation`: the simulation to read the status for.
* `folder::String="."`: the folder in which to search for the status file.
* `verbose::Bool=true`: show simulation status in console.
"""
end
function readSimulationStatus(sim::Simulation;
                              folder::String=".",
                              verbose::Bool=true)
    readSimulationStatus(sim.id, folder=folder, verbose=verbose)
end

export status
"""
Shows the status of all simulations with output files written under the 
specified `folder`, which is the current working directory by default.
"""
function status(folder::String=".";
                loop::Bool=false,
                t_int::Int=10,
                colored_output::Bool=true,
                write_header::Bool=true)

    if colored_output
        id_color_complete = :green
        id_color_in_progress = :yellow
        time_color = :white
        percentage_color = :blue
        lastfile_color = :cyan
    else
        id_color_complete = :default
        id_color_in_progress = :default
        time_color = :default
        percentage_color = :default
        lastfile_color = :default
    end

    repeat = true
    while repeat

        status_files = String[]
        println(Dates.format(DateTime(now()), Dates.RFC1123Format))

        for (root, dirs, files) in walkdir(folder, follow_symlinks=false)

            for file in files
                if contains(file, ".status.txt")
                    push!(status_files, joinpath(root, file))
                end
            end
        end

        if length(status_files) > 0
            if write_header
                println("--------------------------------------" * 
                        "--------------------------------------")
                print_with_color(:default, "simulation folder \t")
                print_with_color(time_color, "      time \t")
                print_with_color(percentage_color, "      completed  ")
                print_with_color(lastfile_color, "last file \n")
                println("--------------------------------------" * 
                        "--------------------------------------")
            end

            for file in status_files
                data = readdlm(file)
                id = replace(file, ".status.txt", "")
                id = replace(id, "./", "")
                id = replace(id, r".*/", "")
                time_s = @sprintf "%6.2fs" data[1]
                time_h = @sprintf "%5.1fh" data[1]/(60. * 60.)
                percentage = @sprintf "%3.0f%%" data[2]
                lastfile = @sprintf "%5d" data[3]
                if data[2] < 99.
                    print_with_color(id_color_in_progress, "$id \t")
                else
                    print_with_color(id_color_complete, "$id \t")
                end
                print_with_color(time_color, "$time_s ($time_h) \t")
                print_with_color(percentage_color, "$percentage \t")
                print_with_color(lastfile_color, "$lastfile \n")
            end
            if write_header
                println("--------------------------------------" * 
                        "--------------------------------------")
            end
        else
            warn("no simulations found in $(pwd())/$folder")
        end

        if loop && t_int > 0
            sleep(t_int)
        end
        if !loop
            repeat = false
        end
    end
    nothing
end

export writeVTK
"""
Write a VTK file to disk containing all grains in the `simulation` in an 
unstructured mesh (file type `.vtu`).  These files can be read by ParaView and 
can be visualized by applying a *Glyph* filter.

If the simulation contains an `Ocean` data structure, it's contents will be 
written to separate `.vtu` files.  This can be disabled by setting the argument 
`ocean=false`.  The same is true for the atmosphere.

The VTK files will be saved in a subfolder named after the simulation.
"""
function writeVTK(simulation::Simulation;
                  folder::String=".",
                  verbose::Bool=true,
                  ocean::Bool=true,
                  atmosphere::Bool=true)

    simulation.file_number += 1
    folder = folder * "/" * simulation.id
    mkpath(folder)

    filename = string(folder, "/", simulation.id, ".grains.", 
                      simulation.file_number)
    writeGrainVTK(simulation, filename, verbose=verbose)

    filename = string(folder, "/", simulation.id, ".grain-interaction.", 
                      simulation.file_number)
    writeGrainInteractionVTK(simulation, filename, verbose=verbose)

    if typeof(simulation.ocean.input_file) != Bool && ocean
        filename = string(folder, "/", simulation.id, ".ocean.", 
                        simulation.file_number)
        writeGridVTK(simulation.ocean, filename, verbose=verbose)
    end

    if typeof(simulation.atmosphere.input_file) != Bool && atmosphere
        filename = string(folder, "/", simulation.id, ".atmosphere.", 
                        simulation.file_number)
        writeGridVTK(simulation.atmosphere, filename, verbose=verbose)
    end
    nothing
end

export writeGrainVTK
"""
Write a VTK file to disk containing all grains in the `simulation` in an 
unstructured mesh (file type `.vtu`).  These files can be read by ParaView and 
can be visualized by applying a *Glyph* filter.  This function is called by 
`writeVTK()`.
"""
function writeGrainVTK(simulation::Simulation,
                         filename::String;
                         verbose::Bool=false)

    ifarr = convertGrainDataToArrays(simulation)
    
    # add arrays to VTK file
    vtkfile = WriteVTK.vtk_grid(filename, ifarr.lin_pos, WriteVTK.MeshCell[])

    WriteVTK.vtk_point_data(vtkfile, ifarr.density, "Density [kg m^-3]")

    WriteVTK.vtk_point_data(vtkfile, ifarr.thickness, "Thickness [m]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.contact_radius*2.,
                            "Diameter (contact) [m]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.areal_radius*2.,
                            "Diameter (areal) [m]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.circumreference,
                            "Circumreference  [m]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.horizontal_surface_area,
                            "Horizontal surface area [m^2]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.side_surface_area,
                            "Side surface area [m^2]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.volume, "Volume [m^3]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.mass, "Mass [kg]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.moment_of_inertia,
                            "Moment of inertia [kg m^2]")

    WriteVTK.vtk_point_data(vtkfile, ifarr.lin_vel, "Linear velocity [m s^-1]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.lin_acc,
                            "Linear acceleration [m s^-2]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.force, "Sum of forces [N]")

    WriteVTK.vtk_point_data(vtkfile, ifarr.ang_pos, "Angular position [rad]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.ang_vel,
                            "Angular velocity [rad s^-1]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.ang_acc,
                            "Angular acceleration [rad s^-2]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.torque, "Sum of torques [N*m]")

    WriteVTK.vtk_point_data(vtkfile, ifarr.fixed, "Fixed in space [-]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.rotating, "Free to rotate [-]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.enabled, "Enabled [-]")

    WriteVTK.vtk_point_data(vtkfile, ifarr.contact_stiffness_normal,
                            "Contact stiffness (normal) [N m^-1]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.contact_stiffness_tangential,
                            "Contact stiffness (tangential) [N m^-1]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.contact_viscosity_normal,
                            "Contact viscosity (normal) [N m^-1 s]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.contact_viscosity_tangential,
                            "Contact viscosity (tangential) [N m^-1 s]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.contact_static_friction,
                            "Contact friction (static) [-]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.contact_dynamic_friction,
                            "Contact friction (dynamic) [-]")

    WriteVTK.vtk_point_data(vtkfile, ifarr.youngs_modulus,
                            "Young's modulus [Pa]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.poissons_ratio,
                            "Poisson's ratio [-]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.tensile_strength,
                            "Tensile strength [Pa]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.compressive_strength_prefactor,
                            "Compressive strength prefactor [m^0.5 Pa]")

    WriteVTK.vtk_point_data(vtkfile, ifarr.ocean_drag_coeff_vert,
                            "Ocean drag coefficient (vertical) [-]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.ocean_drag_coeff_horiz,
                            "Ocean drag coefficient (horizontal) [-]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.atmosphere_drag_coeff_vert,
                            "Atmosphere drag coefficient (vertical) [-]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.atmosphere_drag_coeff_horiz,
                            "Atmosphere drag coefficient (horizontal) [-]")

    WriteVTK.vtk_point_data(vtkfile, ifarr.pressure,
                            "Contact pressure [Pa]")

    WriteVTK.vtk_point_data(vtkfile, ifarr.n_contacts,
                            "Number of contacts [-]")

    WriteVTK.vtk_point_data(vtkfile, ifarr.granular_stress,
                            "Granular stress [Pa]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.ocean_stress,
                            "Ocean stress [Pa]")
    WriteVTK.vtk_point_data(vtkfile, ifarr.atmosphere_stress,
                            "Atmosphere stress [Pa]")

    deleteGrainArrays!(ifarr)
    ifarr = 0
    gc()

    outfiles = WriteVTK.vtk_save(vtkfile)
    if verbose
        info("Output file: " * outfiles[1])
    end
    nothing
end

export writeGrainInteractionVTK
"""
    writeGrainInteractionVTK(simulation::Simulation,
                               filename::String;
                               verbose::Bool=false)

Saves grain interactions to `.vtp` files for visualization with VTK, for 
example in Paraview.  Convert Cell Data to Point Data and use with Tube filter.
"""
function writeGrainInteractionVTK(simulation::Simulation,
                                    filename::String;
                                    verbose::Bool=false)

    i1 = Int64[]
    i2 = Int64[]
    inter_particle_vector = Vector{Float64}[]
    force = Float64[]
    effective_radius = Float64[]
    contact_area = Float64[]
    contact_stiffness = Float64[]
    tensile_stress = Float64[]
    shear_displacement = Vector{Float64}[]
    contact_age = Float64[]
    for i=1:length(simulation.grains)
        for ic=1:simulation.Nc_max
            if simulation.grains[i].contacts[ic] > 0
                j = simulation.grains[i].contacts[ic]

                if !simulation.grains[i].enabled ||
                    !simulation.grains[j].enabled
                    continue
                end

                p = simulation.grains[i].lin_pos -
                    simulation.grains[j].lin_pos
                dist = norm(p)

                r_i = simulation.grains[i].contact_radius
                r_j = simulation.grains[j].contact_radius
                δ_n = dist - (r_i + r_j)
                R_ij = harmonicMean(r_i, r_j)

                if simulation.grains[i].youngs_modulus > 0. &&
                    simulation.grains[j].youngs_modulus > 0.
                    E_ij = harmonicMean(simulation.grains[i].
                                        youngs_modulus,
                                        simulation.grains[j].
                                        youngs_modulus)
                    A_ij = R_ij*min(simulation.grains[i].thickness, 
                                    simulation.grains[j].thickness)
                    k_n = E_ij*A_ij/R_ij
                else
                    k_n = harmonicMean(simulation.grains[i].
                                       contact_stiffness_normal,
                                       simulation.grains[j].
                                       contact_stiffness_normal)
                end

                
                push!(i1, i)
                push!(i2, j)
                push!(inter_particle_vector, p)

                push!(force, k_n*δ_n)
                push!(effective_radius, R_ij)
                push!(contact_area, A_ij)
                push!(contact_stiffness, k_n)
                push!(tensile_stress, k_n*δ_n/A_ij)

                push!(shear_displacement, simulation.grains[i].
                      contact_parallel_displacement[ic])

                push!(contact_age, simulation.grains[i].contact_age[ic])
            end
        end
    end

    # Insert a piece for each grain interaction using grain positions as 
    # coordinates and connect them with lines by referencing their indexes.
    open(filename * ".vtp", "w") do f
        write(f, "<?xml version=\"1.0\"?>\n")
        write(f, "<VTKFile type=\"PolyData\" version=\"0.1\" " *
              "byte_order=\"LittleEndian\">\n")
        write(f, "  <PolyData>\n")
        write(f, "    <Piece " *
              "NumberOfPoints=\"$(length(simulation.grains))\" " *
              "NumberOfVerts=\"0\" " *
              "NumberOfLines=\"$(length(i1))\" " *
              "NumberOfStrips=\"0\" " *
              "NumberOfPolys=\"0\">\n")
        write(f, "      <PointData>\n")
        write(f, "      </PointData>\n")
        write(f, "      <CellData>\n")

        # Write values associated to each line
        write(f, "        <DataArray type=\"Float32\" " *
              "Name=\"Inter-particle vector [m]\" " *
              "NumberOfComponents=\"3\" format=\"ascii\">\n")
        for i=1:length(i1)
            write(f, "$(inter_particle_vector[i][1]) ")
            write(f, "$(inter_particle_vector[i][2]) ")
            write(f, "0.0 ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")
        write(f, "        <DataArray type=\"Float32\" " *
              "Name=\"Shear displacement [m]\" " *
              "NumberOfComponents=\"3\" format=\"ascii\">\n")
        for i=1:length(i1)
            @inbounds write(f, "$(shear_displacement[i][1]) ")
            @inbounds write(f, "$(shear_displacement[i][2]) ")
            write(f, "0.0 ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")

        write(f, "        <DataArray type=\"Float32\" Name=\"Force [N]\" " *
              "NumberOfComponents=\"1\" format=\"ascii\">\n")
        for i=1:length(i1)
            @inbounds write(f, "$(force[i]) ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")

        write(f, "        <DataArray type=\"Float32\" " *
              "Name=\"Effective radius [m]\" " *
              "NumberOfComponents=\"1\" format=\"ascii\">\n")
        for i=1:length(i1)
            @inbounds write(f, "$(effective_radius[i]) ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")

        write(f, "        <DataArray type=\"Float32\" " *
              "Name=\"Contact area [m^2]\" " *
              "NumberOfComponents=\"1\" format=\"ascii\">\n")
        for i=1:length(i1)
            @inbounds write(f, "$(contact_area[i]) ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")

        write(f, "        <DataArray type=\"Float32\" " *
              "Name=\"Contact stiffness [N/m]\" " *
              "NumberOfComponents=\"1\" format=\"ascii\">\n")
        for i=1:length(i1)
            @inbounds write(f, "$(contact_stiffness[i]) ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")

        write(f, "        <DataArray type=\"Float32\" " *
              "Name=\"Tensile stress [Pa]\" " *
              "NumberOfComponents=\"1\" format=\"ascii\">\n")
        for i=1:length(i1)
            @inbounds write(f, "$(tensile_stress[i]) ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")

        write(f, "        <DataArray type=\"Float32\" " *
              "Name=\"Contact age [s]\" NumberOfComponents=\"1\" 
        format=\"ascii\">\n")
        for i=1:length(i1)
            @inbounds write(f, "$(contact_age[i]) ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")

        write(f, "      </CellData>\n")
        write(f, "      <Points>\n")

        # Write line endpoints (grain centers)
        #write(f, "        <DataArray Name=\"Position [m]\" type=\"Float32\" " *
        write(f, "        <DataArray type=\"Float32\" Name=\"Points\" " *
              "NumberOfComponents=\"3\" format=\"ascii\">\n")
        for i in simulation.grains
            @inbounds write(f, "$(i.lin_pos[1]) $(i.lin_pos[2]) 0.0 ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")
        write(f, "      </Points>\n")
        write(f, "      <Verts>\n")
        write(f, "        <DataArray type=\"Int64\" Name=\"connectivity\" " *
              "format=\"ascii\">\n")
        write(f, "        </DataArray>\n")
        write(f, "        <DataArray type=\"Int64\" Name=\"offsets\" " *
              "format=\"ascii\">\n")
        write(f, "        </DataArray>\n")
        write(f, "      </Verts>\n")
        write(f, "      <Lines>\n")

        # Write contact connectivity by referring to point indexes
        write(f, "        <DataArray type=\"Int64\" Name=\"connectivity\" " *
              "format=\"ascii\">\n")
        for i=1:length(i1)
            @inbounds write(f, "$(i1[i] - 1) $(i2[i] - 1) ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")
        
        # Write 0-indexed offset for the connectivity array for the end of each 
        # cell
        write(f, "        <DataArray type=\"Int64\" Name=\"offsets\" " *
              "format=\"ascii\">\n")
        for i=1:length(i1)
            @inbounds write(f, "$((i - 1)*2 + 2) ")
        end
        write(f, "\n")
        write(f, "        </DataArray>\n")

        write(f, "      </Lines>\n")
        write(f, "      <Strips>\n")
        write(f, "        <DataArray type=\"Int64\" Name=\"connectivity\" " *
              "format=\"ascii\">\n")
        write(f, "        </DataArray>\n")
        write(f, "        <DataArray type=\"Int64\" Name=\"offsets\" " *
              "format=\"ascii\">\n")
        write(f, "        </DataArray>\n")
        write(f, "      </Strips>\n")
        write(f, "      <Polys>\n")
        write(f, "        <DataArray type=\"Int64\" Name=\"connectivity\" " *
              "format=\"ascii\">\n")
        write(f, "        </DataArray>\n")
        write(f, "        <DataArray type=\"Int64\" Name=\"offsets\" " *
              "format=\"ascii\">\n")
        write(f, "        </DataArray>\n")
        write(f, "      </Polys>\n")
        write(f, "    </Piece>\n")
        write(f, "  </PolyData>\n")
        write(f, "</VTKFile>\n")
    end
    nothing
end

export writeOceanVTK
"""
Write a VTK file to disk containing all ocean data in the `simulation` in a 
structured grid (file type `.vts`).  These files can be read by ParaView and can 
be visualized by applying a *Glyph* filter.  This function is called by 
`writeVTK()`.
"""
function writeGridVTK(grid::Any,
                      filename::String;
                      verbose::Bool=false)
    
    # make each coordinate array three-dimensional
    xq = similar(grid.u[:,:,:,1])
    yq = similar(grid.u[:,:,:,1])
    zq = similar(grid.u[:,:,:,1])

    for iz=1:size(xq, 3)
        @inbounds xq[:,:,iz] = grid.xq
        @inbounds yq[:,:,iz] = grid.yq
    end
    for ix=1:size(xq, 1)
        for iy=1:size(xq, 2)
            @inbounds zq[ix,iy,:] = grid.zl
        end
    end

    # add arrays to VTK file
    vtkfile = WriteVTK.vtk_grid(filename, xq, yq, zq)

    WriteVTK.vtk_point_data(vtkfile, grid.u[:, :, :, 1],
                            "u: Zonal velocity [m/s]")
    WriteVTK.vtk_point_data(vtkfile, grid.v[:, :, :, 1],
                            "v: Meridional velocity [m/s]")
    # write velocities as 3d vector
    vel = zeros(3, size(xq, 1), size(xq, 2), size(xq, 3))
    for ix=1:size(xq, 1)
        for iy=1:size(xq, 2)
            for iz=1:size(xq, 3)
                @inbounds vel[1, ix, iy, iz] = grid.u[ix, iy, iz, 1]
                @inbounds vel[2, ix, iy, iz] = grid.v[ix, iy, iz, 1]
            end
        end
    end
    
    WriteVTK.vtk_point_data(vtkfile, vel, "Velocity vector [m/s]")

    if typeof(grid) == Ocean
        WriteVTK.vtk_point_data(vtkfile, grid.h[:, :, :, 1],
                                "h: Layer thickness [m]")
        WriteVTK.vtk_point_data(vtkfile, grid.e[:, :, :, 1],
                                "e: Relative interface height [m]")
    end

    outfiles = WriteVTK.vtk_save(vtkfile)
    if verbose
        info("Output file: " * outfiles[1])
    end
    nothing
end

export writeParaviewPythonScript
"""
function writeParaviewPythonScript(simulation,
                                   [filename, folder, vtk_folder, verbose])

Create a `".py"` script for visualizing the simulation VTK files in Paraview.
The script can be run from the command line with `pvpython` (bundled with
Paraview), or from the interactive Python shell inside Paraview.

# Arguments
* `simulation::Simulation`: input simulation file containing the data.
* `filename::String`: output file name for the Python script. At its default
    (blank) value, the script is named after the simulation id (`simulation.id`).
* `folder::String`: output directory, current directory the default.
* `vtk_folder::String`: directory containing the VTK output files, by default
    points to the full system path equivalent to `"./<simulation.id>/"`.
* `save_animation::Bool`: make the generated script immediately save a rendered
    animation to disk when the `".py"` script is called.
* `verbose::Bool`: show diagnostic information during
function call, on by
    default.
"""
function writeParaviewPythonScript(simulation::Simulation;
                                   filename::String="",
                                   folder::String=".",
                                   vtk_folder::String="",
                                   save_animation::Bool=true,
                                   save_images::Bool=false,
                                   width::Integer=1920,
                                   height::Integer=1080,
                                   framerate::Integer=10,
                                   grains_color_scheme::String="X Ray",
                                   verbose::Bool=true)
    if filename == ""
        folder = string(folder, "/", simulation.id)
        mkpath(folder)
        filename = string(folder, "/", simulation.id, ".py")
    end
    if vtk_folder == ""
        vtk_folder = string(pwd(), "/", simulation.id)
    end

    if simulation.file_number == 0
        simulation.file_number = readSimulationStatus(simulation)
    end

    open(filename, "w") do f
        write(f, """from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()
FileName=[""")
        for i=1:simulation.file_number
            write(f, "'$(vtk_folder)/$(simulation.id).grains.$(i).vtu', ")
        end
        write(f, """]
imagegrains = XMLUnstructuredGridReader(FileName=FileName)

imagegrains.PointArrayStatus = [
'Density [kg m^-3]',
'Thickness [m]',
'Diameter (contact) [m]',
'Diameter (areal) [m]',
'Circumreference  [m]',
'Horizontal surface area [m^2]',
'Side surface area [m^2]',
'Volume [m^3]',
'Mass [kg]',
'Moment of inertia [kg m^2]',
'Linear velocity [m s^-1]',
'Linear acceleration [m s^-2]',
'Sum of forces [N]',
'Angular position [rad]',
'Angular velocity [rad s^-1]',
'Angular acceleration [rad s^-2]',
'Sum of torques [N*m]',
'Fixed in space [-]',
'Free to rotate [-]',
'Enabled [-]',
'Contact stiffness (normal) [N m^-1]',
'Contact stiffness (tangential) [N m^-1]',
'Contact viscosity (normal) [N m^-1 s]',
'Contact viscosity (tangential) [N m^-1 s]',
'Contact friction (static) [-]',
'Contact friction (dynamic) [-]',
"Young's modulus [Pa]",
"Poisson's ratio [-]",
'Tensile strength [Pa]'
'Compressive strength prefactor [m^0.5 Pa]',
'Ocean drag coefficient (vertical) [-]',
'Ocean drag coefficient (horizontal) [-]',
'Atmosphere drag coefficient (vertical) [-]',
'Atmosphere drag coefficient (horizontal) [-]',
'Contact pressure [Pa]',
'Number of contacts [-]',
'Granular stress [Pa]',
'Ocean stress [Pa]',
'Atmosphere stress [Pa]']

animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [2478, 1570]

# show data in view
imagegrainsDisplay = Show(imagegrains, renderView1)
# trace defaults for the display properties.
imagegrainsDisplay.Representation = 'Surface'
imagegrainsDisplay.AmbientColor = [0.0, 0.0, 0.0]
imagegrainsDisplay.ColorArrayName = [None, '']
imagegrainsDisplay.OSPRayScaleArray = 'Angular acceleration [rad s^-2]'
imagegrainsDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
imagegrainsDisplay.SelectOrientationVectors = 'Angular acceleration [rad s^-2]'
imagegrainsDisplay.ScaleFactor = 6.050000000000001
imagegrainsDisplay.SelectScaleArray = 'Angular acceleration [rad s^-2]'
imagegrainsDisplay.GlyphType = 'Arrow'
imagegrainsDisplay.GlyphTableIndexArray = 'Angular acceleration [rad s^-2]'
imagegrainsDisplay.DataAxesGrid = 'GridAxesRepresentation'
imagegrainsDisplay.PolarAxes = 'PolarAxesRepresentation'
imagegrainsDisplay.ScalarOpacityUnitDistance = 64.20669746996803
imagegrainsDisplay.GaussianRadius = 3.0250000000000004
imagegrainsDisplay.SetScaleArray = ['POINTS', 'Atmosphere drag coefficient (horizontal) [-]']
imagegrainsDisplay.ScaleTransferFunction = 'PiecewiseFunction'
imagegrainsDisplay.OpacityArray = ['POINTS', 'Atmosphere drag coefficient (horizontal) [-]']
imagegrainsDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
imagegrainsDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
imagegrainsDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
imagegrainsDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
imagegrainsDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
imagegrainsDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
imagegrainsDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
imagegrainsDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
imagegrainsDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
imagegrainsDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
imagegrainsDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
imagegrainsDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Glyph'
glyph1 = Glyph(Input=imagegrains,
    GlyphType='Arrow')
glyph1.Scalars = ['POINTS', 'Atmosphere drag coefficient (horizontal) [-]']
glyph1.Vectors = ['POINTS', 'Angular acceleration [rad s^-2]']
glyph1.ScaleFactor = 6.050000000000001
glyph1.GlyphTransform = 'Transform2'

# Properties modified on glyph1
glyph1.Scalars = ['POINTS', 'Diameter (areal) [m]']
glyph1.Vectors = ['POINTS', 'Angular position [rad]']
glyph1.ScaleMode = 'scalar'
glyph1.ScaleFactor = 1.0
glyph1.GlyphMode = 'All Points'

# get color transfer function/color map for 'Diameterarealm'
diameterarealmLUT = GetColorTransferFunction('Diameterarealm')

# show data in view
glyph1Display = Show(glyph1, renderView1)
# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'
glyph1Display.AmbientColor = [0.0, 0.0, 0.0]
glyph1Display.ColorArrayName = ['POINTS', 'Diameter (areal) [m]']
glyph1Display.LookupTable = diameterarealmLUT
glyph1Display.OSPRayScaleArray = 'Diameter (areal) [m]'
glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph1Display.SelectOrientationVectors = 'GlyphVector'
glyph1Display.ScaleFactor = 6.1000000000000005
glyph1Display.SelectScaleArray = 'Diameter (areal) [m]'
glyph1Display.GlyphType = 'Arrow'
glyph1Display.GlyphTableIndexArray = 'Diameter (areal) [m]'
glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
glyph1Display.PolarAxes = 'PolarAxesRepresentation'
glyph1Display.GaussianRadius = 3.0500000000000003
glyph1Display.SetScaleArray = ['POINTS', 'Diameter (areal) [m]']
glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph1Display.OpacityArray = ['POINTS', 'Diameter (areal) [m]']
glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
glyph1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
glyph1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
glyph1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
glyph1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
glyph1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
glyph1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
glyph1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
glyph1Display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
glyph1Display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
glyph1Display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
glyph1Display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# reset view to fit data
renderView1.ResetCamera()

# Properties modified on glyph1
glyph1.GlyphType = 'Sphere'

# update the view to ensure updated data information
renderView1.Update()

# hide color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, False)

# rescale color and/or opacity maps used to exactly fit the current data range
glyph1Display.RescaleTransferFunctionToDataRange(False, True)

# get opacity transfer function/opacity map for 'Diameterarealm'
diameterarealmPWF = GetOpacityTransferFunction('Diameterarealm')

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
diameterarealmLUT.ApplyPreset('$(grains_color_scheme)', True)

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
""")
        if save_animation
            write(f, """
SaveAnimation('$(vtk_folder)/$(simulation.id).avi', renderView1,
ImageResolution=[$(width), $(height)],
FrameRate=$(framerate),
FrameWindow=[0, $(simulation.file_number)])
""")
        end

        if save_images
            write(f, """
SaveAnimation('$(folder)/$(simulation.id).png', renderView1,
ImageResolution=[$(width), $(height)],
FrameRate=$(framerate),
FrameWindow=[0, $(simulation.file_number)])
""")
        end
    end
    if verbose
        info("$(filename) written, execute with " *
             "'pvpython $(vtk_folder)/$(simulation.id).py'")
    end
end

export render
"""
    render(simulation[, pvpython, images, animation])

Wrapper function which calls `writeParaviewPythonScript(...)` and executes it
from the shell using the supplied `pvpython` argument.

# Arguments
* `simulation::Simulation`: simulation object containing the grain data.
* `pvpython::String`: path to the `pvpython` executable to use.  By default, the
    script uses the pvpython in the system PATH.
* `images::Bool`: render images to disk (default: true)
* `gif::Bool`: merge images as GIF and save to disk (default: false, requires
    `images=true`)
* `animation::Bool`: render animation as AVI to disk (default: false)
* `trim::Bool`: trim images in animated sequence (default: true)
* `reverse::Bool`: if `images=true` additionally render reverse-animated gif
    (default: false)
"""
function render(simulation::Simulation; pvpython::String="pvpython",
                images::Bool=true,
                gif::Bool=false,
                animation::Bool=false,
                trim::Bool=true,
                reverse::Bool=false)

    writeParaviewPythonScript(simulation, save_animation=animation,
                              save_images=images, verbose=false)
    try
        run(`$(pvpython) $(simulation.id)/$(simulation.id).py`)

        # if available, use imagemagick to create gif from images
        if images && gif
            try
                trim_string = ""
                if trim
                    trim_string = "-trim"
                end

                # use ImageMagick installed with Homebrew.jl if available,
                # otherwise search for convert in $PATH
                convert = "convert"

                run(`$convert $trim_string +repage -delay 10 
                    -transparent-color white 
                    -loop 0 $(simulation.id)/$(simulation.id).'*'.png 
                    $(simulation.id)/$(simulation.id).gif`)
                if reverse
                    run(`$convert -trim +repage -delay 10 -transparent-color white 
                        -loop 0 -reverse
                        $(simulation.id)/$(simulation.id).'*'.png 
                        $(simulation.id)/$(simulation.id)-reverse.gif`)
                end
            catch return_signal
                if isa(return_signal, Base.UVError)
                    info("Skipping gif merge since `$convert` was not found.")
                end
            end
        end
    catch return_signal
        if isa(return_signal, Base.UVError)
            error("`pvpython` was not found.")
        end
    end
end

export removeSimulationFiles
"""
    removeSimulationFiles(simulation[, folder])

Remove all simulation output files from the specified folder.
"""
function removeSimulationFiles(simulation::Simulation; folder::String=".")
    folder = folder * "/" * simulation.id
    run(`bash -c "rm -rf $(folder)/$(simulation.id).*.vtu"`)
    run(`bash -c "rm -rf $(folder)/$(simulation.id).*.vtp"`)
    run(`bash -c "rm -rf $(folder)/$(simulation.id).*.vts"`)
    run(`bash -c "rm -rf $(folder)/$(simulation.id).status.txt"`)
    run(`bash -c "rm -rf $(folder)/$(simulation.id).*.jld"`)
    run(`bash -c "rm -rf $(folder)/$(simulation.id).py"`)
    run(`bash -c "rm -rf $(folder)/$(simulation.id).avi"`)
    run(`bash -c "rm -rf $(folder)/$(simulation.id).*.png"`)
    run(`bash -c "rm -rf $(folder)"`)
    nothing
end

export plotGrainSizeDistribution
"""
    plotGrainSizeDistribution(simulation, [filename_postfix, nbins,
                                size_type, filetype, gnuplot_terminal,
                                skip_fixed, log_y, verbose)

Plot the grain size distribution as a histogram and save it to the disk.  The 
plot is saved accoring to the simulation id, the optional `filename_postfix` 
string, and the `filetype`, and is written to the current folder.

# Arguments
* `simulation::Simulation`: the simulation object containing the grains.
* `filename_postfix::String`: optional string for the output filename.
* `nbins::Int`: number of bins in the histogram (default = 12).
* `size_type::String`: specify whether to use the `contact` or `areal` radius 
    for the grain size.  The default is `contact`.
* `filetype::String`: the output file type (default = "png").
* `gnuplot_terminal::String`: the gnuplot output terminal to use (default =
    "png").
* `skip_fixed::Bool`: ommit grains that are fixed in space from the size 
    distribution (default = true).
* `log_y::Bool`: plot y-axis in log scale.
* `verbose::String`: show output file as info message in stdout (default = 
    true).
"""
function plotGrainSizeDistribution(simulation::Simulation;
                                     filename_postfix::String = "",
                                     nbins::Int=12,
                                     size_type::String = "contact",
                                     filetype::String = "png",
                                     gnuplot_terminal::String = "png",
                                     skip_fixed::Bool = true,
                                     log_y::Bool = false,
                                     verbose::Bool = true)

    diameters = Float64[]
    for i=1:length(simulation.grains)
        if simulation.grains[i].fixed && skip_fixed
            continue
        end
        if size_type == "contact"
            push!(diameters, simulation.grains[i].contact_radius*2.)
        elseif size_type == "areal"
            push!(diameters, simulation.grains[i].areal_radius*2.)
        else
            error("size_type '$size_type' not understood")
        end
    end

    filename = string(simulation.id * filename_postfix * 
                      "-grain-size-distribution." * filetype)

    # write data to temporary file on disk
    datafile = Base.Filesystem.tempname()
    writedlm(datafile, diameters)
    gnuplotscript = Base.Filesystem.tempname()

    #if maximum(diameters) ≈ minimum(diameters)
        #info("Overriding `nbins = $nbins` -> `nbins = 1`.")
        #nbins = 1
    #end

    open(gnuplotscript, "w") do f

        write(f, """#!/usr/bin/env gnuplot
              set term $gnuplot_terminal
              set out "$(filename)"\n""")
        if log_y
            write(f, "set logscale y\n")
        end
        write(f, """set xlabel "Diameter [m]"
              set ylabel "Count [-]"
              binwidth = $((maximum(diameters) - minimum(diameters)+1e-7)/nbins)
              binstart = $(minimum(diameters))
              set boxwidth 1.0*binwidth
              set style fill solid 0.5
              set key off
              hist = 'u (binwidth*(floor((\$1-binstart)/binwidth)+0.5)+binstart):(1.0) smooth freq w boxes'
              plot "$(datafile)" i 0 @hist ls 1
              """)
    end

    try
        run(`gnuplot $gnuplotscript`)
    catch return_signal
        if isa(return_signal, Base.UVError)
            error("Could not launch external gnuplot process")
        end
    end

    if verbose
        info(filename)
    end
end

export plotGrains
"""
    plotGrains(simulation, [filetype, gnuplot_terminal, verbose])

Plot the grains using Gnuplot and save the figure to disk.

# Arguments
* `simulation::Simulation`: the simulation object containing the grains.
* `filetype::String`: the output file type (default = "png").
* `gnuplot_terminal::String`: the gnuplot output terminal to use (default =
    "png crop size 1200,1200").
* `verbose::String`: show output file as info message in stdout (default = 
    true).
"""
function plotGrains(sim::Simulation;
                    filetype::String = "png",
                    gnuplot_terminal::String = "png crop size 1200,1200",
                    show_figure::Bool = true,
                    verbose::Bool = true)

    mkpath(sim.id)
    filename = string(sim.id, "/", sim.id, ".grains.", sim.file_number, ".",
                      filetype)

    x = Float64[]
    y = Float64[]
    r = Float64[]
    for grain in sim.grains
        push!(x, grain.lin_pos[1])
        push!(y, grain.lin_pos[2])
        push!(r, grain.contact_radius)
    end

    # write data to temporary file on disk
    datafile = Base.Filesystem.tempname()
    writedlm(datafile, [x y r])
    gnuplotscript = Base.Filesystem.tempname()

    open(gnuplotscript, "w") do f

        write(f, """#!/usr/bin/env gnuplot
              set term $(gnuplot_terminal)
              set out "$(filename)"
              set palette defined (0 "gray", 1 "white")
              set xlabel "x [m]"
              set ylabel "y [m]"\n""")
        if typeof(sim.ocean.input_file) != Bool
            write(f, "set xrange ")
            write(f, "[$(sim.ocean.xq[1,1]):$(sim.ocean.xq[end,end])]\n")
            write(f, "set yrange ")
            write(f, "[$(sim.ocean.yq[1,1]):$(sim.ocean.yq[end,end])]\n")
        else
            write(f, "set xrange [$(minimum(x - r)):$(maximum(x + r))]\n")
            write(f, "set yrange [$(minimum(y - r)):$(maximum(y + r))]\n")
        end
        write(f, """set cblabel "Diameter [m]"
              set size ratio -1
              set key off
              plot "$(datafile)" with circles lt 1 lc rgb "black" t "Particle"
              """)
    end

    try
        run(`gnuplot $gnuplotscript`)
    catch return_signal
        if isa(return_signal, Base.UVError)
            error("Could not launch external gnuplot process")
        end
    end

    if verbose
        info(filename)
    end

    if show_figure
        if is_apple()
            run(`open $(filename)`)
        elseif is_linux()
            run(`xdg-open $(filename)`)
        end
    end
end
