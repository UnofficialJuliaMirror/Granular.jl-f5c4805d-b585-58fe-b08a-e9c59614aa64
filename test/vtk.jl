#!/usr/bin/env julia
import Compat

# Check the contact search and geometry of a two-particle interaction

Compat.@info "#### $(basename(@__FILE__)) ####"

Compat.@info "Writing simple simulation to VTK file"
sim = Granular.createSimulation(id="test")
Granular.addGrainCylindrical!(sim, [ 0., 0.], 10., 1., verbose=false)
Granular.addGrainCylindrical!(sim, [18., 0.], 10., 1., verbose=false)
sim.ocean = Granular.createRegularOceanGrid([10, 20, 5], [10., 25., 2.])  
Granular.findContacts!(sim, method="all to all")
Granular.writeVTK(sim, verbose=false)

cmd_post = ""
if Compat.Sys.islinux()
    cmd = "sha256sum"
elseif Compat.Sys.isapple()
    cmd = ["shasum", "-a", "256"]
elseif Compat.Sys.iswindows()
    Compat.@info "checksum verification not yet implemented on Windows"
    exit()
    cmd = ["powershell", "-Command", "\"Get-FileHash", "-Algorithm", "SHA256"]
    cmd_post = "\""
else
    error("checksum verification of VTK file not supported on this platform")
end

grainpath = "test/test.grains.1.vtu"
grainchecksum = 
"eff130f14975dd5da06186c5b61c0c1a9d679f2188d03019c8144898ddc902b6  " *
grainpath * "\n"

graininteractionpath = "test/test.grain-interaction.1.vtp"
graininteractionchecksum = 
"c61f314a997c4405eaa98e6a4e3f81368ab2e1c60d06d9a0d998b69d7c8fb1bf  " *
graininteractionpath * "\n"

oceanpath = "test/test.ocean.1.vts"
oceanchecksum =
"b65f00942f1cbef7335921948c9eb73d137574eb806c33dea8b0e9b638665f3b  " *
oceanpath * "\n"

@test read(`$(cmd) $(grainpath)$(cmd_post)`, String) == grainchecksum
@test read(`$(cmd) $(graininteractionpath)$(cmd_post)`, String) == 
    graininteractionchecksum
@test read(`$(cmd) $(oceanpath)$(cmd_post)`, String) == oceanchecksum

Granular.removeSimulationFiles(sim)

Compat.@info "Testing VTK write during run!()"
Granular.setOutputFileInterval!(sim, 1e-9)
Granular.setTotalTime!(sim, 1.5)
Granular.setTimeStep!(sim)
sim.file_number = 0
Granular.run!(sim, single_step=true)
@test Granular.readSimulationStatus(sim.id) == 1
Granular.setOutputFileInterval!(sim, 0.1)
Granular.run!(sim)

Granular.status()

Compat.@info "Testing generation of Paraview Python script"
Granular.writeParaviewPythonScript(sim,
                                 save_animation=true,
                                 save_images=false)
@test isfile("$(sim.id)/$(sim.id).py") && filesize("$(sim.id)/$(sim.id).py") > 0

Compat.@info "Testing Paraview rendering if `pvpython` is present"
try
    run(`pvpython $(sim.id)/$(sim.id).py`)
catch return_signal
    if !isa(return_signal, Base.UVError)
        @test isfile("$(sim.id)/$(sim.id).avi")
    end
end

Granular.writeParaviewPythonScript(sim,
                                 save_animation=false,
                                 save_images=true)
try
    run(`pvpython $(sim.id)/$(sim.id).py`)
catch return_signal
    if !isa(return_signal, Base.UVError)
        @test isfile("$(sim.id)/$(sim.id).0000.png")
        @test isfile("$(sim.id)/$(sim.id).0014.png")
        Granular.render(sim)
        @test isfile("$(sim.id)/$(sim.id).0001.png")
    end
end

@test read(`$(cmd) $(grainpath)$(cmd_post)`, String) == grainchecksum
@test read(`$(cmd) $(graininteractionpath)$(cmd_post)`, String) == 
    graininteractionchecksum
@test read(`$(cmd) $(oceanpath)$(cmd_post)`, String) == oceanchecksum

Granular.removeSimulationFiles(sim)
