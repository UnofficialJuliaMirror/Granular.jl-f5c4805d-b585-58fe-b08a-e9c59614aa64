#!/usr/bin/env julia
import Compat

Compat.@info "Determining if JLD is installed"
if typeof(Compat.Pkg.installed("JLD")) == VersionNumber
    Compat.@info "JLD found, proceeding with JLD-specific tests"

    Compat.@info "Writing simple simulation to JLD file"
    sim = Granular.createSimulation(id="test")
    Granular.addGrainCylindrical!(sim, [ 0., 0.], 10., 1., verbose=false)
    Granular.addGrainCylindrical!(sim, [18., 0.], 10., 1., verbose=false)
    sim.ocean = Granular.createRegularOceanGrid([10, 20, 5], [10., 25., 2.])  
    Granular.findContacts!(sim, method="all to all")
    Granular.writeVTK(sim, verbose=false)

    Granular.writeSimulation(sim)
    Granular.writeSimulationStatus(sim)

    Compat.@info "Reading from JLD file by specifying the input file name"
    sim2 = Granular.readSimulation("./test/test.1.jld")
    Granular.compareSimulations(sim, sim2)

    Compat.@info "Reading and overwriting from JLD file by simulation id"
    sim3 = Granular.createSimulation("test")
    @test 1 == Granular.readSimulationStatus(sim3)
    sim3 = Granular.readSimulation(sim3)
    Granular.compareSimulations(sim, sim3)

    rm("./test/test.1.jld")
end
