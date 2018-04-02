#!/usr/bin/env julia

# Check the basic icefloe functionality

Compat.@info "#### $(basename(@__FILE__)) ####"

Compat.@info "Writing simple simulation to VTK file"
sim = Granular.createSimulation(id="test")
Granular.addGrainCylindrical!(sim, [ 0., 0.], 10., 1., verbose=false)
Granular.printGrainInfo(sim.grains[1])

Compat.@info "Testing grain value checks "
@test_throws ErrorException Granular.addGrainCylindrical!(sim, [.1, .1, .1],
                                                          10., 1.)
@test_throws ErrorException Granular.addGrainCylindrical!(sim, [.1, .1],
                                                          10., 1., 
                                                          lin_vel=[.2,.2,.2])
@test_throws ErrorException Granular.addGrainCylindrical!(sim, [.1, .1],
                                                          10., 1., 
                                                          lin_acc=[.2,.2,.2])
@test_throws ErrorException Granular.addGrainCylindrical!(sim, [.1, .1],
                                                          0., 1.)
@test_throws ErrorException Granular.addGrainCylindrical!(sim, [.1, .1],
                                                          10., 1., density=-2.)
@test_throws ErrorException Granular.disableGrain!(sim, 0)

Compat.@info "Testing grain comparison "
sim = Granular.createSimulation(id="test")
Granular.addGrainCylindrical!(sim, [ 0., 0.], 10., 1., verbose=false)
Granular.addGrainCylindrical!(sim, [ 0., 0.], 10., 1., verbose=false)
Granular.compareGrains(sim.grains[1], sim.grains[2])

gnuplot = true
try
    run(`gnuplot --version`)
catch return_signal
    if isa(return_signal, Base.UVError)
        Compat.@warn "Skipping plotting routines: Could not launch gnuplot process"
        gnuplot = false
    end
end
if gnuplot
    Compat.@info "Testing GSD plotting "
    Granular.plotGrainSizeDistribution(sim)
    @test isfile("test-grain-size-distribution.png")
    rm("test-grain-size-distribution.png")
    Granular.plotGrainSizeDistribution(sim, skip_fixed=false)
    @test isfile("test-grain-size-distribution.png")
    rm("test-grain-size-distribution.png")
    Granular.plotGrainSizeDistribution(sim, size_type="areal")
    @test isfile("test-grain-size-distribution.png")
    rm("test-grain-size-distribution.png")
    @test_throws ErrorException Granular.plotGrainSizeDistribution(sim, size_type="asdf")
end

Compat.@info "Testing external body force routines"
sim = Granular.createSimulation(id="test")
Granular.addGrainCylindrical!(sim, [ 0., 0.], 10., 1., verbose=false)
Granular.setBodyForce!(sim.grains[1], [1., 2.])
@test sim.grains[1].external_body_force ≈ [1., 2.]
Granular.addBodyForce!(sim.grains[1], [1., 2.])
@test sim.grains[1].external_body_force ≈ [2., 4.]

Compat.@info "Testing zeroKinematics!()"
sim.grains[1].force .= ones(2)
sim.grains[1].lin_acc .= ones(2)
sim.grains[1].lin_vel .= ones(2)
sim.grains[1].torque = 1.
sim.grains[1].ang_acc = 1.
sim.grains[1].ang_vel = 1.
Granular.zeroKinematics!(sim)
@test Granular.totalGrainKineticTranslationalEnergy(sim) ≈ 0.
@test Granular.totalGrainKineticRotationalEnergy(sim) ≈ 0.
