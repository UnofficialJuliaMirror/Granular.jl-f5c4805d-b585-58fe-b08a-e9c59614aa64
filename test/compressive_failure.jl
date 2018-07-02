#!/usr/bin/env julia
using Compat.Test
import Granular

verbose = true

Compat.@info "Testing compressive failure"
sim = Granular.createSimulation()
Granular.addGrainCylindrical!(sim, [0.,0.], 1., 1., compressive_strength=1.,
                              lin_vel=[1., 0.], fixed=true)
Granular.addGrainCylindrical!(sim, [2.,0.], 1., 1., compressive_strength=1.,
                              fixed=true)
@test count(x->x==true, sim.grains[1].compressive_failure) == 0
Granular.setTimeStep!(sim, verbose=verbose)
Granular.setTotalTime!(sim, 1.0)
Granular.run!(sim, verbose=verbose)
@test sim.grains[1].compressive_failure[1] == true
@test count(x->x==true, sim.grains[1].compressive_failure) == 1


