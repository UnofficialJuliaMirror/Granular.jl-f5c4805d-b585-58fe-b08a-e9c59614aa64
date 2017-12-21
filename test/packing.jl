#!/usr/bin/env julia
using Compat.Test
import Granular

verbose = true

info("#### $(basename(@__FILE__)) ####")

info("Testing regular packing generation (power law GSD)")
sim = Granular.createSimulation()
Granular.regularPacking!(sim, [2, 2], 1., 1., size_distribution="powerlaw")
@test 4 == length(sim.grains)
for grain in sim.grains
    @test grain.contact_radius ≈ 1.
end

sim = Granular.createSimulation()
Granular.regularPacking!(sim, [10, 10], 1., 10., size_distribution="powerlaw")
@test 100 == length(sim.grains)
for grain in sim.grains
    @test grain.contact_radius >= 1.
    @test grain.contact_radius <= 10.
end

info("Testing regular packing generation (uniform GSD)")
sim = Granular.createSimulation()
Granular.regularPacking!(sim, [2, 2], 1., 1., size_distribution="uniform")
@test 4 == length(sim.grains)
for grain in sim.grains
    @test grain.contact_radius ≈ 1.
end

sim = Granular.createSimulation()
Granular.regularPacking!(sim, [10, 10], 1., 10., size_distribution="uniform")
@test 100 == length(sim.grains)
for grain in sim.grains
    @test grain.contact_radius >= 1.
    @test grain.contact_radius <= 10.
end


info("Testing irregular (Poisson-disk) packing generation (monodisperse size)")
sim = Granular.createSimulation("poisson1-monodisperse-nopadding")
sim.ocean = Granular.createRegularOceanGrid([1, 1, 1], [1., 1., 1.])
Granular.irregularPacking!(sim,
                           radius_max=.1,
                           radius_min=.1,
                           padding_factor=0.,
                           verbose=true)

info("Testing irregular (Poisson-disk) packing generation (wide PSD)")
sim = Granular.createSimulation("poisson2-wide-nopadding")
sim.ocean = Granular.createRegularOceanGrid([1, 1, 1], [1., 1., 1.])
Granular.irregularPacking!(sim,
                           radius_max=.1,
                           radius_min=.001,
                           padding_factor=0.,
                           verbose=true)
sim = Granular.createSimulation("poisson3-wide-padding")
sim.ocean = Granular.createRegularOceanGrid([1, 1, 1], [1., 1., 1.])
Granular.irregularPacking!(sim,
                           radius_max=.1,
                           radius_min=.001,
                           padding_factor=2.,
                           verbose=true)
