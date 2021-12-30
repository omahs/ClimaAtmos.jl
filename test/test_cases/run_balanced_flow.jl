include("initial_conditions/dry_shallow_baroclinic_wave.jl")

# Set up parameters
using CLIMAParameters
struct BalancedFlowParameters <: CLIMAParameters.AbstractEarthParameterSet end

function run_balanced_flow(
    ::Type{FT};
    stepper = SSPRK33(),
    nelements = (6, 10),
    npolynomial = 3,
    case = :default,
    dt = 0.02,
    callbacks = (),
    mode = :regression,
) where {FT}
    params = BalancedFlowParameters()

    domain = SphericalShell(
        FT,
        radius = CLIMAParameters.Planet.planet_radius(params),
        height = FT(30.0e3),
        nelements = nelements,
        npolynomial = npolynomial,
    )

    model = Nonhydrostatic3DModel(
        domain = domain,
        boundary_conditions = nothing,
        parameters = params,
        hyperdiffusivity = FT(100),
    )

    # execute differently depending on testing mode
    if mode == :integration
        simulation = Simulation(model, stepper, dt = dt, tspan = (0.0, 1.0))
        @test simulation isa Simulation

        # test set function
        @unpack ρ, uh, w, ρe_tot =
            init_dry_shallow_baroclinic_wave(FT, params, isbalanced = true)
        set!(simulation, ρ = ρ, uh = uh, w = w, ρe_tot = ρe_tot)

        # test successful integration
        @test step!(simulation) isa Nothing # either error or integration runs
    elseif mode == :regression
        # TODO!: Implement meaningful(!) regression test
    elseif mode == :validation
        # TODO!: Implement the rest plots and analyses
        # 1. sort out saveat kwarg for Simulation
        # 2. create animation for a rising bubble; timeseries of total energy
    else
        throw(ArgumentError("$mode incompatible with test case."))
    end

    nothing
end
