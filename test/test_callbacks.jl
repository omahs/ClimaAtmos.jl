# Generate simple model using shallow water equations for bickley jet problem
include("test_cases/initial_conditions/bickley_jet_2d_plane.jl");
FT = Float64;
npolynomial = 4;
nelements = (2,2);
params = map(FT, (
               g = 9.8,  # gravitational constant
               D₄ = 1e-4,  # hyperdiffusion constant
               ϵ = 0.1,  # perturbation size for initial condition
               l = 0.5,  # Gaussian width
               k = 0.5,  # sinusoidal wavenumber
               h₀ = 1.0,  # reference density
           ));

@unpack h, u, c = init_bickley_jet_2d_plane(params);
domain = PeriodicPlane(
               FT,
               xlim = (-2π, 2π),
               ylim = (-2π, 2π),
               nelements = nelements,
               npolynomial = npolynomial,
           );
model = ShallowWaterModel(domain = domain, parameters = params);

# Begin Tests
cb_1 = JLD2Output(model, "TempTestDir1", "TestFilename1", 1);
cb_2 = JLD2Output(model, "TempTestDir2", "TestFilename2", 2);

@test generate_callback(cb_1) isa DiffEqBase.DiscreteCallback
@test generate_callback(cb_2) isa DiffEqBase.DiscreteCallback
@test DiffEqBase.CallbackSet(cb1,cb2) isa DiffEqBase.CallbackSet
@test isfile(joinpath(@__DIR__, cb_1.filedir, cb_1.filename*".jl")) == false
@test isfile(joinpath(@__DIR__, cb_2.filedir, cb_2.filename*".jl")) == false

# Generate simple simulation data for test
simulation = Simulation(model, stepper, dt = 0.01, tspan = (0.0,0.04))

include("initial_conditions/bickley_jet_2d_plane.jl")

