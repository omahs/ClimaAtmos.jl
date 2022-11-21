ca_dir = joinpath(dirname(@__DIR__));
include(joinpath(ca_dir, "examples", "hybrid", "cli_options.jl"));

ENV["CI_PERF_SKIP_INIT"] = true # we only need haskey(ENV, "CI_PERF_SKIP_INIT") == true

filename = joinpath(ca_dir, "examples", "hybrid", "driver.jl")
dict = parsed_args_per_job_id(; trigger = "driver.jl")
parsed_args_prescribed = parsed_args_from_ARGS(ARGS)

# Start with performance target, but override anything provided in ARGS
parsed_args_target = dict["compressible_edmf_bomex_jfnk"];
parsed_args = merge(parsed_args_target, parsed_args_prescribed);
parsed_args["enable_threading"] = false;

# The callbacks flame graph is very expensive, so only do 2 steps.
const n_samples = occursin("callbacks", parsed_args["job_id"]) ? 2 : 20

try # capture integrator
    include(filename)
catch err
    if err.error !== :exit_profile_init
        rethrow(err.error)
    end
end

function do_work!(integrator_args, integrator_kwargs)
    for _ in 1:2
        integrator = get_integrator(integrator_args, integrator_kwargs)
    end
    return nothing
end

do_work!(integrator_args, integrator_kwargs) # compile first
import Profile
Profile.clear_malloc_data()
Profile.clear()
prof = Profile.@profile begin
    do_work!(integrator_args, integrator_kwargs)
end

(; output_dir, job_id) = simulation

import ProfileCanvas

if haskey(ENV, "BUILDKITE_COMMIT") || haskey(ENV, "BUILDKITE_BRANCH")
    output_dir = job_id
    mkpath(output_dir)
    ProfileCanvas.html_file(joinpath(output_dir, "flame.html"))
else
    ProfileCanvas.view(Profile.fetch())
end
