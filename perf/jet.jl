ca_dir = joinpath(dirname(@__DIR__));
include(joinpath(ca_dir, "examples", "hybrid", "cli_options.jl"));

ENV["CI_PERF_SKIP_RUN"] = true # we only need haskey(ENV, "CI_PERF_SKIP_RUN") == true

import RRTMGP
import ClimaComms
import CUDA

RRTMGP.Device.array_type() = CUDA.has_cuda_gpu() ? CUDA.CuArray : Array
ClimaCore.Device.device() = ClimaComms.CPU()

filename = joinpath(ca_dir, "examples", "hybrid", "driver.jl")
dict = parsed_args_per_job_id(; trigger = "benchmark.jl")
parsed_args_prescribed = parsed_args_from_ARGS(ARGS)

# Start with performance target, but override anything provided in ARGS
parsed_args_target = dict["perf_target_unthreaded"];
parsed_args = merge(parsed_args_target, parsed_args_prescribed);

try # capture integrator
    include(filename)
catch err
    if err.error !== :exit_profile
        rethrow(err.error)
    end
end

import JET

OrdinaryDiffEq.step!(integrator) # Make sure no errors
# JET.@test_opt OrdinaryDiffEq.step!(integrator)
JET.@report_opt OrdinaryDiffEq.step!(integrator)
