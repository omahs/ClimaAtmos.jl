import Random
Random.seed!(1234)
import ClimaAtmos as CA

import JET

config = CA.AtmosConfig(;
    parsed_args = CA.AtmosTargetParsedArgs(; target_job = "edmfx_adv_test_box"),
) # make sure no errors

JET.@test_opt CA.AtmosConfig(;
    parsed_args = CA.AtmosTargetParsedArgs(; target_job = "edmfx_adv_test_box"),
)

integrator = CA.get_integrator(config); # make sure no errors

JET.@test_opt CA.get_integrator(config)

