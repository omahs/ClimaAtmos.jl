const ca_dir = joinpath(@__DIR__, "..", "..")
include(joinpath(ca_dir, "src", "solver", "yaml_helper.jl"))
using PrettyTables

function print_repl_script(config)
    ib = """"""
    ib *= """\n"""
    ib *= """using Revise; import ClimaAtmos as CA;\n"""
    ib *= """\n"""
    ib *= """config_dict = Dict();\n"""
    for (flag, val) in config
        if val isa AbstractString
            ib *= "config_dict[\"$flag\"] = \"$val\";\n"
        else
            ib *= "config_dict[\"$flag\"] = $val;\n"
        end
    end
    ib *= """\n"""
    ib *= """config = CA.AtmosConfig(; config_dict);\n"""
    ib *= """\n"""
    ib *= """include("examples/hybrid/driver.jl")\n"""
    println(ib)
end

configs = configs_per_job_id()
@assert length(configs) > 0

for (job_id, config) in configs
    println("### Buildkite job `$job_id`")
    print_repl_script(config)
    println()
end
