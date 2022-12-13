#=
A simple script for adding a specific dependency
to all specified environments.

Example:
```
julia .dev/add_pkg.jl ClimaCore#MyBranch
```
=#

@show ARGS
if isempty(ARGS)
    error("add_pkg.jl requires a package to be added, e.g., `ClimaCore#MyBranch`")
end
br = ARGS[1]
pkgname=first(split(br, "#"))
rev=last(split(br, "#"))

root = dirname(@__DIR__)
dirs = (
    root,
    joinpath(root, "test"),
    # joinpath(root, ".dev"),
    joinpath(root, "perf"),
    # joinpath(root, "docs"),
    joinpath(root, "examples"),
    joinpath(root, "ode_compat_examples"),
)
@info "Adding $pkgname, version $rev to environments $dirs"

cd(root) do
    for dir in dirs
        reldir = relpath(dir, root)
        @info "Adding package version $br to environment `$reldir`"
        cmd = if dir == root
            `$(Base.julia_cmd()) --project -e """import Pkg; Pkg.add(name=\"$pkgname\",rev=\"$rev\")"""`
        elseif dir == joinpath(root, ".dev")
            `$(Base.julia_cmd()) --project=$reldir -e """import Pkg; Pkg.add(name=\"$pkgname\",rev=\"$rev\")"""`
        else
            `$(Base.julia_cmd()) --project=$reldir -e """import Pkg; Pkg.develop(;path=\".\"); Pkg.add(name=\"$pkgname\",rev=\"$rev\")"""`
        end
        run(cmd)
    end
end

# https://github.com/JuliaLang/Pkg.jl/issues/3014
for dir in dirs
    cd(dir) do
        rm("LocalPreferences.toml"; force = true)
    end
end
