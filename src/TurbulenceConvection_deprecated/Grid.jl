import ClimaCore.Spaces as Spaces
import ClimaCore.Fields as Fields
import ClimaCore.Domains as Domains
import ClimaCore.Meshes as Meshes
import ClimaCore.Utilities as Utilities
struct Grid{NZ, CS, FS, SC, SF}
    cs::CS
    fs::FS
    zc::SC
    zf::SF
    function Grid(space::Spaces.CenterFiniteDifferenceSpace)
        nz = length(space)
        cs = space
        fs = Spaces.FaceFiniteDifferenceSpace(cs)
        zc = Fields.coordinate_field(cs)
        zf = Fields.coordinate_field(fs)
        CS = typeof(cs)
        FS = typeof(fs)
        SC = typeof(zc)
        SF = typeof(zf)
        return new{nz, CS, FS, SC, SF}(cs, fs, zc, zf)
    end
end

Grid(mesh::Meshes.IntervalMesh) = Grid(Spaces.CenterFiniteDifferenceSpace(mesh))

function Grid(Δz::FT, nz::Int) where {FT <: AbstractFloat}
    z₀, z₁ = FT(0), FT(nz * Δz)

    domain = Domains.IntervalDomain(
        CC.Geometry.ZPoint{FT}(z₀),
        CC.Geometry.ZPoint{FT}(z₁),
        boundary_tags = (:bottom, :top),
    )

    mesh = Meshes.IntervalMesh(domain, nelems = nz)
    return Grid(mesh)
end

fc_index(
    i,
    ::Union{
        Spaces.FaceExtrudedFiniteDifferenceSpace,
        Spaces.FaceFiniteDifferenceSpace,
    },
) = Utilities.PlusHalf(i)

fc_index(
    i,
    ::Union{
        Spaces.CenterExtrudedFiniteDifferenceSpace,
        Spaces.CenterFiniteDifferenceSpace,
    },
) = i

# Index of the first interior cell above the surface
surf(f::Fields.Field) = Spaces.level(f, fc_index(1, axes(f)))[]
toa(f::Fields.Field) =
    Spaces.level(f, fc_index(Spaces.nlevels(axes(f)), axes(f)))[]
