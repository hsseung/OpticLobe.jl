"""
    OpticLobe is for exploring the neurons and connections of the Drosophila optic lobe.
    It is based on data files that are downloadable from the FlyWire Codex.

    See https://codex.flywire.ai for more information

    Basic usage:
    ```julia
    using OpticLobe
    Wtt["Tm1", "Dm3p"] # number of synapses from Tm1 cells to Dm3p cells
    ```
"""
# BEWARE: information is exported as global variables

module OpticLobe

const PKG_ROOT = pkgdir(OpticLobe)
const DATADIR = joinpath(PKG_ROOT, "data")

ENV["DATADEPS_ALWAYS_ACCEPT"] = true

#include("githubdependencies.jl")
#include("neurontable.jl")

include("codexdependencies.jl")

include("cellids.jl")
export ind2id, id2ind

include("celltypes.jl")
export intrinsictypes, boundarytypes, centraltypes, ind2type, visualtypes, alltypes

include("weightmatrix.jl")
export W, A

include("inoutaverages.jl")
export Wtt, Wct, Wtc, infraction, outfraction, inmean, outmean

#println("reading corrections")
#include("corrections.jl")

W = NamedArray(W, names = (ind2id, ind2id), dimnames = ("cellid", "cellid"))
A = NamedArray(A, names = (ind2id, alltypes), dimnames = ("cellid", "celltype"))
Ai = A[:, intrinsictypes]
export Ai

println("build class-family hierarchy above types")
include("hierarchy.jl")
export class2families, family2types

println("general utilities")
include("utils.jl")
export toppre, toppost   # functions for printing top input and output types
export showall, strings2ticks, type2ids, convert2arrows

println("codex and neuroglancer utilities")
include("flywire.jl")
export codex_open, ng_open, ng_hyper
end
