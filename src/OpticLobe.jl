"""
    OpticLobe is for exploring the neurons and connections of the Drosophila optic lobe.
    It is based on data files that are downloadable from the FlyWire Codex.

    See https://codex.flywire.ai for more information

    Basic usage:
    ```julia
    using OpticLobe
    Wtt["Tm1", "Dm3p"] # number of synapses from Tm1 cells to Dm3p cells
    ```

    ## Synapse Versions
    
    The package supports two synapse detection methods:
    - **Princeton** (default): Newer synapse predictions
    - **Buhmann**: Original synapse predictions from Buhmann et al.
    
    ### Single synapse version mode (default)
    By default, only Princeton synapses are loaded into `W`:
    ```julia
    using OpticLobe
    W[cellid1, cellid2]  # Uses Princeton synapses by default
    ```
    
    To switch to Buhmann synapses as default:
    ```julia
    set_default_synapses("Buhmann")
    # Restart Julia for changes to take effect
    ```
    
    ### Both synapse versions mode
    To load both synapse versions simultaneously:
    ```julia
    enable_both_synapses(true)
    # Restart Julia for changes to take effect
    ```
    
    After enabling, access both versions:
    ```julia
    W_Princeton[cellid1, cellid2]  # Princeton synapses
    W_Buhmann[cellid1, cellid2]    # Buhmann synapses
    W[cellid1, cellid2]             # Default version (Princeton unless changed)
    ```
"""
# BEWARE: information is exported as global variables

module OpticLobe

const PKG_ROOT = pkgdir(@__MODULE__)
const DATADIR = joinpath(PKG_ROOT, "data")

ENV["DATADEPS_ALWAYS_ACCEPT"] = true

include("codexdependencies.jl")

include("cellids.jl")
export ind2id, id2ind

include("cellstats.jl")
export cell_length, cell_area, cell_volume

include("celltypes.jl")
export ind2category, ind2side, intrinsictypes, boundarytypes, othertypes, ind2type, visualtypes, alltypes, A

include("classification.jl")
export ind2superclass, ind2class, ind2subclass

include("neurotransmitters.jl")
export ind2nt, type2nt

using Preferences

function set_default_synapses(new_synapses::String)
    if !(new_synapses in ("Buhmann", "Princeton"))
        throw(ArgumentError("Invalid synapses version: \"$(new_synapses)\""))
    end
    @set_preferences!("default_synapses" => new_synapses)
    @info("Default synapses set to $(new_synapses); restart your Julia session for this change to take effect!")
end
export set_default_synapses

function enable_both_synapses(enable::Bool=true)
    @set_preferences!("load_both_synapses" => enable)
    @info("Load both synapses set to $(enable); restart your Julia session for this change to take effect!")
end
export enable_both_synapses

include("weightmatrix.jl")
export W

# Export additional weight matrices if both were loaded
if @load_preference("load_both_synapses", false)
    export W_Buhmann, W_Princeton
end

include("inoutaverages.jl")
export Wtt, Wct, Wtc, infraction, outfraction, inmean, outmean

W = NamedArray(W, names = (ind2id, ind2id), dimnames = ("cellid", "cellid"))
if @load_preference("load_both_synapses", false)
    global W_Buhmann = NamedArray(W_Buhmann, names = (ind2id, ind2id), dimnames = ("cellid", "cellid"))
    global W_Princeton = NamedArray(W_Princeton, names = (ind2id, ind2id), dimnames = ("cellid", "cellid"))
end
A = NamedArray(A, names = (ind2id, alltypes), dimnames = ("cellid", "celltype"))
Ai = A[:, intrinsictypes]
export Ai
export Name

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

println("pq coordinates of columns")
include("columncoordinates.jl")
include("columncell.jl")
#export column2pq, pq2column, columns_df
export pq2column

include("cellcoordinates.jl")
export id2pq

include("hexgraphics.jl")
export rect2hex, square2hex, crop, montage
export eyehot
export triad, eyetriad, typetriad, celltriad
export ellipsesummary, drawellipse
export hexproject, drawpqaxes, hexannulus
export HexagonEye

include("maptrace.jl")
export tracebacktypes, tracetypes, inmaps, outmaps, scorepath
export preimage, prepreimage, preprepreimage # legacy code superseded by `inmaps`
#export postimage

include("spatial.jl")
export findcenter, convcluster
export seven

end
