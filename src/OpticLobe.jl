# BEWARE: information is exported as global variables

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
export Wtt, Wct, Wtc, infraction, outfraction, inmean, outmean, importance, inrank, outrank

"""
    W::NamedArray{Int32, 2}

Main synaptic connectivity matrix where `W[i, j]` is the number of synapses from neuron `i` to neuron `j`.

This is the core data structure containing synaptic connections between all neurons in the dataset.
The matrix is sparse for memory efficiency and indexed by FlyWire root IDs using NamedArrays.

# Indexing
- By position: `W[1, 2]` - synapses from first cell to second cell  
- By cell ID: `W[Name(720575940599333574), Name(720575940620875399)]` - synapses between specific cells
- Row/column slicing: `W[cellid, :]` - all outputs from a cell, `W[:, cellid]` - all inputs to a cell

# Examples
```julia
# Check connectivity between specific cells
W[Name(720575940599333574), Name(720575940620875399)]

# Get all postsynaptic partners of a cell (with synapse counts)
outputs = W[Name(720575940599333574), :]
top_targets = sort(outputs, rev=true)[1:10]

# Get all presynaptic partners of a cell  
inputs = W[:, Name(720575940620875399)]
top_sources = sort(inputs, rev=true)[1:10]

# Total output/input synapse counts for a cell
total_out = sum(W[Name(720575940599333574), :])
total_in = sum(W[:, Name(720575940620875399)])
```

# Notes
- Sparse matrix (~130K × 130K cells, ~77M synapses for Princeton version)
- Autapses (self-connections) are set to zero  
- T1 cell outputs are zeroed based on prior biological knowledge
- Synapse version depends on package configuration (Princeton or Buhmann)
- Use `W_Princeton` and `W_Buhmann` if both versions are loaded
"""
W = NamedArray(W, names = (ind2id, ind2id), dimnames = ("cellid", "cellid"))
if @load_preference("load_both_synapses", false)
    global W_Buhmann = NamedArray(W_Buhmann, names = (ind2id, ind2id), dimnames = ("cellid", "cellid"))
    global W_Princeton = NamedArray(W_Princeton, names = (ind2id, ind2id), dimnames = ("cellid", "cellid"))
end
A = NamedArray(A, names = (ind2id, alltypes), dimnames = ("cellid", "celltype"))

"""
    Ai::NamedArray{Bool, 2}

Intrinsic cell assignment matrix - subset of `A` containing only intrinsic cell types.

Boolean matrix where `Ai[cellid, intrinsic_type]` is `true` if the cell belongs to that intrinsic type.
This is a filtered version of the full assignment matrix `A` containing only the 230 intrinsic types.

# Examples
```julia
Ai[Name(720575940599333574), Name("Tm1")]  # Check if specific cell is Tm1
Ai[:, Name("Tm1")]                         # All Tm1 cells (boolean vector)
findall(Ai[:, Name("Tm1")])                # Indices of all Tm1 cells
```

# Notes
- Subset of `A` restricted to `intrinsictypes`: `A[:, intrinsictypes]`
- Row indices: cell IDs from `ind2id`
- Column indices: intrinsic type names from `intrinsictypes` 
- Dimensions: ~130K cells × 230 intrinsic types
- Used for analyses focused on optic lobe intrinsic neurons
"""
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

println("pq coordinates of columns and cells")
#include("columncoordinates.jl")
#include("columncell.jl")
#export column2pq, pq2column, columns_df
#export pq2column

#include("cellcoordinates.jl")
#export id2pq

include("columnassignment.jl")
export id2pq, pq2column

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
