# %% [markdown]
# ## neurotransmitters for cells and types

# The `neurons` table should have been previously read by `cellids.jl`

# %%
using Missings, MissingsAsFalse
using NamedArrays
using StatsBase

# %%
# For standalone execution, include dependencies
if (@__MODULE__) == Main
    include("celltypes.jl")
end

# %%
function modethreshold(a::AbstractArray, threshold = 0.3)
    winner = mode(a)
    if ismissing(winner)
        return missing
    else
        @mfalse winningcount = sum(a .== winner)
        return winningcount/length(a) > threshold ? winner : missing
    end
end

# %%
"""
    ind2nt::Vector{Union{Missing, String}}

Maps cell indices to predicted neurotransmitter type.

Neurotransmitter predictions for individual cells:
- Examples: "gaba", "acetylcholine", "glutamate", "dopamine"
- Based on computational predictions from FlyWire analysis
- `missing`: Cells without neurotransmitter predictions

# Examples
```julia
ind2nt[1]                          # Neurotransmitter of first cell
ind2nt[id2ind[720575940599333574]] # Neurotransmitter of specific cell
findall(ind2nt .== "gaba")         # All GABAergic cell indices
```

# Notes
- Based on computational predictions, not direct measurements
- Individual cells may have uncertain or missing predictions
- Use `type2nt` for more reliable type-level neurotransmitter assignments
"""
ind2nt = passmissing(String).(neurons.nt_type)

# %% [markdown]

# %%
"""
    type2nt::NamedArray{Union{Missing, String}, 1}

Maps cell type names to consensus neurotransmitter type.

Neurotransmitter assignment for each cell type based on majority voting:
- Uses `modethreshold()` to require >30% consensus among cells of each type
- More reliable than individual cell predictions (`ind2nt`)
- `missing`: Types without sufficient neurotransmitter consensus

# Examples
```julia
type2nt[Name("Tm1")]              # Neurotransmitter of Tm1 cells
type2nt[Name("Dm3v")]             # Neurotransmitter of Dm3v cells
findall(type2nt .== "gaba")       # All GABAergic cell type names
```

# Notes
- Requires >30% of cells in a type to agree on neurotransmitter
- Based on consensus from `ind2nt` predictions
- Indexed by cell type names from `alltypes`
"""
@mfalse type2nt = [modethreshold(ind2nt[ind2type .== celltype]) for celltype in alltypes]

# %%
type2nt = NamedArray(type2nt, alltypes)

