# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:percent
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: Julia 1.10.3
#     language: julia
#     name: julia-1.10
# ---

# %% [markdown]
# # Input and Output Averages
# These are all Named SparseArrays
#
# Unnormalized synapse counts
# - `Wtt[pre, post]` - number of synapses from `pre` type to `post` type
#
# Normalized synapse counts
# - `infraction[pre, post]` - fraction of input synapses to `post` type coming from `pre` type
# - `outfraction[pre, post]` - fraction of output synapses from `pre` type going to `post` type
#
# mean and standard deviation over cells in a type
# - `inmean`, `instd` - number of input synapses to `post` cell coming from `pre` type
# - `outmean`, `outstd` - number of output synapses from `pre` cell going to `post` type
#
# quartiles over cells in a type
# - `in25`, `inmedian`, `in75`
# - `out25`, `outmedian` `out75`

# %%
using SparseArrays, NamedArrays
using StatsBase

# %%
using MissingsAsFalse
using ProgressMeter

# %%
# For standalone execution, include dependencies and set defaults
if (@__MODULE__) == Main
    include("weightmatrix.jl")
end

# %%
ntypes = length(alltypes)

# %%
"""
    Wtt::NamedArray{Int32, 2}

Type-to-type connectivity matrix showing synapse counts between cell types.

Matrix where `Wtt[pretype, posttype]` gives the total number of synapses 
from all cells of `pretype` to all cells of `posttype`.

# Examples
```julia
Wtt["Tm1", "Dm3v"]        # Total synapses from Tm1 cells to Dm3v cells
Wtt["Tm1", :]             # All outputs from Tm1 cells (by target type)
Wtt[:, "Dm3v"]            # All inputs to Dm3v cells (by source type)
```

# Notes
- Computed as `A'*W*A` where `A` is the assignment matrix and `W` is cell-to-cell connectivity
- Indices are cell type names from `alltypes`
- Foundation for computing normalized connectivity measures
"""
Wtt = A'*W*A           # type to type

"""
    Wtc::NamedArray{Int32, 2}

Type-to-cell connectivity matrix showing synapse counts from cell types to individual cells.

Matrix where `Wtc[pretype, cellid]` gives the total number of synapses 
from all cells of `pretype` to the specific cell `cellid`.

# Examples
```julia
Wtc["Tm1", Name(720575940599333574)]  # Synapses from Tm1 cells to specific cell
Wtc["Tm1", :]                         # All Tm1 outputs (by target cell)
Wtc[:, Name(720575940599333574)]      # All inputs to specific cell (by source type)
```

# Notes
- Computed as `A'*W` where `A` is the assignment matrix and `W` is cell-to-cell connectivity
- Row indices: cell type names from `alltypes`
- Column indices: cell IDs from `ind2id`
"""
Wtc = A'*W             # type to cell

"""
    Wct::NamedArray{Int32, 2}

Cell-to-type connectivity matrix showing synapse counts from individual cells to cell types.

Matrix where `Wct[cellid, posttype]` gives the total number of synapses 
from the specific cell `cellid` to all cells of `posttype`.

# Examples
```julia
Wct[Name(720575940599333574), "Dm3v"]  # Synapses from specific cell to Dm3v cells
Wct[Name(720575940599333574), :]       # All outputs from specific cell (by target type)
Wct[:, "Dm3v"]                         # All inputs to Dm3v cells (by source cell)
```

# Notes
- Computed as `W*A` where `W` is cell-to-cell connectivity and `A` is the assignment matrix
- Row indices: cell IDs from `ind2id`  
- Column indices: cell type names from `alltypes`
"""
Wct = W*A              # cell to type

normalization = 1 ./ sum(Wct, dims=1)  # this is contorted, but seems necessary to preserve sparsity
normalization[isinf.(normalization)] .= 0 # needed because AN_GNG_124 has zero inputs

"""
    infraction::NamedArray{Float32, 2}

Input fraction matrix showing the fraction of input synapses to each postsynaptic type from each presynaptic type.

Matrix where `infraction[pretype, posttype]` gives the fraction of all input synapses 
to `posttype` that come from `pretype`. Values sum to 1.0 across each column (postsynaptic type).

# Examples
```julia
infraction["Tm1", "Dm3v"]     # Fraction of Dm3v inputs from Tm1 cells
infraction[:, "Dm3v"]         # All input fractions to Dm3v (sums to ~1.0)
sum(infraction[:, "Dm3v"])    # Should be ~1.0 (allowing for rounding)
```

# Notes
- Computed by normalizing `Wtt` columns: `Wtt ./ sum(Wtt, dims=1)`
- Values range from 0.0 to 1.0
- Each column sums to 1.0 (fraction of total inputs)
- Types with zero inputs have fraction 0.0 for all presynaptic partners
"""
infraction = Wtt .* normalization

normalization = 1 ./ sum(Wtc, dims=2)
normalization[isinf.(normalization)] .= 0    # setting 0/0 = 0  This takes care of T1, and also DNge061 and DNge074

"""
    outfraction::NamedArray{Float32, 2}

Output fraction matrix showing the fraction of output synapses from each presynaptic type to each postsynaptic type.

Matrix where `outfraction[pretype, posttype]` gives the fraction of all output synapses 
from `pretype` that go to `posttype`. Values sum to 1.0 across each row (presynaptic type).

# Examples
```julia
outfraction["Tm1", "Dm3v"]    # Fraction of Tm1 outputs going to Dm3v cells
outfraction["Tm1", :]         # All output fractions from Tm1 (sums to ~1.0)
sum(outfraction["Tm1", :])    # Should be ~1.0 (allowing for rounding)
```

# Notes
- Computed by normalizing `Wtt` rows: `Wtt ./ sum(Wtt, dims=2)`
- Values range from 0.0 to 1.0
- Each row sums to 1.0 (fraction of total outputs)
- Types with zero outputs (like T1) have fraction 0.0 for all postsynaptic partners
"""
outfraction = Wtt .* normalization

normalization = 1 ./ sum(A, dims = 1)

"""
    inmean::NamedArray{Float32, 2}

Mean input synapse count per postsynaptic cell for each type-to-type connection.

Matrix where `inmean[pretype, posttype]` gives the mean number of input synapses 
per `posttype` cell that come from `pretype` cells.

# Examples
```julia
inmean["Tm1", "Dm3v"]         # Mean synapses per Dm3v cell from Tm1 cells
inmean[:, "Dm3v"]             # Mean inputs to each Dm3v cell (by presynaptic type)
sum(inmean[:, "Dm3v"])        # Total mean inputs per Dm3v cell
```

# Notes
- Computed as `Wtt / number_of_postsynaptic_cells`
- Represents average connectivity strength per postsynaptic cell
- Useful for comparing connection strengths independent of cell type sizes
- Units: synapses per cell
"""
inmean = Wtt .* normalization

normalization = 1 ./ sum(A, dims = 1)'

"""
    outmean::NamedArray{Float32, 2}

Mean output synapse count per presynaptic cell for each type-to-type connection.

Matrix where `outmean[pretype, posttype]` gives the mean number of output synapses 
per `pretype` cell that go to `posttype` cells.

# Examples
```julia
outmean["Tm1", "Dm3v"]        # Mean synapses per Tm1 cell to Dm3v cells
outmean["Tm1", :]             # Mean outputs from each Tm1 cell (by postsynaptic type)
sum(outmean["Tm1", :])        # Total mean outputs per Tm1 cell
```

# Notes
- Computed as `Wtt / number_of_presynaptic_cells`
- Represents average connectivity strength per presynaptic cell
- Useful for comparing connection strengths independent of cell type sizes
- Units: synapses per cell
"""
outmean = Wtt .* normalization

# %% [markdown]
# ## check for NaNs
# %%
@assert sum(isnan.(Wtt)) == 0
@assert sum(isnan.(infraction)) == 0
@assert sum(isnan.(outfraction)) == 0
@assert sum(isnan.(inmean)) == 0
@assert sum(isnan.(outmean)) == 0

# %% [markdown]
# ## function for converting to Named sparse array with `names = (alltypes, alltypes)``

# %%
namedsparse(x) = NamedArray(SparseMatrixCSC{Float32, Int32}(x), names = (alltypes, alltypes), dimnames = ("celltype", "celltype"))

# %% [markdown]
# ## convert to Named sparse arrays
println("convert to Named sparse arrays")

# %%
inmean, outmean = namedsparse.([inmean, outmean])
infraction, outfraction = namedsparse.([infraction, outfraction])

Wtt = NamedArray(Wtt, names = (alltypes, alltypes), dimnames = ("celltype", "celltype"))      # type to type
Wct = NamedArray(Wct, names = (ind2id, alltypes), dimnames = ("cellid", "celltype"))      # type to type
Wtc = NamedArray(Wtc, names = (alltypes, ind2id), dimnames = ("celltype", "cellid"))      # type to type

# compute importance after this correction
"""
    importance::NamedArray{Float32, 2}

Connection importance from either presynaptic or postsynaptic perspective.

Matrix where `importance[pretype, posttype]` gives the maximum of:
- `infraction`: fraction of postsynaptic type's inputs (important to postsynaptic type)
- `outfraction`: fraction of presynaptic type's outputs (important to presynaptic type)

This identifies connections that are significant to at least one of the partner types.

# Examples
```julia
importance["Tm1", "Dm3v"]     # Importance of Tm1→Dm3v connection
importance["Tm1", :]          # Importance of all Tm1 output connections
importance[:, "Dm3v"]         # Importance of all Dm3v input connections
```

# Notes
- Computed as `max.(infraction, outfraction)`
- Values range from 0.0 to 1.0
- High values indicate the connection is a major pathway for at least one partner
- A connection can be important even if asymmetric (e.g., major output for pre, minor input for post)
"""
importance = max.(infraction, outfraction)

# %%
"""
    inrank::NamedArray{Int64, 2}

Rank ordering of presynaptic inputs for each postsynaptic intrinsic type.

Matrix where `inrank[pretype, posttype]` gives the rank of `pretype` among all 
presynaptic partners of `posttype`, based on mean synapse count per postsynaptic cell.
Rank 1 indicates the strongest input partner.

# Examples
```julia
inrank["Tm1", "Dm3v"]         # Rank of Tm1 among Dm3v's input partners
inrank[:, "Dm3v"]             # All input ranks for Dm3v (1 = strongest)
findall(inrank[:, "Dm3v"] .== 1)  # Strongest input partner(s) to Dm3v
```

# Notes
- Only includes intrinsic types (not boundary/central/visual types)
- Based on `inmean` values (mean synapses per postsynaptic cell)
- Lower rank numbers indicate stronger connections
- Computed using `ordinalrank` with `rev=true` (highest mean gets rank 1)
"""
inrank = stack(ordinalrank.(eachcol(collect(inmean[intrinsictypes, intrinsictypes])), rev=true), dims=2)
inrank = NamedArray(inrank, (intrinsictypes, intrinsictypes))

"""
    outrank::NamedArray{Int64, 2}

Rank ordering of postsynaptic outputs for each presynaptic intrinsic type.

Matrix where `outrank[pretype, posttype]` gives the rank of `posttype` among all 
postsynaptic partners of `pretype`, based on mean synapse count per presynaptic cell.
Rank 1 indicates the strongest output partner.

# Examples
```julia
outrank["Tm1", "Dm3v"]        # Rank of Dm3v among Tm1's output partners
outrank["Tm1", :]             # All output ranks for Tm1 (1 = strongest)
findall(outrank["Tm1", :] .== 1)  # Strongest output partner(s) of Tm1
```

# Notes
- Only includes intrinsic types (not boundary/central/visual types)
- Based on `outmean` values (mean synapses per presynaptic cell)
- Lower rank numbers indicate stronger connections
- Computed using `ordinalrank` with `rev=true` (highest mean gets rank 1)
"""
outrank = stack(ordinalrank.(eachrow(collect(outmean[intrinsictypes, intrinsictypes])), rev=true), dims=1)
outrank = NamedArray(outrank, (intrinsictypes, intrinsictypes))
