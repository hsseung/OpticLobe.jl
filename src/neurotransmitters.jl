# %% [markdown]
# ## neurotransmitters for cells and types

# The `neurons` table should have been previously read by `cellids.jl`

# %%
using Missings, MissingsAsFalse
using NamedArrays
using StatsBase

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
ind2nt = passmissing(String).(neurons.nt_type)

# %% [markdown]

# %%
@mfalse type2nt = [modethreshold(ind2nt[ind2type .== celltype]) for celltype in alltypes]

# %%
type2nt = NamedArray(type2nt, alltypes)

