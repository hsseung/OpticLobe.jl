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
ntypes = length(alltypes)

# %%
Wtt = A'*W*A           # type to type
Wtc = A'*W             # type to cell
Wct = W*A              # cell to type

normalization = 1 ./ sum(Wct, dims=1)  # this is contorted, but seems necessary to preserve sparsity
infraction = Wtt .* normalization

normalization = 1 ./ sum(Wtc, dims=2)
normalization[isinf.(normalization)] .= 0    # setting 0/0 = 0  This takes care of T1, and also DNge061 and DNge074
outfraction = Wtt .* normalization

normalization = 1 ./ sum(A, dims = 1)
inmean = Wtt .* normalization

normalization = 1 ./ sum(A, dims = 1)'
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
importance = max.(infraction, outfraction)
