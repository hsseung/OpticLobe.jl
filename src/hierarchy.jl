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
# # Class-Family-Type Hierarchy for intrinsic types

# %%
#using JLD2

# %%
#using OpticLobe

# %%
using OrderedCollections
using NaturalSort

# %% [markdown]
# ## five classes at top level

# %%
# each class consists of neuropil-defined families specified by ordered dict.
class2families = OrderedDict{String, Vector{String}}()

# %%
class2families["receptor"] = ["R1-6", "R7-8"]
class2families["columnar"] = ["L", "C", "Lawf", "T1", "Mi", "Tm", "TmY", "Y", "T2", "T3", "T4", "T5", "Tlp",]
class2families["interneuron"] = ["Lai", "Dm", "Pm", "Sm", "Li", "LPi"]
# cross-neuropil tangential
class2families["x-tangential"] = ["Lat", "MLt", "LMt", "LLPt", "PDt"]
# cross-neuropil amacrine
class2families["x-amacrine"] = ["LMa", "MLLPa"]

# %%
for (k, v) in class2families
    println(k, ": ", v)
end

# %% [markdown]
# ## create dictionary of type families
# filter `intrinsictypes` using prefix corresponding to each family

# %%
# each family consists of neuropil-defined types
family2types = OrderedDict{String, Vector{String}}()

# %% [markdown]
# ### receptor
family2types["R1-6"] = ["R1-6"]
family2types["R7-8"] = filter(startswith(r"R[78]"), intrinsictypes)  # includes also R7-DRA and R8-DRA

# %% [markdown]
# ### columnar

# %%
columnarprefixes = OrderedDict{String, Union{Regex, String}}(class2families["columnar"] .=> class2families["columnar"])
columnarprefixes["L"] = r"L\d+"
columnarprefixes["C"] = r"C\d+"
columnarprefixes["Tm"] = r"Tm[^Y]"
for (f, prefix) in columnarprefixes
    family2types[f] = filter(startswith(prefix), intrinsictypes)
end

# %% [markdown]
# ### interneuron (intra-neuropil non-columnar)

# %%
for prefix in class2families["interneuron"]
    family2types[prefix] = filter(startswith(prefix), intrinsictypes)
end

# %% [markdown]
# ### tangential (cross-neuropil)

# %%
for prefix in class2families["x-tangential"]
    family2types[prefix] = filter(startswith(prefix), intrinsictypes)
end

# %% [markdown]
# ### amacrine (cross-neuropil)

# %%
family2types["LMa"] = vcat(filter(startswith("LMa"), intrinsictypes), ["CT1"])
family2types["MLLPa"] = ["Am1"]

# %%
for (k, v) in family2types
    family2types[k] = sort(v, lt=natural)
end
## %%
for (k, v) in family2types
    println(k, ": ", v)
end

println(sum([length(v) for v in values(family2types)]), " intrinsic types + photoreceptors")

@assert issetequal(intrinsictypes, vcat(values(family2types)...))

# %%
# @save "hierarchyclassfamilytype.jld2" family2types class2families
