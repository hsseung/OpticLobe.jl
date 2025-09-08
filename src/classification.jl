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
# # Classification/Hierarchical Annotations
#
# if you want to run as a standalone script, first run these scripts:
#   include("codexdependencies.jl")
#   include("cellids.jl")
#
# and define this function
#   (id2ind::Dict)(k) = id2ind[k]
#
# read Codex downloads
# - `classification.csv.gz` for all hierarchical annotations

# %%
using DataDeps, CSV, DataFrames

# %%
# For standalone execution, include dependencies
if (@__MODULE__) == Main
    include("cellids.jl")
end

# %%
classification = sort(CSV.read(datadep"Codex classification/classification.csv.gz", DataFrame))

# Verify that rows of `classification` have the canonical order, the same order as `neurons`
@assert classification.root_id == ind2id

# %%
"""
    ind2superclass::Vector{Union{Missing, String}}

Maps cell indices to superclass hierarchical classification.

Highest level of the three-tier hierarchical classification system:
- Examples: "visual", "sensory", "motor", "central"  
- Provides broad functional/anatomical groupings
- `missing`: Cells without hierarchical classification

# Examples
```julia
ind2superclass[1]                     # Superclass of first cell
ind2superclass[id2ind[720575940599333574]]  # Superclass of specific cell
findall(ind2superclass .== "visual") # All visual superclass cell indices
```

# Notes
- Part of hierarchical classification: superclass → class → subclass
- Based on FlyWire Codex hierarchical annotations
- Used for high-level functional analysis and cell grouping
"""
ind2superclass = Vector{Union{Missing, String}}(missing, length(ind2id))
ind2superclass[id2ind.(classification.root_id)] .= classification.super_class

# %%
"""
    ind2class::Vector{Union{Missing, String}}

Maps cell indices to class hierarchical classification.

Mid-level of the three-tier hierarchical classification system:
- Examples: "optic lobe", "central complex", "mushroom body"
- More specific than superclass, broader than subclass
- `missing`: Cells without hierarchical classification

# Examples
```julia
ind2class[1]                        # Class of first cell
ind2class[id2ind[720575940599333574]]  # Class of specific cell
findall(ind2class .== "optic lobe") # All optic lobe class cell indices
```

# Notes
- Part of hierarchical classification: superclass → class → subclass
- Based on FlyWire Codex hierarchical annotations
- Provides intermediate-level functional/anatomical groupings
"""
ind2class = Vector{Union{Missing, String}}(missing, length(ind2id))
ind2class[id2ind.(classification.root_id)] .= classification.class

# %%
"""
    ind2subclass::Vector{Union{Missing, String}}

Maps cell indices to subclass hierarchical classification.

Most specific level of the three-tier hierarchical classification system:
- Examples: "medulla intrinsic", "lamina monopolar", "projection neurons"
- Most detailed functional/anatomical categorization
- `missing`: Cells without hierarchical classification

# Examples
```julia
ind2subclass[1]                           # Subclass of first cell
ind2subclass[id2ind[720575940599333574]]  # Subclass of specific cell
findall(ind2subclass .== "medulla intrinsic")  # All medulla intrinsic cell indices
```

# Notes
- Part of hierarchical classification: superclass → class → subclass
- Based on FlyWire Codex hierarchical annotations
- Provides finest-grained functional/anatomical categories
"""
ind2subclass = Vector{Union{Missing, String}}(missing, length(ind2id))
ind2subclass[id2ind.(classification.root_id)] .= classification.sub_class
