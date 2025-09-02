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
classification = sort(CSV.read(datadep"Codex classification/classification.csv.gz", DataFrame))

# Verify that rows of `classification` have the canonical order, the same order as `neurons`
@assert classification.root_id == ind2id

# %%
ind2superclass = Vector{Union{Missing, String}}(missing, length(ind2id))
ind2superclass[id2ind.(classification.root_id)] .= classification.super_class

# %%
ind2class = Vector{Union{Missing, String}}(missing, length(ind2id))
ind2class[id2ind.(classification.root_id)] .= classification.class

# %%
ind2subclass = Vector{Union{Missing, String}}(missing, length(ind2id))
ind2subclass[id2ind.(classification.root_id)] .= classification.sub_class
