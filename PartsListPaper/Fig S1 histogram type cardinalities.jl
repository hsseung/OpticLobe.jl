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
# # Figure S1. histogram of type cardinalities, nonuniform bucketing

# Needs updating to not depend on CSV file for central brain cardinalities

# %%
using OpticLobe

# %% jupyter={"outputs_hidden": false}
using NamedArrays
using StatsBase
using Plots
using Measures

# %% [markdown]
# ## families types and cells per class

# %%
cellnumbers = NamedArray(sum(Ai.array, dims=1)[:], intrinsictypes)

# %% [markdown]
# ## central brain number of cells per type

# %%
using CSV, DataFrames

# %%
cb_df = CSV.read("../v630/central_brain_type_cardinalities.csv", DataFrame)

# %%
maximum(cb_df[:, "Num Cells"])

# %%
binedges = [1, 2, 3, 11, 26, 51, 101, 201, 700, 2500]

# %%
h_cb = fit(Histogram, cb_df[:, "Num Cells"], binedges, closed=:left)

# %%
h_optic = fit(Histogram, cellnumbers, binedges, closed=:left)

# %%
cmap = palette(:default)
default(
    fontfamily="Helvetica",
    label = "",
    titlefontsize = 16,
    tickfontsize = 14,
    guidefontsize = 18,
    xrotation = 45,
    bottom_margin = 10mm,
    left_margin = 10mm,
    top_margin = 5mm,
    xlabel = "number of cells",
    ylabel = "number of types",
    xticks = (1:9, [1, 2, "3-10", "11-25", "26-50", "51-100", "101-200", "201-700", "701-"]),
)

# %%
plot(
    bar(
        h_optic.weights,
        title = "right optic lobe"
    ),
    bar(
        h_cb.weights,
        title = "central brain"
    ),
    layout = (1, 2), size = (1200, 450)
)

# %%
#savefig("HistogramNumberOfCellsPerTypeBuckets.svg")
savefig("HistogramTypeCardinalitiesBuckets.svg")
