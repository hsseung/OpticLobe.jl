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
#     display_name: Julia 1.10.4
#     language: julia
#     name: julia-1.10
# ---

# %% [markdown]
# # Hexagonal grid coordinates for *Drosophila* eye and optic lobe
#
# read CSV file to create $37\times 34$ array `pq2column` and $796\times 2$ array `column2pq`:
#
#   - `pq2column[p, q]` yields
#     - the column ID of the $(p, q)$ location on the hexagonal grid
#     - `missing` if there is no column at that location
#   - `column2pq[i, :]` returns `[p, q]` for column ID = `i`
#
# Column IDs run from 1 to 796. The $p-q$ coordinate system is defined below.
#
# Note: To use the coordinate system, we also need `columns_Mi1.csv`, which maps column IDs to cell IDs (64 bit integers), and is handled elsewhere. Column ID $i$ is occupied by the $i$th cell ID in the list. 

# %% [markdown]
# ## p-q coordinates for columns
#
# For the right eye, we have:
# ```
#          dorsal
#            ^
#            |
# anterior <-+-> posterior
#            |
#            v
#         ventral
# ```
#
# Axes for the hexagonal grid are defined in Fig. 2 of Arthur Zhao et al. 2022.
#
# $+p$ is anterodorsal, and $+q$ is posterodorsal.
# $+v$ is dorsal, and $+h$ is posterior.
#    
# ```
# +p   +v   +q
#  \   |   /
#   \  |  /
#    \ | /
#     \|/
#      +----- +h
# ```
#
# The $+p$ and $+q$ directions would be $120^\circ$ apart on a perfect hexagonal grid. Zhao et al. approximate them as $90^\circ$ apart.

# %%
using CSV, DataFrames

# %% [markdown]
# ## spreadsheet
# read CSV exported from
# https://docs.google.com/spreadsheets/d/1ol_8nJ1K1bwJ_gcY_esZ-S0Dqk-NKYVHp1kJ3S7cFl8/edit?usp=sharing
#
# This spreadsheet contains column IDs entered into a grid.\
# The grid is rotated clockwise so that the top right is the top of the eye.\
# Columns are aligned along the $p$ axis, and rows are aligned along the $q$ axis.
#
# ```
# p
# ^
# |
# |
# |
# |-----> q
# ```
#
# The assignments to $p-q$ coordinates were performed manually.
#
# In the spreadsheet, column 628 is the origin at (19, 17). This is just a rough guess though, and doesn't matter for analyses.

# %%
RightEyeGrid = CSV.read(joinpath(DATADIR, "RightEyeGrid.csv"), DataFrame)

# %% [markdown]
# We need to reverse the order of the rows if we want `RightEyeGrid[p, q]` to contain the $(p, q)$ column.
# That puts the top of the eye at the bottom right.

# %%
reverse!(RightEyeGrid)

# %%
pq2column = Array(RightEyeGrid[:, 2:end])

# %%
nright = maximum(skipmissing(pq2column))

# %%
right = 1:nright

# %%
# check that all right eye columns are represented
issetequal(right, collect(skipmissing(pq2column)))

# %% [markdown]
# ## define inverse map that goes from column ID to (p, q) coordinates

# %%
pmax, qmax = size(pq2column)

# %%
column2pq = zeros(Int32, nright, 2)
for p = 1:pmax
    for q = 1:qmax
        id = pq2column[p, q]
        if ~ismissing(id)
            column2pq[id, :] = [p q]
        end
    end
end

# %% [markdown]
# ### how to turn inverse map into a DataFrame

# %%
DataFrame(column2pq, [:p, :q])

# %%
#using JLD2
#@save "ColumnCoordinates.jld2" pq2column column2pq
