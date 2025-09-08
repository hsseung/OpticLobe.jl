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
# Reads visual neuron columns from codex download and creates:
#
# - `id2pq` - dictionary mapping cell IDs to hexagonal (p,q) coordinates  
# - `pq2column` - matrix mapping (p,q) coordinates to column IDs for right hemisphere
#
# This module supersedes columncoordinates.jl, columncell.jl, and cellcoordinates.jl
#

# %%
using CSV, DataFrames

# %%
# For standalone execution, include dependencies
if (@__MODULE__) == Main
    println("Running in standalone mode")

    include("codexdependencies.jl")
end

# %%
visualcolumns = CSV.read(datadep"Codex visual neuron columns/column_assignment.csv.gz", DataFrame)

# id2pq is defined for both eyes
# The offset is backward compatible with right eye, and makes 1 the minimum value of p and q for the right eye
# But 0 is the minimum value for the left eye. This should be fixed later.
pq = [visualcolumns.p .+ 19 visualcolumns.q .+ 17]

"""
    id2pq::Dict{Int, Vector{Int}}

Maps FlyWire cell IDs to hexagonal (p,q) column coordinates in the optic lobe.

Dictionary mapping visual neuron cell IDs to their spatial coordinates in the hexagonal
column lattice of the optic lobe. Each coordinate is a 2-element vector [p, q] representing
the position in the hexagonal grid system.

# Examples
```julia
id2pq[720575940599333574]          # Get (p,q) coordinates for specific cell
coords = id2pq[720575940599333574]  # Returns [p, q] as Vector{Int}
p, q = id2pq[720575940599333574]    # Destructure coordinates

# Find all cells at a specific location
cells_at_pq = [id for (id, coord) in id2pq if coord == [10, 15]]
```

# Notes
- Only available for visual neurons that have been assigned to columns
- Coordinates use offset system: right eye p,q ≥ 1, left eye p,q ≥ 0
- Based on FlyWire Codex column assignment data
- Used with `pq2column` for spatial analysis and visualization
- Coordinates correspond to the hexagonal lattice structure of optic lobe columns
"""
id2pq = Dict(visualcolumns.root_id .=> map(collect, eachrow(pq)))

# Note that pq2column is only defined at present for the right eye
# Create OffsetArray for right hemisphere columns indexed by p and q
using OffsetArrays
right_rows = visualcolumns[visualcolumns.hemisphere .== "right", :]
p_range = minimum(right_rows.p):maximum(right_rows.p)
q_range = minimum(right_rows.q):maximum(right_rows.q)
pq2column = OffsetArray(Matrix{Union{Missing, Int}}(missing, length(p_range), length(q_range)), p_range, q_range)

for row in eachrow(right_rows)
    pq2column[row.p, row.q] = row.column_id
end

"""
    pq2column::Matrix{Union{Missing, Int}}

Maps hexagonal (p,q) coordinates to column IDs for the right optic lobe hemisphere.

Matrix where `pq2column[p, q]` returns the column ID at hexagonal coordinates (p, q),
or `missing` if no column exists at that location. Currently only defined for the
right hemisphere of the optic lobe.

# Examples
```julia
pq2column[10, 15]                  # Get column ID at position (10, 15)
column_id = pq2column[p, q]        # Returns Int or missing

# Find all valid column positions
valid_positions = findall(!ismissing, pq2column)
column_ids = pq2column[valid_positions]
```

# Notes
- Only covers right optic lobe hemisphere (left hemisphere not implemented)
- Uses raw (p, q) coordinates from the original data (not the offset coordinates in `id2pq`)
- Returns `missing` for positions without columns (gaps in hexagonal lattice)
- Based on FlyWire Codex column assignment data
- Used with `id2pq` for spatial mapping and visualization
- Matrix indices correspond to the natural hexagonal coordinate system
"""
pq2column = parent(pq2column)
