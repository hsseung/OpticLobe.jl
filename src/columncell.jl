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
# # column assignments

# %%
# using OpticLobe

# %%
using CSV, DataFrames

# %%
using NamedArrays

# %%
columns_df = CSV.read(joinpath(DATADIR, "columns_Mi1.csv"), DataFrame)

right = 1:nrow(columns_df)
ncolumn = length(right)

# %% [markdown]
# ## connectivity-based assignment

# %%
using Hungarian

# %%
function hungarianmax(W)
    hungarian(maximum(W) .- W)
end

# %%
"""
For default (`dims=1`), for each cell in `posttype`, find the strongest presynaptic partner in `pretype`
The strength must be at least `threshold`. Otherwise the answer is `missing`.
"""
function findstrongest(pretype::Vector{<:Union{Missing, Int64}}, posttype::Vector{<:Union{Missing, Int64}}; dims = 1, threshold = 1)
    npre = length(pretype)
    npost = length(posttype)
    
    if dims == 1
        out = Array{Union{Missing, Int64}}(missing, 1, npost)
        ind = .~ismissing.(posttype)
        ass, ignore = hungarianmax(W[id2ind.(pretype), id2ind.(posttype[ind])]')
        for (i, j) in enumerate(findall(ind))
            if ass[i] > 0
                if W[id2ind(pretype[ass[i]]), id2ind(posttype[j])] < threshold
                    ass[i] = 0
                end
            end
        end
        out[ind] = map(x-> x==0 ? missing : pretype[x], ass)
    else
        out = Array{Union{Missing, Int64}}(missing, npre, 1)
        ind = .~ismissing.(pretype)
        ass, ignore = hungarianmax(W[id2ind.(pretype[ind]), id2ind.(posttype)])
        for (i, j) in enumerate(findall(ind))
            if ass[i] > 0
                if W[id2ind(pretype[j]), id2ind(posttype[ass[i]])] < threshold
                    ass[i] = 0
                end
            end
        end
        out[ind] = map(x-> x==0 ? missing : posttype[x], ass)
    end
    return out
end

# %%
"""
For default (`dims=1`), for each cell in `posttype`, find the strongest presynaptic partner in `pretype`
The strength must be above `threshold`. Otherwise the answer is `missing`.
"""
function findstrongest0(pretype::Vector{<:Union{Missing, Int64}}, posttype::Vector{<:Union{Missing, Int64}}; dims = 1, threshold = 1)
    npre = length(pretype)
    npost = length(posttype)
    
    if dims == 1
        out = Array{Union{Missing, Int64}}(missing, 1, npost)
        for (i, id) = enumerate(posttype)
            if ~ismissing(id)
                val, ind = findmax(W[id2ind.(pretype), id2ind.(id)][:])
                if val > threshold
                    out[i] = pretype[ind]
                else
                    out[i] = missing
                end
            end
        end
    else
        out = Array{Union{Missing, Int64}}(missing, npre, 1)
        for (i, id) = enumerate(pretype)
            if ~ismissing(id)
                val, ind = findmax(W[id2ind.(id), id2ind.(posttype)][:])
                if val > threshold
                    out[i] = posttype[ind]
                else
                    out[i] = missing
                end
            end
        end
    end
    return out
end

# %%
function mergepredictions(pred1, pred2)
    n1 = length(pred1)
    n2 = length(pred2)
    @assert n1 == n2
    out = Array{Union{Missing, Int64}}(missing, n1)
    for i = 1:n1
        if ismissing(pred1[i])
            if ismissing(pred2[i])
                out[i] = missing
            else
                out[i] = pred2[i]
            end
        else
            if ismissing(pred2[i])
                out[i] = pred1[i]
            else
                if pred1[i] == pred2[i]
                    out[i] = pred1[i]
                else
                    out[i] = missing
                end
            end
        end
    end
    return out
end             

# %%
mergepredictions([3, missing, 4], [3, 2, 5])

# %% [markdown]
# ## modular

# %% [markdown]
# ### L1
# all assigned to proper coordinates, based on visualization, except for a few at border of eye
#
# Hungarian assignment of L1 and Mi1 are same as naive, except for four L1 that lack postsynaptic Mi1, and one L1 that synapses onto two Mi1 that are assigned to other L1

# %%
# L1 -> Mi1
columns_df.L1 = findstrongest(type2ids("L1"), columns_df.Mi1, dims = 1, threshold = 1)[:]
println(sum(ismissing.(columns_df.L1)), " columns have no L1 assigned")

# %%
# leftover L1 cells that were not assigned to columns
setdiff(type2ids("L1"), columns_df.L1)

# %% [markdown]
# #### manual additions
# 720575940625354644 and 720575940628173394 are colocated with Mi1 cells, but no synapses\
# fix manually

# %%
columns_df.L1[columns_df.Mi1 .== 720575940636037983] .= 720575940625354644
columns_df.L1[columns_df.Mi1 .== 720575940624959927] .= 720575940628173394

# %%
W[Name.(720575940625354644), Name.(720575940636037983)]  # zero synapses but this is wrong

# %%
W[Name.(720575940628173394), Name.(720575940624959927)]

# %% [markdown]
# #### other assignment failures
# 720575940605282848 has no Mi1 target. It should be at (18, 30), a location outside the current Mi1 map.\
# 720575940630884087 synapses onto two Mi1 cells, each of which is assigned to an L1 cell\
# 720575940636953445 has no Mi1 target

# %% [markdown]
# ### L5

# %%
# L5 -> Mi1 (preferred)
columns_df.L5 = findstrongest(type2ids("L5"), columns_df.Mi1, dims = 1)[:]
ismissing.(columns_df.L5) |> sum

# %%
# L1 -> L5
# @mfalse columns_df.L5 = findstrongest(columns_df.L1, ind2id[ind2type .== "L5"], dims = 2)[:]
# ismissing.(columns_df.L5) |> sum

# %%
# L5 -> Mi1, L1 -> L5
# @mfalse columns_df.L5 = mergepredictions(findstrongest(ind2id[ind2type .== "L5"], columns_df.Mi1, dims = 1)[:], findstrongest(columns_df.L1, ind2id[ind2type .== "L5"], dims = 2)[:])
# ismissing.(columns_df.L5) |> sum

# %% [markdown]
# ### L2

# %%
# L2 -> L5 (preferred)
columns_df.L2 = findstrongest(type2ids("L2"), columns_df.L5, dims = 1)[:]
ismissing.(columns_df.L2) |> sum

# %%
# C3 -> L2 
# columns_df.L2 = findstrongest(columns_df.C3, L2, dims = 2)[:]
# ismissing.(columns_df.L2) |> sum

# %%
# C3 -> L2, L2 -> L5
# columns_df.L2 = mergepredictions(findstrongest(columns_df.C3, L2, dims = 2)[:], findstrongest(L2, columns_df.L5, dims = 1)[:])
# ismissing.(columns_df.L2) |> sum

# %% [markdown]
# ### L3

# %%
# L3 -> Mi1
columns_df.L3 = findstrongest(type2ids("L3"), columns_df.Mi1, dims = 1)[:]
ismissing.(columns_df.L3) |> sum

# %% [markdown]
# ### C3

# %%
# L1 -> C3 (preferred)
columns_df.C3 = findstrongest(columns_df.L1, type2ids("C3"), dims = 2)[:]
ismissing.(columns_df.C3) |> sum

# %%
# L5 -> C3 
# some of these must be disagreeing with L1 -> C3
# columns_df.C3 = findstrongest(columns_df.L5, C3, dims = 2)[:]
# ismissing.(columns_df.C3) |> sum

# %%
# L1 -> C3, L5 -> C3
# columns_df.C3 = mergepredictions(findstrongest(columns_df.L1, C3, dims = 2)[:], findstrongest(columns_df.L5, C3, dims = 2)[:])
# ismissing.(columns_df.C3) |> sum

# %%
# L3 -> C3
# columns_df.C3 = findstrongest(columns_df.L3, C3, dims = 2)[:]
# ismissing.(columns_df.C3) |> sum

# %%
# L1 -> C3, L3 -> C3
# columns_df.C3 = mergepredictions(findstrongest(columns_df.L1, C3, dims = 2)[:], findstrongest(columns_df.L3, C3, dims = 2)[:])
# ismissing.(columns_df.C3) |> sum

# %% [markdown]
# ### C2
# unicolumnar input from L1, but L5 and Mi1 input includes nearest neighbors

# %%
# L1 -> C2 (preferred, looks better in ng)
columns_df.C2 = findstrongest(columns_df.L1, type2ids("C2"), dims = 2)[:]
ismissing.(columns_df.C2) |> sum

# %%
# L5 -> C2 
# columns_df.C2 = findstrongest(columns_df.L5, C2, dims = 2)[:]
# ismissing.(columns_df.C2) |> sum

# %%
# L5 -> C2, L1 -> C2
# lots of disagreement
# columns_df.C2 = mergepredictions(findstrongest(columns_df.L5, C2, dims = 2)[:], findstrongest(columns_df.L1, C2, dims = 2)[:])
# ismissing.(columns_df.C2) |> sum

# %%
# C2 -> L2
# columns_df.C2 = findstrongest(C2, columns_df.L2, dims = 1)[:]
# ismissing.(columns_df.C2) |> sum

# %%
# L1 -> C2, C2 -> L2
# columns_df.C2 = mergepredictions(findstrongest(columns_df.L1, C2, dims = 2)[:], findstrongest(C2, columns_df.L2, dims = 1)[:])
# ismissing.(columns_df.C2) |> sum

# %% [markdown]
# ### Tm1
# L2 is top input of Tm1

# %%
columns_df.Tm1 = findstrongest(columns_df.L2, type2ids("Tm1"), dims = 2)[:]
ismissing.(columns_df.Tm1) |> sum

# %% [markdown]
# ### Tm2
# L2 is top input of Tm2

# %%
columns_df.Tm2 = findstrongest(columns_df.L2, type2ids("Tm2"), dims = 2)[:]
ismissing.(columns_df.Tm2) |> sum

# %% [markdown]
# ### Mi4

# %%
columns_df.Mi4 = findstrongest(columns_df.L5, type2ids("Mi4"), dims = 2)[:]
ismissing.(columns_df.Mi4) |> sum

# %% [markdown]
# ### Mi9

# %%
columns_df.Mi9 = findstrongest(columns_df.Mi4, type2ids("Mi9"), dims = 2)[:]
ismissing.(columns_df.Mi9) |> sum

# %% [markdown]
# ### Tm9
# L3 and Mi4 are top inputs

# %%
# L3 -> Tm9
columns_df.Tm9 = findstrongest(columns_df.L3, type2ids("Tm9"), dims = 2)[:]
ismissing.(columns_df.Tm9) |> sum

# %%
# Mi4 -> Tm9 (preferred)
# @mfalse columns_df.Tm9 = findstrongest(columns_df.Mi4, ind2id[ind2type .== "Tm9"], dims = 2)[:]
# ismissing.(columns_df.Tm9) |> sum

# %% [markdown]
# ### Tm20
# L3 and Mi4 are top inputs

# %%
# Mi4 -> Tm20 (preferred)
columns_df.Tm20 = findstrongest(columns_df.Mi4, type2ids("Tm20"), dims = 2)[:]
ismissing.(columns_df.Tm20) |> sum

# %%
# L3 -> Tm20
# @mfalse columns_df.Tm20 = findstrongest(columns_df.L3, ind2id[ind2type .== "Tm20"], dims = 2)[:]
# ismissing.(columns_df.Tm20) |> sum

# %% [markdown]
# ### L4
# This is numerous, but is it modular? Some connectivity with neighboring columns.

# %%
# preferred
columns_df.L4 = findstrongest(type2ids("L4"), columns_df.Tm2, dims = 1)[:]
ismissing.(columns_df.L4) |> sum

# %%
# @mfalse columns_df.L4 = findstrongest(ind2id[ind2type .== "L4"], columns_df.L2, dims = 1)[:]
# ismissing.(columns_df.L4) |> sum

# %%
# C2 -> L4 
# @mfalse columns_df.L4 = findstrongest(columns_df.C2, ind2id[ind2type .== "L4"], dims = 2)[:]
# ismissing.(columns_df.L4) |> sum

# %%
# L4 -> L2, C2 -> L4
# clearly inconsistencies
# columns_df.L4 = mergepredictions(findstrongest(L4, columns_df.L2, dims = 1)[:], findstrongest(columns_df.C2, L4, dims = 2)[:])
# ismissing.(columns_df.L4) |> sum

# %% [markdown]
# ### T1

# %%
# L2 -> T1 (preferred)
columns_df.T1 = findstrongest(columns_df.L2, type2ids("T1"), dims = 2)[:]
ismissing.(columns_df.T1) |> sum

# %%
columns_df

# %%
# @save "ColumnAssignments.jld2" columns_df

# %% [markdown]
# ## not modular

# %%
# @load "ColumnAssignments.jld2"

# %% [markdown]
# ### Tm3

# %%
columns_df.Tm3 = findstrongest(columns_df.Mi1, type2ids("Tm3"), dims = 2)[:]
ismissing.(columns_df.Tm3) |> sum

# %% [markdown]
# ### Tm4

# %%
columns_df.Tm4 = findstrongest(columns_df.L2, type2ids("Tm4"), dims = 2)[:]
ismissing.(columns_df.Tm4) |> sum

# %% [markdown]
# ### T2

# %%
columns_df.T2 = findstrongest(columns_df.L5, type2ids("T2"), dims = 2)[:]
ismissing.(columns_df.T2) |> sum

# %% [markdown]
# ### T2a

# %%
columns_df.T2a = findstrongest(columns_df.Mi1, type2ids("T2a"), dims = 2)[:]
ismissing.(columns_df.T2a) |> sum

# %% [markdown]
# ### T3

# %%
columns_df.T3 = findstrongest(columns_df.Mi1, type2ids("T3"), dims = 2)[:]
ismissing.(columns_df.T3) |> sum

# %%
# @save "ColumnAssignmentsExtended.jld2" columns_df

# %% [markdown]
# ### T4
# Mi1 is top input

# %%
for celltype in ["T4a", "T4b", "T4c", "T4d"]
    columns_df[!, celltype] = findstrongest(columns_df.Mi1, type2ids(celltype), dims = 2)[:]
    println(celltype, " ", sum(ismissing.(columns_df[!, celltype])))
end

# %% [markdown]
# ### T5

# %%
for celltype in ["T5a", "T5b", "T5c", "T5d"]
    columns_df[!, celltype] = findstrongest(columns_df.Tm1, type2ids(celltype), dims = 2)[:]
    println(celltype, " ", sum(ismissing.(columns_df[!, celltype])))
end

# %% [markdown]
# ### TmY5a

# %%
columns_df.TmY5a = findstrongest(columns_df.Tm4, type2ids("TmY5a"), dims = 2)[:]
ismissing.(columns_df.TmY5a) |> sum

# %% [markdown]
# ### Y3

# %%
# Tm1 -> Y3
columns_df.Y3 = findstrongest(columns_df.Tm1, type2ids("Y3"), dims = 2)[:]
ismissing.(columns_df.Y3) |> sum
