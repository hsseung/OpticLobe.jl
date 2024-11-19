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
#     display_name: Julia 1.10.4
#     language: julia
#     name: julia-1.10
# ---

# %% [markdown]
# # Discriminating with 2D projections

# %%
using OpticLobe

# %%
using NamedArrays, SparseArrays

# %%
using Plots

# %%
using Measures

# %%
using Printf

# %%
using ProgressMeter

# %%
using StatsPlots, StatsBase

# %%
using FreqTables

# %%
using Distances, Clustering

# %%
using MissingsAsFalse

# %%
using OrderedCollections

# %%
#TARGETDIR = "/Users/sseung/Active/v783/Discriminators2D"
TARGETDIR = "."

# %%
Plots.default(fontfamily = "Helvetica")

# %% [markdown]
# ## projection vectors for cells

# %%
Av = A[:, visualtypes]

# %%
predicates = vcat("in-" .* visualtypes, "out-" .* visualtypes)
normalization = 1 ./ sum(W, dims = 1).array
infrac = Wtc[visualtypes, :].array .* normalization
normalization = 1 ./ sum(W, dims = 2).array
outfrac = Wct[:, visualtypes].array .* normalization
inoutfraction = NamedArray(vcat(infrac, outfrac'), names = (predicates, ind2id))

# %% [markdown]
# ## function definitions

# %%
# scatter plot of input/output fractions for cell inds
function scatterfeature(featureX, featureY, neuropil_class, targettype = nothing)
    inds = getindex.(findall(Ao[:, cellclasses[neuropil_class]] .== 1), 1)
    if targettype == nothing
        scatter(inoutfraction[featureX, inds], inoutfraction[featureY, inds], legend=:none, xlabel = featureX, ylabel=featureY)
    else
        mcolors = fill(palette(:default)[1], length(inds))
        mcolors[ind2type[inds] .== targettype] .= palette(:default)[2]
        scatter(inoutfraction[featureX, inds], inoutfraction[featureY, inds], legend=:none, xlabel = featureX, ylabel=featureY, markercolor = mcolors, markerstrokecolor = mcolors)
    end
end

# %%
# scatter plot of input/output fractions for cell inds
function scatterfeature2(featureX, featureY, neuropil_class, targettype = nothing)
    inds = getindex.(findall(Ao[:, cellclasses[neuropil_class]] .== 1), 1)
    if targettype == nothing
        scatter(inoutfraction[featureX, inds], inoutfraction[featureY, inds], legend=:none, xlabel = featureX, ylabel=featureY)
    else
        istarget = ind2type[inds] .== targettype        
        scatter(inoutfraction[featureX, inds[.~istarget]], inoutfraction[featureY, inds[.~istarget]], legend=:none, xlabel = featureX, ylabel=featureY, msc=:auto)
        scatter!(inoutfraction[featureX, inds[istarget]], inoutfraction[featureY, inds[istarget]], legend=:none, xlabel = featureX, ylabel=featureY, msc=:auto)
    end
end

# %%
# scatter plot of input/output fractions for cell inds
function scatterfeature3(featureX, featureY, backgroundtypes, targettype, thresholdpair = nothing, fscore = nothing)
    inds = getindex.(findall(Av[:, backgroundtypes] .== 1), 1)
    istarget = ind2type[inds] .== targettype        
    scatter(inoutfraction[featureX, inds[.~istarget]], inoutfraction[featureY, inds[.~istarget]], 
        legend = :none, 
        xlabel = featureX, ylabel=featureY, 
        msc=:auto,
        markersize = 2,
    )
    if thresholdpair != nothing
        vline!([thresholdpair[1]])
        hline!([thresholdpair[2]])
    end
    if fscore != nothing
        annotate!((0.75, 0.9), text(@sprintf("F-score = %3.3f", fscore), "Helvetica"))
    end
    scatter!(inoutfraction[featureX, inds[istarget]], inoutfraction[featureY, inds[istarget]], 
        legend=:none, 
        xlabel = featureX, ylabel=featureY, 
        msc=:auto,
        markersize = 2
    )
end

# %%
function histfeature(feature, backgroundtypes, targettype)
    inds = getindex.(findall(Av[:, backgroundtypes] .== 1), 1)
    istarget = ind2type[inds] .== targettype  
    groupedhist(inoutfraction[feature, inds], group=istarget, bar_position=:stack, ylim=(0, 20), nbins = 100, legend=:none, xlabel = feature, ylabel = "number of cells")
end

# %%
function findpredicates(targettype, family, precthreshold, recallthreshold)
    prec = NamedArray(zeros(length(predicates), length(thresholds)), names=(predicates, 1:length(thresholds)))
    recall = NamedArray(zeros(length(predicates), length(thresholds)), names=(predicates, 1:length(thresholds)))
    @mfalse targetinds = ind2type .== targettype
    truenum = sum(Ao[:, targettype])
    familyinds = getindex.(findall(Ao[:, cellclasses[family]] .== 1), 1)
    for (i, thres) in enumerate(thresholds)
        tpos = sum(inoutfraction[predicates, targetinds] .> thres, dims=2)[:]
        pos = sum(inoutfraction[predicates, familyinds] .> thres, dims=2)[:]
        recall[:, i] = tpos./truenum
        prec[:, i] = tpos./(eps(Float32) .+ pos)
    end
    candidates = findall(recall .> recallthreshold .&& prec .> precthreshold)
    return unique(predicates[getindex.(candidates, 1)]), prec, recall
end

# %%
"""
evaluate precision and recall for discrimination of `targetids` using conjunctions of two thresholded features

    `backgroundids` should be superset of `targetids`

    `prec[pred1, pred2, i, j]` - precision 
    `recall[pred1, pred2, i, j]` - recall 
"""
function tryconjunctions(targetids, backgroundids, totrythreshold = 0.01, nlevel = 25)
    targetinds = id2ind.(targetids)
    backgroundinds = id2ind.(backgroundids)

    ntarget = length(targetinds)
    levels = (nlevel:-1:1)/nlevel

    # only try predicates for which inoutfraction averaged over target cells exceeds totrythreshold
    totry = filter(x -> mean(inoutfraction[x, targetinds]) .> totrythreshold, predicates)   
    ntry = length(totry)
    
    prec = zeros(ntry, ntry, nlevel, nlevel)
    recall = zeros(ntry, ntry, nlevel, nlevel)

    targetfrac = Matrix(inoutfraction[totry, targetinds])
    backgroundfrac = Matrix(inoutfraction[totry, backgroundinds])
    
    maxthreshold = mean(targetfrac, dims=2)   # maximum value of threshold in conjunctions

    @showprogress for k = 1:ntry
        thresholds1 = maxthreshold[k]*levels
        for l = 1:k
            thresholds2 = maxthreshold[l]*levels
            for (i, thres1) in enumerate(thresholds1)
                targetthres1 = targetfrac[k, :] .> thres1
                backgroundthres1 = backgroundfrac[k, :] .> thres1
                for (j, thres2) in enumerate(thresholds2)
                    tpos = sum(targetfrac[l, targetthres1] .> thres2)
                    pos = sum(backgroundfrac[l, backgroundthres1] .> thres2)
                    recall[k, l, i, j] = tpos/ntarget
                    prec[k, l, i, j] = pos > 0 ? tpos/pos : 0
                end
            end
        end
    end
    return prec, recall, totry, maxthreshold
end

# %%
"""
evaluate precision and recall for discrimination of `targetids` using conjunctions of two thresholded features

    `backgroundids` should be superset of `targetids`

    `prec[pred1, pred2, i, j]` - precision 
    `recall[pred1, pred2, i, j]` - recall 
"""
# slower version but easier to understand?
function tryconjunctions2(targetids, backgroundids, totrythreshold = 0.01, nlevel = 25)
    targetinds = id2ind.(targetids)
    backgroundinds = id2ind.(backgroundids)

    ntarget = length(targetinds)
    levels = (nlevel:-1:1)/nlevel

    # only try predicates for which inoutfraction averaged over target cells exceeds totrythreshold
    totry = filter(x -> mean(inoutfraction[x, targetinds]) .> totrythreshold, predicates)   
    ntry = length(totry)
    
    prec = zeros(ntry, ntry, nlevel, nlevel)
    recall = zeros(size(prec))
    thresholds1 = zeros(size(prec))
    thresholds2 = zeros(size(prec))

    targetfrac = Matrix(inoutfraction[totry, targetinds])
    backgroundfrac = Matrix(inoutfraction[totry, backgroundinds])
    
    maxthreshold = mean(targetfrac, dims=2)   # maximum value of threshold in conjunctions

    @showprogress for i = 1:nlevel
        for k = 1:ntry
            thres1 = maxthreshold[k]*levels[i]
            targetthres1 = targetfrac[k, :] .> thres1
            backgroundthres1 = backgroundfrac[k, :] .> thres1
            for j = 1:nlevel
                for l = 1:ntry
                    thres2 = maxthreshold[l]*levels[j]
                    thresholds1[k, l, i, j] = thres1
                    thresholds2[k, l, i, j] = thres2
                    tpos = sum(targetfrac[l, targetthres1] .> thres2)
                    pos = sum(backgroundfrac[l, backgroundthres1] .> thres2)
                    recall[k, l, i, j] = tpos/ntarget
                    prec[k, l, i, j] = pos > 0 ? tpos/pos : 0
                end
            end
        end
    end
    return prec, recall, totry, maxthreshold
end

# %% [markdown]
# ## example in Figure S3c

# %%
# example in paper
plot(
    scatterfeature3( "out-TmY3", "in-C3", family2types["Pm"], "Pm04"),
#    fontfamily = "Helvetica",
    title = "Pm04",
    fontsize = 14,
    tickfontsize = 12,
    size = (400, 400)
    )

# %%
savefig("Fig S3c example 2D discriminator.pdf")

# %% [markdown]
# ## 2D discriminators for all families that contain more than one type, excluding photoreceptors

# %%
for (familyname, typenames) in family2types
    if length(typenames) > 1 && !startswith(familyname, "R")
        println(familyname, " ", length(typenames))
    end
end

# %%
nlevel = 25
npanel = 25
for (familyname, typenames) in family2types
    if length(typenames) > 1 && !startswith(familyname, "R")
        familyids = vcat(type2ids.(typenames)...)

        for targettype in typenames
            println(targettype)        
            @mfalse targetids = type2ids(targettype)
            prec, recall, totry, maxthreshold = tryconjunctions(targetids, familyids, 0.01, nlevel)
            
            Fscore = 2*prec.*recall
            valid = Fscore .> 0    
            Fscore[valid] = Fscore[valid]./(prec[valid] .+ recall[valid])
            vals, inds = findmax(Fscore, dims=(3, 4))
            vals = vals[:, :, 1, 1]
            inds = inds[:, :, 1, 1]
            
            rankedlist = sort([[vals[i, j], [i, j]] for i=1:length(totry) for j=1:i], rev=true)
            
            desiredlist = rankedlist[1:npanel]
            
            thresholdindexpair = [collect(Tuple(inds[i, j])[3:4]) for i=1:length(totry), j=1:length(totry)]
            
            thresholdvaluepair = [maxthreshold[item[2]].*(nlevel .+ 1 .- thresholdindexpair[item[2]...])/nlevel for item in desiredlist]
            
            plot(
                [scatterfeature3(totry[desiredlist[i][2]]..., typenames, targettype, thresholdvaluepair[i], desiredlist[i][1]) for i=1:npanel]...,
                size=(2200, 2000),
                plot_title = targettype,
                left_margin = 10mm,
                right_margin = 10mm,
                )
            savefig(joinpath(TARGETDIR, "$targettype.pdf"))
        end
    end
end
