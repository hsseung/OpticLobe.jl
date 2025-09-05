# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:percent
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.3
#   kernelspec:
#     display_name: Julia 1.11.6
#     language: julia
#     name: julia-1.11
# ---

# %% [markdown]
# # Form vision paper (spatial figures)
#
# For reproducing figures in Seung, [Predicting visual function by interpreting a neuronal wiring diagram](https://doi.org/10.1038/s41586-024-07953-5), Nature 634:113-123 (2024).
#
# This notebook is restricted to figures about spatial analysis.

# %%
using OpticLobe

include("config.jl")

# %%
using Missings, MissingsAsFalse

# %%
using StatsBase
using StatsPlots # boxplots
using Measures

# %%
using NamedArrays, SparseArrays, LinearAlgebra

# %%
using OrderedCollections

# %%
#using LaTeXStrings, MathTeXEngine

# %%
using ProgressMeter

# %%
using Printf

# %%
using Luxor

# %%
using ColorSchemes

# %% [markdown]
# ## cell types

# %%
Dm3types = ["Dm3v", "Dm3p", "Dm3q"]
TmYtypes = ["TmY4", "TmY9q", "TmY9q⊥"]

# %%
hexel = ["L1", "L2", "L3", "L4", "L5", "Tm1", "Tm2", "Tm9", "Mi1", "Mi4", "Mi9"]
nothexel = setdiff(intrinsictypes, hexel)

# %% [markdown]
# ## helper functions

# %%
"""
for orientation boxplots
"""
function wraparound(theta; center = 0)
    if theta > center + pi/2 
        return theta - pi
    end
    if theta < center - pi/2
        return theta + pi
    end
    return theta
end

# %%
# Deprecated because LaTeX fonts don't play well with Illustrator
# had to give up on superscript ⊥
#function type2label(s::String)
#    return latexstring(replace(s, "⊥" => "^⊥"))
#end

#function path2label(v::Vector{String})
#    return latexstring(replace(join(v, "→"), "⊥" => "^⊥"))
#end

# %%
function type2label(s::String)
    return s
end

function path2label(v::Vector{String})
    return join(v, "→")
end

# %% [markdown]
# ## add Dm3 and TmY locations to `id2pq`

# %%
# missing means that there is zero connectivity
for posttype in vcat(Dm3types, TmYtypes, "Y3")
    rfs = preimage("Tm1", posttype);
    @mfalse for (i, id) in enumerate(ind2id[(ind2type .== posttype) .& (ind2side .== "right")])
        center = findcenter(float(rfs[i]), seven)
        if ~ismissing(center)
            id2pq[id] = center
        end
    end
end

# special case that has no Tm1 inputs (bad proofreading)
id2pq[720575940618964225] = findcenter(prepreimage("Mi1", "Y3", 720575940618964225), seven)

# %% [markdown]
# ## Fig 1 Tm1-Dm3 connectivity maps

# %% [markdown]
# ### Fig 1bcd average Tm1-Dm3 maps

# %%
#Drawing(800, 800, "Tm1Dm3Averages.pdf")
Drawing(400, 450)
origin()
kernels = typetriad("Tm1", Dm3types, radius = 3, hexelsize = 16, normalize = :common)
finish()
preview()

# %% [markdown]
# ### Fig 1f Tm1-Dm3v example cells on eye maps

# %%
rfs = preimage("Tm1", "Dm3v");
maximum.(rfs[6:8])

# %%
#Drawing(600, 700, "Tm1Dm3vExamples.pdf")
Drawing(600, 700)
origin()
eyetriad("Tm1", type2ids("Dm3v")[6:8])
finish()
preview()

# %% [markdown]
# ### Fig 1g Tm1-Dm3v example cells, aligned and cropped

# %%
# Drawing(450, 500, "Tm1Dm3vExamplesAlignedCropped.pdf")
Drawing(450, 500)
origin()
celltriad("Tm1", type2ids("Dm3v")[6:8], radius = 5, hexelsize = 12, ellipse = true)
finish()
preview()

# %% [markdown]
# ### Fig 1h Tm1-Dm3 ellipses

# %%
#Drawing(600, 600, "Dm3Ellipses.pdf")
Drawing(600, 600)
origin()
scale(50)
setline(0.5)

for (celltype, ellipsecolor) in zip(Dm3types, ["red", "green", "blue"])
    setcolor(ellipsecolor)
    ellipses = ellipsesummary.(preimage("Tm1", celltype))
    for e in sample(ellipses, 50, replace=false)
        if !ismissing(e)
            Luxor.rotate(-pi/6 + e.theta)
            Luxor.ellipse(0, 0, e.minor, e.major, :stroke)   # sqrt(3) for conversion to units of latticeconstant
            Luxor.rotate(pi/6 - e.theta)
        end
    end
end
sethue("black")
setline(4)
drawpqaxes(1)
finish()
preview()

# %% [markdown]
# ### Fig 1ij Tm1-Dm3 angles and aspect ratio boxplots

# %%
default(; # Plots defaults
    fontfamily="Arial Unicode",
    label="", # only explicit legend entries
    guidefontsize = 12,
    tickfontsize = 12,
    legendfontsize = 10,
#    linewidth = 2,
#    ylim = (0, Inf), 
    widen = true,
#    xrotation = 30
    )

# %%
ellipses = [ellipsesummary.(preimage("Tm1", Dm3type)) for Dm3type in Dm3types]
aspectratios = [[e.major/e.minor for e in skipmissing(ellipse)] for ellipse in ellipses]
angles = [[(wraparound(e.theta - pi/6)) for e in skipmissing(ellipse)] for ellipse in ellipses]
diamsmajor = [[e.major for e in skipmissing(ellipse)] for ellipse in ellipses]
diamsminor = [[e.minor for e in skipmissing(ellipse)] for ellipse in ellipses]
#diamsmajor = [[passmissing(x -> x.major)(e) for e in ellipse] for ellipse in ellipses]   # need to skip because quantiles don't work with missing values
#diamsminor = [[passmissing(x -> x.major)(e) for e in ellipse] for ellipse in ellipses]
Dm3ellipses = OrderedDict(Dm3types .=> ellipses)
Dm3aspectratios = OrderedDict(Dm3types .=> aspectratios)
Dm3angles = OrderedDict(Dm3types .=> angles)
celltypes = vcat([fill(i, length(v)) for (i, v) in enumerate(values(Dm3aspectratios))]...)

haspect = boxplot(celltypes, vcat(values(Dm3aspectratios)...), ylabel = "aspect ratio", 
    xtick = strings2ticks(Dm3types), ylim = (1, Inf)
)

hangle = boxplot(celltypes, vcat(values(Dm3angles)...)*180/pi, ylabel = "orientation (deg)", 
    xtick = strings2ticks(Dm3types), yticks = [-60, 0, 60]
)

plot(hangle, haspect, size = (800, 300), left_margin = 5Measures.mm, bottom_margin = 3Measures.mm)

# %%
savefig(joinpath(TARGETDIR, "Tm1Dm3anglesaspects.pdf"))

# %% [markdown]
# ### Fig 1klm Dm3 projections

# %%
paths = [["Tm1", Dm3type] for Dm3type in Dm3types]
rfs = inmaps.(tracetypes.(paths))
radius = 4
rfscrop = [passmissing(crop).(rf.array, id2pq.(type2ids(Dm3type)), radius, 1) for (rf, Dm3type) in zip(rfs, Dm3types)];

# %%
long = [hexproject.(skipmissing(ims), ax) for (ims, ax) in zip(rfscrop, [:v, :p, :q])]
trans = [hexproject.(skipmissing(ims), ax) for (ims, ax) in zip(rfscrop, [:h, :p⊥, :q⊥])];

# %%
Plots.reset_defaults()

# %%
default(; # Plots defaults
    fontfamily="Arial Unicode",
    label="", # only explicit legend entries
    guidefontsize = 12,
    tickfontsize = 12,
    legendfontsize = 10,
#    linewidth = 2,
#    ylim = (0, Inf), 
    widen = true,
#    xrotation = 30
    )

# %%
h = [groupedbar(-radius:radius, [a c], yerror = [b d], 
    linewidth = 1, 
    label = ["longitudinal" "transverse"], 
    ylabel = "synapses",
    xlabel = "displacement",
        annotations = (-3, 25, Plots.text(e, "Helvetica")),
        ylim = (0, 30)
) for (a, b, c, d, e) in zip(mean.(long), std.(long), mean.(trans), std.(trans), Dm3types)]

# %%
plot(h..., layout = (1, 3), size = (1000, 200), bottom_margin = 10Measures.mm, left_margin = 5Measures.mm)

# %%
savefig(joinpath(TARGETDIR, "Dm3projections.pdf"))

# %% [markdown]
# ## Fig 2 Tm1-TmY connectivity maps

# %% [markdown]
# ### Fig 2def average Tm1-TmY maps

# %%
#Drawing(800, 800, "Tm1TmYAverages.pdf")
Drawing(400, 450)
origin()
kernels = typetriad("Tm1", TmYtypes, radius = 3, hexelsize = 16, normalize = :common)
finish()
preview()

# %%
#Drawing(800, 800, "Tm1TmYAverages.pdf")
Drawing(400, 450)
origin()
kernels = typetriad("Tm1", TmYtypes, radius = 3, hexelsize = 16, normalize = :separate)
finish()
preview()

# %% [markdown]
# ### Fig 2g example TmY4 cells

# %%
ids = type2ids("TmY4")
# Drawing(800, 800, "Tm1TmY4ExamplesAligned.pdf")
@drawsvg rfs = celltriad("Tm1", ids[13:15], radius = 5, hexelsize = 12) 430 500

# %% [markdown]
# ### Fig 2h example TmY9q cells

# %%
@mfalse ids = type2ids("TmY9q")
#Drawing(800, 800, "Tm1TmY9qExamplesAligned.pdf")
Drawing(800, 800)
@drawsvg celltriad("Tm1", ids[1:3], radius = 5, hexelsize = 12) 430 500

# %% [markdown]
# ### Fig 2i example TmY9q⊥ cells

# %%
@mfalse ids = ind2id[ind2type .== "TmY9q⊥"]
#Drawing(800, 800, "Tm1TmY9qperpExamplesAligned.pdf")
Drawing(800, 800)
@drawsvg celltriad("Tm1", ids[1:3], radius = 5, hexelsize = 12) 430 500

# %% [markdown]
# ### Fig 2j Tm1-TmY ellipses

# %%
#Drawing(600, 600, "TmYEllipses.pdf")
Drawing(600, 600)
origin()
scale(50)
setline(0.5)

for (celltype, ellipsecolor) in zip(TmYtypes, ["red", "green", "blue"])
    setcolor(ellipsecolor)
    ellipses = ellipsesummary.(preimage("Tm1", celltype))
    for e in sample(ellipses, 50, replace=false)
        if ~ismissing(e)
            Luxor.rotate(-pi/6 + e.theta)
            Luxor.ellipse(0, 0, e.minor, e.major, :stroke)   # sqrt(3) for conversion to units of latticeconstant
            Luxor.rotate(pi/6 - e.theta)
        end
    end
end
sethue("black")
setline(4)
drawpqaxes(1)
finish()
preview()

# %% [markdown]
# ### Fig 2kl Tm1-TmY angles and aspect ratios

# %%
ellipses = [ellipsesummary.(preimage("Tm1", TmYtype)) for TmYtype in TmYtypes]
aspectratios = [[e.major/e.minor for e in skipmissing(ellipse)] for ellipse in ellipses]
angles = [[wraparound(e.theta - pi/6, center = pi/2) for e in skipmissing(ellipse)] for ellipse in ellipses]
diamsmajor = [[e.major for e in skipmissing(ellipse)] for ellipse in ellipses]
diamsminor = [[e.minor for e in skipmissing(ellipse)] for ellipse in ellipses]
TmYellipses = OrderedDict(TmYtypes .=> ellipses)
TmYaspectratios = OrderedDict(TmYtypes .=> aspectratios)
TmYangles = OrderedDict(TmYtypes .=> angles)
celltypes = vcat([fill(i, length(v)) for (i, v) in enumerate(values(TmYaspectratios))]...)

# %%
default(; # Plots defaults
    fontfamily="Arial Unicode",
    label="", # only explicit legend entries
    guidefontsize = 12,
    tickfontsize = 12,
    legendfontsize = 10,
#    linewidth = 2,
#    ylim = (0, Inf), 
    widen = true,
#    xrotation = 30
    )

# %%
haspect = boxplot(celltypes, vcat(values(TmYaspectratios)...), ylabel = "aspect ratio", 
    xtick = strings2ticks(TmYtypes), ylim = (1, Inf)
)

# %%
celltypes = vcat([fill(i, length(v)) for (i, v) in enumerate(values(TmYangles))]...)
hangle = boxplot(celltypes, vcat(values(TmYangles)...)*180/pi, ylabel = "orientation (deg)", 
    xtick = strings2ticks(TmYtypes), yticks = [60, 90, 150]
)

# %%
plot(hangle, haspect, size = (800, 300), left_margin = 5Measures.mm, bottom_margin = 3Measures.mm)

# %%
savefig(joinpath(TARGETDIR, "Tm1TmYanglesaspects.pdf"))

# %% [markdown]
# ### Fig 2mno TmY projections

# %%
paths = [["Tm1", TmYtype] for TmYtype in TmYtypes]
rfs = inmaps.(tracetypes.(paths))
radius = 4
rfscrop = [passmissing(crop).(rf.array, id2pq.(type2ids(TmYtype)), radius, 1) for (rf, TmYtype) in zip(rfs, TmYtypes)];

# %%
trans = [hexproject.(skipmissing(ims), ax) for (ims, ax) in zip(rfscrop, [:v, :q⊥, :q])]
long = [hexproject.(skipmissing(ims), ax) for (ims, ax) in zip(rfscrop, [:h, :q, :q⊥])];

# %%
Plots.reset_defaults()

# %%
default(; # Plots defaults
    fontfamily="Arial Unicode",
    label="", # only explicit legend entries
    guidefontsize = 12,
    tickfontsize = 12,
    legendfontsize = 10,
#    linewidth = 2,
#    ylim = (0, Inf), 
    widen = true,
#    xrotation = 30
    )

# %%
h = [groupedbar(-radius:radius, [a c], yerror = [b d], 
    linewidth = 1, 
    label = ["longitudinal" "transverse"], 
    ylabel = e == "TmY4" ? "synapses" : "",
    xlabel = "displacement",
        annotations = (-3, 15, Plots.text(e, "Arial Unicode")),
        ylim = (0, 20)
) for (a, b, c, d, e) in zip(mean.(long), std.(long), mean.(trans), std.(trans), TmYtypes)]

# %%
plot(h..., layout = (1, 3), size = (1000, 200), bottom_margin = 10Measures.mm, left_margin = 5Measures.mm)

# %%
savefig(joinpath(TARGETDIR, "TmYprojections.pdf"))

# %% [markdown]
# ## Fig 4

# %% [markdown]
# ### Fig 4abcd Dm3 average ERFs

# %%
# threshold eliminated. now just the six strongest disynaptic pathways for each target
hexelsize = 7

monosynaptic = [["Tm1", Dm3type] for Dm3type in Dm3types]
crfs = inmaps.(tracetypes.(monosynaptic))
crfsellipses = [ellipsesummary.(z.array) for z in crfs]
crfscrop = [passmissing(crop).(ims.array, id2pq.(ids), 5) for (ims, ids) in zip(crfs, type2ids.(Dm3types))]
crfsavg = mean.(skipmissing.(crfscrop))
crfsavgellipses = NamedArray(ellipsesummary.(crfsavg), Dm3types)  # only the ellipse is used here

corner = 260
Drawing(1200, 1200, joinpath(TARGETDIR, "Dm3ERF.pdf"))
#Drawing(1200, 1200)

target = "Dm3p"

paths = [["Tm1", t, target] for t in intrinsictypes]
p = sortperm(scorepath.(paths), rev = true)
paths = paths[p[1:6]]
intermediaries = getindex.(paths, 2)

rfs = inmaps.(tracebacktypes.(paths))
rfscrop = [passmissing(crop).(rf.array, id2pq.(type2ids(target)), 7, 1) for rf in rfs]
rfsavg = [mean(skipmissing(rf)) for rf in rfscrop]
cseries = map(x -> x ? "green" : "magenta", type2nt[intermediaries].array .== "ACH")

s = maximum.(rfsavg)
s[cseries .== "green"] ./= maximum(s[cseries .== "green"])
s[cseries .!= "green"] ./= maximum(s[cseries .!= "green"])

origin()
Luxor.translate(-corner, -corner)
setline(2)
hexannulus(rfsavg, orbitscale = 1.5, text = type2label.(intermediaries), ellipses = ellipsesummary.(rfsavg), ellipsecolors = cseries, maxvals = true)
Luxor.scale(1.5)
Luxor.rotate(-pi/6)
Luxor.translate(-crfsavgellipses[target].center)
Luxor.rotate(pi/6)
setcolor("blue")
fontsize(12)
Luxor.text(target, Point(0, -35), halign = :center)
setline(3)
setdash("dotted")
drawellipse(crfsavgellipses[target])
setdash("solid")
for (i, e) in enumerate(ellipsesummary.(rfsavg))
    setopacity(s[i])
    drawellipse(e, color = cseries[i])
end
Luxor.rotate(-pi/6)
setline(1)
setcolor("black")
Luxor.poly(hexcenter.([HexagonEye(1, 0, hexelsize), HexagonEye(0, 0, hexelsize), HexagonEye(0, 1, hexelsize)]), :stroke)
Luxor.rotate(pi/6)


####################

target = "Dm3q"

paths = [["Tm1", t, target] for t in intrinsictypes]
p = sortperm(scorepath.(paths), rev = true)
paths = paths[p[1:6]]
intermediaries = getindex.(paths, 2)

rfs = inmaps.(tracebacktypes.(paths))
rfscrop = [passmissing(crop).(rf.array, id2pq.(type2ids(target)), 7, 1) for rf in rfs]
rfsavg = [mean(skipmissing(rf)) for rf in rfscrop]
cseries = map(x -> x ? "green" : "magenta", type2nt[intermediaries].array .== "ACH")

s = maximum.(rfsavg)
s[cseries .== "green"] ./= maximum(s[cseries .== "green"])
s[cseries .!= "green"] ./= maximum(s[cseries .!= "green"])

origin()
Luxor.translate(corner, -corner)
setline(2)
hexannulus(rfsavg, orbitscale = 1.5, text = type2label.(intermediaries), ellipses = ellipsesummary.(rfsavg), ellipsecolors = cseries, maxvals = true)
Luxor.scale(1.5)
Luxor.rotate(-pi/6)
Luxor.translate(-crfsavgellipses[target].center)
Luxor.rotate(pi/6)
setcolor("blue")
fontsize(12)
Luxor.text(target, Point(0, -35), halign = :center)
setline(3)
setdash("dotted")
drawellipse(crfsavgellipses[target])
setdash("solid")
for (i, e) in enumerate(ellipsesummary.(rfsavg))
    setopacity(s[i])
    drawellipse(e, color = cseries[i])
end
Luxor.rotate(-pi/6)
setline(1)
setcolor("black")
Luxor.poly(hexcenter.([HexagonEye(1, 0, hexelsize), HexagonEye(0, 0, hexelsize), HexagonEye(0, 1, hexelsize)]), :stroke)
Luxor.rotate(pi/6)


################
target = "Dm3v"

paths = [["Tm1", t, target] for t in intrinsictypes]
p = sortperm(scorepath.(paths), rev = true)
paths = paths[p[1:6]]
intermediaries = getindex.(paths, 2)

rfs = inmaps.(tracebacktypes.(paths))
rfscrop = [passmissing(crop).(rf.array, id2pq.(type2ids(target)), 7, 1) for rf in rfs]
rfsavg = [mean(skipmissing(rf)) for rf in rfscrop]
cseries = map(x -> x ? "green" : "magenta", type2nt[intermediaries].array .== "ACH")

s = maximum.(rfsavg)
s[cseries .== "green"] ./= maximum(s[cseries .== "green"])
s[cseries .!= "green"] ./= maximum(s[cseries .!= "green"])

origin()
Luxor.translate(-corner, corner)
setline(2)
hexannulus(rfsavg, orbitscale = 1.5, text = type2label.(intermediaries), ellipses = ellipsesummary.(rfsavg), ellipsecolors = cseries, maxvals = true)
Luxor.scale(1.5)
Luxor.rotate(-pi/6)
Luxor.translate(-crfsavgellipses[target].center)
Luxor.rotate(pi/6)
setcolor("blue")
fontsize(12)
Luxor.text(target, Point(0, -35), halign = :center)
setline(3)
setdash("dotted")
drawellipse(crfsavgellipses[target])
setdash("solid")
for (i, e) in enumerate(ellipsesummary.(rfsavg))
    setopacity(s[i])
    drawellipse(e, color = cseries[i])
end
Luxor.rotate(-pi/6)
setline(1)
setcolor("black")
Luxor.poly(hexcenter.([HexagonEye(1, 0, hexelsize), HexagonEye(0, 0, hexelsize), HexagonEye(0, 1, hexelsize)]), :stroke)
Luxor.rotate(pi/6)


###### summary
origin()

Luxor.translate(corner, corner)

ellipses = [ellipsesummary.(rfs.array) for rfs in inmaps.(tracetypes.([["Tm1", "T2a", Dm3type] for Dm3type in Dm3types]))]
crfellipses = [ellipsesummary.(rfs.array) for rfs in inmaps.(tracetypes.([["Tm1", Dm3type] for Dm3type in Dm3types]))];
println(length.(ellipses))

#Luxor.translate(0, corner/6)
scale(3)
setline(0.5)
setopacity(1)

cseries = distinguishable_colors(3, [RGB(1,1,1), RGB(0,0,0), colorant"magenta", colorant"green"], dropseed=true)

for (ellipse, crfellipse, ellipsecolor) in zip(ellipses, crfellipses, cseries)
    sethue(ellipsecolor)
    for icell in sample(1:length(ellipse), 50, replace=false)
        if ~ismissing(ellipse[icell]) && ~ismissing(crfellipse[icell])
            Luxor.rotate(-pi/6)
            Luxor.translate(-hexelsize.*crfellipse[icell].center)
            Luxor.rotate(pi/6)
            drawellipse(ellipse[icell], color = ellipsecolor, hexelsize = hexelsize)
            Luxor.rotate(-pi/6)
            Luxor.translate(hexelsize.*crfellipse[icell].center)
            Luxor.rotate(pi/6)
        end
    end
end
sethue("black")
setline(6)
drawpqaxes(hexelsize)


# origin()

# Luxor.translate(corner, corner)

# paths = [["Tm1", "TmY9q⊥", "Dm3v"], ["Tm1", "TmY9q", "Dm3p"], ["Tm1", "TmY9q⊥", "Dm3q"]]
# ellipses = [ellipsesummary.(rfs.array) for rfs in inmaps.(tracetypes.(paths))];

# Luxor.translate(corner/2, -corner/2)
# scale(1.5)
# setline(0.5)
# setopacity(1)

# hexelsize = 7

# for (ellipse, crfellipse, ellipsecolor) in zip(ellipses, crfellipses, ["red", "green", "blue"])
#     sethue(ellipsecolor)
#     for icell in sample(1:length(ellipse), 50, replace=false)
#         if ~ismissing(ellipse[icell]) && ~ismissing(crfellipse[icell])
#             Luxor.rotate(-pi/6)
#             Luxor.translate(-hexelsize.*crfellipse[icell].center)
#             Luxor.rotate(pi/6)
#             drawellipse(ellipse[icell], color = ellipsecolor, hexelsize = hexelsize)
#             Luxor.rotate(-pi/6)
#             Luxor.translate(hexelsize.*crfellipse[icell].center)
#             Luxor.rotate(pi/6)
#         end
#     end
# end
# sethue("black")
# setline(6)
# Luxor.line(hexelsize.*Point(-6, 0), hexelsize.*Point(-6, sqrt(3)), :stroke)

finish()

preview()

# %%
paths = [["Tm1", "T2a", Dm3type] for Dm3type in Dm3types]
ellipses = [ellipsesummary.(rfs.array) for rfs in inmaps.(tracetypes.(paths))];

# %%
crfellipses = [ellipsesummary.(rfs.array) for rfs in inmaps.(tracetypes.([["Tm1", Dm3type] for Dm3type in Dm3types]))];

# %%
cseries = distinguishable_colors(3, [RGB(1,1,1), RGB(0,0,0), colorant"magenta", colorant"green"], dropseed=true)

# %%
#Drawing(600, 600, "Dm3Ellipses.pdf")
Drawing(900, 900)
origin()
scale(5)
setline(0.5)

hexelsize = 7

for (ellipse, crfellipse, ellipsecolor) in zip(ellipses, crfellipses, cseries)
    sethue(ellipsecolor)
    for icell in sample(1:length(ellipse), 50, replace=false)
        if ~ismissing(ellipse[icell]) && ~ismissing(crfellipse[icell])
            Luxor.rotate(-pi/6)
            Luxor.translate(-hexelsize.*crfellipse[icell].center)
            Luxor.rotate(pi/6)
            drawellipse(ellipse[icell], color = ellipsecolor, hexelsize = hexelsize)
            Luxor.rotate(-pi/6)
            Luxor.translate(hexelsize.*crfellipse[icell].center)
            Luxor.rotate(pi/6)
        end
    end
end
sethue("black")
setline(6)
drawpqaxes(hexelsize)
finish()
preview()

# %% [markdown]
# ### Fig 4efg T2a-Dm3 projections

# %%
erfpaths = [["Tm1", "T2a", Dm3type] for Dm3type in Dm3types]
crfpaths = [["Tm1", Dm3type] for Dm3type in Dm3types]

# %%
erfs = inmaps.(tracebacktypes.(erfpaths))
radius = 7
erfscrop = [passmissing(crop).(erf.array, id2pq.(type2ids(Dm3type)), radius, 1) for (erf, Dm3type) in zip(erfs, Dm3types)];

# %%
crfs = inmaps.(tracebacktypes.(crfpaths))
radius = 7
crfscrop = [passmissing(crop).(crf.array, id2pq.(type2ids(Dm3type)), radius, 1) for (crf, Dm3type) in zip(crfs, Dm3types)];

# %%
erfprojs = [hexproject.(skipmissing(ims), ax) for (ims, ax) in zip(erfscrop, [:v, :p, :q])]
crfprojs = [hexproject.(skipmissing(ims), ax) for (ims, ax) in zip(crfscrop, [:v, :p, :q])];

# %%
Plots.reset_defaults()

# %%
default(; # Plots defaults
    fontfamily="Arial Unicode",
    label="", # only explicit legend entries
    guidefontsize = 12,
    tickfontsize = 12,
    legendfontsize = 10,
    linewidth = 2,
#    ylim = (0, Inf), 
    widen = true,
#    xrotation = 30
    )

# %%
h = []
for (a, b, c, d, e) in zip(mean.(erfprojs), std.(erfprojs), mean.(crfprojs), std.(crfprojs), Dm3types)
plot(-radius:radius, 100*a, ribbon = 100*b, 
    xlabel = "displacement",
    xflip = e == "Dm3p",
    xticks = (-6:2:6),
    ylim = (0, 0.3),
    label = "T2a",
        ylabel = e == "Dm3v" ? "input fraction (%)" : "",
    legend = :topleft,
    color = 1,
    annotations = (:topcenter, Plots.text(e, "Helvetica", :center)),
)
push!(h, plot!(twinx(), -radius:radius, 
    100*c,
    ribbon = 100*d, 
    ylim = (0, 10),
    yticks = (0:2:10),
    label = "direct",
    legend = :topright,
    color = 2
))
end

# %%
plot(h..., layout = (1, 3), size = (1000, 250), bottom_margin = 10Measures.mm, left_margin = 5Measures.mm)

# %%
savefig(joinpath(TARGETDIR, "T2aDm3ERFprojections.pdf"))

# %% [markdown]
# ## Fig 5

# %% [markdown]
# ### Fig 5abcd TmY average ERFs

# %%
thres = 0.0025

# %%
monosynaptic = [["Tm1", TmYtype] for TmYtype in TmYtypes]
crfs = inmaps.(tracetypes.(monosynaptic))
crfsellipses = [ellipsesummary.(z.array) for z in crfs]
crfscrop = [passmissing(crop).(ims.array, id2pq.(ids), 5) for (ims, ids) in zip(crfs, type2ids.(TmYtypes))]
crfsavg = mean.(skipmissing.(crfscrop))
crfsavgellipses = NamedArray(ellipsesummary.(crfsavg), TmYtypes)  # only the ellipse is used here

# %%
corner = 260
hexelsize = 7

Drawing(1200, 1200, joinpath(TARGETDIR, "TmYERF.pdf"))
#Drawing(1200, 1200)

target = "TmY4"

intermediaries = intrinsictypes[[scorepath(["Tm1", celltype, target]) for celltype in intrinsictypes] .> thres]
paths = [["Tm1", t, target] for t in intermediaries]

p = sortperm(scorepath.(paths), rev = true)
rfs = inmaps.(tracebacktypes.(paths[p]))

rfscrop = [passmissing(crop).(rf.array, id2pq.(type2ids(target)), 7, 1) for rf in rfs]
rfsavg = [mean(skipmissing(rf)) for rf in rfscrop]
cseries = map(x -> x ? "green" : "magenta", type2nt[intermediaries[p]].array .== "ACH")

s = maximum.(rfsavg)
s[cseries .== "green"] ./= maximum(s[cseries .== "green"])
s[cseries .!= "green"] ./= maximum(s[cseries .!= "green"])

origin()
Luxor.translate(-corner, -corner)
setline(2)
hexannulus(rfsavg, orbitscale = 1.5, text = type2label.(intermediaries[p]), ellipses = ellipsesummary.(rfsavg), ellipsecolors = cseries, maxvals = true)
Luxor.scale(1.5)
Luxor.rotate(-pi/6)
Luxor.translate(-crfsavgellipses[target].center)
Luxor.rotate(pi/6)
setcolor("blue")
fontsize(12)
Luxor.text(type2label(target), Point(0, -35), halign = :center)
setline(3)
setdash("dotted")
drawellipse(crfsavgellipses[target])
setdash("solid")
for (i, e) in enumerate(ellipsesummary.(rfsavg))
    setopacity(s[i])
    drawellipse(e, color = cseries[i])
end
Luxor.rotate(-pi/6)
setline(1)
setcolor("black")
Luxor.poly(hexcenter.([HexagonEye(1, 0, hexelsize), HexagonEye(0, 0, hexelsize), HexagonEye(0, 1, hexelsize)]), :stroke)
Luxor.rotate(pi/6)

####################

target = "TmY9q"
intermediaries = intrinsictypes[[scorepath(["Tm1", celltype, target]) for celltype in intrinsictypes] .> thres]
paths = [["Tm1", t, target] for t in intermediaries]

p = sortperm(scorepath.(paths), rev = true)
rfs = inmaps.(tracebacktypes.(paths[p]))

rfscrop = [passmissing(crop).(rf.array, id2pq.(type2ids(target)), 7, 1) for rf in rfs]
rfsavg = [mean(skipmissing(rf)) for rf in rfscrop]
cseries = map(x -> x ? "green" : "magenta", type2nt[intermediaries[p]].array .== "ACH")

s = maximum.(rfsavg)
s[cseries .== "green"] ./= maximum(s[cseries .== "green"])
s[cseries .!= "green"] ./= maximum(s[cseries .!= "green"])

origin()
Luxor.translate(corner, -corner)
setline(2)
hexannulus(rfsavg, orbitscale = 1.5, text = type2label.(intermediaries[p]), ellipses = ellipsesummary.(rfsavg), ellipsecolors = cseries, maxvals = true)
Luxor.scale(1.5)
Luxor.rotate(-pi/6)
Luxor.translate(-crfsavgellipses[target].center)
Luxor.rotate(pi/6)
setcolor("blue")
fontsize(12)
Luxor.text(type2label(target), Point(0, -35), halign = :center)
setline(3)
setdash("dotted")
drawellipse(crfsavgellipses[target])
setdash("solid")
for (i, e) in enumerate(ellipsesummary.(rfsavg))
    setopacity(s[i])
    drawellipse(e, color = cseries[i])
end
Luxor.rotate(-pi/6)
setline(1)
setcolor("black")
Luxor.poly(hexcenter.([HexagonEye(1, 0, hexelsize), HexagonEye(0, 0, hexelsize), HexagonEye(0, 1, hexelsize)]), :stroke)
Luxor.rotate(pi/6)

################
target = "TmY9q⊥"
intermediaries = intrinsictypes[[scorepath(["Tm1", celltype, target]) for celltype in intrinsictypes] .> thres]
paths = [["Tm1", t, target] for t in intermediaries]

p = sortperm(scorepath.(paths), rev = true)
rfs = inmaps.(tracebacktypes.(paths[p]))

rfscrop = [passmissing(crop).(rf.array, id2pq.(type2ids(target)), 7, 1) for rf in rfs]
rfsavg = [mean(skipmissing(rf)) for rf in rfscrop]
cseries = map(x -> x ? "green" : "magenta", type2nt[intermediaries[p]].array .== "ACH")

s = maximum.(rfsavg)
s[cseries .== "green"] ./= maximum(s[cseries .== "green"])
s[cseries .!= "green"] ./= maximum(s[cseries .!= "green"])

origin()
Luxor.translate(-corner, corner)
setline(2)
hexannulus(rfsavg, orbitscale = 1.5, text = type2label.(intermediaries[p]), ellipses = ellipsesummary.(rfsavg), ellipsecolors = cseries, maxvals = true)
Luxor.scale(1.5)
Luxor.rotate(-pi/6)
Luxor.translate(-crfsavgellipses[target].center)
Luxor.rotate(pi/6)
setcolor("blue")
fontsize(12)
Luxor.text(type2label(target), Point(0, -35), halign = :center)
setline(3)
setdash("dotted")
drawellipse(crfsavgellipses[target])
setdash("solid")
for (i, e) in enumerate(ellipsesummary.(rfsavg))
    setopacity(s[i])
    drawellipse(e, color = cseries[i])
end
Luxor.rotate(-pi/6)
setline(1)
setcolor("black")
Luxor.poly(hexcenter.([HexagonEye(1, 0, hexelsize), HexagonEye(0, 0, hexelsize), HexagonEye(0, 1, hexelsize)]), :stroke)
Luxor.rotate(pi/6)

###### summary
origin()

Luxor.translate(corner, corner)

paths = [["Tm1", "TmY4", "TmY4"], ["Tm1", "TmY9q", "TmY9q"], ["Tm1", "TmY9q⊥", "TmY9q⊥"]]
ellipses = [ellipsesummary.(rfs.array) for rfs in inmaps.(tracetypes.(paths))];

crfellipses = [ellipsesummary.(rfs.array) for rfs in inmaps.(tracetypes.([["Tm1", TmYtype] for TmYtype in TmYtypes]))];

#Luxor.translate(-corner/2, -corner/2)
scale(3)
setline(0.5)
setopacity(1)

cseries = distinguishable_colors(3, [RGB(1,1,1), RGB(0,0,0), colorant"magenta", colorant"green"], dropseed=true)

hexelsize = 7

for (ellipse, crfellipse, ellipsecolor) in zip(ellipses, crfellipses, cseries)
    sethue(ellipsecolor)
    for icell in sample(1:length(ellipse), 50, replace=false)
        if ~ismissing(ellipse[icell]) && ~ismissing(crfellipse[icell])
            Luxor.rotate(-pi/6)
            Luxor.translate(-hexelsize.*crfellipse[icell].center)
            Luxor.rotate(pi/6)
            drawellipse(ellipse[icell], color = ellipsecolor, hexelsize = hexelsize)
            Luxor.rotate(-pi/6)
            Luxor.translate(hexelsize.*crfellipse[icell].center)
            Luxor.rotate(pi/6)
        end
    end
end
sethue("black")
setline(6)
drawpqaxes(hexelsize)

finish()
preview()

# %% [markdown]
# ### Fig 5efg TmY-TmY projections

# %%
erfpaths = [["Tm1", TmYtype, TmYtype] for TmYtype in TmYtypes]
crfpaths = [["Tm1", TmYtype] for TmYtype in TmYtypes]

# %%
radius = 6
erfs = inmaps.(tracebacktypes.(erfpaths))
erfscrop = [passmissing(crop).(erf.array, id2pq.(type2ids(TmYtype)), radius, 1) for (erf, TmYtype) in zip(erfs, TmYtypes)];

# %%
crfs = inmaps.(tracebacktypes.(crfpaths))
crfscrop = [passmissing(crop).(crf.array, id2pq.(type2ids(TmYtype)), radius, 1) for (crf, TmYtype) in zip(crfs, TmYtypes)];

# %%
erfprojs = [hexproject.(skipmissing(ims), ax) for (ims, ax) in zip(erfscrop, [:v, :p, :q])]
crfprojs = [hexproject.(skipmissing(ims), ax) for (ims, ax) in zip(crfscrop, [:v, :p, :q])];

# %%
Plots.reset_defaults()

# %%
default(; # Plots defaults
    fontfamily="Helvetica",
    label="", # only explicit legend entries
    guidefontsize = 12,
    tickfontsize = 12,
    legendfontsize = 10,
    linewidth = 2,
#    ylim = (0, Inf), 
    widen = true,
#    xrotation = 30
    )

# %%
h = []
for (a, b, c, d, e) in zip(mean.(erfprojs), std.(erfprojs), mean.(crfprojs), std.(crfprojs), TmYtypes)
plot(-radius:radius, 100*a, ribbon = 100*b, 
    xlabel = "displacement",
    xticks = (-6:2:6),
    yticks = (0:0.1:0.2),
    ylim = (0, 0.2),
    label = "TmY",
    ylabel = e == "TmY4" ? "input fraction (%)" : "",
    legend = :topleft,
    annotations = (:topcenter, Plots.text(e, "Helvetica", :center)),
        color = 1
)
push!(h, 
        plot!(twinx(), -radius:radius, 100*c, ribbon = 100*d, 
        ylim = (0, 8),
        yticks = (0:2:8),
        label = "direct",
        legend = :topright,
        color = 2
))
end

# %%
plot(h..., layout = (1, 3), size = (1000, 250), bottom_margin = 10Measures.mm, left_margin = 5Measures.mm)

# %%
savefig(joinpath(TARGETDIR, "TmYTmYERFprojections.pdf"))

# %% [markdown]
# ## Fig 6 LC15

# %% [markdown]
# ### LC15 locations

# %%
# anchor LC15 cells on Mi1→T3→LC15 maps
# missing means that there is zero connectivity
rfs = prepreimage("Mi1", "T3", "LC15")
for (i, id) in enumerate(type2ids("LC15"))
    center = findcenter(float(rfs[i]), seven)
    if ~ismissing(center)
        id2pq[id] = center
    end
end

# %% [markdown]
# ### Fig 6b LC15 average disynaptic maps

# %%
scores = [scorepath([c1, c2, "LC15"]) for c1 in hexel, c2 in nothexel]
scores = NamedArray(scores, names = (hexel, nothexel))
scoresums = NamedArray(sum(scores.array, dims = 1)[:], nothexel)
scoremaxs = NamedArray(maximum(scores.array, dims = 1)[:], nothexel)

# no need to exclude inhibitory intermediaries here, as top six are predicted excitatory
# contrast with LC10ev
#thres = 0.01
#intermediaries = nothexel[scoresums.array .> thres]
intermediaries = nothexel[sortperm(scoresums.array, rev=true)[1:6]]   # top six

p = sortperm(scoresums[intermediaries], rev = true)

scoresums[intermediaries[p]]

# %%
# source hexel types
ignore, winners = findmax(scores[:, intermediaries], dims=1)
sources = hexel[getindex.(winners, 1)][:]
sources[p]

# %%
paths = [[source, intermediary, "LC15"] for (source, intermediary) in zip(sources[p], intermediaries[p])]
rfs = inmaps.(tracebacktypes.(paths));

# %%
cseries = distinguishable_colors(6, ColorSchemes.hot[:], dropseed=true)

# %%
rfscrop = [passmissing(crop).(rf.array, id2pq.(type2ids("LC15")), 8, 1) for rf in rfs]
rfsavg = [mean(skipmissing(rf)) for rf in rfscrop]

#Drawing(700, 600)
Drawing(700, 600, joinpath(TARGETDIR, "LC15.pdf"))
    origin()
    ellipses = ellipsesummary.(rfsavg)
    hexannulus(rfsavg, orbitscale = 1.5, text = type2label.(intermediaries[p]), ellipses = ellipses, ellipsecolors = cseries, maxvals = true)
    Luxor.scale(1.5)
    drawpqaxes(7)
    Luxor.rotate(-pi/6)
    Luxor.translate(-ellipses[1].center)   # place T3 ellipse at center (doesn't matter so much for this plot)
    Luxor.rotate(pi/6)
    setline(3)
    for (i, e) in enumerate(ellipses)
        if in(i, [3 4 5])       # make TmY ellipses more visible
            setopacity(1)
        else
            setopacity(0.3)
        end
        drawellipse(e, color = cseries[i])
    end
finish()
preview()

# %% [markdown]
# ### Fig 6c LC15 ellipses

# %%
ellipses = [ellipsesummary.(rf.array) for rf in rfs]  # for individual cells, rfs computed above
hexelsize = 7

Drawing(500, 600, joinpath(TARGETDIR, "Tm1LC15Ellipses.pdf"))
#Drawing(500, 600)
origin()
scale(3)
setline(1)
for (e1, e2, e3, e4) in zip(ellipses[1], ellipses[3], ellipses[4], ellipses[5])
    Luxor.rotate(-pi/6)
    Luxor.translate(-hexelsize * e1.center)  # anchor on T3 ellipse, which isn't shown
    Luxor.rotate(pi/6)
    drawellipse(e2, color = cseries[3], hexelsize = hexelsize)
    drawellipse(e3, color = cseries[4], hexelsize = hexelsize)
    drawellipse(e4, color = cseries[5], hexelsize = hexelsize)
    Luxor.rotate(-pi/6)
    Luxor.translate(hexelsize*e1.center)
    Luxor.rotate(pi/6)
end
setline(5)
setcolor("black")
drawpqaxes(hexelsize)

finish()
preview()

# %% [markdown]
# ## Fig 6 LC10e ventral

# %%
### new clustering
LC10ed = [720575940603823136, 720575940604935089, 720575940610075029, 720575940614496735, 720575940614746178, 720575940616179050, 720575940618360741, 720575940620416176, 720575940621175979, 720575940622744414, 720575940622805854, 720575940623979815, 720575940624787504, 720575940625143396, 720575940625981080, 720575940628023296, 720575940628359783, 720575940628615752, 720575940629385943, 720575940629398492, 720575940629958764, 720575940631576546, 720575940633027615, 720575940637401321, 720575940637706574, 720575940637731955, 720575940638701375, 720575940644393070]
LC10ev = [720575940604125920, 720575940606274592, 720575940608977988, 720575940609733387, 720575940609890402, 720575940612305250, 720575940612332401, 720575940614897662, 720575940621233157, 720575940621602047, 720575940622486057, 720575940622602992, 720575940623601193, 720575940624050407, 720575940624595893, 720575940628632064, 720575940629076042, 720575940629756556, 720575940630438651, 720575940631305467, 720575940632034385, 720575940634088602, 720575940640974555]

# %%
issetequal(type2ids("LC10e"), union(LC10ed, LC10ev))

# %% [markdown]
# ### ventral locations in pq coordinates

# %%
# anchor LC10 cells on Tm1→TmY9q→LC10 maps
# missing means that there is zero connectivity
rfs = prepreimage("Tm1", "TmY9q", LC10ev);
# other possible anchors
#    rfs = prepreimage("Tm1", "TmY9q⊥", LC10ev);
#    rfs = prepreimage("Mi9", "Tm8a", LC10ev);
@mfalse for (i, id) in enumerate(LC10ev)
    center = findcenter(float(rfs[i]), seven)
    if !ismissing(center)
        id2pq[id] = center
    end
end

# %% [markdown]
# ### Fig 6e LC10ev average disynaptic maps

# %%
infractionLC10ev = sum(A'*W[:, Name.(LC10ev)], dims = 2)/sum(W[:, Name.(LC10ev)])
infractionLC10ev = NamedArray(infractionLC10ev[:], alltypes)

thres = 0.003
scores = [infraction[c1, c2]*infractionLC10ev[c2] for c1 in hexel, c2 in nothexel]
scores = NamedArray(scores, names = (hexel, nothexel))
scoresums = NamedArray(sum(scores.array, dims = 1)[:], nothexel)
scoremaxs = NamedArray(maximum(scores.array, dims = 1)[:], nothexel)

intermediaries = nothexel[scoresums.array .> thres .&& type2nt[nothexel].array .== "ACH"]
#intermediaries = nothexel[scoresums.array .> thres]
p = sortperm(scoresums[intermediaries], rev = true)

scoresums[intermediaries[p]]

# %%
ignore, winners = findmax(scores[:, intermediaries], dims=1)
sources = hexel[getindex.(winners, 1)][:]
sources[p]

# %%
paths = Vector{Union{String, Vector{<:Integer}}}[[source, intermediary, LC10ev] for (source, intermediary) in zip(sources[p], intermediaries[p])]
rfs = inmaps.(tracebacktypes.(paths))
pathtexts = [[source, intermediary, "LC10ev"] for (source, intermediary) in zip(sources[p], intermediaries[p])]

# %%
cseries = distinguishable_colors(6, ColorSchemes.hot[:], dropseed=true)

# %%
rfscrop = [passmissing(crop).(rf.array, id2pq.(LC10ev), 8, 1) for rf in rfs]
rfsavg = [mean(skipmissing(rf)) for rf in rfscrop]

Drawing(700, 600, joinpath(TARGETDIR, "LC10ev.pdf"))
#Drawing(700, 600)
    origin()
    ellipses = ellipsesummary.(rfsavg)
    hexannulus(rfsavg, orbitscale = 1.5, text = type2label.(intermediaries[p]), ellipses = ellipses, ellipsecolors = cseries, maxvals = true)
    Luxor.scale(1.5)
    drawpqaxes(7)
    Luxor.rotate(-pi/6)
    Luxor.translate(-ellipses[1].center)
    Luxor.rotate(pi/6)
    setline(3)
    for (i, e) in enumerate(ellipses)
        if i in [1, 2]
            setopacity(1)
        else
            setopacity(0.3)
        end
        drawellipse(e, color = cseries[i])
    end
finish()
preview()

# %% [markdown]
# ### Fig 6f ventral ellipses

# %%
# overload function (hexgraphics.jl) to handle LC10ev as list of IDs, rather than a string indicating a type
function OpticLobe.tracebacktypes(celltypes::Vector{Any})
    @mfalse inds = [isa(celltype, String) ? ind2type .== celltype : id2ind.(celltype) for celltype in celltypes]
    return *([Wback[inds[i], inds[i+1]] for i = 1:length(inds)-1]...)
end

# %%
hexelsize = 7

#paths = [["Mi9", "Tm8a", LC10ev], ["Tm1", "TmY9q", LC10ev], ["Tm1", "TmY9q⊥", LC10ev]]
paths = Vector{Union{String, Vector{<:Integer}}}[["Tm1", "TmY9q", LC10ev], ["Tm1", "TmY9q⊥", LC10ev]]

rfs = inmaps.(tracebacktypes.(paths))
ellipses = [ellipsesummary.(rf.array) for rf in rfs]

Drawing(600, 800, joinpath(TARGETDIR, "Tm1LC10evEllipses.pdf"))
#Drawing(600, 800)
origin()
scale(3)
setline(1)
for (e1, e2) in zip(ellipses[1], ellipses[2])
    Luxor.rotate(-pi/6)
    Luxor.translate(-hexelsize * e1.center)
    Luxor.rotate(pi/6)
    drawellipse(e1, color = cseries[1], hexelsize = hexelsize)
    drawellipse(e2, color = cseries[2], hexelsize = hexelsize)
    Luxor.rotate(-pi/6)
    Luxor.translate(hexelsize*e1.center)
    Luxor.rotate(pi/6)
end

setline(5)
setcolor("black")
drawpqaxes(hexelsize)
finish()
preview()

# %% [markdown]
# ### Fig 6f inset displacement vectors

# %%
rot30 = [cos(pi/6) sin(pi/6); -sin(pi/6) cos(pi/6)]

# %%
rfs = [prepreimage("Mi9", "Tm8a", LC10ev), prepreimage("Tm1", "TmY9q", LC10ev), prepreimage("Tm1", "TmY9q⊥", LC10ev)]
ellipses = [ellipsesummary.(ims) for ims in rfs]

d = [collect(e.center)/sqrt(3) for e in ellipses[3]] - [collect(e.center)/sqrt(3) for e in ellipses[2]]

# %%
hexelsize = 50
Drawing(600, 600, joinpath(TARGETDIR, "LC10eDisplacements.pdf"))
#Drawing(600, 600)
origin()
setcolor("gray")
    for x in d
        Luxor.arrow(Point(0, 0), Point(hexelsize*rot30*x...), arrowheadlength=20, linewidth = 5)
    end
Luxor.rotate(-pi/6)
setline(10)
setcolor("black")
Luxor.poly(hexcenter.([HexagonEye(1, 0, hexelsize), HexagonEye(0, 0, hexelsize), HexagonEye(0, 1, hexelsize)]), :stroke)
Luxor.rotate(pi/6)
finish()
preview()

# %% [markdown] toc-hr-collapsed=true
# ## Fig S3 Tm1Mi4Tm2L3-Dm3

# %%
using Format

# %%
function valuesandtext(values, textstring, fmtstring = "%3d")
    setcolor("black")
    fontsize(20)
    fontface("Arial Unicode MS")
    Luxor.text(cfmt(fmtstring, values[1]), Point(-170, 115))
    Luxor.text(cfmt(fmtstring, values[2]), Point(-190, 75))
    Luxor.text(cfmt(fmtstring, values[3]), Point(-150, 75))
    Luxor.text(textstring, Point(-150, 155), halign = :center)
end    

# %%
Drawing(900, 800, joinpath(TARGETDIR, "Mi4Tm2L3Tm9-Dm3.pdf"))
#Drawing(900, 800)
origin()
Luxor.translate(-200, -200)
rfs = typetriad("Mi4", Dm3types, radius = 3, hexelsize = 15, normalize = :separate)
valuesandtext(maximum.(rfs), "Mi4→Dm3", "%3.1f")

origin()
Luxor.translate(200, -200)
rfs = typetriad("Tm2", Dm3types, radius = 3, hexelsize = 15, normalize = :separate)
valuesandtext(maximum.(rfs), "Tm2→Dm3", "%3.1f")

origin()
Luxor.translate(-200, 200)
rfs = typetriad("L3", Dm3types, radius = 3, hexelsize = 15, normalize = :separate)
valuesandtext(maximum.(rfs), "L3→Dm3", "%3.1f")

origin()
Luxor.translate(200, 200)
rfs = typetriad("Tm9", Dm3types, radius = 3, hexelsize = 15, normalize = :separate)
valuesandtext(maximum.(rfs), "Tm9→Dm3", "%3.1f")

finish()
preview()

# %% [markdown] toc-hr-collapsed=true
# ## Fig S6 lateral interactions Dm3-TmY circuit

# %%
function maparray(pretypes, posttypes; radius = 6, spacing = 140)
    rfscrop = [passmissing(crop).(preimage(pretype, posttype), id2pq.(type2ids(posttype)), radius, 1) for pretype in pretypes, posttype in posttypes]
    rfscropavg = sum.(skipmissing.(rfscrop))./length.(rfscrop)
    immax = maximum(hcat(rfscropavg[:]...))
        
    for pre = 1:length(pretypes)
        for post = 1:length(posttypes)
            x, y = (post-(length(posttypes)+1)/2)*spacing, (pre-(length(pretypes)+1)/2)*spacing
            Luxor.translate(x, y)
            rect2hex(square2hex(get(ColorSchemes.hot, rfscropavg[pre, post]/immax)))
            Luxor.translate(-x, -y)
        end
    end
    return immax
end

# %%
Drawing(1000, 1000, joinpath(TARGETDIR, "Dm3TmYAllAverages.pdf"))
spacing = 250
#Drawing(1000, 1000)
origin()

Luxor.translate(-spacing, -spacing)
maxvalue = maparray(Dm3types, Dm3types)
fontsize(16)
Luxor.text(@sprintf("maxvalue = %3.1f", maxvalue), Point(0, -220), halign = :center)

origin()
Luxor.translate(spacing, -spacing)
maxvalue = maparray(Dm3types, TmYtypes)
fontsize(16)
Luxor.text(@sprintf("maxvalue = %3.1f", maxvalue), Point(0, -220), halign = :center)

origin()
Luxor.translate(-spacing, spacing)
maxvalue = maparray(TmYtypes, Dm3types)
fontsize(16)
Luxor.text(@sprintf("maxvalue = %3.1f", maxvalue), Point(0, -220), halign = :center)

origin()
Luxor.translate(spacing, spacing)
maxvalue = maparray(TmYtypes, TmYtypes)
fontsize(16)
Luxor.text(@sprintf("maxvalue = %3.1f", maxvalue), Point(0, -220), halign = :center)

finish()
preview()

# %% [markdown]
# ## Data S3 Dm3 individual cells
# TODO: update to `distinguishable_colors`
#
# function generate_target_montages(target::String, n_monosynaptic::Int=5; extra_paths::Vector=[])
#     # Create target directory
#     target_dir = joinpath(TARGETDIR, target)
#     if !isdir(target_dir)
#         mkpath(target_dir)
#     end
#     
#     # Find monosynaptic connections
#     monosynaptic = names(first(sort(infraction[hexel, target], rev=true), n_monosynaptic), 1)
#     
#     # Calculate disynaptic scores
#     scores = [infraction[c1, c2]*infraction[c2, target] for c1 in hexel, c2 in nothexel]
#     scores = NamedArray(scores, names = (hexel, nothexel))
#     ignore, winners = findmax(scores, dims=1)
#     sources = hexel[getindex.(winners, 1)][:]
#     
#     scoresums = NamedArray(sum(scores.array, dims = 1)[:], nothexel)
#     p = sortperm(scoresums, rev = true)
#     disynaptic = [sources[p] nothexel[p]][1:10, :]
#     
#     # Create paths
#     paths = vcat(
#         [[t, target] for t in monosynaptic],
#         [[t..., target] for t in eachrow(disynaptic)],
#         extra_paths
#     )
#
#     # Compute receptive fields
#     rfs = inmaps.(tracebacktypes.(paths))
#     
#     # Set up colors
#     tocolor = length.(paths) .>= 3
#     tocolor[1] = true
#     ellipsecolors = Vector{Any}(nothing, length(paths))
#     ellipsecolors[tocolor] .= distinguishable_colors(sum(tocolor))
#     
#     # Generate montages
#     @showprogress for id in type2ids(target)
#         ims = [rf[Name(id)] for rf in rfs]
#         montage(ims, fname = joinpath(target_dir, "$id.pdf"), 
#                labels = path2label.(paths), hexelsize = 6, ellipses = true, 
#                ellipsecolors = ellipsecolors, maxvals = true, 
#                centers = repeat([id2pq[id]], length(ims)), summary = 1)
#     end
# end

# %% [markdown]
# ### Dm3p

# %%
generate_target_montages("Dm3p", 5)

# %% [markdown]
# ### Dm3q

# %%
generate_target_montages("Dm3q", 5)

# %% [markdown]
# ### Dm3v

# %%
generate_target_montages("Dm3v", 5)

# %% [markdown]
# ## Data S4 TmY individual cells (supplementary data)

# %% [markdown]
# ### TmY4

# %%
generate_target_montages("TmY4", 4, extra_paths=[["Tm1", "TmY4", "Dm3v", "TmY4"]])

# %% [markdown]
# ### TmY9q

# %%
generate_target_montages("TmY9q", 4, extra_paths=[["Tm1", "TmY9q", "Dm3p", "TmY9q"]])

# %% [markdown]
# ### TmY9q⊥

# %%
generate_target_montages("TmY9q⊥", 4, extra_paths=[["Tm1", "TmY9q⊥", "Dm3q", "TmY9q⊥"]])

# %% [markdown]
# ## Data S5 LC10ev individual cells

# %%
infractionLC10ev = sum(A'*W[:, Name.(LC10ev)], dims = 2)/sum(W[:, Name.(LC10ev)])
infractionLC10ev = NamedArray(infractionLC10ev[:].array, alltypes)

thres = 0.001
scores = [infraction[c1, c2]*infractionLC10ev[c2] for c1 in hexel, c2 in nothexel]
scores = NamedArray(scores, names = (hexel, nothexel))
scoresums = NamedArray(sum(scores.array, dims = 1)[:], nothexel)
scoremaxs = NamedArray(maximum(scores.array, dims = 1)[:], nothexel)

intermediaries = nothexel[scoresums.array .> thres .&& type2nt[nothexel].array .== "ACH"]
p = sortperm(scoresums[intermediaries], rev = true)

scoresums[intermediaries[p]]

# %%
ignore, winners = findmax(scores[:, intermediaries[p]], dims=1)
sources = hexel[getindex.(winners, 1)][:]

# %%
paths = Vector{Union{String, Vector{<:Integer}}}[[source, intermediary, LC10ev] for (source, intermediary) in zip(sources, intermediaries[p])]

rfs = inmaps.(tracebacktypes.(paths))

pathtexts = [[source, intermediary, "LC10ev"] for (source, intermediary) in zip(sources, intermediaries[p])]

# %%
cseries = distinguishable_colors(length(paths), ColorSchemes.hot[:], dropseed=true)

# %%
hexelsize = 6

# Create LC10ev subdirectory if it doesn't exist
lc10ev_dir = joinpath(TARGETDIR, "LC10ev")
if !isdir(lc10ev_dir)
    mkpath(lc10ev_dir)
end

@showprogress for id in LC10ev
    ims = [rf[Name(id)] for rf in rfs]
    montage(ims, fname = joinpath(lc10ev_dir, "$id.pdf"), 
    labels = path2label.(pathtexts), hexelsize = hexelsize, ellipses = true, ellipsecolors = cseries, 
    maxvals = true, centers = repeat([id2pq[id]], length(ims)), summary = 1)
end

# %% [markdown]
# ## Data S5 LC15 individual cells

# %%
#thres = 0.0025 
thres = 0.005
target = "LC15"
scores = [scorepath([c1, c2, target]) for c1 in hexel, c2 in nothexel]
scores = NamedArray(scores, names = (hexel, nothexel))
scoresums = NamedArray(sum(scores.array, dims = 1)[:], nothexel)
scoremaxs = NamedArray(maximum(scores.array, dims = 1)[:], nothexel)

# %%
intermediaries = nothexel[scoresums.array .> thres .&& type2nt[nothexel].array .== "ACH"]
p = sortperm(scoresums[intermediaries], rev = true)

# %%
ignore, winners = findmax(scores[:, intermediaries], dims=1)
sources = hexel[getindex.(winners, 1)][:]

# %%
sources[p]

# %%
paths = [[source, intermediary, "LC15"] for (source, intermediary) in zip(sources[p], intermediaries[p])]

rfs = inmaps.(tracebacktypes.(paths))

# %%
cseries = distinguishable_colors(length(paths), ColorSchemes.hot[:], dropseed=true)

# %%
hexelsize = 6

# Create LC15 subdirectory if it doesn't exist
lc15_dir = joinpath(TARGETDIR, "LC15")
if !isdir(lc15_dir)
    mkpath(lc15_dir)
end

@showprogress for id in type2ids("LC15")
    ims = [rf[Name(id)] for rf in rfs]
    montage(ims, fname = joinpath(lc15_dir, "$id.pdf"), 
    labels = path2label.(paths), hexelsize = hexelsize, ellipses = true, ellipsecolors = cseries, 
    maxvals = true, centers = repeat([id2pq[id]], length(ims)), summary = 1)
end

# %%
