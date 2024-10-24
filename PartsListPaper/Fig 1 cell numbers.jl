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
# # Figure 1. Cell numbers

# %%
using OpticLobe

# %% jupyter={"outputs_hidden": false}
using NamedArrays
using Plots, Measures
using OrderedCollections

# %%
const TARGETDIR = "Fig 1"
# const TARGETDIR = "~/sseung@princeton.edu/OpticLobeCellTypesPaper/panels"

# %% [markdown]
# ## right panel

# %% [markdown]
# ### right top: families types and cells per class

# %%
classes = collect(keys(class2families))
class2types = OrderedDict(keys(class2families) .=> [vcat([family2types[f] for f in v]...) for v in values(class2families)])

# %%
cmap = palette(:default)
default(
    fontfamily="Helvetica",
    label = "",
    )

# %%
numfamiliesperclass = length.(values(class2families))

familiesperclassbar = bar(
        numfamiliesperclass,
        xticks = strings2ticks(classes),
        yticks = [0, 5, 10],
        ylabel = "# families",
        permute = (:y, :x),
        xflip = true,
        xlim = (0.5, length(classes) + 0.5),
        left_margin = 5mm,
        c = cmap[1:5]
    )

# %%
numtypesperclass = length.(values(class2types))
numtypesperclass[1] = 8  # hack for R1-6
typesperclassbar = bar(
        numtypesperclass,
        ylabel = "# types",
        permute = (:y, :x),
        xlim = (0.5, length(classes) + 0.5),
        yticks = [0, 50, 100],
        widen = true,
        xflip = true,
        c = cmap[1:5]
    )
plot!(typesperclassbar, yformatter = :none)

# %%
cellnumbers = NamedArray(sum(Ai.array, dims=1)[:], intrinsictypes)
numcellsperclass = [sum(cellnumbers[typelist]) for typelist in values(class2types)]
cellsperclassbar = bar(
        numcellsperclass/10000,
        ylabel = "# cells/10000",
        permute = (:y, :x),
        xlim = (0.5, length(classes) + 0.5),
        ylim = (-Inf, 3.2),
        widen = true,
        xflip = true,
        c = cmap[1:5]
    )
plot!(cellsperclassbar, yformatter = :none)

# %%
plot(familiesperclassbar, typesperclassbar, cellsperclassbar, 
    layout = (1, 3),
    size = (600, 200),
    tickfontsize = 14,
    bottom_margin = 7mm,
    guidefontsize = 16,)

# %%
savefig(joinpath(TARGETDIR, "FamiliesTypesCellsPerClass.pdf"))

# %% [markdown]
# ### right middle: number of cells and types per family

# %%
cmap = palette(:default)
default(
    fontfamily="Helvetica",
    label = ""
    )

# %%
families = collect(keys(family2types))
familycardinalities = length.(values(family2types))
familycardinalities[1] = 6  # hack for R1-6
classcardinalities = length.(values(class2families))
familycolors = cmap[vcat([fill(i, classcardinalities[i]) for i in 1:length(class2families)]...)]
b1 = bar(
        familycardinalities,
        xticks = strings2ticks(families),
        ylabel = "# cell types",
        xlabel = "type family",
        permute = (:y, :x),
        xlim = (0.5, length(families) + 0.5),
        xflip = true,
        color = familycolors,
    )

# %%
b2 = bar(
    sum.([cellnumbers[family2types[f]] for f in families])/1000,
    xticks = (1:length(families)),
    ylabel = "# cells/1000",
    permute = (:y, :x),
    xlim = (0.5, length(families) + 0.5),
    ylim = (-Inf, 9),
    xflip = true,
    color = familycolors,
    )
plot!(b2, yformatter = :none)

# %%
# standalone version with family names and legend
temp = -ones(1, 5)
plot(b2, temp, color=cmap[1:5]', 
    linewidth = 3, 
    legend = :bottomright, 
    yticks = strings2ticks(families),
    ) 

# %%
plot(b1, b2,
    tickfontsize = 9,
    guidefontsize = 10,
)

# %%
savefig(joinpath(TARGETDIR, "CellTypesAndCellsPerFamily.pdf"))

# %% [markdown]
# ### right bottom: histogram of cells per type

# %%
cellnumbers = NamedArray(sum(Ai.array, dims=1)[:], intrinsictypes)
#cellnumbers["R1-6"] = round(cellnumbers["R1-6"]/6)
cellnumbers = cellnumbers[intrinsictypes .!= "R1-6"]   # leave out, because it has 6x the number of one type

# %%
#default()
default(
    fontfamily="Helvetica",
    label = ""
    )

# %%
cellspertypehist = histogram(
    log10.(cellnumbers), 
    nbins = 50, 
    xticks = (0:3, [1, 10, 100, 1000]),
    yticks = [0, 10, 20],
    xlabel = "# cells",
    ylabel = "# types",
    legend = false,
)

# %%
savefig(joinpath(TARGETDIR, "DistributionLogCellNumberPerType.pdf"))

# %% [markdown]
# ### combined right panel

# %%
hright = plot( 
    familiesperclassbar, typesperclassbar, cellsperclassbar, 
    b1, b2, 
    cellspertypehist, 
    layout = @layout([
            grid(1,3){0.15h}
            grid(1,2)
            c{0.15h}
        ]), 
    tickfontsize = 9,
    guidefontsize = 11,
)

# %%
# this needs some padding adjustments to show up properly
plot(hright, size = (400, 800))
savefig(joinpath(TARGETDIR, "ClassFamilyTypeCell.pdf"))

# %% [markdown]
# ## left panel: number of cells per type

# %%
cellnumbers = NamedArray(sum(Ai.array, dims=1)[:], intrinsictypes)
#cellnumbers["R1-6"] = round(cellnumbers["R1-6"]/6)
cellnumbers = cellnumbers[intrinsictypes .!= "R1-6"]   # leave out, because it has 6x the number of one type

# %%
first(sort(cellnumbers, rev=true), 30) |> showall

# %%
default(; # Plots defaults
    fontfamily="Helvetica",
    label="" # only explicit legend entries
    )

# %%
# five column layout
46*[1 2 3 4 5]

# %%
# four column layout
57*[1 2 3 4]

# %%
perm = sortperm(cellnumbers)
noptictypes = length(cellnumbers)

r = [1:57, 58:114, 115:171, 172:noptictypes]  # four column layout
#r = [1:46, 47:92, 93:138, 139:184, 185:noptictypes]   # five column layout
ncol = length(r)
    
columns = [plot(
        cellnumbers[perm[r[i]]].array, 
        permute = (:y, :x),
        xlim = (0.5, length(perm[r[i]])+0.5),
        xticks=((1:length(perm[r[i]])), names(cellnumbers)[1][perm[r[i]]]),
        ylabel = "# cells",
        seriestype = :stepmid,
        yrotation = 45,
    ) for i = ncol:-1:1]

hleft = plot(columns..., 
#    size = (200*ncol, 800*4/ncol), 
    layout = (1, ncol),
    ytickfontsize = 9,
    xtickfontsize = 8,
    guidefontsize = 11,
)

# %%
savefig(joinpath(TARGETDIR, "CellNumberVsTypeFourColumns.pdf"))
#savefig(joinpath(TARGETDIR, "CellNumberVsTypeFiveColumns.svg"))

# %% [markdown]
# ## left and right panels combined

# %%
plot(hleft, hright, layout = @layout([a{0.66w} b{0.34w}]), size = (1200, 800), bottom_margin = 3mm)

# %%
savefig(joinpath(TARGETDIR, "Fig 1 bottom.pdf"))
