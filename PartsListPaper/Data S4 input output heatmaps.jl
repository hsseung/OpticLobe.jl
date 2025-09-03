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
# # Data S4. in out fraction heatmaps
#
# note that unicode may not show up correctly in Apple Quick Look, but it renders OK in browser

# %%
using OpticLobe

include("config.jl")

# %% jupyter={"outputs_hidden": false}
# using NamedArrays
using Plots
using PlotlyKaleido   # use because Plots.savefig is broken for PlotlyJS, not paying attention to size attribute
using Measures
using ProgressMeter

# %%
plotlyjs()    # use instead of GR, because GR produces too much whitespace

# %%
PlotlyKaleido.start()

# %% [markdown]
# ## visualize infraction and outfraction for types of each family

# %%
"""
find partners with infraction or outfraction exceeding threshold for any type in `typestrings`
"""
function findpartners(typestrings, thres)
    preind = getindex.(findall(infraction[:, typestrings] .> thres), 1)
    postind = getindex.(findall(outfraction[typestrings, :] .> thres), 2)
    return sort(unique(vcat(preind, postind)))
end

# %%
function vizclasspartners(family::String, thres)
    typestrings = family2types[family]
    winners = findpartners(typestrings, thres)
    return inoutheat(typestrings, winners)
end

# %%
function vizclasspartners(families::Vector{String}, thres)
    typestrings = vcat([family2types[f] for f in families]...)
    winners = findpartners(typestrings, thres)
    return inoutheat(typestrings, winners)
end

# %%
# default colormap is `cgrad(:inferno)`
# Note that the `scale` option doesn't work in `cgrad` in `plotly`. It worked in `GR`, but I'm not sure what it was doing anyway.
# Now a `tanh` saturating nonlinearity is applied manually
function inoutheat(typestrings, winners)
    pixelsize = 20 # heatmap pixel size in screen pixels
    textwidth = 15*maximum(length.(alltypes[winners])) + 40 # estimated width of text in screen pixels + white space
    height = pixelsize*length(winners)
    width = 2*pixelsize*length(typestrings) + textwidth
    width = width*10/9   # since heatmaps take up 90% of width in layout
    f(x) = tanh(x/0.1)  # compressive nonlinearity applied before heatmap to make small fractions more visible
    nbar = 100
    cmax = 0.2 # maximum value on colorbar
    colorbar = reshape(collect(1:nbar)*cmax/nbar, nbar, 1)
    default(
        fontfamily = "Helvetica",
        label = "",
        xrotation = 60,
        cbar = false, 
        yflip = true,
        xticks = strings2ticks(typestrings),
        clim = (0, 1),
        tickfontsize = 12,
        guidefontsize = 10,
        titlefontsize = 14,
    )
    hh = plot(
        heatmap(
            f.(infraction[winners, typestrings]),                       
            yticks = strings2ticks(alltypes[winners]),
            title = "in fraction",
        ),
        heatmap(
            f.(outfraction[typestrings, winners]'),
            yticks = :none,
            title = "out fraction",
        ),
        heatmap(f.(colorbar),
            xaxis = false,
            yflip = false,
            xlim = (0.5, 1.5),
            yticks = (0:nbar/4:nbar, (0:cmax/4:cmax)),
            ylim = (0, nbar),
            left_margin = 5mm,
            aspect_ratio = 0.1,
            tick_dir = :out
            ),
        layout = @layout([a{0.45w} b{0.45w} [_; c{0.1h}; _]]),
        size = (width, height)
    )
    return hh
end

# %%
# vizclasspartners( ["Lat", "LLPt", "MLt", "LMt", "PDt", "LMa", "MLLPa"], 0.02)

# %%
cgrad(:inferno)

# %%
# TARGETDIR = "~/sseung@princeton.edu/OpticLobeCellTypesPaper/DataS4"
TARGETDIR = "Data S4"

# %%
@showprogress for family in [
    ["L", "C", "Lawf", "Lai", "T1", "T2", "T3"],
    "Mi",
    "Tm",
    "TmY",
    ["T4", "T5"],
    ["Y", "Tlp"],
    ["Lat", "LLPt", "MLt", "LMt", "PDt", "LMa", "MLLPa"],
    "Dm", "Pm", "Sm", "Li", "LPi",
    ]
    hh = vizclasspartners(family, 0.02)
    width, height = hh.attr[:size]
    display(hh) # populate the hh.o field, which is a SyncPlot object
    PlotlyKaleido.savefig(joinpath(TARGETDIR, "$(join(family)).pdf"), hh.o, height = height, width = width)
end

# %%
PlotlyKaleido.kill_kaleido()
