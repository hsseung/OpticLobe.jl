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
# # Dendrogram of type centers

# %%
using OpticLobe

include("config.jl")

# %%
using Clustering, Distances

# %%
using NewickTree, Phylo

# %%
using NamedArrays, OffsetArrays

# %%
using StatsPlots, Measures

# %%
using ColorSchemes, Colors

# %%
using Printf

# %%
using Preferences

# Configuration for different synapse versions
# Load preference from OpticLobe package since this is a standalone script
const SYNAPSE_VERSION = Preferences.load_preference(OpticLobe, "default_synapses", "Princeton")

# Synapse-specific thresholds and settings
const SYNAPSE_CONFIG = Dict(
    "Princeton" => Dict(
        :thres_coarse => 0.89,
        :thres_fine => 0.885, 
        :thres_finest => 0.86,
        :colormap => :auto        # Use automatic colormap
    ),
    "Buhmann" => Dict(
        :thres_coarse => 0.9,    # May need adjustment
        :thres_fine => 0.885,     # May need adjustment
        :thres_finest => 0.86,    # May need adjustment
        :colormap => :manual      # Use manual colormap from Arie for coarse threshold
    )
)

println("Using $(SYNAPSE_VERSION) synapses with configuration: ", SYNAPSE_CONFIG[SYNAPSE_VERSION])

# %%
# remove photoreceptors from list of types
OpticLobe.intrinsictypes = filter(!startswith("R"), intrinsictypes)

# %% [markdown]
# ### average linkage clustering

# %%
# feature vectors and pairwise distances
X = vcat(infraction[intrinsictypes, intrinsictypes], outfraction[intrinsictypes, intrinsictypes]')
dist = pairwise(jaccard, X.array)

# %%
# hierarchical clustering
c = hclust(dist, linkage = :average)

# %%
# flat clustering, ordered by cluster size
# This code remains as it might be convenient for future explorations.
# For the figures, the flat clusterings will be recomputed later using Phylo
config = SYNAPSE_CONFIG[SYNAPSE_VERSION]
thres_absolute = config[:thres_coarse]
ass = cutree(c, h=thres_absolute)
order = sortperm(counts(ass), rev=true)

for i in order
    intrinsictypes[ass .== i] |> println
end

# %% [markdown]
# ## convert to Newick format, and use cell types as node names
# Conversion is necessary to use `Phylo` for plotting dendrograms

# %%
# code modified from https://github.com/JuliaStats/Clustering.jl/issues/209
using NewickTree: setdistance!

# Global constant for minimum distance baseline
# BEWARE: all distances relative to this baseline
const MINDIST = 0.4

function extract_clusters(tree, threshold_absolute, dfsnames, mindist=MINDIST)
    """Extract clusters from tree using given threshold"""
    forest = deepcopy(tree)
    h = nodeheights(forest)
    thres = maximum(h) - (threshold_absolute - mindist)
    
    for node in h.axes[1][h .< thres]
        deletenode!(forest, node)
    end
    
    roots = getroots(forest)
    clusters = Vector{String}[]
    for node in roots
        d = filter(e -> !startswith(e, "Node"), getdescendants(forest, node.name))
        push!(clusters, isempty(d) ? [node.name] : d)
    end
    
    rootnames = [r.name for r in roots]
    p = sortperm(indexin(rootnames, dfsnames))
    rootnames = rootnames[p]
    clusters = reverse.(clusters[p])
    
    return clusters, rootnames
end

function color_clusters(tree, rootnames)
    """Assign colors to tree branches based on cluster roots"""
    coloring = []
    for (i, root) in enumerate(rootnames)
        push!(coloring, colorchildren(tree, root, i))
    end
    colorings = +(coloring...)
    
    unique_colors = unique(colorings)
    reassign = Dict(unique_colors .=> 0:(length(unique_colors)-1))
    for i = 1:length(colorings)
        if colorings[i] > 0
            colorings[i] = reassign[colorings[i]]
        end
    end
    
    ncluster = length(rootnames)
    return colorings, ncluster
end

function print_cluster_stats(clusters)
    """Print statistics about cluster sizes"""
    println(length(clusters), " clusters total, including ", 
            sum(length.(clusters) .== 1), " singletons")
    println("cluster sizes ", sort(length.(clusters), rev=true))
end

function get_tree(hc, labels)
   nodes = [NewickTree.Node(i, n=string(n), d=0.) for (i,n) in zip(hc.order,labels)]
   n = length(nodes)
   idfun(x) = x > 0 ? x + n : abs(x)
   for i=1:size(hc.merges, 1)
       nid = n + i
       j, k = idfun.(hc.merges[i,:])
       a = nodes[j]
       b = nodes[k]
       h = hc.heights[i] - MINDIST
#       h = hc.heights[i]
       newnode = NewickTree.Node(nid, n="$nid", d=h)
       setdistance!(a, h - NewickTree.distance(a))
       setdistance!(b, h - NewickTree.distance(b))
       push!(newnode, a)
       push!(newnode, b)
       push!(nodes, newnode)
   end
   setdistance!(nodes[end], 0.)
   return nodes[end]
end

# %%
# convert to Newick format written as a string
tr = get_tree(c, intrinsictypes)
io = IOBuffer()
writenw(io, tr)
s = String(take!(io))

# %% [markdown]
# ## use `Phylo` to plot radial dendrogram, a.k.a. polar or circular dendrogram

# %%
unsortedtree = parsenewick(s)

# %%
tree = deepcopy(unsortedtree)
sort!(tree, rev=true)   # descendants from each node are sorted in order of their size. This is called ladderize in some other packages.

# %%
# node names in the order in which they will be plotted
dfsnames = [n.name for n in traversal(tree, inorder)]

# %% [markdown]
# ### threshold tree to obtain a forest

# %%
clusters, rootnames = extract_clusters(tree, thres_absolute, dfsnames)

# %%
print_cluster_stats(clusters)

# %%
clusternames = Vector{String}(undef, length(clusters))
issingleton = length.(clusters) .== 1
ncluster = sum(.!issingleton)   # number of clusters that are not singletons
clusternames[.!issingleton] .= ["Cluster$i" for i = 1:ncluster]
clusternames[issingleton] .= vcat(clusters[issingleton]...)

# %% [markdown]
# ### colormap configuration
# Use manual colormap for Buhmann, automatic for Princeton

# %%
if config[:colormap] == :manual
    # Manual colormap from Arie (for Buhmann)
    cmap = OffsetArray(
        parse.(Colorant,
            ["#000000",
             "#ffff00",
             "#008000",
             "#ffffaf",
             "#ffd700",
             "#00ffff",
             "#ff00ff",
             "#ee82ee",
             "#4169E1",
             "#7fff00",
             "#1e90ff",
             "#ff4500",
             "#00fa9a",
             "#F0E68C",
             "#F5DEB3",
             "#BDB76B",
             "#b8860b",
             ]
        ),
        0:16)
else
    # Automatic distinguishable colors (for Princeton)
    cmap = OffsetArray(vcat([RGB(0,0,0)], distinguishable_colors(ncluster, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)), 0:ncluster)
end

# %%
function colorchildren(tree, start, i)
    map_depthfirst((val, node) -> node == start ? i : val, 0, tree)
end

# %%
# label descendants of each root
colorings, ncluster = color_clusters(tree, rootnames)

# %%
default(; # Plots defaults
    fontfamily="Helvetica",
    )

# %%
plot(tree, treetype = :fan, linecolor = cmap[colorings], linewidth=2)
temp = -ones(1, 16)
plot!(temp, color = permutedims(cmap[1:end]), linewidth = 3, legend = :bottomright, 
    label = permutedims(["Cluster$i" for i=1:ncluster]),
)

# %%
println("saving dendrogram")
savefig(joinpath(TARGETDIR, "Fig 2c TypePolarDendrogram.pdf"))

# %% [markdown]
# ## heatmap of connectivity between clusters

# %%
Wclustertocluster = [sum(Wtt[c1, c2]) for c1 in clusters, c2 in clusters]

# %%
outfractioncluster = Wclustertocluster ./ [sum(Wtc[c, :]) for c in clusters]   # out fraction
infractioncluster = Wclustertocluster ./ [sum(Wct[:, c]) for c in clusters]'   # in fraction
hout = heatmap(outfractioncluster, yflip = :true, 
    xticks = strings2ticks(clusternames), xrot = 60, 
    yticks = strings2ticks(clusternames), 
    legend = :none
    )   # out fraction
hin = heatmap(infractioncluster, yflip = true, 
    xticks = strings2ticks(clusternames), xrot = 60, 
    yticks = strings2ticks(clusternames), 
    legend = :none
    )   # in fraction
@printf("maximum infraction and outfraction = %3.2f %3.2f\n", maximum(infractioncluster), maximum(outfractioncluster))

# %%
default(; # Plots defaults
    fontfamily="Helvetica",
    )
plot(hin, hout, size = (800, 400), bottom_margin = 11mm)

# %%
println("saving input and output fractions for connections between flat clusters")
savefig(joinpath(TARGETDIR, "Fig S11b inoutfractionclusters.svg"))  # pdf doesn't render properly

# %% [markdown]
# ## generate separate colorbar

# %%
io = open("Fig S11b colorbar.svg", "w")
show(io, "image/svg+xml", ColorSchemes.inferno)
close(io)

# %% [markdown]
# ## dendrogram with finer clusters

# %%
thres_absolute = config[:thres_fine]
clusters, rootnames = extract_clusters(tree, thres_absolute, dfsnames)

# %%
print_cluster_stats(clusters)

# %%
issingleton = length.(clusters) .== 1
ncluster = sum(.!issingleton)   # number of clusters that are not singletons

# %%
# label descendants of each root
colorings, ncluster = color_clusters(tree, rootnames)

# %%
default(; # Plots defaults
    fontfamily="Helvetica",
    )

# %%
# c12 = distinguishable_colors(12, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
c12 = collect(palette(:mk_12)) # https://mk.bcgsc.ca/colorblind/

# %%
cmap = OffsetArray(vcat([RGB(0, 0, 0)], c12, c12, c12), 0:36)

# %%
plot(tree, treetype = :fan, linecolor = cmap[colorings], linewidth=2)

# %%
println("saving dendrogram with finer flat clustering")
savefig(joinpath(TARGETDIR, "Fig S10a TypePolarDendrogramFiner.pdf"))

# %% [markdown]
# ## dendrogram with finest clusters

# %%
thres_absolute = config[:thres_finest]
clusters, rootnames = extract_clusters(tree, thres_absolute, dfsnames)

# %%
print_cluster_stats(clusters)

# %%
issingleton = length.(clusters) .== 1
ncluster = sum(.!issingleton)   # number of clusters that are not singletons

# %%
# label descendants of each root
coloring = []
for (i, root) in enumerate(rootnames)
    push!(coloring, colorchildren(tree, root, i))
end
colorings = +(coloring...)

# %%
reassign = Dict(unique(colorings) .=> 0:ncluster)
for i = 1:length(colorings)
    if colorings[i] > 0
        colorings[i] = reassign[colorings[i]]
    end
end

# %%
default(; # Plots defaults
    fontfamily="Helvetica",
    )

# %%
plot(tree, treetype = :fan, linecolor = cmap[colorings], linewidth=2)

# %%
println("saving dendrogram with finest flat clustering")
savefig(joinpath(TARGETDIR, "Fig S10b TypePolarDendrogramFinest.pdf"))
