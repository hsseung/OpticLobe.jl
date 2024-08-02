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

# %%
using Clustering, Distances

# %%
using NewickTree, Phylo

# %%
using NamedArrays, OffsetArrays

# %%
using StatsPlots, Measures

# %%
using ColorSchemes

# %%
using Printf

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
ass = cutree(c, h=0.9)
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

function get_tree(hc, labels)
   mindist = 0.4   # lowest possible distance. BEWARE: all distances relative to this baseline
   nodes = [NewickTree.Node(i, n=string(n), d=0.) for (i,n) in zip(hc.order,labels)]
   n = length(nodes)
   idfun(x) = x > 0 ? x + n : abs(x)
   for i=1:size(hc.merges, 1)
       nid = n + i
       j, k = idfun.(hc.merges[i,:])
       a = nodes[j]
       b = nodes[k]
       h = hc.heights[i] - mindist
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
forest = deepcopy(tree)
h = nodeheights(forest)
thres = maximum(h) - 0.5   # 0.5 is relative to `mindist` 0.4, so real threshold is 0.9

# %%
for node in h.axes[1][h .< thres]
    deletenode!(forest, node)
end

roots = getroots(forest)

clusters = Vector{String}[]
for node in roots
    d = filter( e -> !startswith(e, "Node"), getdescendants(forest, node.name))
    push!(clusters, isempty(d) ? [node.name] : d)
end

# %%
rootnames = [r.name for r in roots]
p = sortperm(indexin(rootnames, dfsnames))
rootnames = rootnames[p]
clusters = reverse.(clusters[p])

# %%
println(length(clusters), " clusters total, including ", sum(length.(clusters) .== 1), " singletons")
println("cluster sizes ", sort(length.(clusters), rev = true))

# %%
clusternames = Vector{String}(undef, length(clusters))
issingleton = length.(clusters) .== 1
ncluster = sum(.!issingleton)   # number of clusters that are not singletons
clusternames[.!issingleton] .= ["Cluster$i" for i = 1:ncluster]
clusternames[issingleton] .= vcat(clusters[issingleton]...)

# %% [markdown]
# ### colormap from Arie
# OffsetArray is convenient so 0th element can be black

# %%
cmap = OffsetArray(
    color.(
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

# %%
# use if manual palette not available
# cmap = OffsetArray(vcat([RGB(0,0,0)], distinguishable_colors(ncluster, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)), 0:ncluster)

# %%
function colorchildren(tree, start, i)
    map_depthfirst((val, node) -> node == start ? i : val, 0, tree)
end

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
temp = -ones(1, 16)
plot!(temp, color = permutedims(cmap[1:end]), linewidth = 3, legend = :bottomright, 
    label = permutedims(["Cluster$i" for i=1:ncluster]),
)

# %%
println("saving dendrogram")
StatsPlots.savefig("Fig 2c TypePolarDendrogram.svg")

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
savefig("Fig S11b inoutfractionclusters.svg")

# %% [markdown]
# ## generate separate colorbar

# %%
io = open("Fig S11b colorbar.svg", "w")
show(io, "image/svg+xml", ColorSchemes.inferno)
close(io)

# %% [markdown]
# ## dendrogram with finer clusters

# %%
forest = deepcopy(tree)
h = nodeheights(forest)
thres = maximum(h) - 0.485

# %%
for node in h.axes[1][h .< thres]
    deletenode!(forest, node)
end

roots = getroots(forest)

clusters = Vector{String}[]
for node in roots
    d = filter( e -> !startswith(e, "Node"), getdescendants(forest, node.name))
    push!(clusters, isempty(d) ? [node.name] : d)
end

# %%
rootnames = [r.name for r in roots]
p = sortperm(indexin(rootnames, dfsnames))
rootnames = rootnames[p]
clusters = reverse.(clusters[p])

# %%
println(length(clusters), " clusters total, including ", sum(length.(clusters) .== 1), " singletons")
println("cluster sizes ", sort(length.(clusters), rev = true))

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
c12 = distinguishable_colors(12, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)

# %%
cmap = OffsetArray(vcat([RGB(0, 0, 0)], c12, c12, c12), 0:36)

# %%
plot(tree, treetype = :fan, linecolor = cmap[colorings], linewidth=2)

# %%
println("saving dendrogram with finer flat clustering")
StatsPlots.savefig("Fig S10a TypePolarDendrogramFiner.svg")

# %% [markdown]
# ## dendrogram with finest clusters

# %%
forest = deepcopy(tree)
h = nodeheights(forest)
thres = maximum(h) - 0.46

# %%
for node in h.axes[1][h .< thres]
    deletenode!(forest, node)
end

roots = getroots(forest)

clusters = Vector{String}[]
for node in roots
    d = filter( e -> !startswith(e, "Node"), getdescendants(forest, node.name))
    push!(clusters, isempty(d) ? [node.name] : d)
end

# %%
rootnames = [r.name for r in roots]
p = sortperm(indexin(rootnames, dfsnames))
rootnames = rootnames[p]
clusters = reverse.(clusters[p])

# %%
println(length(clusters), " clusters total, including ", sum(length.(clusters) .== 1), " singletons")
println("cluster sizes ", sort(length.(clusters), rev = true))

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
StatsPlots.savefig("Fig S10b TypePolarDendrogramFinest.svg")
