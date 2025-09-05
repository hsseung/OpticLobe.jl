using MissingsAsFalse

using PrettyTables

using StatsBase

showall = x->show(stdout, "text/plain", x)

# BEWARE first version returns fraction of synapses, second version returns number of synapses

function toppost(celltype::String; sort=:out, nresults = 15)
    backend = isdefined(Main, :IJulia) && Main.IJulia.inited ? Val(:html) : Val(:text)
    if sort == :in
        perm = sortperm(infraction[celltype, :], rev=true)
    else
        perm = sortperm(outfraction[celltype, :], rev=true)
    end
    data = [outfraction[celltype, perm[1:nresults]] infraction[celltype, perm[1:nresults]]]
    pretty_table(hcat(["out", "in"], 100*data');
        header = vcat([celltype*"-post"], names(outfraction)[1][perm[1:nresults]]), 
        backend = backend, 
        formatters = ft_printf("%2d")
    )
end

function toppre(celltype::String; sort=:in, nresults = 15)
    backend = isdefined(Main, :IJulia) && Main.IJulia.inited ? Val(:html) : Val(:text)
    if sort == :out
        perm = sortperm(outfraction[:, celltype], rev=true)
    else
        perm = sortperm(infraction[:, celltype], rev=true)
    end
    data = [infraction[perm[1:nresults], celltype] outfraction[perm[1:nresults], celltype]]
    pretty_table(hcat(["in", "out"], 100*data');
        header = vcat([celltype*"-pre"], names(infraction)[1][perm[1:nresults]]), 
        backend = backend, 
        formatters = ft_printf("%2d")
    )
end

function toppost(ids::Vector{Int64}; nresults = 15, normalize=false)
    backend = isdefined(Main, :IJulia) && Main.IJulia.inited ? Val(:html) : Val(:text)
    if normalize
        post = sum(W[id2ind.(ids), :]*A, dims=1)[:].array
        perm = sortperm(post, rev = true)
        posttotal = sum(W[id2ind.(ids), :])
        data = [post[perm[1:nresults]]/posttotal post[perm[1:nresults]]./sum(W*A, dims=1)[perm[1:nresults]]]
        pretty_table(hcat(["out", "in"], 100*data'); 
            header = vcat(["-post"], alltypes[perm[1:nresults]]),
            backend = backend, 
            formatters = ft_printf("%2d")
        )
    else
        post = mean(W[id2ind.(ids), :]*A, dims=1)[:].array
        perm = sortperm(post, rev = true)
        pretty_table(hcat(["avg syn #"], post[perm[1:nresults]]'); 
            header = vcat(["-post"], alltypes[perm[1:nresults]]),
            backend = backend, 
            formatters = ft_printf("%2d")
        )
    end
end

function toppre(ids::Vector{Int64}; nresults = 15, normalize=false)
    backend = isdefined(Main, :IJulia) && Main.IJulia.inited ? Val(:html) : Val(:text)
    if normalize
        pre = sum(A'*W[:, id2ind.(ids)], dims=2)[:].array
        perm = sortperm(pre, rev = true)
        pretotal = sum(W[:, id2ind.(ids)])
        data = [pre[perm[1:nresults]]/pretotal pre[perm[1:nresults]]./sum(A'*W, dims=2)[perm[1:nresults]]]
        pretty_table(hcat(["% in", "% out"], 100*data'); 
            header = vcat(["-pre"], alltypes[perm[1:nresults]]),
            backend = backend, 
            formatters = ft_printf("%2d")
        )
    else
        pre = mean(A'*W[:, id2ind.(ids)], dims=2)[:].array
        perm = sortperm(pre, rev = true)
        pretty_table(hcat(["avg syn #"], pre[perm[1:nresults]]'); 
            header = vcat(["-pre"], alltypes[perm[1:nresults]]),
            backend = backend, 
            formatters = ft_printf("%2d")
        )
    end
end

function strings2ticks(sarray::Vector{String})
    return (1:length(sarray), sarray)
end

"""
    type2ids(celltype::String; side::String="right") -> Vector{Int64}

Get cell IDs for all cells of a specified type, with optional side filtering for visual types.

For visual cell types (those in `visualtypes`), applies side filtering to return only cells 
from the specified hemisphere. For non-visual cell types (central brain, etc.), returns all 
cells regardless of side since ind2side doesn't yet contain information about non-visual types.

# Arguments
- `celltype::String`: Name of the cell type
- `side::String`: Side filter for visual types ("left", "right"). Default "right". 
  Ignored for non-visual cell types.

# Returns
- `Vector{Int64}`: Cell IDs of the specified type and side (for visual types)

# Examples
```julia
# Visual cell type - applies side filtering
tm1_right = type2ids("Tm1")                    # right side (default)
tm1_left = type2ids("Tm1", side="left")        # left side

# Non-visual cell type - ignores side parameter  
central_cells = type2ids("SomeCentralType")     # returns all cells
```

# Notes
- Uses global variables `visualtypes`, `ind2id`, `ind2type`, `ind2side`
- Visual vs non-visual determination based on membership in `visualtypes`
- TODO: Extend side filtering to work for non-visual types that have side information
"""
function type2ids(celltype::String; side::String="right")
    if celltype in visualtypes
        @mfalse ind2id[(ind2type .== celltype) .& (ind2side .== side)]
    else
        @mfalse ind2id[ind2type .== celltype]
    end
end

"""
    convert2arrows(s)

convert string to use doublearrows
formerly used when we had composite names
"""
function convert2arrows(s)
    replace.(s, "__from__" => "⇐", "__to__" => "⇒")
end
