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

function toppost(ids; nresults = 15, normalize=false)
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

function toppre(ids; nresults = 15, normalize=false)
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

function type2ids(celltype::String)
    @mfalse ind2id[ind2type .== celltype]
end

"""
    convert2arrows(s)

convert string to use doublearrows
formerly used when we had composite names
"""
function convert2arrows(s)
    replace.(s, "__from__" => "⇐", "__to__" => "⇒")
end
