## spatial maps of connectivity
## `inmaps` is the central function, and supersedes legacy code (preimage, postimage, ...)

"""
normalize the columns of W to sum to one
generally relevant if W is nonnegative sparse matrix
(in Julia, division does not preserve sparsity)
"""
function normalize_columns(W::SparseMatrixCSC)
    # https://stackoverflow.com/questions/24296856/in-julia-how-can-i-column-normalize-a-sparse-matrix
    # get the column sums of A
    S = float(vec(sum(W, dims=1)))
    S[S .== 0] .= Inf    # if the column is all zeros, it will remain unchanged (rather than divide-by-zero yielding NaN)

    # get the nonzero entries in A. ei is row index, ej is col index, ev is the value in A
    ei, ej, ev = findnz(W)

    # get the number or rows and columns in A
    m, n = size(W)

    # create a new normalized matrix. For each nonzero index (ei,ej), its new value will be
    # the old value divided by the sum of that column, which can be obtained by S[ej]
    return sparse(ei, ej, Float32.(ev)./S[ej], m, n)
end

## define backward Markov process on cells
Wback = normalize_columns(W.array)
Wback = NamedArray(Wback, names = (ind2id, ind2id))

"""
    tracetypes(celltypes)

sum over all cellular paths through the series of `celltypes`

global variable W assumed to be a NamedArray with ids as names
"""
function tracetypes(celltypes::Vector{String})
    @mfalse inds = [ind2type .== celltype for celltype in celltypes]
    return *([W[inds[i], inds[i+1]] for i = 1:length(inds)-1]...)
end

"""
    tracebacktypes(celltypes)

weighted sum over all cellular paths through the series of `celltypes`

global variable W assumed to be a NamedArray with ids as names
"""
function tracebacktypes(celltypes::Vector{String})
    @mfalse inds = [ind2type .== celltype for celltype in celltypes]
    return *([Wback[inds[i], inds[i+1]] for i = 1:length(inds)-1]...)
end

function tracebacktypes(celltypes::Vector{Union{String, Vector{<:Integer}}})
    @mfalse inds = [isa(celltype, String) ? ind2type .== celltype : id2ind.(celltype) for celltype in celltypes]
    return *([Wback[inds[i], inds[i+1]] for i = 1:length(inds)-1]...)
end

"""
    scorepath(celltypes)

scores a type path rather than summing over cellular paths
"""
function scorepath(celltypes::Vector{String})
    return *([infraction[celltypes[i], celltypes[i+1]] for i = 1:length(celltypes)-1]...)
end

"""
    inmaps(paths)

`paths` - NamedArray with ids as names

for each postid, map all preids to eye coordinates
return NamedArray where postids are names, and elements are eye matrices.
"""
function inmaps(paths::NamedArray)
    preids, postids = names(paths)
#    rfs = OrderedDict{Integer, Matrix}()   # change to Int64?
#    rfs = NamedArray([zeros(Int64, size(pq2column)) for i = 1:length(postids)], postids)
    rfs = NamedArray([zeros(Float32, size(pq2column)) for i = 1:length(postids)], postids)
    for postid in postids
        for preid in preids
            if haskey(id2pq, preid)
                if ~ismissing(id2pq[preid])
                    rfs[Name.(postid)][id2pq[preid]...] = paths[Name.(preid), Name.(postid)]
                end
            end
        end
    end
    return rfs
end

#### legacy code that has been superseded by `inmaps` and other functions above
#### deprecated

"""
create image showing cells in `pretype` presynaptic to cell `idpost`
`pretype` - list of possible presynaptic cell IDs, usually those belonging to one type
`im` - image showing locations of cells that are presynaptic to cell `idpost`. pixel values are synapse number.
"""
function preimage(pretype::String, idpost::Integer)
    @mfalse preinds = findall(ind2type .== pretype)
    im = zeros(Int64, size(pq2column))
    for ind in preinds
        id = ind2id[ind]
        if haskey(id2pq, id)
            if ~ismissing(id2pq[id])
                im[id2pq[id]...] = W[ind, id2ind[idpost]]
            end
        end
    end
    return im
end

function preimage(pretype::String, posttype::String)
    @mfalse postids = ind2id[ind2type .== posttype]
    return [preimage(pretype, idpost) for idpost in postids]
end

function prepreimage(prepretype::String, pretype::String, idpost::Integer)
    @mfalse preinds = findall(ind2type .== pretype)
    @mfalse prepreinds = findall(ind2type .== prepretype)
    
    im = zeros(Int64, size(pq2column))
    
    W2 = W[prepreinds, preinds]*W[preinds, id2ind[idpost]]
    for (i, ind) in enumerate(prepreinds)
        id = ind2id[ind]
        if haskey(id2pq, id)
            im[id2pq[id]...] = W2[i]
        end
    end
    return im
end

function prepreimage(prepretype::String, pretype::String, idspost::Vector{<: Integer})
    @mfalse postinds = id2ind.(idspost)
    @mfalse preinds = findall(ind2type .== pretype)
    @mfalse prepreinds = findall(ind2type .== prepretype)
    
    ims = [zeros(Int64, size(pq2column)) for _ in 1:length(postinds)]

    W2 = W[prepreinds, preinds]*W[preinds, postinds]
    for j = 1:length(postinds)
        for (i, ind) in enumerate(prepreinds)
            id = ind2id[ind]
            if haskey(id2pq, id)
                ims[j][id2pq[id]...] = W2[i, j]
            end
        end
    end
    return ims
end

function prepreimage(prepretype::String, pretype::String, posttype::String)
    @mfalse postinds = findall(ind2type .== posttype)
    @mfalse preinds = findall(ind2type .== pretype)
    @mfalse prepreinds = findall(ind2type .== prepretype)
    ims = [zeros(Int64, size(pq2column)) for _ in 1:length(postinds)]
    W2 = W[prepreinds, preinds]*W[preinds, postinds]
    for j in 1:length(postinds)
        for (i, ind) in enumerate(prepreinds)
            id = ind2id[ind]
            if haskey(id2pq, id)
                ims[j][id2pq[id]...] = W2[i, j]
            end
        end
    end
    return ims
end

function preprepreimage(preprepretype::String, prepretype::String, pretype::String, idpost::Integer)
    @mfalse preinds = findall(ind2type .== pretype)
    @mfalse prepreinds = findall(ind2type .== prepretype)
    @mfalse preprepreinds = findall(ind2type .== preprepretype)
    
    im = zeros(Int64, size(pq2column))
    
    W3 = W[preprepreinds, prepreinds]*W[prepreinds, preinds]*W[preinds, id2ind[idpost]]
    for (i, ind) in enumerate(preprepreinds)
        id = ind2id[ind]
        if haskey(id2pq, id)
            im[id2pq[id]...] = W3[i]
        end
    end
    return im
end

function preprepreimage(preprepretype::String, prepretype::String, pretype::String, posttype::String)
    @mfalse postinds = findall(ind2type .== posttype)
    @mfalse preinds = findall(ind2type .== pretype)
    @mfalse prepreinds = findall(ind2type .== prepretype)
    @mfalse preprepreinds = findall(ind2type .== preprepretype)
    
    ims = [zeros(Int64, size(pq2column)...) for _ in 1:length(postinds)]
    
    W3 = W[preprepreinds, prepreinds]*W[prepreinds, preinds]*W[preinds, postinds]
    for (j, postind) in enumerate(postinds)
        for (i, ind) in enumerate(preprepreinds)
            id = ind2id[ind]
            if haskey(id2pq, id)
                ims[j][id2pq[id]...] = W3[i, j]
            end
        end
    end
    return ims
end


######## WARNING the post functions need to be revised to produce vector of images

"""
create image showing cells in `pretype` presynaptic to cell `idpost`
`pretype` - list of possible presynaptic cell IDs, usually those belonging to one type
`im` - image showing locations of cells that are presynaptic to cell `idpost`. pixel values are synapse number.
"""
function postimage(idpre::Integer, posttype::String)
    @mfalse postinds = findall(ind2type .== posttype)
    postids = ind2id[postinds]
    im = zeros(Int64, size(pq2column))
    if ~isempty(postinds)
        for ind in postinds
            id = ind2id[ind]
            if haskey(id2pq, id)
                im[id2pq[id]...] = W[id2ind[idpre], ind]
            end
        end
    end
    return im
end

function postimage(idspre::Vector{<:Integer}, posttype::String)
    im = zeros(Int64, size(pq2column)..., length(idspre))
    for (i, idpre) in enumerate(idspre)
        im[:, :, i] = postimage(idpre, posttype)
    end
    return im
end

function postimage(pretype::String, posttype::String)
    @mfalse preids = ind2id[ind2type .== pretype]
    im = zeros(Int64, size(pq2column)..., length(preids))
    for (i, idpre) in enumerate(preids)
        im[:, :, i] = postimage(idpre, posttype)
    end
    return im
end
