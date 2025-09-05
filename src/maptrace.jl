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
    tracetypes(celltypes::Vector{String}) -> NamedArray

Compute multi-step connectivity matrix for synaptic pathways traversing a sequence of cell types.

Given a sequence of cell types [T₀, T₁, ..., Tₙ], this computes the matrix product:
W[T₀, T₁] × W[T₁, T₂] × ... × W[Tₙ₋₁, Tₙ]

where W[A, B] represents the submatrix of the global connectivity matrix W restricted to 
presynaptic cells of type A and postsynaptic cells of type B.

# Arguments
- `celltypes::Vector{String}`: Ordered sequence of cell type names [T₀, T₁, ..., Tₙ] defining 
  the pathway. Must contain at least 2 cell types.
- `side::String`: Restrict to cells on given side ("left", "right"). Default "right".

# Returns  
- `NamedArray`: Multi-step connectivity matrix with cell IDs as row/column names.
  Entry (i,j) represents the summed synaptic strength over all possible n-step paths from 
  cell i (type T₀) to cell j (type Tₙ), where intermediate steps are constrained to traverse
  cells of types T₁, T₂, ..., Tₙ₋₁ in order.

# Examples
```julia
# Direct connectivity (1-step pathway)
direct = tracetypes(["Tm1", "Dm3v"])

# Disynaptic pathway (2-step: Tm1→T2a→Dm3v)
disynaptic = tracetypes(["Tm1", "T2a", "Dm3v"])  

# Trisynaptic pathway (3-step)
trisynaptic = tracetypes(["L1", "Tm1", "T2a", "Dm3v"])
```

# Notes
- Requires global connectivity matrix `W` as NamedArray with cell IDs as names
- For n+1 cell types, computes n matrix multiplications (n synaptic steps)
- Result dimensions: (# cells in T₀) × (# cells in Tₙ)
- Each entry sums over all intermediate cellular pathways satisfying the type constraints
"""
function tracetypes(celltypes::Vector{String}; side::String="right")
    @mfalse inds = [(ind2type .== celltype) .& (ind2side .== side) for celltype in celltypes]
    return *([W[inds[i], inds[i+1]] for i = 1:length(inds)-1]...)
end

"""
    tracebacktypes(celltypes::Vector{String}) -> NamedArray
    tracebacktypes(celltypes::Vector{Union{String, Vector{<:Integer}}}) -> NamedArray

Like `tracetypes`, but uses column-normalized connectivity matrix `Wback` instead of raw `W`.

`Wback` represents a backward Markov process where each column sums to 1, giving transition 
probabilities rather than absolute synapse counts. This is useful for analyzing relative 
input strengths and pathway probabilities.

# Arguments
- First method: sequence of cell type names, with optional `side` parameter
- Second method: allows mixing cell type names with vectors of specific cell IDs, with optional `side` parameter (applies only to String cell types)

# Examples
```julia
# Compare raw vs normalized connectivity
raw = tracetypes(["Tm1", "T2a", "Dm3v"])
normalized = tracebacktypes(["Tm1", "T2a", "Dm3v"])

# Use specific cell IDs for final target
lc10ev_cells = [720575940604125920, 720575940606274592]
mixed = tracebacktypes(["Tm1", "TmY9q", lc10ev_cells])
```
"""
function tracebacktypes(celltypes::Vector{String}; side::String="right")
    @mfalse inds = [(ind2type .== celltype) .& (ind2side .== side) for celltype in celltypes]
    return *([Wback[inds[i], inds[i+1]] for i = 1:length(inds)-1]...)
end

function tracebacktypes(celltypes::Vector{Union{String, Vector{<:Integer}}}; side::String="right")
    @mfalse inds = [isa(celltype, String) ? (ind2type .== celltype) .& (ind2side .== side) : id2ind.(celltype) for celltype in celltypes]
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
    inmaps(paths::NamedArray) -> NamedArray

Convert pathway connectivity matrix to spatial input maps (receptive fields).

Takes a connectivity matrix (typically from `tracetypes` or `tracebacktypes`) and creates 
spatial maps showing the locations of presynaptic inputs for each postsynaptic cell.

# Arguments
- `paths::NamedArray`: Connectivity matrix with presynaptic cell IDs as row names and 
  postsynaptic cell IDs as column names. Values represent connection strengths.

# Returns
- `NamedArray`: Collection of spatial input maps indexed by postsynaptic cell IDs. Each 
  element is a 2D matrix matching `pq2column` dimensions, where pixel values represent 
  the summed connectivity from presynaptic cells at those spatial locations.

# Examples
```julia
# Create receptive fields for Dm3v cells from Tm1 inputs
tm1_dm3v = tracetypes(["Tm1", "Dm3v"])
dm3v_rfs = inmaps(tm1_dm3v)

# Access specific cell's receptive field
cell_rf = dm3v_rfs[Name(some_dm3v_id)]

# Create ERF maps from multi-step pathways  
erf_paths = tracebacktypes(["Tm1", "T2a", "Dm3v"])
erf_maps = inmaps(erf_paths)
```

# Notes
- Requires global variables `id2pq` (cell spatial coordinates) and `pq2column` (spatial grid)
- Only maps cells that have known spatial coordinates in `id2pq`
- Result matrices use eye/column coordinate system for spatial visualization
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

"""
    outmaps(paths::NamedArray) -> NamedArray

Transpose of `inmaps`: create spatial output maps (projective fields) for each presynaptic cell.

While `inmaps` shows where each postsynaptic cell receives inputs from, `outmaps` shows 
where each presynaptic cell sends outputs to.

# Returns
- `NamedArray` indexed by presynaptic cell IDs, with spatial maps showing postsynaptic target locations

# Examples
```julia
# Compare input vs output spatial patterns
paths = tracetypes(["Tm1", "Dm3v"])
receptive_fields = inmaps(paths)   # what each Dm3v cell receives from
projective_fields = outmaps(paths) # where each Tm1 cell projects to
```
"""
function outmaps(paths::NamedArray)
    preids, postids = names(paths)
    pfs = NamedArray([zeros(Float32, size(pq2column)) for i = 1:length(preids)], preids)
    for preid in preids
        for postid in postids
            if haskey(id2pq, postid)
                if ~ismissing(id2pq[postid])
                    pfs[Name.(preid)][id2pq[postid]...] = paths[Name.(preid), Name.(postid)]
                end
            end
        end
    end
    return pfs  # projective fields
end

#### The following is deprecated legacy code.
#### It has been superseded by `inmaps` and other functions above

"""
    preimage(pretype, idpost) -> Matrix
    preimage(pretype, posttype; side="right") -> Vector{Matrix}

Create image showing cells in `pretype` presynaptic to specified postsynaptic cell or cell type.

For the single-cell method, automatically restricts presynaptic cells to the same side 
as the postsynaptic cell. For the cell-type method, allows manual side specification.

# Arguments
- `pretype::String`: Presynaptic cell type name
- `idpost::Integer`: Postsynaptic cell ID (first method)
- `posttype::String`: Postsynaptic cell type name (second method)
- `side::String`: For second method, restrict to cells on given side ("left", "right")

# Returns
- First method: Matrix showing presynaptic cell locations with synapse counts as pixel values
- Second method: Vector of matrices, one for each postsynaptic cell of the specified type and side
"""
function preimage(pretype::String, idpost::Integer)
    side = ind2side[id2ind[idpost]]     # could be a problem if side = missing
    @mfalse preinds = findall((ind2type .== pretype) .& (ind2side .== side))
    
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

function preimage(pretype::String, posttype::String; side::String="right")
    @mfalse postids = ind2id[(ind2type .== posttype) .& (ind2side .== side)]
    return [preimage(pretype, idpost) for idpost in postids]
end

function prepreimage(prepretype::String, pretype::String, idpost::Integer)
    side = ind2side[id2ind[idpost]]
    @assert !ismissing(side) "Postsynaptic cell must have a valid side"
    
    @mfalse preinds = findall((ind2type .== pretype) .& (ind2side .== side))
    @mfalse prepreinds = findall((ind2type .== prepretype) .& (ind2side .== side))
    
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
    
    # Assert all posts have valid sides and are on the same side
    sides = ind2side[postinds]
    @assert all(.!ismissing.(sides)) "All postsynaptic cells must have a valid side"
    @assert length(unique(sides)) == 1 "All postsynaptic cells must be on the same side"
    
    # Get the common side
    side = first(sides)
    
    # Apply side filtering
    @mfalse preinds = findall((ind2type .== pretype) .& (ind2side .== side))
    @mfalse prepreinds = findall((ind2type .== prepretype) .& (ind2side .== side))
    
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

function prepreimage(prepretype::String, pretype::String, posttype::String; side::String="right")
    @mfalse postinds = findall((ind2type .== posttype) .& (ind2side .== side))
    @mfalse preinds = findall((ind2type .== pretype) .& (ind2side .== side))
    @mfalse prepreinds = findall((ind2type .== prepretype) .& (ind2side .== side))
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
