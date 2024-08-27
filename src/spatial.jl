function centroid(im; whole = false)
    if sum(im) > 0
        mu = [sum((1:size(im, 1)) .* im)/sum(im), sum(im .* (1:size(im, 2))')/sum(im)]
        return whole ? round.(Int64, mu) : mu
    else
        return missing    # special case of zero image: centroid not defined
    end
end

"""
find center by computing centroid or finding maximum
"""
function findcenter0(im; centroid = true)   # alternative is center defined by location of maximum
    if sum(im) > 0
        if centroid
            return round.(Int64, [sum((1:size(im, 1)) .* im)/sum(im), sum(im .* (1:size(im, 2))')/sum(im)])
        else
            p, q = Tuple(argmax(im))
            return [p, q]
        end
    else
        return missing    # special case of zero image: center not defined
    end
end


"""
find center by convolving `im` with `filt` and finding maximum
"""
function findcenter(im, filt) 
    if sum(im) > 0
        p, q = Tuple(argmax(conv2(im, filt)))
        return [p, q]
    else
        return missing    # special case of zero image: center not defined
    end
end

using OffsetArrays, PaddedViews

# align and average images in rfs
# centroid alignment
function kernel(rfs, maxsep = 4)
    s = OffsetArray(zeros(2*maxsep+1, 2*maxsep+1), -maxsep:maxsep, -maxsep:maxsep)
    pmax, qmax = size(rfs, 1), size(rfs, 2)
    for z = 1:size(rfs, 3)
        center = findcenter(rfs[:, :, z], centroid = true)
        if ~ismissing(center) 
            im = PaddedView(0, rfs[:, :, z], (-maxsep:pmax+maxsep+1, -maxsep:qmax+maxsep+1))
            x, y = center
            s[-maxsep:maxsep, -maxsep:maxsep] += im[x-maxsep:x+maxsep, y-maxsep:y+maxsep]
        end
    end
    return s/size(rfs, 3)
end

using NNlib

const newaxis = [CartesianIndex()]

# note that NNlib.jl convolution is really convolution, not cross-correlation
# convolution for 2D image and kernel (no channels or batch)
function conv2(x, w)    
    m, n = size(w)
    if isodd(m)
        padm = floor(Int64, m/2)
    else
        padm = (div(m, 2) - 1, div(m, 2))   
    end
    if isodd(n)
        padn = floor(Int64, n/2)
    else
        padn = (div(n, 2) - 1, div(n, 2))  # right side has more padding for even-sized kernel
    end
    return NNlib.conv(x[:, :, newaxis, newaxis], w[:, :, newaxis, newaxis], pad = (padm..., padn...))[:, :, 1, 1]
end


"""
convolutional clustering
"""
function convcluster(inmaps, nout = 2; init = nothing, radius = 1, niter = 20)
    k = 2*radius + 1
    nin = size(inmaps, 3)
    if isnothing(init)
        w = randn(k, k, nin, nout)
        w ./= sqrt.(sum(w.^2, dims = (1, 2, 3)))
    else
        w = init
    end
    inmapspadded = PaddedView(0, inmaps, (1-radius:size(inmaps,1) + radius, 1-radius:size(inmaps,2) + radius, 1:size(inmaps, 3), 1:size(inmaps, 4)))

    overlap = 0
    nwinners = 0
    winners = zeros(Int64, size(inmaps)[end])
    for iter = 1:niter
        outmaps = NNlib.conv(inmaps, w, pad=radius, flipped=true);

        winnerinfo = argmax(outmaps, dims=(1, 2, 3))[:]
        winners = getindex.(winnerinfo, 3)

        wnew = zeros(k, k, nin, nout)
        overlap = 0
        for i = 1:size(outmaps, 4)
            p, q, whichout, whichim = Tuple(winnerinfo[i])
            inmapscropped = inmapspadded[p-radius:p+radius, q-radius:q+radius, :, whichim] 
            wnew[:, :, :, whichout] += inmapscropped
            overlap += sum(w[:, :, :, whichout].*inmapscropped)
        end
#        nwinners = freqtable(getindex.(winners, 3))
#        println(overlap, nwinners)
        w = wnew./sqrt.(sum(wnew.^2, dims = (1, 2, 3)))
    end
    return w, overlap, winners
end

## various hand-designed filters
diamond = [1 1;
    1 1]

seven = 
[1 1 0;
1 1.1 1;
0 1 1]

### bar filters for Dm3

threev = 
[1 0 0;
0 1.1 0;
0 0 1]

threep = 
[0 1 0;
0 1.1 0;
0 1 0]

threeq = 
[0 0 0;
1 1.1 1;
0 0 0]

### bar filters for TmY

threeq⊥ = 
[0.5 0.5 0;
0 1.1 0;
0 0.5 0.5]

threeh = 
[0 0 1;
0 1.1 0;
1 0 0]
