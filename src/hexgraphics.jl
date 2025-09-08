# eye map, hexagonal lattice

using Luxor, ColorSchemes
using Printf
using LinearAlgebra
using MissingsAsFalse   # for triad functions

# eye map (depends on ColumnAssignments.jld2)
# zero at ommatidia/columns
# `missing` outside eye
eye = map(x -> ismissing(x) ? missing : 0, pq2column)

struct EyeMap
    im::Matrix
    function EyeMap(im)
        eyeim = allowmissing(im)
        eyeim[ismissing.(eye)] .= missing
        return new(eyeim)
    end
end


"""
    HexagonEye(p, q, hexelsize)

`p`, `q` - Zhao's eye coordinates

returns a Luxor hexagon

The conversion to axial coordinates is (p, q) -> (q, -p), followed by 30 degree rotation
"""
HexagonEye(p, q, hexelsize) = HexagonAxial(q, -p, hexelsize)


"""
    rect2hex(rect; hexelsize=6, pointannotation=nothing)

Display rectangular array data on a hexagonal lattice visualization.

Renders data stored in a rectangular array as hexagonal tiles (hexels) arranged in the 
natural hexagonal pattern of the optic lobe. Each array element becomes a colored hexagon
positioned according to its (p,q) coordinate indices.

# Arguments
- `rect`: Rectangular array containing color data for hexagonal tiles
- `hexelsize::Int`: Radius of hexagons from center to vertex (default: 6)
- `pointannotation::Tuple`: Optional (p,q) location to mark with cyan star

# Examples
```julia
# Display eye map data
data = rand(RGB, size(pq2column))
rect2hex(data, hexelsize=8)

# Mark a specific location
rect2hex(data, pointannotation=(15, 20))

# Typical usage with eye data
heatmap = eyehot(some_data)
rect2hex(heatmap)
```

# Notes
- Uses Luxor.jl for hexagonal rendering with 30° rotation
- Array center is placed at origin, with (1,1) direction pointing upward
- `missing` values in the array are skipped (not rendered)
- Coordinate system matches optic lobe (p,q) conventions
- Commonly used with `eyehot()` for heatmap visualization
"""
function rect2hex(rect; hexelsize = 6, pointannotation = nothing)
    Luxor.rotate(-pi/6)
    center = ceil.(Int64, size(rect)./2)  # approximate center of square, to place at origin
 
    for k in keys(rect)
        if ~ismissing(rect[k])
            setcolor(rect[k])
            p, q = Tuple(k) .- Tuple(center)
            poly(hextile(HexagonEye(p, q, hexelsize)), :fill)
        end
    end
    if ~isnothing(pointannotation) && ~ismissing(pointannotation)
        p, q = Tuple(pointannotation) .- Tuple(center)
        c = hexcenter(HexagonEye(p, q, hexelsize))
        setcolor("cyan")
        star(c, hexelsize, 6, 0.2, 0.0, :fill)
    end
    Luxor.rotate(pi/6)
end



"""
convert to RGB values using the `hot` colormap

`x` - scalar or array
"""
function hot(x)
    get(ColorSchemes.hot, x)
end

"""
    eyehot(im[, sat])

Convert one channel to heatmap on eye map

Default is to normalize output so that maximum value is one 
(or approximately one due to presence of `eps` in denominator).

Call with `sat = 1` if no normalization of image is desired

`im` - array that is same size as the eye map/mask

convert to image using `hot` 

"""
function eyehot(im::AbstractArray, sat::Real = maximum(im))
    passmissing(hot).((im + eye)/(sat + eps()))
end

"""
    eyeheat(im[, sat])

Convert one channel to heatmap on eye map

Default is to normalize output so that maximum value is one 
(or approximately one due to presence of `eps` in denominator).

Call with `sat = 1` if no normalization of image is desired

`im` - array that is same size as the eye map/mask

convert to image using `hot` 

"""
function eyeheat(im::AbstractArray, sat::Real = maximum(im); cmap = ColorSchemes.inferno)
    passmissing(x -> get(cmap, x)).((im + eye)/(sat + eps()))
end

"""
    eyergb(r, g[, b, sat])

Convert three channels to RGB image on eye map.

Default is to normalize output so that maximum value of each channel is one 
(or approximately one due to presence of `eps` in denominator).

Call with `sat = [1, 1, 1]` if no normalization of image is desired

If called with two channels, blue channel is set to zero.

`r`, `g`, `b` - arrays containing RGB values, same size as the eye map/mask
`sat` - array containing saturation values for channels
"""
function eyergb(r::AbstractArray, g::AbstractArray, b::AbstractArray = zeros(size(r)), sat = [maximum(r), maximum(g), maximum(b)])
    rgbstack = stack([r/(sat[1] + eps()), g/(sat[2] + eps()), b/(sat[3] + eps())], dims=1)
    rgbcombined = allowmissing(colorview(RGB, rgbstack))
    rgbcombined[ismissing.(eye)] .= missing
    return rgbcombined
end

"""
    square2hex(square::AbstractArray)

Convert square array to hexagonal shape by masking corner regions.

Transforms a square array into a hexagonal shape by setting corner regions to `missing`,
creating a hexagonal cutout suitable for display on hexagonal lattices. The resulting
shape approximates a hexagon by removing pixels where the row-column index difference
exceeds half the array size.

# Arguments
- `square::AbstractArray`: Square input array to convert to hexagonal shape

# Returns
- Array of same type with corners set to `missing` to form hexagonal shape

# Examples
```julia
# Convert square kernel to hexagonal shape
kernel = rand(11, 11)
hex_kernel = square2hex(kernel)

# Typical usage in visualization pipeline
data = crop(some_image, center, radius)
hex_data = square2hex(data)
rect2hex(hex_data)
```

# Notes
- Used primarily for displaying kernels and cropped regions on hexagonal grids
- Corner regions where `abs(i - j) > size(square,1)/2` are set to `missing`
- Commonly used with `crop()` function for spatial analysis visualization
- Essential step before `rect2hex()` for proper hexagonal display
- Creates approximate hexagon - not perfect but visually effective
"""
function square2hex(square::AbstractArray) # (used to display kernel, or hexagonal cutout from eyemap)
    out = allowmissing(similar(square))
    m = floor(size(square, 1)/2)
    for i = 1:size(square, 1)
        for j = 1:size(square, 2)
            if abs(i - j) > m
                out[i, j] = missing
            else
                out[i, j] = square[i, j]
            end
        end
    end
    return out
end


"""
    crop(im::Matrix, center::Vector, radius, sat=maximum(im))

Extract a square region from an image centered at specified coordinates.

Crops a square patch of size `(2*radius+1) × (2*radius+1)` from the input image,
centered at the given (p,q) coordinates. Values outside the original image bounds
are padded with zeros. The result is normalized by the saturation value.

# Arguments
- `im::Matrix`: Input image to crop from
- `center::Vector`: [p, q] coordinates for crop center
- `radius`: Half-width of the square crop region
- `sat`: Saturation value for normalization (default: maximum of image)

# Returns
- Cropped and normalized square array of size `(2*radius+1, 2*radius+1)`

# Examples
```julia
# Crop 11×11 region around a cell location
cell_center = id2pq[some_cell_id]
cropped = crop(connectivity_map, cell_center, 5)

# Crop with custom normalization
cropped = crop(data, [15, 20], 3, sat=1.0)

# Typical workflow for hexagonal display
patch = crop(heatmap, center, radius)
hex_patch = square2hex(patch) 
rect2hex(hex_patch)
```

# Notes
- Uses `PaddedView` to handle out-of-bounds regions with zero padding
- Normalization: divides by `sat + eps()` to prevent division by zero
- Commonly used for extracting local connectivity patterns around cells
- Output is suitable for `square2hex()` conversion for hexagonal visualization
- Coordinates follow (p,q) convention matching optic lobe spatial system
"""
function crop(im::Matrix, center::Vector, radius, sat = maximum(im))
    p, q = center[1], center[2]
    return collect(PaddedView(0, im/(sat + eps()), (p-radius:p+radius, q-radius:q+radius)))
end

function crophot(im::Matrix, center::Vector, radius, sat = maximum(im))
    p, q = center[1], center[2]
    return square2hex(get(ColorSchemes.hot, collect(PaddedView(0, im/(sat + eps()), (p-radius:p+radius, q-radius:q+radius)))))
end



"""
    montage(ims; hexelsize=4, fname=nothing, labels=nothing, ellipses=false, ellipsecolors=nothing, maxvals=false, fontsize=2.5*hexelsize, montagesize=nothing, centers=nothing, major="row", summary=nothing)

Create a grid layout of hexagonal eye maps for comparative visualization.

Arranges multiple eye map images in a montage layout with hexagonal rendering,
suitable for comparing connectivity patterns, cell types, or other spatial data
across the optic lobe. Each image is displayed as a hexagonal lattice.

# Arguments
- `ims`: Vector of images or 3D/4D array to arrange in montage
- `hexelsize::Int`: Size of hexagonal tiles (default: 4)
- `fname::String`: Optional filename to save the montage
- `labels::Vector`: Text labels for each image panel
- `ellipses::Bool`: Whether to overlay ellipse fits on images
- `ellipsecolors`: Colors for ellipse overlays
- `maxvals::Bool`: Whether to display maximum values as text
- `fontsize::Real`: Size of text labels (default: 2.5 × hexelsize)
- `montagesize::Tuple`: Override automatic grid size (rows, cols)
- `centers`: Center coordinates for cropping each image
- `major::String`: Layout order - "row" (default) or "column" major
- `summary`: Optional summary panel to include

# Examples
```julia
# Basic montage of connectivity maps
maps = [outmaps("Tm1", target) for target in ["Dm1", "Dm2", "Dm3"]]
montage(maps, labels=["Tm1→Dm1", "Tm1→Dm2", "Tm1→Dm3"])

# Save montage with custom layout
montage(data, fname="comparison.png", montagesize=(2, 3), hexelsize=6)

# Include ellipse fits and maximum values
montage(maps, ellipses=true, maxvals=true, fontsize=12)
```

# Notes
- Automatically determines grid layout to be approximately square if not specified
- Uses hexagonal rendering via `rect2hex()` for each panel
- Supports both row-major (default) and column-major ordering
- Can handle 3D arrays (WHN) or vectors of 2D arrays
- Essential tool for comparative analysis of spatial connectivity patterns
- Integrates with Luxor.jl for high-quality vector graphics output
"""
function montage(ims; hexelsize = 4, fname = nothing, labels = nothing, ellipses = false, ellipsecolors = nothing, maxvals = false, fontsize = 2.5 * hexelsize, montagesize = nothing, centers = nothing, major = "row", summary = nothing)
    m, n = size(ims[1])
    spacing = ceil(min(m, n)*hexelsize*1.7)
    
    nim = length(ims)
    
    if ndims(ims) == 2
        nchannel, ncluster = size(ims)
    elseif ndims(ims) == 1
        if isnothing(montagesize)
            nchannel = ~isnothing(summary) ? ceil(Int64, sqrt(nim + 1)) : ceil(Int64, sqrt(nim))
            ncluster = nchannel
        else
            nchannel, ncluster = montagesize
        end
    end
    
    if major == "row" || major == :row
        ncol, nrow = nchannel, ncluster
    else
        nrow, ncol = nchannel, ncluster
    end

    if isnothing(fname)
        Drawing(round(ncol*spacing*1.1), round(nrow*spacing*1.1))
    else
        Drawing(ncol*spacing*1.1, nrow*spacing*1.1, fname)
    end
    
    origin()

    Luxor.translate(-spacing*(ncol-1)/2, -spacing*(nrow-1)/2)

    for i = 1:nim
        if ~isnothing(centers)
            rect2hex(eyehot(ims[i]), hexelsize = hexelsize, pointannotation = centers[i])
        else
            rect2hex(eyehot(ims[i]), hexelsize = hexelsize)
        end
        if ~isnothing(labels)
            setcolor("black")
            Luxor.fontsize(fontsize)
            Luxor.text(labels[i], Point(-0.97*spacing/2, 0.17*spacing/2), angle = 0.9*pi/3, valign = :middle)
        end
        if ellipses
            e = ellipsesummary(ims[i])
            if ~isnothing(ellipsecolors)
                drawellipse(e, color = ellipsecolors[i], hexelsize = hexelsize)
            else
                drawellipse(e, hexelsize = hexelsize)
            end
        end
        if maxvals
            setcolor("black")
            Luxor.fontsize(fontsize)
            Luxor.text(@sprintf("max %d, sum %d", 10000*maximum(ims[i]), 10000*sum(ims[i])), Point(-spacing/2, -0.23*spacing/2), angle = -pi/3, valign = :middle)
        end
        Luxor.translate(spacing, 0)
        if rem(i, ncol) == 0
            Luxor.translate(-spacing*ncol, spacing)
        end
    end
    # draw ellipses together at 3x scale
    if ~isnothing(summary)
        anchor = ellipsesummary(ims[summary])
        if ~ismissing(anchor)
            Luxor.rotate(-pi/6)
            Luxor.translate(-3*hexelsize*anchor.center)
            Luxor.rotate(pi/6)
        end
        for (i, e) in enumerate(ellipsesummary.(ims))
            if i == summary
                setdash("dashed")
            else
                setdash("solid")
            end
            drawellipse(e, color = ellipsecolors[i], hexelsize = 3*hexelsize)
        end
        if ~ismissing(anchor)
            Luxor.rotate(-pi/6)
            Luxor.translate(3*hexelsize*anchor.center)
            Luxor.rotate(pi/6)
        end
        sethue("black")
        drawpqaxes(3*hexelsize)
    end
    finish()
    if isnothing(fname)
        preview()
    end
end


"""
    triad(imv, imp, imq; kwargs...)

display three hexel images "orbiting" the origin, clockwise from bottom
hexels must be colorants
"""
function triad(imv, imp, imq; hexelsize = 7, orbitscale = 1, textorbitscale = 1, ellipses = nothing, ellipsecolors = nothing, text = nothing, textcolors = nothing, fontsize = 20)
    imsize = size(imv, 1)
    orbit = orbitscale*hexelsize*imsize

    # directions in Luxor's coordinate system (x increases from left to right, y from up to down)
    v = (0, 1)  # -v axis
    p = (-sqrt(3)/2, -0.5)  # +p axis
    q = (sqrt(3)/2, -0.5)  # +q axis
    
#    origin()
#    Luxor.translate(-0.2*orbit.*v...)
    for (i, location, im) in zip(1:3, [v, p, q], [imv, imp, imq])
        Luxor.translate(orbit.*location...) 
        rect2hex(im, hexelsize = hexelsize)   
        if ~isnothing(ellipses)
            c = ~isnothing(ellipsecolors) ? ellipsecolors[i] : "blue"
            if ~isnothing(c)
                drawellipse(ellipses[i], color = c, hexelsize = hexelsize)
            end
        end
        if ~isnothing(text)
            c = ~isnothing(textcolors) ? textcolors[i] : "blue"
            setcolor(c)
            Luxor.fontsize(fontsize)
#            Luxor.text(text[i], Point(0, orbit/1.25), halign = :center)
            if i == 1
                Luxor.text(text[i], Point(0, textorbitscale*orbit), halign = :center)
            else
                Luxor.text(text[i], Point(0, -textorbitscale*orbit), halign = :center)
            end
        end
        Luxor.translate(-orbit.*location...)  
    end
end

function eyetriad(pretype::String, postids::Vector; hexelsize = 5, cmap = ColorSchemes.hot, ellipse = false, ellipsecolors = nothing, text = nothing, textcolors = nothing)
    @mfalse inputmaps = [preimage(pretype, id) for id in postids]

    immax = maximum(hcat(inputmaps...))
    println(immax)
    
    if ellipse
        triad(eyehot.(inputmaps)..., hexelsize = hexelsize, ellipses = ellipsesummary.(inputmaps[1:3]), ellipsecolors = ellipsecolors, text = text, textcolors = textcolors)
    else
        triad(eyehot.(inputmaps)..., hexelsize = hexelsize, text = text, textcolors = textcolors)
    end
end 

function celltriad(pretype::String, postids::Vector; radius = 5, hexelsize = 7, cmap = ColorSchemes.hot, ellipse = false)
    @mfalse inputmaps = [passmissing(crop)(preimage(pretype, id), id2pq(id), radius, 1) for id in postids]
    immax = maximum.(inputmaps)
    println(immax)
    
    if ellipse
        triad([square2hex(get(cmap, inputmaps[i]/immax[i])) for i = 1:3]..., hexelsize = hexelsize, ellipses = ellipsesummary.(inputmaps[1:3]))
    else
        triad([square2hex(get(cmap, inputmaps[i]/immax[i])) for i = 1:3]..., hexelsize = hexelsize)
    end
    return inputmaps
end

function typetriad(pretype::String, posttypes::Vector{String}; radius = 6, hexelsize = 7, cmap = ColorSchemes.hot, ellipse = false, normalize = :common, side = "right")
    @mfalse inputmaps = [passmissing(crop).(preimage(pretype, posttype; side=side), id2pq.(ind2id[(ind2type .== posttype) .& (ind2side .== side)]), radius, 1) for posttype in posttypes]
    inputmaps = mean.(skipmissing.(inputmaps))

    immax = maximum.(inputmaps)
    if normalize == :common
        immax .= maximum(immax)
    end
    println(immax)

    if ellipse
        triad([square2hex(get(cmap, inputmaps[i]/immax[i])) for i = 1:3]..., hexelsize = hexelsize, ellipses = ellipsesummary.(inputmaps[1:3]))
    else
        triad([square2hex(get(cmap, inputmaps[i]/immax[i])) for i = 1:3]..., hexelsize = hexelsize)
    end
    return inputmaps
end

# function typetriad(pretype::String, posttypes::Vector{String}; radius = 6, hexelsize = 7, cmap = ColorSchemes.hot, ellipse = false)
#     @mfalse inputmaps = [passmissing(crop).(preimage(pretype, posttype), id2pq.(ind2id[ind2type .== posttype]), radius, 1) for posttype in posttypes]
#     inputmaps = mean.(skipmissing.(inputmaps))

#     immax = maximum.(inputmaps)
#     println(immax)

#     if ellipse
#         triad([square2hex(get(cmap, inputmaps[i]/immax[i])) for i = 1:3]..., hexelsize = hexelsize, ellipses = ellipsesummary.(inputmaps[1:3]))
#     else
#         triad([square2hex(get(cmap, inputmaps[i]/immax[i])) for i = 1:3]..., hexelsize = hexelsize)
#     end
#     return inputmaps
# end

function typetriad(prepretype::String, pretype::String, posttypes::Vector{String}; radius = 6, hexelsize = 7, cmap = ColorSchemes.hot, ellipse = false, side = "right")
    @mfalse inputmaps = [passmissing(crop).(prepreimage(prepretype, pretype, posttype), id2pq.(ind2id[(ind2type .== posttype) .& (ind2side .== side)]), radius, 1) for posttype in posttypes]
    inputmaps = mean.(skipmissing.(inputmaps))

    immax = maximum.(inputmaps)
    
    if ellipse
        triad([square2hex(get(cmap, inputmaps[i]/immax[i])) for i = 1:3]..., hexelsize = hexelsize, ellipses = ellipsesummary.(inputmaps[1:3]))
    else
        triad([square2hex(get(cmap, inputmaps[i]/immax[i])) for i = 1:3]..., hexelsize = hexelsize)
    end
    return inputmaps
end


"""
    tetrad(im1, im2, im3, im4; kwargs...)

display four hexel images
hexels must be colorants
"""
function tetrad(im1, im2, im3, im4; hexelsize = 7, orbitscale = 1, textorbitscale = 1, ellipses = nothing, ellipsecolors = nothing, text = nothing, textcolors = nothing, fontsize = 20)
    imsize = size(im1, 1)
    orbit = orbitscale*hexelsize*imsize

    # directions in Luxor's coordinate system (x increases from left to right, y from up to down)
    v1 = (-1, -2/sqrt(3))
    v2 = (1, -2/sqrt(3))
    v3 = (-1, 2/sqrt(3))
    v4 = (1, 2/sqrt(3))
    
    v = [(-1, -1), (1, -1), (-1, 1), (1, 1)]
    
#    Luxor.translate(-0.2*orbit.*v...)
    for (i, location, im) in zip(1:4, v, [im1, im2, im3, im4])
        Luxor.translate(orbit.*location...) 
        rect2hex(im, hexelsize = hexelsize)   
        if ~isnothing(ellipses)
            c = ~isnothing(ellipsecolors) ? ellipsecolors[i] : "blue"
            if ~isnothing(c)
                drawellipse(ellipses[i], color = c, hexelsize = hexelsize)
            end
        end
        if ~isnothing(text)
            c = ~isnothing(textcolors) ? textcolors[i] : "blue"
            setcolor(c)
            Luxor.fontsize(fontsize)
#            Luxor.text(text[i], Point(0, orbit/1.25), halign = :center)
            if i == 3 || i == 4
                Luxor.text(text[i], Point(0, textorbitscale*orbit), halign = :center, valign = :middle)
            else
                Luxor.text(text[i], Point(0, -textorbitscale*orbit), halign = :center, valign = :middle)
            end
        end
        Luxor.translate(-orbit.*location...)  
    end
end


######

"""
    hexannulus(ims; kwargs...)

display hexel images arranged to form an annulus
hexels must be colorants
"""
function hexannulus(ims; hexelsize = 7, orbitscale = 1, textorbitscale = 1, ellipses = nothing, ellipsecolors = nothing, text = nothing, textcolors = nothing, maxvals = false, fontsize = 16)
    imsize = size(ims[1], 1)
    orbit = orbitscale*hexelsize*imsize
    
    nims = length(ims)
    ndirs = 6

    # directions in Luxor's coordinate system (x increases from left to right, y from up to down)
    v = [Point(cos(theta), sin(theta)) for theta in (0:nims-1).*2*pi/ndirs]
    vtext = [Point(cos(theta), sin(theta)) for theta in ((0:nims-1) .+ 0.9).*2*pi/ndirs]
    
    s = maximum.(ims)
#    s[ellipsecolors .== "green"] .= maximum(s[ellipsecolors .== "green"])
#   s[ellipsecolors .!= "green"] .= maximum(s[ellipsecolors .!= "green"])
    
#    Luxor.translate(-0.2*orbit.*v...)
    for (i, location, im) in zip(1:nims, v, ims)
        Luxor.translate(orbit.*location...) 
        rect2hex(square2hex(hot(im/s[i])), hexelsize = hexelsize)   
        if ~isnothing(ellipses)
            c = ~isnothing(ellipsecolors) ? ellipsecolors[i] : "blue"
            if ~isnothing(c)
                drawellipse(ellipses[i], color = c, hexelsize = hexelsize)
            end
        end
        if ~isnothing(text)
            c = ~isnothing(textcolors) ? textcolors[i] : "black"
            setcolor(c)
            Luxor.fontsize(fontsize)
            Luxor.text(text[i], vtext[i].*0.7*orbit + (0, -10), halign = :center, valign = :middle)
        end
        if maxvals
            setcolor("black")
            Luxor.fontsize(16)
            Luxor.text(@sprintf("%3d %3d", 10000*maximum(im), 10000*sum(im)), vtext[i].*0.7*orbit + (0, 10), valign = :middle, halign = :center)
#            Luxor.text(@sprintf("sum %3d", 10000*sum(im)), location.*0.75*orbit + (0, 20), valign = :middle, halign = :center)
        end
        Luxor.translate(-orbit.*location...)  
    end
end

########## ellipses

struct Ellipse
    center::Point
    theta::AbstractFloat
    major::AbstractFloat   # major diameter
    minor::AbstractFloat   # minor diameter
end

# im is hexel image with pq coordinates
# WARNING ellipse parameters are for hexelsize = 1. divide by sqrt(3) to convert to units of latticeconstant
# WARNING XY coordinate system is tilted by 30 degrees because of Luxor's hex coordinates
function ellipsesummary(im::Matrix) 
    hexelsize = 1
    if sum(im) > 0
        m, n = ceil.(Int64, size(im)./2)
        XY = [hexcenter(HexagonAxial(ind[2] - n, - ind[1] + m, hexelsize)) for ind in keys(im)]    # conversion from pq to xy coordinates
        mu = sum(XY.*im)/sum(im)
        covmatrix = sum(map(a -> collect(a)*collect(a)', XY .- mu).*im)/sum(im) .+ 5/24*I(2)
        F = eigen(covmatrix)
        theta = atan(F.vectors[1, 2]/F.vectors[1, 1])
        sigmamin, sigmamax = sqrt.(abs.(sort(F.values)))
        return Ellipse(Point(mu...), theta, 2*sigmamax, 2*sigmamin)
    else
        return missing
    end
end

function ellipsesummary(im::Array{<:Real, 3})
    params = []
    for i = 1:size(im, 3)
        push!(params, ellipsesummary(im[:, :, i]))
    end
    return params
end

# draw ellipse for eye map
function drawellipse(params; color = "blue", hexelsize = 7, action = :stroke)
    if ~isnothing(params) && ~isnothing(color) && ~ismissing(params)
        Luxor.rotate(-pi/6)
            Luxor.translate(hexelsize*params.center)
                Luxor.rotate(params.theta)
                    sethue(color)
                    Luxor.ellipse(O, hexelsize*params.minor, hexelsize*params.major, action = action)
                Luxor.rotate(-params.theta)
            Luxor.translate(-hexelsize*params.center)
        Luxor.rotate(pi/6)
    end
end 

function drawline(pt1, pt2; color = "blue", hexelsize = 6)
    Luxor.rotate(-pi/6)
    setcolor(color)
    Luxor.line(hexelsize*pt1, hexelsize*pt2, action = :stroke)
    Luxor.rotate(pi/6)
end 

"""
    hexproject(im, axis)

    `im` - square matrix of size m = 2k + 1
    output is always vector of length m
    locations are:
        -k:k for cardinal orientations
        (-k:k)*sqrt(3)/2 for orthogonal orientations
"""
function hexproject(im::Matrix, axis)
    @assert size(im, 1) == size(im, 2)
    m = size(im, 1)
    @assert isodd(m)
    w = [0.5, 1, 0.5]
    k = m ÷ 2
    P = zeros(Int32, m, m) .+ collect(-k:k)
    Q = zeros(Int32, m, m) .+ collect(-k:k)'
# cardinal orientations, lattice spacing
    if axis == :v    # P + Q runs from -2k:2k, output will run from -k:k
        proj = counts(P + Q, weights(im))
        proj = dropdims(conv(proj[:, newaxis, :], w[:, newaxis, newaxis], pad = 1, stride = 2), dims = 2)
        return proj
    elseif axis == :p  # P + Q runs from -3k:3k. crop to -2k:2k
        proj = counts(2*P - Q, -2*k:2*k, weights(im))
        proj = dropdims(conv(proj[:, newaxis, :], w[:, newaxis, newaxis], pad = 1, stride = 2), dims = 2)
        return proj
    elseif axis == :q
        proj = counts(2*Q - P, -2*k:2*k, weights(im))
        proj = dropdims(conv(proj[:, newaxis, :], w[:, newaxis, newaxis], pad = 1, stride = 2), dims = 2)
        return proj
# orthogonal orientations, sqrt(3)/2 times lattice spacing
    elseif axis == :h
        proj = counts(Q - P, -k:k, weights(im))
        return proj
    elseif axis == :p⊥
        proj = counts(Q, weights(im))
        return proj
    elseif axis == :q⊥
        proj = counts(P, weights(im))
        return proj
    end
end

# original version of code from notebook...should be the same as above but keeping just in case
# function hexproject(im::Matrix, axis)
#     @assert size(im, 1) == size(im, 2)
#     m = size(im, 1)
#     @assert isodd(m)
#     w = [0.5, 1, 0.5]
#     k = m ÷ 2
#     P = zeros(Int32, m, m) .+ collect(-k:k)
#     Q = zeros(Int32, m, m) .+ collect(-k:k)'
# # cardinal orientations, lattice spacing
#     if axis == :v    # P + Q runs from -2k:2k, output will run from -k:k
#         proj = counts(P + Q, weights(im))
#         proj = dropdims(conv(proj[:, newaxis, :], w[:, newaxis, newaxis], pad = 1, stride = 2), dims = 2)
#         return proj
#     elseif axis == :p  # P + Q runs from -3k:3k. crop to -2k:2k
#         proj = counts(2*P - Q, -2*k:2*k, weights(im))
#         proj = dropdims(conv(proj[:, newaxis, :], w[:, newaxis, newaxis], pad = 1, stride = 2), dims = 2)
#         return proj
#     elseif axis == :q
#         proj = counts(2*Q - P, -2*k:2*k, weights(im))
#         proj = dropdims(conv(proj[:, newaxis, :], w[:, newaxis, newaxis], pad = 1, stride = 2), dims = 2)
#         return proj
# # orthogonal orientations, sqrt(3)/2 times lattice spacing
#     elseif axis == :h
#         proj = counts(Q - P, -k:k, weights(im))
#         return proj
#     elseif axis == :p⊥
#         proj = counts(Q, weights(im))
#         return proj
#     elseif axis == :q⊥
#         proj = counts(P, weights(im))
#         return proj
#     end
# end

function drawpqaxes(hexelsize)
    Luxor.rotate(-pi/6)
    Luxor.poly(hexcenter.([HexagonEye(1, 0, hexelsize), HexagonEye(0, 0, hexelsize), HexagonEye(0, 1, hexelsize)]), :stroke)
    Luxor.rotate(pi/6)
end

