# Visual stimuli for hexagonal lattice

using ColorSchemes

"""
    sine_grating(wavevector, phase; amplitude=1.0, offset=0.0)

Generate a sine grating visual stimulus on the hexagonal lattice.

Creates a sine grating of the form: `amplitude * sin(wavevector ⋅ [p, q] + phase) + offset`
where (p,q) are the hexagonal coordinates of each lattice point and wavevector is [kp, kq].

# Arguments
- `wavevector::Vector{<:Real}`: 2-element wave vector [kp, kq] specifying spatial frequency and direction
- `phase::Real`: Phase offset of the grating in radians
- `amplitude::Real`: Amplitude of the sine modulation (default: 1.0)
- `offset::Real`: DC offset added to the grating (default: 0.0)

# Returns
- Matrix with same dimensions as `pq2column`, containing grating values at each (p,q) position

# Examples
```julia
# Grating along p+q direction (equivalent to previous spatialfrequency=0.1)
grating = sine_grating([0.1, 0.1], 0)

# Grating along p direction only
grating_p = sine_grating([0.1, 0], 0)

# Grating along q direction only
grating_q = sine_grating([0, 0.1], 0)

# Diagonal grating
grating_diag = sine_grating([0.1, -0.05], π/4)

# Display as heatmap
rect2hex(eyeheat(grating))

# High contrast grating
high_contrast = sine_grating([0.15, 0.15], 0, amplitude=1.0, offset=0.0)
rect2hex(eyeheat(high_contrast, cmap=ColorSchemes.grays))
```

# Notes
- Uses the existing `pq2column` matrix to define the hexagonal lattice structure
- Grating values are only computed for valid hexagonal positions (non-missing entries)
- Wave vector [kp, kq] allows arbitrary orientation and spatial frequency
- For visualization, typically used with `eyeheat()` or `eyehot()` followed by `rect2hex()`
- `missing` values are preserved in regions outside the defined hexagonal lattice
"""
function sine_grating(wavevector::Vector{<:Real}, phase::Real; amplitude::Real=1.0, offset::Real=0.0)
    if length(wavevector) != 2
        throw(ArgumentError("wavevector must be a 2-element vector [kp, kq]"))
    end

    kp, kq = wavevector
    grating = similar(pq2column, Union{Float64, Missing})

    for (idx, column_id) in pairs(pq2column)
        if !ismissing(column_id)
            p, q = Tuple(idx)
            grating[idx] = amplitude * sin(kp * p + kq * q + phase) + offset
        else
            grating[idx] = missing
        end
    end

    return grating
end