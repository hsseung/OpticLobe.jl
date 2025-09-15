using MissingsAsFalse

using PrettyTables

using StatsBase

"""
    showall(x)

Display the complete contents of a data structure using plain text formatting.

Convenient function for showing all elements of vectors, matrices, or other data structures
without truncation. Particularly useful for named arrays and large data structures that 
would normally be abbreviated in the REPL.

# Examples
```julia
showall(type2ids("Tm1"))          # Show all Tm1 cell IDs
showall(Wtt["Tm1", :])            # Show all Tm1 outputs (full vector)
showall(intrinsictypes)           # Show all intrinsic type names
```

# Notes
- Equivalent to `show(stdout, "text/plain", x)`
- Bypasses REPL display limits for large arrays
- Useful for inspecting complete contents of connectivity vectors
"""
showall = x->show(stdout, "text/plain", x)

# BEWARE first version returns fraction of synapses, second version returns number of synapses

"""
    toppost(celltype::String; sort=:out, nresults=15, values=:fraction)

Display top postsynaptic partners of a cell type with connectivity fractions or mean synapse counts.

Shows the strongest postsynaptic connections from a given cell type, displaying either:
- Connectivity fractions: output/input fractions as percentages
- Mean synapse counts: average number of synapses per cell

# Arguments
- `celltype::String`: Name of the presynaptic cell type
- `sort::Symbol`: Sort by `:out` (output) or `:in` (input). Default `:out`.
- `nresults::Int`: Number of top connections to display. Default 15.
- `values::Symbol`: Display `:fraction` (percentages) or `:mean` (synapse counts). Default `:fraction`.

# Examples
```julia
toppost("Tm1")                           # Top 15 Tm1 postsynaptic partners (output fractions)
toppost("Tm1", sort=:in)                 # Sorted by input fraction instead
toppost("Dm3v", nresults=10)             # Show only top 10 connections
toppost("Tm1", values=:mean)             # Show mean synapse counts instead of fractions
toppost("Tm1", values=:mean, sort=:in)  # Sort by input mean
```

# Output
Pretty-printed table with:
- Row headers: postsynaptic cell type names
- When `values=:fraction`:
  - "out" column: percentage of presynaptic type's outputs going to each target
  - "in" column: percentage of each target's inputs coming from presynaptic type
- When `values=:mean`:
  - "out" column: mean synapses per presynaptic cell to each target type
  - "in" column: mean synapses per postsynaptic cell from presynaptic type

# Notes
- Uses `outfraction`/`infraction` matrices for fractions, `outmean`/`inmean` for means
- Output adapts to environment: HTML tables in IJulia, text tables in REPL or Emacs
- Fractions displayed as percentages (0-100), means as floating point values
"""
function toppost(celltype::String; sort=:out, nresults = 15, values=:fraction)
    if values == :mean
        # Use mean matrices
        if sort == :in
            perm = sortperm(inmean[celltype, :], rev=true)
        else
            perm = sortperm(outmean[celltype, :], rev=true)
        end
        data = [outmean[celltype, perm[1:nresults]] inmean[celltype, perm[1:nresults]]]
        formatter = ft_printf("%.1f")
        col_headers = ["out #", "in #"]
    else
        # Use fraction matrices (default)
        if sort == :in
            perm = sortperm(infraction[celltype, :], rev=true)
        else
            perm = sortperm(outfraction[celltype, :], rev=true)
        end
        data = [outfraction[celltype, perm[1:nresults]] infraction[celltype, perm[1:nresults]]]
        data = 100*data  # Convert to percentages
        formatter = ft_printf("%2d")
        col_headers = ["out %", "in %"]
    end
    
    # Choose backend based on environment
    if detect_emacs_jupyter()
        backend = Val(:text)
    elseif isdefined(Main, :IJulia) && Main.IJulia.inited
        backend = Val(:html)
    else
        backend = Val(:text)
    end
    
    pretty_table(hcat(col_headers, data');
        header = vcat([celltype*"-post"], names(outfraction)[1][perm[1:nresults]]), 
        backend = backend, 
        formatters = formatter
    )
end

"""
    toppre(celltype::String; sort=:in, nresults=15, values=:fraction)

Display top presynaptic partners of a cell type with connectivity fractions or mean synapse counts.

Shows the strongest presynaptic connections to a given cell type, displaying either:
- Connectivity fractions: input/output fractions as percentages
- Mean synapse counts: average number of synapses per cell

# Arguments
- `celltype::String`: Name of the postsynaptic cell type
- `sort::Symbol`: Sort by `:in` (input) or `:out` (output). Default `:in`.
- `nresults::Int`: Number of top connections to display. Default 15.
- `values::Symbol`: Display `:fraction` (percentages) or `:mean` (synapse counts). Default `:fraction`.

# Examples
```julia
toppre("Dm3v")                          # Top 15 Dm3v presynaptic partners (input fractions)
toppre("Dm3v", sort=:out)               # Sorted by output fraction instead
toppre("T2a", nresults=10)              # Show only top 10 connections
toppre("Dm3v", values=:mean)            # Show mean synapse counts instead of fractions
toppre("Dm3v", values=:mean, sort=:out) # Sort by output mean
```

# Output
Pretty-printed table with:
- Row headers: presynaptic cell type names
- When `values=:fraction`:
  - "in" column: percentage of postsynaptic type's inputs coming from each source
  - "out" column: percentage of each source's outputs going to postsynaptic type
- When `values=:mean`:
  - "in" column: mean synapses per postsynaptic cell from each presynaptic type
  - "out" column: mean synapses per presynaptic cell to postsynaptic type

# Notes
- Uses `infraction`/`outfraction` matrices for fractions, `inmean`/`outmean` for means
- Output adapts to environment: HTML tables in IJulia, text tables in REPL or Emacs
- Fractions displayed as percentages (0-100), means as floating point values
"""
function detect_emacs_jupyter()
    # Check if parent process contains "emacs"
    try
        if Sys.islinux() || Sys.isapple()
            # Get parent process info
            ppid = @ccall getppid()::Cint
            parent_cmd = read(`ps -p $ppid -o comm=`, String) |> strip
            return occursin("emacs", lowercase(parent_cmd))
        end
    catch
        # Fallback if ps command fails
    end
    return false
end

function toppre(celltype::String; sort=:in, nresults = 15, values=:fraction)
    if values == :mean
        # Use mean matrices
        if sort == :out
            perm = sortperm(outmean[:, celltype], rev=true)
        else
            perm = sortperm(inmean[:, celltype], rev=true)
        end
        data = [inmean[perm[1:nresults], celltype] outmean[perm[1:nresults], celltype]]
        formatter = ft_printf("%.1f")
        col_headers = ["in #", "out #"]
    else
        # Use fraction matrices (default)
        if sort == :out
            perm = sortperm(outfraction[:, celltype], rev=true)
        else
            perm = sortperm(infraction[:, celltype], rev=true)
        end
        data = [infraction[perm[1:nresults], celltype] outfraction[perm[1:nresults], celltype]]
        data = 100*data  # Convert to percentages
        formatter = ft_printf("%2d")
        col_headers = ["in %", "out %"]
    end
    
    # Choose backend based on environment
    if detect_emacs_jupyter()
        backend = Val(:text)
    elseif isdefined(Main, :IJulia) && Main.IJulia.inited
        backend = Val(:html)
    else
        backend = Val(:text)
    end
    
    pretty_table(hcat(col_headers, data');
        header = vcat([celltype*"-pre"], names(infraction)[1][perm[1:nresults]]), 
        backend = backend, 
        formatters = formatter
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

"""
    strings2ticks(sarray::Vector{String}) -> Tuple{UnitRange{Int}, Vector{String}}

Convert a vector of strings to tick positions and labels for plotting.

Returns a tuple of (positions, labels) suitable for use with plotting packages
like Plots.jl for setting custom tick marks on axes.

# Arguments
- `sarray::Vector{String}`: Vector of strings to use as tick labels

# Returns
- `Tuple`: (1:length(sarray), sarray) - positions and corresponding labels

# Examples
```julia
ticks = strings2ticks(["Tm1", "Tm2", "Tm3"])     # ([1, 2, 3], ["Tm1", "Tm2", "Tm3"])
ticks = strings2ticks(intrinsictypes[1:10])      # First 10 intrinsic types as ticks

# Usage in plotting:
using Plots
xticks = strings2ticks(["Type1", "Type2", "Type3"])
plot(data, xticks=xticks)
```

# Notes
- Positions start at 1 and increment by 1
- Commonly used for categorical axis labels in connectivity plots
- Output format matches plotting package expectations for custom ticks
"""
function strings2ticks(sarray::Vector{String})
    return (1:length(sarray), sarray)
end

"""
    convert2arrows(s)

convert string to use doublearrows
formerly used when we had composite names
"""
function convert2arrows(s)
    replace.(s, "__from__" => "⇐", "__to__" => "⇒")
end
