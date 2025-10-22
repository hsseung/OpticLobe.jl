# interface to FlyWire resources: Codex, neuroglancer

using DataStructures, JSON, URIs, MissingsAsFalse

"""
    codex_open(idlist::Vector{Int}, version=783)

Open FlyWire Codex in a new browser window to search and view specified cells.

Launches the FlyWire Codex web interface with a pre-filled search query for the provided 
cell IDs, allowing you to explore detailed cell information, morphology, and connectivity data.

# Arguments
- `idlist::Vector{Int}`: Vector of FlyWire root IDs to search for in Codex
- `version::Int`: FlyWire data version (default: 783 for v783 proofreading)

# Examples
```julia
codex_open([720575940599333574])                    # Open single cell in Codex
codex_open(type2ids("Tm1")[1:5])                    # Open first 5 Tm1 cells
codex_open([720575940599333574], 630)               # Use v630 data instead
```

# Notes
- Opens a new browser tab/window with the Codex interface
- Works in both Jupyter notebooks (JavaScript) and REPL/command line (system browser)
- Codex provides detailed cell information including morphology, connectivity, and annotations
- Default uses v783 proofreading data
"""
function codex_open(idlist::Vector{Int}, version = 783)
    query = join(idlist, "+")
    url = "https://codex.flywire.ai/app/search?filter_string=$query&sort_by=&page_size=10&data_version=$version"
    
    # Check if we're in a notebook environment that can render JavaScript
    if isdefined(Main, :IJulia) && Main.IJulia.inited
        urlquoted = "\"" * url * "\""
        display("text/javascript", """window.open($urlquoted)""")
    else
        # REPL or command line - open in system browser
        try
            if Sys.isapple()
                run(`open $url`)
            elseif Sys.islinux()
                run(`xdg-open $url`)
            elseif Sys.iswindows()
                run(`start $url`)
            else
                println("Cannot open browser automatically. Please open this URL manually:")
                println(url)
            end
        catch
            println("Cannot open browser automatically. Please open this URL manually:")
            println(url)
        end
    end
end


"""
    ng_url(idlist[, version; darkmode=false])

return a string containing a neuroglancer URL

# Arguments
- `idlist::Vector{Int}`: Vector of FlyWire root IDs to visualize
- `version::Int`: FlyWire data version (default: 783)
- `darkmode::Bool`: Use dark background in 3D view (default: false)
"""
function ng_url(idlist::Vector{Int}, version = 783; darkmode::Bool = false)
    if version == 783
        seg = "flywire_v141_m783"
    else
        seg = "flywire_v141_m630"
    end
    
    s=
"""
{
  "dimensions": {
    "x": [
      1.6e-8,
      "m"
    ],
    "y": [
      1.6e-8,
      "m"
    ],
    "z": [
      4e-8,
      "m"
    ]
  },
  "position": [
    32736.037109375,
    14171.48828125,
    4156.3671875
  ],
  "crossSectionScale": 1,
  "projectionScale": 30000,
  "layers": [
    {
      "type": "image",
      "source": "precomputed://https://bossdb-open-data.s3.amazonaws.com/flywire/fafbv14",
      "tab": "source",
      "name": "EM"
    },
    {
      "type": "segmentation",
      "source": "precomputed://gs://flywire_neuropil_meshes/whole_neuropil/brain_mesh_v3",
      "tab": "source",
      "objectAlpha": 0.05,
      "hideSegmentZero": false,
      "segments": [
        "1"
      ],
      "segmentColors": {
        "1": "#b5b5b5"
      },
      "name": "brain_mesh_v3"
    },
    {
      "type": "segmentation",
      "source": "precomputed://gs://$seg",
      "tab": "segments",
      "name": "$seg"
    }
  ],
  "showDefaultAnnotations": false,
  "showSlices": false,
  "selectedLayer": {
    "visible": true,
    "layer": "$seg"
  },
  "projectionBackgroundColor": "#ffffff",
  "layout": "3d"
}
"""
    j = JSON.parse(s, dicttype=DataStructures.OrderedDict)   # regular Dict messes up order of dimensions
#    j["layers"][3]["segmentQuery"] = join(idlist, ' ')
    j["layers"][3]["segments"] = string.(idlist)

    # Add dark mode background if requested
    if darkmode
        j["perspectiveViewBackgroundColor"] = "#000000"
    end

    jencoded = URIs.escapeuri(json(j))
#    return "https://neuroglancer-demo.appspot.com/#!$jencoded"
    return "https://ngl.cave-explorer.org/#!$jencoded"
end

"""
    ng_url(idlists::Vector{Vector{Int}}, version=783; names=nothing, darkmode=false)

return a string containing a neuroglancer URL with multiple segmentation layers

# Arguments
- `idlists::Vector{Vector{Int}}`: Vector of cell ID lists, each shown in a separate layer
- `version::Int`: FlyWire data version (default: 783)
- `names::Vector{String}`: Optional custom names for each layer (default: auto-generated)
- `darkmode::Bool`: Use dark background in 3D view (default: false)
"""
function ng_url(idlists::Vector{Vector{Int}}, version = 783; names::Union{Nothing,Vector{String}} = nothing, darkmode::Bool = false)
    if version == 783
        seg = "flywire_v141_m783"
    else
        seg = "flywire_v141_m630"
    end

    s=
"""
{
  "dimensions": {
    "x": [
      1.6e-8,
      "m"
    ],
    "y": [
      1.6e-8,
      "m"
    ],
    "z": [
      4e-8,
      "m"
    ]
  },
  "position": [
    32736.037109375,
    14171.48828125,
    4156.3671875
  ],
  "crossSectionScale": 1,
  "projectionScale": 30000,
  "layers": [
    {
      "type": "image",
      "source": "precomputed://https://bossdb-open-data.s3.amazonaws.com/flywire/fafbv14",
      "tab": "source",
      "name": "EM"
    },
    {
      "type": "segmentation",
      "source": "precomputed://gs://flywire_neuropil_meshes/whole_neuropil/brain_mesh_v3",
      "tab": "source",
      "objectAlpha": 0.05,
      "hideSegmentZero": false,
      "segments": [
        "1"
      ],
      "segmentColors": {
        "1": "#b5b5b5"
      },
      "name": "brain_mesh_v3"
    }
  ],
  "showDefaultAnnotations": false,
  "showSlices": false,
  "projectionBackgroundColor": "#ffffff",
  "layout": "3d"
}
"""
    j = JSON.parse(s, dicttype=DataStructures.OrderedDict)

    # Add a segmentation layer for each idlist
    for (i, idlist) in enumerate(idlists)
        layer_name = names !== nothing && i <= length(names) ? names[i] : "$seg-$i"
        layer = DataStructures.OrderedDict(
            "type" => "segmentation",
            "source" => "precomputed://gs://$seg",
            "tab" => "segments",
            "segments" => string.(idlist),
            "name" => layer_name
        )
        push!(j["layers"], layer)
    end

    # Set selected layer to the first segmentation layer
    first_layer_name = names !== nothing && length(names) >= 1 ? names[1] : "$(seg)-1"
    j["selectedLayer"] = DataStructures.OrderedDict(
        "visible" => true,
        "layer" => first_layer_name
    )

    # Add dark mode background if requested
    if darkmode
        j["perspectiveViewBackgroundColor"] = "#000000"
    end

    jencoded = URIs.escapeuri(json(j))
    return "https://ngl.cave-explorer.org/#!$jencoded"
end

######### open neuroglancer window

"""
    ng_open(idlist::Vector{Int}, version=783; darkmode=false)
    ng_open(celltype::String, version=783; darkmode=false)

Open Neuroglancer 3D viewer in a new browser window to visualize specified cells.

Launches the Neuroglancer web interface with a pre-configured 3D view showing the
specified cells overlaid on the FlyWire brain volume. Useful for exploring cell
morphology, spatial relationships, and anatomical context.

# Arguments
- `idlist::Vector{Int}`: Vector of FlyWire root IDs to visualize
- `celltype::String`: Cell type name (shows all cells of that type)
- `version::Int`: FlyWire data version (default: 783 for v783 proofreading)
- `darkmode::Bool`: Use dark background in 3D view (default: false)

# Examples
```julia
ng_open([720575940599333574])                     # Open single cell in 3D
ng_open(type2ids("Tm1")[1:10])                    # Open first 10 Tm1 cells
ng_open("Dm3v")                                   # Open all Dm3v cells
ng_open([720575940599333574], 630)                # Use v630 data instead
ng_open([720575940599333574], darkmode=true)      # Open with dark background
```

# Notes
- Opens a new browser tab/window with Neuroglancer interface
- Works in both Jupyter notebooks (JavaScript) and REPL/command line (system browser)
- Shows cells in 3D context with brain mesh and EM data
- Default uses v783 proofreading data with appropriate segmentation layer
- Useful for morphological analysis and spatial relationship exploration
"""
function ng_open(idlist::Vector{Int}, version = 783; darkmode::Bool = false)
    url = ng_url(idlist, version, darkmode=darkmode)
    
    # Check if we're in a notebook environment that can render JavaScript
    if isdefined(Main, :IJulia) && Main.IJulia.inited
        urlquoted = "\"" * url * "\""
        display("text/javascript", """window.open($urlquoted)""")
    else
        # REPL or command line - open in system browser
        try
            if Sys.isapple()
                run(`open $url`)
            elseif Sys.islinux()
                run(`xdg-open $url`)
            elseif Sys.iswindows()
                run(`start $url`)
            else
                println("Cannot open browser automatically. Please open this URL manually:")
                println(url)
            end
        catch
            println("Cannot open browser automatically. Please open this URL manually:")
            println(url)
        end
    end
end

function ng_open(celltype::String, version = 783; darkmode::Bool = false)
    @mfalse ng_open(ind2id[ind2type .== celltype], version, darkmode=darkmode)
end

"""
    ng_open(idlists::Vector{Vector{Int}}, version=783; names=nothing, darkmode=false)

Open Neuroglancer with multiple segmentation layers, one for each list of cell IDs.

Opens a single browser window with multiple segmentation layers in one view,
allowing easy comparison between different cell groups.

# Arguments
- `idlists::Vector{Vector{Int}}`: Vector of cell ID lists, each shown in a separate segmentation layer
- `version::Int`: FlyWire data version (default: 783)
- `names::Vector{String}`: Optional custom names for each layer (default: auto-generated)
- `darkmode::Bool`: Use dark background in 3D view (default: false)

# Examples
```julia
# Open with three segmentation layers
ng_open([[id1, id2], [id3, id4], [id5, id6]])

# Compare different cell types in separate layers with custom names
ng_open([type2ids("Tm1")[1:5], type2ids("Tm2")[1:5]], names=["Tm1", "Tm2"])
```
"""
function ng_open(idlists::Vector{Vector{Int}}, version = 783; names::Union{Nothing,Vector{String}} = nothing, darkmode::Bool = false)
    url = ng_url(idlists, version, names=names, darkmode=darkmode)

    # Check if we're in a notebook environment that can render JavaScript
    if isdefined(Main, :IJulia) && Main.IJulia.inited
        urlquoted = "\"" * url * "\""
        display("text/javascript", """window.open($urlquoted)""")
    else
        # REPL or command line - open in system browser
        try
            if Sys.isapple()
                run(`open $url`)
            elseif Sys.islinux()
                run(`xdg-open $url`)
            elseif Sys.iswindows()
                run(`start $url`)
            else
                println("Cannot open browser automatically. Please open this URL manually:")
                println(url)
            end
        catch
            println("Cannot open browser automatically. Please open this URL manually:")
            println(url)
        end
    end
end

######### create hyperlink

"""
    ng_hyper(idlist::Vector{Int}; anchor=join(idlist, ", "), version=783, darkmode=false)
    ng_hyper(celltype::String; version=783, darkmode=false)

Create a clickable hyperlink or URL to view specified cells in Neuroglancer.

Generates either a clickable HTML hyperlink (in Jupyter notebooks) or prints a URL
(in REPL) that opens Neuroglancer with the specified cells pre-loaded for 3D visualization.

# Arguments
- `idlist::Vector{Int}`: Vector of FlyWire root IDs to visualize
- `celltype::String`: Cell type name (shows all cells of that type)
- `anchor::String`: Text to display for the hyperlink (default: comma-separated cell IDs)
- `version::Int`: FlyWire data version (default: 783 for v783 proofreading)
- `darkmode::Bool`: Use dark background in 3D view (default: false)

# Returns
- In Jupyter: HTML hyperlink object that can be displayed
- In REPL: Prints URL and returns `nothing`

# Examples
```julia
# In Jupyter notebooks:
ng_hyper([720575940599333574])                    # Clickable link with cell ID as text
ng_hyper(type2ids("Tm1")[1:5], anchor="5 Tm1 cells")  # Custom link text
ng_hyper("Dm3v")                                  # Link to all Dm3v cells
ng_hyper([720575940599333574], darkmode=true)     # Link with dark background

# In REPL:
ng_hyper([720575940599333574])                    # Prints: "720575940599333574: https://..."
```

# Notes
- Environment detection: shows HTML links in Jupyter, prints URLs in REPL
- Hyperlinks are clickable in Jupyter notebooks for convenient navigation
- URL format compatible with Neuroglancer's state encoding system
- Default uses v783 proofreading data
"""
function ng_hyper(idlist::Vector{Int}; anchor = join(idlist, ", "), version = 783, darkmode::Bool = false)
    url = ng_url(idlist, version, darkmode=darkmode)
    
    # Check if we're in a notebook environment that can render HTML
    if isdefined(Main, :IJulia) && Main.IJulia.inited
        return HTML("<a href=$url>$anchor</a>")
    else
        # REPL or command line - just print the URL
        println("$anchor: $url")
        return nothing
    end
end

function ng_hyper(celltype::String; version = 783, darkmode::Bool = false)
    @mfalse idlist = ind2id[ind2type .== celltype]
    return ng_hyper(idlist, anchor = celltype, version = version, darkmode = darkmode)
end

"""
    ng_hyper(idlists::Vector{Vector{Int}}; anchor="multi-layer view", names=nothing, version=783, darkmode=false)

Create a single hyperlink or URL that opens multiple segmentation layers in Neuroglancer.

# Arguments
- `idlists::Vector{Vector{Int}}`: Vector of cell ID lists, each shown in a separate segmentation layer
- `anchor::String`: Text to display for the hyperlink (default: "multi-layer view")
- `names::Vector{String}`: Optional custom names for each layer (default: auto-generated)
- `version::Int`: FlyWire data version (default: 783)
- `darkmode::Bool`: Use dark background in 3D view (default: false)

# Returns
- In Jupyter: HTML hyperlink object that can be displayed
- In REPL: Prints URL and returns `nothing`

# Examples
```julia
# Generate single link with multiple segmentation layers
ng_hyper([[id1, id2], [id3, id4]], anchor="Compare groups")

# Compare different cell types in layers with custom names
idlists = [type2ids("Tm1")[1:5], type2ids("Tm2")[1:5]]
ng_hyper(idlists, anchor="Tm1 vs Tm2", names=["Tm1", "Tm2"])
```
"""
function ng_hyper(idlists::Vector{Vector{Int}}; anchor = "multi-layer view", names::Union{Nothing,Vector{String}} = nothing, version = 783, darkmode::Bool = false)
    url = ng_url(idlists, version, names=names, darkmode=darkmode)

    # Check if we're in a notebook environment that can render HTML
    if isdefined(Main, :IJulia) && Main.IJulia.inited
        return HTML("<a href=$url>$anchor</a>")
    else
        # REPL or command line - just print the URL
        println("$anchor: $url")
        return nothing
    end
end
