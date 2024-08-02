# interface to FlyWire resources: Codex, neuroglancer

using DataStructures, JSON, URIs, MissingsAsFalse

"""
    codex_open(idlist[, version])

search the FlyWire Codex for the cells in `idlist`

if `version` is unspecified, use v783.
"""
function codex_open(idlist::Vector{Int}, version = 783)
    query = join(idlist, "+")
    url = "https://codex.flywire.ai/app/search?filter_string=$query&sort_by=&page_size=10&data_version=$version"
    urlquoted = "\"" * url * "\""
    display("text/javascript", """window.open($urlquoted)""")
end


"""
    ng_url(idlist[, version])

return a string containing a neuroglancer URL
"""
function ng_url(idlist::Vector{Int}, version = 783)
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
    jencoded = URIs.escapeuri(json(j))
#    return "https://neuroglancer-demo.appspot.com/#!$jencoded"
    return "https://ngl.cave-explorer.org/#!$jencoded"
end

######### open neuroglancer window

"""
    ng_open(idlist[, version])

open neuroglancer window showing cells in `idlist`
"""
function ng_open(idlist::Vector{Int}, version = 783)
    url = "\"" * ng_url(idlist, version) * "\""
    display("text/javascript", """window.open($url)""")
end

"""
    ng_open(celltype[, version])

open neuroglancer window showing all cells in `celltype`
"""
function ng_open(celltype::String, version = 783)
    @mfalse ng_open(ind2id[ind2type .== celltype], version)
end

######### create hyperlink

"""
    ng_hyper(idlist[; version = 783, anchor = join(idlist, ", ")])

create a hyperlink to a neuroglancer view of cells in `idlist`
Anchor text can be specified, and defaults to the cell IDs.
"""
function ng_hyper(idlist::Vector{Int}; anchor = join(idlist, ", "), version = 783)
    url = ng_url(idlist, version)
    return HTML("<a href=$url>$anchor</a>")
end

"""
    ng_hyper(celltype[; version = 783])

create a hyperlink to a neuroglancer view of cells in `celltype`
The anchor text is just the `celltype` string.
"""
function ng_hyper(celltype::String; version = 783)
    @mfalse idlist = ind2id[ind2type .== celltype]
    return ng_hyper(idlist, anchor = celltype, version = version)
end
