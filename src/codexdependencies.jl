"""
    codexdependencies.jl

Registers data dependencies for automatic download via DataDeps.jl. Files are
automatically downloaded on first use, and are cached in the scratchspaces
folder.

This module defines the external data files required by OpticLobe.jl. These
files were originally downloadable from the FlyWire Codex. Here they are being
downloaded from Zenodo, which archives old versions that are useful for reproducing
published papers.

The data includes:
- Cell statistics (morphological measurements)
- Neuron IDs (cell identifiers)
- Visual neuron types (cell type classifications)
- Classification data (cell type hierarchies)
- Synaptic connections (connectivity data, ~202MB)

All data is from FlyWire version 783 proofreading.
Uses MD5 checksums for verification as provided by the data source.
"""

using DataDeps
using MD5    # switching to MD5 checksums, rather than DataDeps default of SHA256

DOWNLOADS = "https://zenodo.org/records/17030629/files"

# this has never changed, so specifying version seems unnecessary
register(DataDep("Codex cell stats",
    """
    Dataset: FlyWire v783
    Author:
    License:
    Website: codex.flywire.ai

    Citation: Matsliah et al.
    """,
                 joinpath(DOWNLOADS, "cell_stats.csv.gz"),
                 (md5, "a79d4875e2b80e51d1269dcb86dd8743"),
                 )
         );


register(DataDep("Codex neuron IDs",
    """
    Dataset: FlyWire v783
    Author:
    License:
    Website: codex.flywire.ai

    Citation: Matsliah et al.
    """,
    joinpath(DOWNLOADS, "neurons.csv.gz"),
                 (md5, "f60333e9e4124160b9b203b1712a6f91"),
));

register(DataDep("Codex visual neuron types",
    """
    Dataset: FlyWire v783
    Author:
    License:
    Website: codex.flywire.ai

    Citation: Matsliah et al.
    """,
    joinpath(DOWNLOADS, "visual_neuron_types.csv.gz"),
                 (md5, "f1870fc7c2bbfe7b99b869fdb83cf7ff"),
));

register(DataDep("Codex classification",
    """
    Dataset: FlyWire v783
    Author:
    License:
    Website: codex.flywire.ai

    Citation: Matsliah et al.
    """,
    joinpath(DOWNLOADS, "classification.csv.gz"),
                 (md5, "6200e3f2016e113a829c92ecdeb4ac9c"),
));

register(DataDep("Codex connections no threshold",
    """
    Dataset: FlyWire v783
    Author:
    License:
    Website: codex.flywire.ai

        Citation: Matsliah et al.

    This file is 202MB
    """,
    joinpath(DOWNLOADS, "connections_no_threshold.csv.gz"),
                 (md5, "8332c0f03acd2419f56f9bbc002a1554"),
    ));
