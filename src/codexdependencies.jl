"""
    codexdependencies.jl

Registers data dependencies for automatic download via DataDeps.jl. Files are
automatically downloaded on first use, and are cached in the scratchspaces
folder.

This module defines the external data files required by OpticLobe.jl. The
files are downloaded from the FlyWire Codex.

The data includes:
- Cell statistics (morphological measurements)
- Neuron IDs (cell identifiers)
- Visual neuron types (cell type annotations)
- Classification data (cell type hierarchies)
- Synaptic connections, both new version from Princeton and original from Buhmann et al.
- Cell type annotations

All data is from FlyWire version 783 proofreading.
MD5 checksums are verified only for the data files that remain unchanged from the original v783 release.
"""

using DataDeps
using MD5    # switching to MD5 checksums, rather than DataDeps default of SHA256

DOWNLOADS = "https://storage.googleapis.com/flywire-data/codex/data/fafb/783"

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
));

register(DataDep("Codex Buhmann connections no threshold",
    """
    Dataset: FlyWire v783
    Author:
    License:
    Website: codex.flywire.ai

        Citation: Matsliah et al.

    This file is 202MB
    """,
    joinpath(DOWNLOADS, "connections_buhmann_no_threshold.csv.gz"),
                 (md5, "8332c0f03acd2419f56f9bbc002a1554"),
    ));

register(DataDep("Codex Princeton connections no threshold",
    """
    Dataset: FlyWire v783
    Author:
    License:
    Website: codex.flywire.ai

        Citation: Matsliah et al.

    This file is 202MB
    """,
    joinpath(DOWNLOADS, "connections_princeton_no_threshold.csv.gz"),
    ));

register(DataDep("Codex cell types",
    """
    Dataset: FlyWire v783
    Author:
    License:
    Website: codex.flywire.ai

        Citation: Matsliah et al.

    This file is 202MB
    """,
    joinpath(DOWNLOADS, "consolidated_cell_types.csv.gz"),
    ));
