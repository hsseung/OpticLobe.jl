using DataDeps

DOWNLOADS = "https://storage.googleapis.com/flywire-data/codex/data/783"

# this has never changed, so specifying version seems unnecessary
register(DataDep("Codex cell stats",
    """
    Dataset: FlyWire v783
    Author:
    License:
    Website: codex.flywire.ai

    Citation: Matsliah et al.
    """,
    joinpath(DOWNLOADS, "cell_stats.csv.gz?generation=1709735415919143"),
    "bd5879e1b5df964bea2f3ca5316348d4276ce2ccaac283f0e36583c04fbd3d8e"
    ));


register(DataDep("Codex neuron IDs",
    """
    Dataset: FlyWire v783
    Author:
    License:
    Website: codex.flywire.ai

    Citation: Matsliah et al.
    """,
    joinpath(DOWNLOADS, "neurons.csv.gz"),
    "6a6b3759e635f0f35a677d169052362131ec61d95f55919298b55c43fce4e719"
));

register(DataDep("Codex visual neuron types",
    """
    Dataset: FlyWire v783
    Author:
    License:
    Website: codex.flywire.ai

    Citation: Matsliah et al.
    """,
    joinpath(DOWNLOADS, "visual_neuron_types.csv.gz?generation=1720619193993712"),
    "4d6aee4af9f0aed3098a4eb91163508a61ded6363a3c441a7d50a5f6fcf5c366"
));

register(DataDep("Codex classification",
    """
    Dataset: FlyWire v783
    Author:
    License:
    Website: codex.flywire.ai

    Citation: Matsliah et al.
    """,
    joinpath(DOWNLOADS, "classification.csv.gz?generation=1721489294315289"),
    "cfd265c5b650df0caaf8b6d02386c9ec53f3c0029102430cf2763d5b40ee6f7d"
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
    "e40e028b47f9bfc9b6b5a8da534be307eaded1f716b9a39e2c4a586c9f7fc67e"
));
