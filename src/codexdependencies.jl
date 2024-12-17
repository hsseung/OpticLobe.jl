using DataDeps
using MD5    # switching as GCS provides human-readable MD5 checksums. DataDeps default is SHA256

DOWNLOADS = "https://storage.googleapis.com/flywire-data/codex/data/fafb/783"

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
    joinpath(DOWNLOADS, "visual_neuron_types.csv.gz?generation=1720619193993712"),
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
    joinpath(DOWNLOADS, "classification.csv.gz?generation=1733428415941889"),
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
