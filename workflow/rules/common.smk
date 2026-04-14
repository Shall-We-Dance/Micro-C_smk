SAMPLES = sorted(list(config["samples"].keys()))
OUTDIR = config.get("output", {}).get("dir", "results")
REF = config.get("reference", {})
THREADS = config.get("threads", {})

PAIR_RESOLUTIONS = config.get("matrix", {}).get("resolutions", [1000, 2000, 5000, 10000])
PAIR_RES_CSV = ",".join(str(r) for r in PAIR_RESOLUTIONS)

FEATURES_CFG = config.get("features", {})
COMPARTMENT_RESOLUTIONS = FEATURES_CFG.get("compartment_resolutions", [5000])
if isinstance(COMPARTMENT_RESOLUTIONS, int):
    COMPARTMENT_RESOLUTIONS = [COMPARTMENT_RESOLUTIONS]
COMPARTMENT_RESOLUTIONS = sorted({int(r) for r in COMPARTMENT_RESOLUTIONS})

FILTER_CFG = config.get("pairs", {}).get("filter", {})
MIN_MAPQ = int(FILTER_CFG.get("min_mapq", 30))
MAX_CIS_DISTANCE_ARTIFACT = FILTER_CFG.get("max_cis_artifact_dist", None)
BLACKLIST_BED = FILTER_CFG.get("blacklist_bed", "")
BLACKLIST_ENABLED = bool(FILTER_CFG.get("enable_blacklist", False) and BLACKLIST_BED)

EXPORT_HIC = bool(config.get("matrix", {}).get("export_hic", False))
