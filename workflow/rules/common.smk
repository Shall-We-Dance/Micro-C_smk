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

FEATURE_RESOLUTION = int(FEATURES_CFG.get("feature_resolution", config.get("matrix", {}).get("base_resolution", 1000)))

# Contigs shorter than this are excluded from expected/insulation/QC analysis.
# dm6 chrom.sizes contains ~1860 tiny random/Un contigs that break cooltools
# (empty-bin IndexError in insulation); a >=1Mb cutoff keeps exactly the 7
# canonical chromosomes.
MIN_CHROM_SIZE = int(config.get("matrix", {}).get("min_chrom_size", 1_000_000))

FILTER_CFG = config.get("pairs", {}).get("filter", {})
MIN_MAPQ = int(FILTER_CFG.get("min_mapq", 30))
MAX_CIS_DISTANCE_ARTIFACT = FILTER_CFG.get("max_cis_artifact_dist", 0)
SAME_FRAGMENT_MAX_DIST = int(FILTER_CFG.get("same_fragment_max_dist", 0))
REQUIRE_UNIQUE = bool(FILTER_CFG.get("require_unique", False))
BLACKLIST_BED = FILTER_CFG.get("blacklist_bed", "")
BLACKLIST_ENABLED = bool(FILTER_CFG.get("enable_blacklist", False) and BLACKLIST_BED)

DEDUP_MAX_MISMATCH = int(config.get("pairs", {}).get("dedup", {}).get("max_mismatch", 1))

# ---- per-sample (protocol-aware) settings -------------------------------
# samples.<name> may declare:
#   protocol: micro-c | hi-c
#   enzyme / cut_sites / min_frag_size / max_frag_size  (Hi-C metadata,
#   informational + used for protocol defaults)
#   base_resolution / resolutions        (matrix overrides)
#   max_cis_artifact_dist                (Hi-C short-cis/self-ligation filter)
#   same_fragment_max_dist               (unligated/self-ligated +- fragments)
# Protocol defaults for Hi-C: 5kb base resolution, 1kb short-cis artifact
# cutoff (MboI/DpnII 4-cutter fragments are mostly <1kb, so cis pairs below
# that are dominated by self-ligation products, cf. HiC-Pro/distiller).
# Micro-C: global base (1kb) and no distance filter by default.

def _sample_cfg(name):
    return config["samples"].get(name, {})


def sample_protocol(name):
    return str(_sample_cfg(name).get("protocol", "micro-c")).lower()


def sample_base_resolution(name):
    if "base_resolution" in _sample_cfg(name):
        return int(_sample_cfg(name)["base_resolution"])
    if sample_protocol(name) == "hi-c":
        return int(config.get("matrix", {}).get("hic_base_resolution", 5000))
    return int(config.get("matrix", {}).get("base_resolution", 1000))


def sample_resolutions(name):
    base = sample_base_resolution(name)
    res = _sample_cfg(name).get(
        "resolutions", config.get("matrix", {}).get("resolutions", [])
    )
    return sorted(
        {int(r) for r in res if int(r) >= base and int(r) % base == 0}
    )


def sample_resolutions_csv(name):
    return ",".join(str(r) for r in sample_resolutions(name))


def sample_max_cis_artifact_dist(name):
    if "max_cis_artifact_dist" in _sample_cfg(name):
        return int(_sample_cfg(name)["max_cis_artifact_dist"])
    if sample_protocol(name) == "hi-c":
        return int(config.get("matrix", {}).get("hic_cis_artifact_dist", 1000))
    return int(MAX_CIS_DISTANCE_ARTIFACT)


def sample_same_fragment_max_dist(name):
    if "same_fragment_max_dist" in _sample_cfg(name):
        return int(_sample_cfg(name)["same_fragment_max_dist"])
    return int(SAME_FRAGMENT_MAX_DIST)


def sample_min_mapq(name):
    return int(_sample_cfg(name).get("min_mapq", MIN_MAPQ))


# Number of read-name buckets for the parallel pairtools parse step
# (split name-sorted bam by read-name hash -> parse buckets in parallel ->
# concatenate). Each bucket runs its own parse process, so total parse
# throughput scales ~linearly with bucket count (subject to --cores).
PARSE_BUCKETS = int(config.get("pairs", {}).get("parse_buckets", 8))

ASSEMBLY = str(config.get("matrix", {}).get("assembly", "")).strip()
ZOOMIFY_THREADS = int(config.get("matrix", {}).get("zoomify_threads", THREADS.get("cooler", 8)))

CONCORDANCE_RESOLUTIONS = sorted({int(r) for r in config.get("qc", {}).get("concordance_resolutions", [100000, 500000])})

APA_FLANK = int(FEATURES_CFG.get("apa_flank", 100000))

GROUPS = config.get("groups", {})
GROUPS = {str(g): list(s) for g, s in GROUPS.items() if s}
if GROUPS:
    flat = [s for g in GROUPS.values() for s in g]
    assert len(flat) == len(set(flat)), "samples must not appear in more than one group"
    assert set(flat) == set(SAMPLES), "groups must cover exactly the configured samples"

EXPORT_HIC = bool(config.get("matrix", {}).get("export_hic", False))
HIC_EXPORT_THREADS = int(config.get("matrix", {}).get("hic_threads", THREADS.get("cooler", 8)))
HIC_EXPORT_TMPDIR = str(config.get("matrix", {}).get("hic_tmpdir", "")).strip()
