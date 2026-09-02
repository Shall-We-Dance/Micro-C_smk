if len(GROUPS) >= 2:

    # Per-group merge rule. Runs once per group AFTER the QC gates exist
    # (qc_gates.tsv is a hard input dependency, P1-8) and decides at runtime
    # whether the group can be merged:
    #   1. protocol guard: all members must share the same protocol (P1-5);
    #   2. QC gate: only eligible samples (qc_gates.tsv status == "pass") are
    #      merged, unless DIFFERENTIAL_INCLUDE_FAILED overrides (warning is
    #      logged);
    #   3. mixed base resolutions: members whose base resolution is finer
    #      than the group common base (max over members) are coarsened first
    #      with `cooler coarsen --factor`, so `cooler merge` only ever sees
    #      identical bins (P1-5);
    #   4. zoomify only the resolutions every eligible member can produce
    #      (existing resolution-intersection logic).
    # The declared output is a per-group status note ({group}.skipped.note,
    # written in every case, content: reason + excluded samples). The merged
    # mcool is written BEFORE the note as a side product and only exists for
    # groups whose note starts with READY; differential_compare re-checks the
    # files and skips groups whose note is not READY.
    for _group in sorted(GROUPS):
        _members = sorted(GROUPS[_group])
        rule:
            input:
                cool=expand(f"{OUTDIR}/matrices/{{sample}}.cool", sample=_members),
                qc_gates=f"{OUTDIR}/report/qc_gates.tsv"
            output:
                note=f"{OUTDIR}/integrative/{_group}.skipped.note"
            params:
                group=_group,
                members=_members,
                member_protocols={s: sample_protocol(s) for s in _members},
                member_bases={s: sample_base_resolution(s) for s in _members},
                member_resolutions={s: sample_resolutions(s) for s in _members},
                include_failed=DIFFERENTIAL_INCLUDE_FAILED,
                outdir=OUTDIR
            conda:
                "envs/cooltools.yaml"
            log:
                f"logs/integrative/merge_{_group}.log"
            threads: ZOOMIFY_THREADS
            shell:
                r"""
                set -euo pipefail
                mkdir -p $(dirname {output.note}) $(dirname {log})
                python - <<'PY' > {log} 2>&1
import os
import shutil
import subprocess
import sys

OUTDIR = {params.outdir!r}
group = {params.group!r}
members = {params.members!r}
member_protocols = {params.member_protocols!r}
member_bases = {params.member_bases!r}
member_resolutions = {params.member_resolutions!r}
include_failed = {params.include_failed}
nproc = {threads}

note_path = "{output.note}"


def write_note(status, reason, excluded):
    with open(note_path, "w") as h:
        h.write(status + "\t" + reason + "\texcluded=" + ",".join(sorted(excluded)) + "\n")
    print(status + ": " + reason)


def read_qc_gates(path):
    eligible = set()
    if not os.path.isfile(path):
        return eligible
    with open(path) as h:
        header = h.readline().rstrip("\n").split("\t")
        si = header.index("status") if "status" in header else -1
        for line in h:
            row = line.rstrip("\n").split("\t")
            if si >= 0 and len(row) > si and row[si] == "pass":
                eligible.add(row[0])
    return eligible


# ---- P1-5 protocol guard: mixed-protocol groups are not allowed ----
protocols = sorted({{member_protocols[s] for s in members}})
if len(protocols) > 1:
    detail = ", ".join(s + "=" + member_protocols[s] for s in members)
    print("[ERROR] mixed-protocol group " + group + " not allowed (members: " + detail + ")", file=sys.stderr)
    sys.exit(1)

# ---- P1-8 QC gate hard dependency: only eligible samples are merged ----
if include_failed:
    eligible = list(members)
    print("[WARN] DIFFERENTIAL_INCLUDE_FAILED=true: ignoring QC gates, merging all group members")
else:
    eligible = [s for s in members if s in read_qc_gates(OUTDIR + "/report/qc_gates.tsv")]
excluded = [s for s in members if s not in eligible]

if len(eligible) < 2:
    write_note("SKIPPED", "fewer than 2 eligible samples in group; group not merged", excluded)
    sys.exit(0)

# ---- existing resolution-intersection logic, applied to eligible members ----
common_base = max(member_bases[s] for s in members)
common_res = None
for s in eligible:
    rs = set(member_resolutions[s])
    common_res = rs if common_res is None else common_res & rs
# after coarsening, the merged matrix has base resolution common_base, so
# zoomify can only build resolutions that are multiples of common_base
# (e.g. 5000 is not derivable from a 2000bp base)
common_res = sorted(r for r in (common_res or []) if r % common_base == 0)
if not common_res:
    write_note("SKIPPED", "no common zoomify resolution across eligible group members; group not merged", excluded)
    sys.exit(0)

# ---- P1-5 mixed-base-resolution merge: coarsen fine-base members ----
coarsen_dir = os.path.join(OUTDIR, "integrative", "tmp", group)
os.makedirs(coarsen_dir, exist_ok=True)
merge_inputs = []
for s in eligible:
    base = member_bases[s]
    cool_path = os.path.join(OUTDIR, "matrices", s + ".cool")
    if base < common_base:
        if common_base % base != 0:
            print("[ERROR] cannot coarsen sample " + s + " base " + str(base) + " to common base " + str(common_base) + "bp: not a multiple", file=sys.stderr)
            sys.exit(1)
        factor = common_base // base
        tmp = os.path.join(coarsen_dir, s + "." + str(common_base) + "bp.cool")
        print("[coarsen] " + s + ": " + str(base) + "bp -> " + str(common_base) + "bp (factor " + str(factor) + ")")
        subprocess.run(["cooler", "coarsen", "-k", str(factor), "-n", str(nproc), cool_path, "-o", tmp], check=True)
        merge_inputs.append(tmp)
    else:
        merge_inputs.append(cool_path)

mcool_path = os.path.join(OUTDIR, "integrative", group + ".mcool")
base_merged = mcool_path + ".base.cool"
print("[merge] merging " + str(len(merge_inputs)) + " cools")
# cooler 0.10.x: OUT_PATH is positional (no -o option), input paths follow.
subprocess.run(["cooler", "merge", base_merged] + merge_inputs, check=True)
print("[zoomify] resolutions: " + ",".join(str(r) for r in common_res))
subprocess.run(
    ["cooler", "zoomify", "--nproc", str(nproc), "--resolutions", ",".join(str(r) for r in common_res),
     "--balance", base_merged, "-o", mcool_path],
    check=True,
)

# cleanup temp coarsened files (and the intermediate merged cool)
try:
    os.removedirs(coarsen_dir)
except OSError:
    shutil.rmtree(coarsen_dir, ignore_errors=True)
if os.path.isfile(base_merged):
    os.remove(base_merged)

write_note("READY", "merged " + str(len(eligible)) + " samples at base " + str(common_base) + "bp; zoomify resolutions " + ",".join(str(r) for r in common_res), excluded)
PY
                """

    rule differential_compare:
        input:
            # per-group status notes (always exist; group mcools are side
            # products and are re-checked inside the python)
            group_notes=expand(f"{OUTDIR}/integrative/{{group}}.skipped.note", group=sorted(GROUPS)),
            qc_gates=f"{OUTDIR}/report/qc_gates.tsv",
            compartments=expand(
                f"{OUTDIR}/features/{{sample}}/compartments_{{res}}bp.bedgraph",
                sample=SAMPLES,
                res=COMPARTMENT_RESOLUTIONS,
            )
        output:
            matrix=f"{OUTDIR}/integrative/matrix_compare.tsv",
            flip=f"{OUTDIR}/integrative/compartment_flip.tsv"
        params:
            resolutions=CONCORDANCE_RESOLUTIONS,
            min_chrom_size=MIN_CHROM_SIZE,
            compartment_resolutions=COMPARTMENT_RESOLUTIONS,
            samples=SAMPLES,
            groups=GROUPS,
            flip_min_e1=FLIP_MIN_E1,
            include_failed=DIFFERENTIAL_INCLUDE_FAILED,
            outdir=OUTDIR
        conda:
            "envs/cooltools.yaml"
        log:
            f"logs/integrative/differential_compare.log"
        threads: 8
        shell:
            r"""
            set -euo pipefail
            mkdir -p $(dirname {output.matrix}) $(dirname {output.flip}) $(dirname {log})
            python - <<'PY' > {log} 2>&1
import os
import itertools
import cooler
import numpy as np
import pandas as pd
from scipy import stats

OUTDIR = {params.outdir!r}
SAMPLES = {params.samples!r}
GROUPS = {params.groups!r}
params_resolutions = {params.resolutions!r}
min_chrom_size = {params.min_chrom_size}
params_compartment_resolutions = {params.compartment_resolutions!r}
flip_min_e1 = {params.flip_min_e1}
include_failed = {params.include_failed}

group_names = sorted(GROUPS.keys())
mcool_dir = OUTDIR + "/integrative/"
gates_path = OUTDIR + "/report/qc_gates.tsv"

# ---- per-group status from the merge notes (READY vs SKIPPED) ----
group_status = dict()
for g in group_names:
    try:
        with open(mcool_dir + g + ".skipped.note") as h:
            group_status[g] = h.readline().split("\t", 1)[0]
    except OSError:
        group_status[g] = "MISSING"

# ---- P1-8 QC gate filtering: only eligible samples contribute ----
if include_failed:
    eligible = set(SAMPLES)
    print("[WARN] DIFFERENTIAL_INCLUDE_FAILED=true: ignoring QC gates")
else:
    eligible = set()
    if os.path.isfile(gates_path):
        with open(gates_path) as h:
            header = h.readline().rstrip("\n").split("\t")
            si = header.index("status") if "status" in header else -1
            for line in h:
                row = line.rstrip("\n").split("\t")
                if si >= 0 and len(row) > si and row[si] == "pass":
                    eligible.add(row[0])
group_members = sorted({{s for g in GROUPS.values() for s in g}})
excluded = [s for s in group_members if s not in eligible]
skipped_groups = [g for g in group_names if group_status.get(g) != "READY"]
print("group status: " + ", ".join(g + "=" + group_status.get(g) for g in group_names))
print("eligible samples: " + ",".join(sorted(eligible)))
print("excluded samples: " + ",".join(sorted(excluded)))

# ---- 1) genome-wide matrix comparison between group aggregates ----
with open("{output.matrix}", "w") as h:
    h.write("group1\tgroup2\tresolution\tpearson\tspearman\tn_bins\n")
for g1, g2 in itertools.combinations(group_names, 2):
    mc1 = mcool_dir + g1 + ".mcool"
    mc2 = mcool_dir + g2 + ".mcool"
    # Group mcools are optional side products: groups whose note is not READY
    # (skipped) produced no mcool, so their pairs are skipped here.
    if group_status.get(g1) != "READY" or group_status.get(g2) != "READY":
        print("skip pair " + g1 + " vs " + g2 + ": group skipped")
        continue
    if not os.path.isfile(mc1) or not os.path.isfile(mc2):
        print("skip pair " + g1 + " vs " + g2 + ": group mcool missing")
        continue
    # Iterate the intersection of resolutions available in BOTH group mcools
    # and the configured CONCORDANCE_RESOLUTIONS (unavailable combos are
    # reported as NA rows below).
    avail1 = {{int(r.rsplit("/", 1)[-1]) for r in cooler.fileops.list_coolers(mc1)}}
    avail2 = {{int(r.rsplit("/", 1)[-1]) for r in cooler.fileops.list_coolers(mc2)}}
    for res in params_resolutions:
        if res not in avail1 or res not in avail2:
            with open("{output.matrix}", "a") as h:
                h.write(g1 + "\t" + g2 + "\t" + str(res) + "\tNA\tNA\t0\n")
            continue
        clr1 = cooler.Cooler(mc1 + "::/resolutions/" + str(res))
        clr2 = cooler.Cooler(mc2 + "::/resolutions/" + str(res))
        chroms = [
            c
            for c, s in zip(clr1.chromnames, clr1.chromsizes)
            if s >= min_chrom_size
        ]
        # Group mcools carry zoomify --balance weights; use them only when the
        # weight column actually exists, otherwise fall back to raw counts.
        use_balance = "weight" in clr1.bins().columns and "weight" in clr2.bins().columns
        v1, v2 = [], []
        for chrom in chroms:
            m1 = clr1.matrix(balance=use_balance).fetch(chrom)
            m2 = clr2.matrix(balance=use_balance).fetch(chrom)
            t1 = m1[np.triu_indices_from(m1, k=1)]
            t2 = m2[np.triu_indices_from(m2, k=1)]
            finite = np.isfinite(t1) & np.isfinite(t2)
            v1.append(t1[finite])
            v2.append(t2[finite])
        if not v1:
            with open("{output.matrix}", "a") as h:
                h.write(g1 + "\t" + g2 + "\t" + str(res) + "\tNA\tNA\t0\n")
            continue
        a = np.concatenate(v1)
        b = np.concatenate(v2)
        if len(a) < 2:  # pearsonr needs >= 2 finite bins
            with open("{output.matrix}", "a") as h:
                h.write(g1 + "\t" + g2 + "\t" + str(res) + "\tNA\tNA\t" + str(len(a)) + "\n")
            continue
        pearson = stats.pearsonr(a, b).statistic
        spearman = stats.spearmanr(a, b).statistic
        with open("{output.matrix}", "a") as h:
            h.write(g1 + "\t" + g2 + "\t" + str(res) + "\t" + format(pearson, ".4f") + "\t" + format(spearman, ".4f") + "\t" + str(len(a)) + "\n")

# ---- 2) per-chromosome compartment (E1) flip between groups ----
with open("{output.flip}", "w") as h:
    h.write("group1\tgroup2\tresolution\tchrom\tpearson\tflip_fraction\tn_bins_total\tn_bins_used\n")
sample_groups = {{s: g for g, ss in GROUPS.items() for s in ss}}
for g1, g2 in itertools.combinations(group_names, 2):
    if group_status.get(g1) != "READY" or group_status.get(g2) != "READY":
        continue
    for res in params_compartment_resolutions:
        tracks = {{}}
        for sample in SAMPLES:
            if sample_groups.get(sample) not in (g1, g2) or sample not in eligible:
                continue
            path = OUTDIR + "/features/" + sample + "/compartments_" + str(res) + "bp.bedgraph"
            if not os.path.isfile(path) or os.path.getsize(path) == 0:
                continue
            df = pd.read_csv(path, sep="\t", header=None, names=["chrom", "start", "end", "E1"])
            df = df.dropna(subset=["E1"])
            tracks[sample] = df.set_index(["chrom", "start"])["E1"]
        if not tracks:
            continue
        samples_a = [s for s in SAMPLES if sample_groups.get(s) == g1 and s in tracks]
        samples_b = [s for s in SAMPLES if sample_groups.get(s) == g2 and s in tracks]
        if not samples_a or not samples_b:
            continue
        mean_a = pd.concat([tracks[s] for s in samples_a], axis=1).mean(axis=1)
        mean_b = pd.concat([tracks[s] for s in samples_b], axis=1).mean(axis=1)
        table = pd.DataFrame({{"A": mean_a, "B": mean_b}}).dropna()
        for chrom, grp in table.groupby(level="chrom"):
            if len(grp) < 10:
                continue
            va = grp["A"].values
            vb = grp["B"].values
            pearson = stats.pearsonr(va, vb).statistic
            # Only bins with a meaningful E1 magnitude (|E1| > flip_min_e1 in
            # at least one group-mean track) can be a compartment flip:
            # near-zero E1 is noise and must not count as a sign change.
            used = (np.abs(va) > flip_min_e1) | (np.abs(vb) > flip_min_e1)
            n_bins_total = int(len(grp))
            n_bins_used = int(np.count_nonzero(used))
            if n_bins_used > 0:
                flip_frac = float(np.mean(np.sign(va[used]) != np.sign(vb[used])))
                flip_frac_str = format(flip_frac, ".4f")
            else:
                flip_frac_str = "NA"
            with open("{output.flip}", "a") as h:
                h.write(g1 + "\t" + g2 + "\t" + str(res) + "\t" + str(chrom) + "\t" + format(pearson, ".4f") + "\t" + flip_frac_str + "\t" + str(n_bins_total) + "\t" + str(n_bins_used) + "\n")

# ---- notes section: excluded samples / skipped groups (P1-8) ----
notes_lines = [
    "#notes\texcluded_samples\t" + ",".join(sorted(excluded)),
    "#notes\tskipped_groups\t" + ",".join(sorted(skipped_groups)),
]
with open("{output.matrix}", "a") as h:
    for ln in notes_lines:
        h.write(ln + "\n")
with open("{output.flip}", "a") as h:
    for ln in notes_lines:
        h.write(ln + "\n")

print("differential compare done for groups: " + ", ".join(group_names))
PY
            """

    rule differential_summary:
        input:
            matrix=f"{OUTDIR}/integrative/matrix_compare.tsv",
            flip=f"{OUTDIR}/integrative/compartment_flip.tsv",
            group_notes=expand(f"{OUTDIR}/integrative/{{group}}.skipped.note", group=sorted(GROUPS)),
            qc_gates=f"{OUTDIR}/report/qc_gates.tsv"
        output:
            summary=f"{OUTDIR}/integrative/differential_summary.tsv"
        run:
            import os

            os.makedirs(os.path.dirname(output.summary), exist_ok=True)

            def _data_rows(path):
                if not os.path.isfile(path) or os.path.getsize(path) == 0:
                    return 0
                n = 0
                for line in open(path):
                    if line.startswith("#"):
                        continue
                    n += 1
                return max(n - 1, 0)  # minus the header row

            n_matrix = _data_rows(input.matrix)
            n_flip = _data_rows(input.flip)

            notes = {}
            for p in input.group_notes:
                g = os.path.basename(p)[: -len(".skipped.note")]
                try:
                    with open(p) as h:
                        notes[g] = h.read().strip()
                except OSError:
                    notes[g] = "MISSING note file"

            if DIFFERENTIAL_INCLUDE_FAILED:
                eligible = set(s for ss in GROUPS.values() for s in ss)
            else:
                eligible = set()
                if os.path.isfile(input.qc_gates):
                    with open(input.qc_gates) as h:
                        header = h.readline().rstrip("\n").split("\t")
                        si = header.index("status") if "status" in header else -1
                        for line in h:
                            row = line.rstrip("\n").split("\t")
                            if si >= 0 and len(row) > si and row[si] == "pass":
                                eligible.add(row[0])
            group_samples = sorted(s for ss in GROUPS.values() for s in ss)
            excluded = [s for s in group_samples if s not in eligible]

            with open(output.summary, "w") as handle:
                handle.write("analysis_module\tstatus\tnote\n")
                handle.write(f"matrix_compare\tready\t{n_matrix} group-pair/resolution combinations; aggregate correlation, not replicate-aware differential test\n")
                handle.write(f"compartment_flip\tready\t{n_flip} group-pair/chromosome combinations; aggregate correlation, not replicate-aware differential test\n")
                for g in sorted(notes):
                    fields = notes[g].split("\t")
                    status = "skipped" if fields[0] == "SKIPPED" else "ready"
                    rest = "\t".join(fields[1:])
                    handle.write(f"group_merge/{g}\t{status}\t{rest}\n")
                handle.write(f"differential_eligible\tready\t{len(eligible)} of {len(group_samples)} group samples pass QC gates (excluded: {','.join(excluded)})\n")
                handle.write("boundaries\tready\tPer-sample insulation boundaries available\n")
                handle.write("loops\tready\tPer-sample dots/loop calls available\n")
                handle.write("integration\tplanned\tJoin with ATAC/RNA/ChIP features in downstream notebooks\n")

else:

    rule differential_summary:
        output:
            summary=f"{OUTDIR}/integrative/differential_summary.tsv"
        run:
            import os

            os.makedirs(os.path.dirname(output.summary), exist_ok=True)
            with open(output.summary, "w") as handle:
                handle.write("analysis_module\tstatus\tnote\n")
                handle.write("compartments\tready\tPer-sample E1 tracks available\n")
                handle.write("boundaries\tready\tPer-sample insulation boundaries available\n")
                handle.write("loops\tready\tPer-sample dots/loop calls available\n")
                handle.write("differential\tskipped\tDefine >=2 groups in config['groups'] to enable group comparisons\n")
                handle.write("integration\tplanned\tJoin with ATAC/RNA/ChIP features in downstream notebooks\n")
