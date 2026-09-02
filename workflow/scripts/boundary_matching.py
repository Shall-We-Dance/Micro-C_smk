"""Shared boundary-matching implementation for the qc_boundary_overlap rule.

Imported both by the production rule (workflow/rules/05_qc_plots.smk) and by
the unit tests (workflow/tests/test_boundary_matching.py) so that the
production algorithm and the tested algorithm can never drift apart.

Matching semantics: greedy one-to-one nearest-neighbor matching. For each a
midpoint in input order, find the nearest UNUSED b midpoint within
[m - tol, m + tol]: locate the first candidate >= m - tol with bisect and
scan forward while <= m + tol, keeping the closest unused candidate. Each b
midpoint is used at most once.
"""

import bisect


def match_boundaries(a_mids, b_mids, tol):
    """Greedy one-to-one nearest-neighbor matching.

    a_mids/b_mids: iterables of boundary midpoints (sorted for the greedy
    order to be stable).
    tol: maximum absolute distance for a match (inclusive).
    Returns the number of matched pairs.
    """
    b_mids = sorted(b_mids)
    used = [False] * len(b_mids)
    matched = 0
    for am in a_mids:
        pos = bisect.bisect_left(b_mids, am - tol)
        best = None
        while pos < len(b_mids) and b_mids[pos] <= am + tol:
            if not used[pos]:
                d = abs(b_mids[pos] - am)
                if best is None or d < best[0]:
                    best = (d, pos)
            pos += 1
        if best is not None:
            used[best[1]] = True
            matched += 1
    return matched


def boundary_metrics(a_mids, b_mids, tol):
    """Return (n_matched, precision, recall, f1, jaccard) for a vs b."""
    m = match_boundaries(a_mids, b_mids, tol)
    n1, n2 = len(a_mids), len(b_mids)
    precision = m / n1 if n1 else 0.0
    recall = m / n2 if n2 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
    denom = n1 + n2 - m
    jaccard = m / denom if denom else 0.0
    assert 0.0 <= jaccard <= 1.0, (
        "impossible jaccard "
        + str(jaccard)
        + " (m="
        + str(m)
        + ", n1="
        + str(n1)
        + ", n2="
        + str(n2)
        + ")"
    )
    return m, precision, recall, f1, jaccard
