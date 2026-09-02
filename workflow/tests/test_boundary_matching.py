#!/usr/bin/env python3
"""Unit tests for the shared boundary-matching module
workflow/scripts/boundary_matching.py, which is imported by the production
rule qc_boundary_overlap in workflow/rules/05_qc_plots.smk so that production
code and tests can never drift apart.

Runnable with plain python:  python3 workflow/tests/test_boundary_matching.py
"""

import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../scripts")
from boundary_matching import boundary_metrics, match_boundaries

TOL = 10000


def test_known_case():
    # a=[100,200,300], b=[105,310], tol=10k: greedy nearest matching gives
    # 100<->105 and 200<->310; 300 has no unused b within tol -> m=2 ->
    # P=2/3, R=1.0, F1=0.8, J=2/3
    m, p, r, f1, j = boundary_metrics([100, 200, 300], [105, 310], TOL)
    assert m == 2, m
    assert abs(p - 2 / 3) < 1e-9, p
    assert abs(r - 1.0) < 1e-9, r
    assert abs(f1 - 0.8) < 1e-9, f1
    assert abs(j - 2 / 3) < 1e-9, j
    # the old broken scenario that produced J>1 (multi-to-one was counted)
    m, p, r, f1, j = boundary_metrics([100, 200, 300], [105, 305], TOL)
    assert 0.0 <= j <= 1.0, j
    print("known-case test OK")


def test_one_to_one_no_overcount():
    # one b boundary cannot be used twice
    m, p, r, f1, j = boundary_metrics([100, 200], [150], TOL)
    assert m == 1 and p == 0.5 and r == 1.0 and abs(j - 1 / 2) < 1e-9, (m, p, r, j)
    m, p, r, f1, j = boundary_metrics([100, 200, 300], [150, 250, 350], TOL)
    assert m == 3 and p == 1.0 and r == 1.0 and j == 1.0, (m, p, r, j)
    print("one-to-one test OK")


def test_crowded_candidates():
    # many b midpoints inside the tolerance: greedy must pick the NEAREST
    # unused one, not merely the first found in scan order. With first-unused
    # semantics this case matches twice (1000->900, 1150->1001); with
    # nearest-unused it matches only once (1000->1001, then 1150 can reach
    # neither 900 (too far) nor the used 1001).
    m = match_boundaries([1000, 1150], [900, 1001], 200)
    assert m == 1, m
    # dense cluster: every a finds a distinct nearest b
    m = match_boundaries([1000, 1001], [900, 1001], 200)
    assert m == 2, m
    m, p, r, f1, j = boundary_metrics([1000], [900, 950, 1000, 1001, 1100], 100)
    assert m == 1 and 0.0 <= j <= 1.0, (m, j)
    print("crowded-candidates test OK")


def test_duplicate_coordinates():
    # duplicate midpoints on either side must still be one-to-one
    m = match_boundaries([100, 100], [100], 100)
    assert m == 1, m
    m = match_boundaries([100, 100], [100, 100], 100)
    assert m == 2, m
    m = match_boundaries([100], [100, 100], 100)
    assert m == 1, m
    m, p, r, f1, j = boundary_metrics([100, 100], [100], 100)
    assert m == 1 and 0.0 <= j <= 1.0, (m, j)
    print("duplicate-coordinates test OK")


def test_tolerance_boundary():
    # m + tol and m - tol are inclusive: exactly at the edge matches, one
    # past the edge does not
    assert match_boundaries([1000], [1200], 200) == 1
    assert match_boundaries([1000], [1201], 200) == 0
    assert match_boundaries([1000], [800], 200) == 1
    assert match_boundaries([1000], [799], 200) == 0
    print("tolerance-boundary test OK")


def test_symmetry_sanity():
    # a->b greedy matching is direction dependent: record both directions but
    # do NOT require them to agree; both Jaccards must stay in [0, 1].
    cases = [
        ([100, 200, 300], [150, 250]),
        ([1000, 1150], [900, 1001]),
        ([10, 20, 30, 40], [15, 45]),
        ([5, 55], [50]),
        ([0, 100, 200], [90]),
    ]
    for a, b in cases:
        m_ab, _, _, _, j_ab = boundary_metrics(a, b, 30)
        m_ba, _, _, _, j_ba = boundary_metrics(b, a, 30)
        assert 0.0 <= j_ab <= 1.0 and 0.0 <= j_ba <= 1.0, (a, b, j_ab, j_ba)
    print("symmetry-sanity test OK")


def test_jaccard_range_random():
    rng = random.Random(42)
    for _ in range(2000):
        n1 = rng.randint(0, 40)
        n2 = rng.randint(0, 40)
        a = sorted(rng.randint(0, 5_000_000) for _ in range(n1))
        b = sorted(rng.randint(0, 5_000_000) for _ in range(n2))
        tol = rng.choice([0, 100, 1000, 10000, 100000])
        m, p, r, f1, j = boundary_metrics(a, b, tol)
        assert 0.0 <= j <= 1.0, (a, b, tol, j)
        assert 0.0 <= p <= 1.0 and 0.0 <= r <= 1.0 and 0.0 <= f1 <= 1.0
    print("random Jaccard-range test OK (2000 cases)")


if __name__ == "__main__":
    test_known_case()
    test_one_to_one_no_overcount()
    test_crowded_candidates()
    test_duplicate_coordinates()
    test_tolerance_boundary()
    test_symmetry_sanity()
    test_jaccard_range_random()
    print("ALL BOUNDARY-MATCHING TESTS PASSED")
