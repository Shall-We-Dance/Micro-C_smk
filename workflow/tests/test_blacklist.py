#!/usr/bin/env python3
"""Unit tests for workflow/scripts/blacklist.py (review P1-4).

Runnable with plain python:  python3 workflow/tests/test_blacklist.py

The production rule workflow/rules/03_pairs.smk (filter_pairs) merges BED
intervals per chromosome with merge_intervals() BEFORE filtering, so nested
or overlapping intervals can no longer shadow each other: the reported bug
was that with [0,100) and [50,60), position 80 was NOT detected as
blacklisted.
"""

import os
import random
import sys

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "scripts"
))
import blacklist


def naive_is_blacklisted(intervals, pos):
    return any(start <= pos < end for start, end in intervals)


def test_overlapping_merge():
    # review P1-4 regression case: [0,100) with nested [50,60); pos 80 must
    # be caught even though 80 lies inside the inner interval only
    merged = blacklist.merge_intervals([(0, 100), (50, 60)])
    assert merged == [(0, 100)], merged
    assert blacklist.is_blacklisted(merged, 80)
    assert blacklist.is_blacklisted(merged, 0)
    assert blacklist.is_blacklisted(merged, 55)
    assert blacklist.is_blacklisted(merged, 99)
    assert not blacklist.is_blacklisted(merged, 100)
    print("overlapping-merge test OK")


def test_nested_merge():
    merged = blacklist.merge_intervals([(100, 200), (120, 150), (130, 140), (180, 250)])
    assert merged == [(100, 250)], merged
    assert blacklist.is_blacklisted(merged, 130)
    assert blacklist.is_blacklisted(merged, 240)
    print("nested-merge test OK")


def test_adjacent_stay_separate():
    # adjacent intervals (end == start) do not overlap under half-open
    # semantics: they must stay separate yet both match correctly
    merged = blacklist.merge_intervals([(10, 20), (20, 30)])
    assert merged == [(10, 20), (20, 30)], merged
    assert blacklist.is_blacklisted(merged, 15)
    assert blacklist.is_blacklisted(merged, 25)
    # pos 20 is the exclusive end of [10,20) but the inclusive start of
    # [20,30): the union covers it, and the adjacent intervals must stay
    # separate for the merge to be correct
    assert blacklist.is_blacklisted(merged, 20)
    assert not blacklist.is_blacklisted(merged, 30)  # past both intervals
    print("adjacent-interval test OK")


def test_boundary_positions():
    merged = blacklist.merge_intervals([(5, 6)])
    assert blacklist.is_blacklisted(merged, 5)    # start inclusive
    assert not blacklist.is_blacklisted(merged, 6)  # end exclusive
    assert not blacklist.is_blacklisted(merged, 4)
    assert not blacklist.is_blacklisted(merged, 7)
    print("boundary test OK")


def test_large_merged_middle():
    # a position deep in the middle of a large merged interval is caught
    merged = blacklist.merge_intervals([(1_000_000, 9_000_000), (2_000_000, 3_000_000)])
    assert merged == [(1_000_000, 9_000_000)], merged
    assert blacklist.is_blacklisted(merged, 4_999_999)
    assert not blacklist.is_blacklisted(merged, 9_000_000)
    print("large-merged-middle test OK")


def test_fuzz_vs_naive():
    rng = random.Random(20260812)
    for case in range(1000):
        intervals = []
        for _ in range(rng.randint(0, 25)):
            start = rng.randint(0, 100_000)
            end = start + rng.randint(1, 5_000)
            intervals.append((start, end))
        merged = blacklist.merge_intervals(intervals)
        # merged output must be sorted and disjoint
        assert all(merged[i][1] <= merged[i + 1][0] for i in range(len(merged) - 1)), merged
        for _ in range(10):
            pos = rng.randint(-10, 110_000)
            got = blacklist.is_blacklisted(merged, pos)
            expected = naive_is_blacklisted(intervals, pos)
            assert got == expected, (case, intervals, merged, pos, got, expected)
    print("fuzz-vs-naive test OK (1000 cases)")


if __name__ == "__main__":
    test_overlapping_merge()
    test_nested_merge()
    test_adjacent_stay_separate()
    test_boundary_positions()
    test_large_merged_middle()
    test_fuzz_vs_naive()
    print("ALL BLACKLIST TESTS PASSED")
