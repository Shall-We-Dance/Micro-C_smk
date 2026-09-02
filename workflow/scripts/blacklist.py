#!/usr/bin/env python3
"""Shared blacklist interval helpers for Micro-C/Hi-C pair filtering
(review P1-4).

Intervals are 0-based half-open [start, end) coordinates, matching BED
fields 2/3 and pairtools' 0-based pair positions.

The old filter_pairs implementation sorted starts/ends separately and only
inspected the last interval with start <= pos, so positions covered by
nested/overlapping intervals were missed (e.g. [0,100) and [50,60): pos 80
was reported as not blacklisted). The production rule now merges intervals
per chromosome with merge_intervals() BEFORE filtering and looks positions
up with is_blacklisted().
"""

import bisect


def merge_intervals(intervals):
    """Merge an iterable of (start, end) 0-based half-open intervals into a
    sorted list of disjoint intervals.

    Overlapping and nested intervals are collapsed into their union;
    adjacent intervals (end == next start) do NOT overlap under half-open
    semantics and are kept separate.
    """
    merged = []
    for start, end in sorted(intervals):
        if merged and start < merged[-1][1]:
            if end > merged[-1][1]:
                merged[-1] = (merged[-1][0], end)
        else:
            merged.append((start, end))
    return merged


def is_blacklisted(merged, pos):
    """Return True if pos falls in any merged [start, end) interval.

    Boundary semantics: start is inclusive (pos == start -> True), end is
    exclusive (pos == end -> False). O(log n) via bisect on the sorted
    interval starts.
    """
    if not merged:
        return False
    starts = [iv[0] for iv in merged]
    i = bisect.bisect_right(starts, pos) - 1
    return i >= 0 and pos < merged[i][1]
