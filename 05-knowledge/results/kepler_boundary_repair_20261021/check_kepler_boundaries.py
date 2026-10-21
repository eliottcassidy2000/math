#!/usr/bin/env python3
# Hardened via AST Assert -> If(not test): raise AssertionError(message).
# Original script SHA-256: df005e84d20f6e5b794f934cb691ed5a7ff712f9cb87a19613fcf565b39f80bb
"""Exact scoped hostile checks for THM-4737; no floating-point geometry."""
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement
from pathlib import Path
import json

def distance(x):
    x %= 1
    return min(x, 1 - x)

def rate_signs(n, x):
    x %= 1
    nx = n * x
    right = nx.numerator // nx.denominator
    left = (right - 1) % n if nx.denominator == 1 else right
    return (1 <= left <= n - 2, 1 <= right <= n - 2)

def value(speeds, t):
    return min((distance(v * t) for v in speeds))

def feasible_intervals(speeds, threshold):
    intervals = [(F(0), F(1))]
    for v in speeds:
        cells = [(F(m, v) + threshold / v, F(m + 1, v) - threshold / v) for m in range(v)]
        intervals = [(max(a, c), min(b, d)) for a, b in intervals for c, d in cells if max(a, c) <= min(b, d)]
    return intervals

def maximum_by_full_pl_envelope(speeds):
    kinks = sorted({F(j, 2 * v) for v in speeds for j in range(2 * v + 1)})
    candidates = set(kinks)
    for a, b in zip(kinks, kinks[1:]):
        mid = (a + b) / 2
        lines = []
        for v in speeds:
            phase = v * mid % 1
            slope = v if phase < F(1, 2) else -v
            lines.append((slope, distance(v * mid) - slope * mid))
        for (s, c), (u, d) in combinations(lines, 2):
            if s != u:
                t = (d - c) / (s - u)
                if a <= t <= b:
                    candidates.add(t)
    best = max((value(speeds, t) for t in candidates))
    return best

def maximum_by_pair_sums(speeds):
    candidates = {F(m, a + b) for a, b in combinations_with_replacement(speeds, 2) for m in range(a + b + 1)}
    return max(((value(speeds, t), t) for t in candidates))
phase_checks = 0
old_mismatches = []
for n in range(3, 33):
    for j in range(2 * n):
        x = F(j, 2 * n)
        left, right = rate_signs(n, x)
        lonely = distance(x) >= F(1, n)
        if not lonely == (left or right):
            raise AssertionError((n, x))
        if lonely != right:
            old_mismatches.append((n, str(x)))
            if not x == F(n - 1, n):
                raise AssertionError()
        phase_checks += 1
if not len(old_mismatches) == 30:
    raise AssertionError()
speeds = (1, 2)
t = F(1, 3)
phases = [v * t % 1 for v in speeds]
if not value(speeds, t) == F(1, 3):
    raise AssertionError()
if not not all((rate_signs(3, x)[1] for x in phases)):
    raise AssertionError()
if not all((any(rate_signs(3, x)) for x in phases)):
    raise AssertionError()
if not feasible_intervals(speeds, F(1, 3)) == [(F(1, 3), F(1, 3)), (F(2, 3), F(2, 3))]:
    raise AssertionError()
if not feasible_intervals(speeds, F(2, 5)) == []:
    raise AssertionError()
row_checks = 0
for k in range(2, 5):
    for speeds in combinations(range(1, 9), k):
        best, witness = maximum_by_pair_sums(speeds)
        if not best == maximum_by_full_pl_envelope(speeds):
            raise AssertionError(speeds)
        if not witness.denominator <= 2 * max(speeds):
            raise AssertionError()
        threshold = F(1, k + 1)
        intervals = feasible_intervals(speeds, threshold)
        if not bool(intervals) == (best >= threshold):
            raise AssertionError()
        if intervals:
            if not value(speeds, intervals[0][0]) >= threshold:
                raise AssertionError()
            if not value(speeds, witness) >= threshold:
                raise AssertionError()
        row_checks += 1
result = {'status': 'PASS (finite exact checks; all-n proof supplied in note)', 'phase_partition': 'n=3..32; every vertex and every open-side midpoint', 'phase_checks': phase_checks, 'old_right_rate_endpoint_mismatches': len(old_mismatches), 'minimal_counterexample': {'n': 3, 'D': [1, 2], 't': '1/3', 'phases': ['1/3', '2/3'], 'lonely': True, 'all_right_rates_positive': False}, 'exact_tight_lonely_set_D_1_2_modulo_1': ['1/3', '2/3'], 'maximizer_checks': row_checks, 'maximizer_universe': 'all distinct positive integer k-sets in 1..8, k=2..4', 'independent_paths': 'complete piecewise-linear envelope versus pair-sum candidate grid', 'nonattainment_control': 'D={1,2}, target 2/5 gives an empty rational closed-time set', 'scope': 'No LRC proof, no general real-flow nonattainment counterexample, no repository changes.'}
print(json.dumps(result, indent=2))
Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2) + '\n')
