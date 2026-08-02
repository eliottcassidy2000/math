#!/usr/bin/env python3
"""Independent audit of the H upper-median two-star finite certificate."""

from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import permutations
from pathlib import Path
from math import inf
import hashlib

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "04-computation/lrc14_j7_reflected_levels_all_q_mass_closure_thm2941.py"
H = (1, 2, 3, 4, 6, 12)
L = 168
J = 90
S = 10**15
LEAVES = (2, 3, 4, 5)
EDGES = ((0, 1), (0, 2), (0, 3), (0, 4), (0, 5),
         (1, 2), (1, 3), (1, 4), (1, 5))


def load():
    spec = spec_from_file_location("H_transverse_raw_interval_engine", BASE)
    if spec is None or spec.loader is None:
        raise RuntimeError(BASE)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def lower(x):
    return x.numerator * S // x.denominator


def upper(x):
    return -((-x.numerator * S) // x.denominator)


def max_matching_dp(levels, pivots, threshold, leaf_gain, debt):
    """Independent maximum-weight matching via level-by-level subset DP."""
    # Entries are (weight, assignments-by-leaf-index), with -1 unassigned.
    state = {0: (0, (-1, -1, -1, -1))}
    for q in levels:
        if q in pivots:
            continue
        nxt = dict(state)
        for mask, (weight, assignment) in state.items():
            for s, leaf in enumerate(LEAVES):
                bit = 1 << s
                if mask & bit or leaf_gain[leaf, q] > threshold:
                    continue
                candidate_assignment = list(assignment)
                candidate_assignment[s] = q
                candidate = (weight + debt[leaf, q], tuple(candidate_assignment))
                old = nxt.get(mask | bit)
                if old is None or candidate > old:
                    nxt[mask | bit] = candidate
        state = nxt
    return state.get(15)


def data(R, D, m):
    levels = tuple(range(m, m + D + 1))
    arcs = {(i, q): R.reflected_level_arcs(L, H[i], q, J)
            for i in range(6) for q in levels}
    exact = {}
    for i, j in EDGES:
        for p in levels:
            for q in levels:
                if p == q:
                    continue
                left, right = arcs[i, p], arcs[j, q]
                exact[i, j, p, q] = (
                    R.interval_mass(left) + R.interval_mass(right)
                    - R.interval_mass(R.merge_intervals(left + right))
                )
    gain = {key: lower(value) for key, value in exact.items()}
    exact_debt = {(i, q): Fraction(H[i], 7 * (q * L - H[i]))
                  for i in range(6) for q in levels}
    debt = {key: upper(value) for key, value in exact_debt.items()}
    return levels, exact, gain, exact_debt, debt


def certify_dp(R, D, m):
    levels, exact, gain, exact_debt, debt = data(R, D, m)
    overall = None
    matching_calls = 0
    feasible = 0
    for q0, q1 in permutations(levels, 2):
        pivot_gain = gain[0, 1, q0, q1]
        leaf_gain = {}
        thresholds = {pivot_gain}
        for leaf in LEAVES:
            for q in levels:
                if q in (q0, q1):
                    continue
                value = max(gain[0, leaf, q0, q], gain[1, leaf, q1, q])
                leaf_gain[leaf, q] = value
                thresholds.add(value)
        base_debt = debt[0, q0] + debt[1, q1]
        relaxed = base_debt + sum(
            max(debt[leaf, q] for q in levels if q not in (q0, q1))
            for leaf in LEAVES
        )
        for threshold in sorted(t for t in thresholds if t >= pivot_gain):
            if threshold > relaxed:
                row = (threshold - relaxed, (q0, q1), threshold, relaxed, "relaxed-tail")
                overall = row if overall is None or row < overall else overall
                break
            matching_calls += 1
            answer = max_matching_dp(levels, (q0, q1), threshold, leaf_gain, debt)
            if answer is None:
                continue
            feasible += 1
            weight, assignment = answer
            row = (threshold - base_debt - weight,
                   (q0, q1, *assignment), threshold, base_debt + weight, "dp-matching")
            overall = row if overall is None or row < overall else overall
    if overall is None or overall[0] <= 0:
        raise RuntimeError((D, m, overall))
    return overall, matching_calls, feasible


def main():
    print("LRC H UPPER-MEDIAN TWO-STAR MATCHING -- INDEPENDENT SUBSET-DP AUDIT")
    R = load()
    ruler, ranges = R.safe_cell_ranges(H)
    if ruler != L or not R.body_cell_is_safe(L, H, J):
        raise RuntimeError((ruler, ranges, J))
    rows = []
    for D, m in ((6, 6), (8, 8), (10, 10), (12, 12), (15, 45), (20, 20)):
        certificate, calls, feasible = certify_dp(R, D, m)
        row = (D, m, certificate, calls, feasible)
        rows.append(row)
        print(f"D={D};m={m};certificate={certificate};matching_calls={calls};feasible={feasible}")
    semantic = hashlib.sha256(repr(tuple(rows)).encode()).hexdigest()
    print(f"semantic_sha256={semantic}")
    print("independent_subset_DP_and_direct_interval_route=PASS")


if __name__ == "__main__":
    main()
