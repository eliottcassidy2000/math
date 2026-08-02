#!/usr/bin/env python3
"""FINITE-EXACT upper-median two-star scout for the reflected body H.

For each (D,m), levels are six distinct integers in [m,m+D].  The selected
cell is the upper median j=90 and the allowed certificate edges are the
two-star incident to labels 0 and 1.  The minimum of

    max_edge intersection_mass - exact_singleton_debt

is computed without enumerating all six-tuples: after fixing the two pivot
levels, thresholding the four leaf gains reduces the hostile search to a
maximum-weight matching with four left vertices.
"""

from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import permutations, product
from pathlib import Path
import argparse
import hashlib

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "04-computation/lrc14_j7_reflected_one_cone_closure_thm2941.py"
H = (1, 2, 3, 4, 6, 12)
L = 168
CELL = 90
SCALE = 10**15
PIVOTS = (0, 1)
LEAVES = (2, 3, 4, 5)
EDGES = tuple((i, j) for i in range(6) for j in range(i + 1, 6)
              if i in PIVOTS or j in PIVOTS)


def load():
    spec = spec_from_file_location("H_transverse_one_cone", SOURCE)
    if spec is None or spec.loader is None:
        raise RuntimeError(SOURCE)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def floor_scaled(value):
    return (value.numerator * SCALE) // value.denominator


def ceil_scaled(value):
    return -((-value.numerator * SCALE) // value.denominator)


def matching_max(levels, pivots, threshold, leaf_gain, debt):
    """Maximum leaf debt and a realizing injective assignment."""
    # Every debt weight decreases strictly with q.  In an optimal four-leaf
    # matching a leaf uses one of its first four allowed levels: otherwise its
    # four better levels cannot all be occupied by the other three leaves.
    choices = tuple(
        tuple(q for q in levels
              if q not in pivots and leaf_gain[leaf, q] <= threshold)[:4]
        for leaf in LEAVES
    )
    if any(not row for row in choices):
        return None
    best = None
    for assignment in product(*choices):
        if len(set(assignment)) != len(LEAVES):
            continue
        row = (
            sum(debt[leaf, q] for leaf, q in zip(LEAVES, assignment)),
            tuple(zip(LEAVES, assignment)),
        )
        if best is None or row > best:
            best = row
    return best


def scan(R, D, m):
    levels = tuple(range(m, m + D + 1))
    arcs = {
        (i, q): R.reflected_level_arcs(L, H[i], q, CELL)
        for i in range(6) for q in levels
    }
    gain_exact = {
        (i, j, p, q): R.interval_mass(arcs[i, p]) + R.interval_mass(arcs[j, q])
        - R.interval_mass(R.merge_intervals(arcs[i, p] + arcs[j, q]))
        for i, j in EDGES for p in levels for q in levels if p != q
    }
    debt_exact = {
        (i, q): F(H[i], 7 * (q * L - H[i]))
        for i in range(6) for q in levels
    }
    gain = {key: floor_scaled(value) for key, value in gain_exact.items()}
    debt = {key: ceil_scaled(value) for key, value in debt_exact.items()}
    certificate_best = None
    representative_best = None
    threshold_checks = 0
    for q0, q1 in permutations(levels, 2):
        base_gain = gain[0, 1, q0, q1]
        leaf_gain = {}
        candidates = {base_gain}
        for leaf in LEAVES:
            for q in levels:
                if q in (q0, q1):
                    continue
                value = max(gain[0, leaf, q0, q], gain[1, leaf, q1, q])
                leaf_gain[leaf, q] = value
                candidates.add(value)
        base_debt = debt[0, q0] + debt[1, q1]
        global_debt_upper = base_debt + sum(
            max(debt[leaf, q] for q in levels if q not in (q0, q1))
            for leaf in LEAVES
        )
        for threshold in sorted(value for value in candidates if value >= base_gain):
            threshold_checks += 1
            if threshold > global_debt_upper:
                row = (
                    threshold - global_debt_upper,
                    (q0, q1),
                    threshold,
                    global_debt_upper,
                    "relaxed-tail",
                )
                if certificate_best is None or row < certificate_best:
                    certificate_best = row
                break
            matching = matching_max(levels, (q0, q1), threshold, leaf_gain, debt)
            if matching is None:
                continue
            leaf_debt, raw_assignment = matching
            assignment = [None] * 6
            assignment[0], assignment[1] = q0, q1
            for leaf, q in raw_assignment:
                assignment[leaf] = q
            assignment = tuple(assignment)
            actual_gain, edge = max(
                (gain[i, j, assignment[i], assignment[j]], (i, j)) for i, j in EDGES
            )
            total_debt = base_debt + leaf_debt
            certificate_row = (
                threshold - total_debt,
                assignment,
                threshold,
                total_debt,
                "matching",
            )
            if certificate_best is None or certificate_row < certificate_best:
                certificate_best = certificate_row
            row = (actual_gain - total_debt, assignment, edge, actual_gain, total_debt)
            if representative_best is None or row < representative_best:
                representative_best = row
    if certificate_best is None:
        raise RuntimeError((D, m, "empty"))
    if certificate_best[0] <= 0:
        raise RuntimeError((D, m, "nonpositive certificate floor", certificate_best))
    # Obtain one exact rational representative.  The certificate minimizer can
    # be the relaxed tail, which names only the pivot levels; in that case use
    # the tightest actual matching encountered.
    assignment = (
        representative_best[1]
        if representative_best is not None
        else tuple(levels[:6])
    )
    actual_gain_exact, actual_edge = max(
        (gain_exact[i, j, assignment[i], assignment[j]], (i, j)) for i, j in EDGES
    )
    actual_debt_exact = sum(
        (debt_exact[i, assignment[i]] for i in range(6)), F(0)
    )
    exact_row = (
        actual_gain_exact - actual_debt_exact,
        assignment,
        actual_edge,
        actual_gain_exact,
        actual_debt_exact,
    )
    return certificate_best, exact_row, threshold_checks


def brute(R, D, m):
    levels = tuple(range(m, m + D + 1))
    arcs = {(i, q): R.reflected_level_arcs(L, H[i], q, CELL)
            for i in range(6) for q in levels}
    gain = {
        (i, j, p, q): R.interval_mass(arcs[i, p]) + R.interval_mass(arcs[j, q])
        - R.interval_mass(R.merge_intervals(arcs[i, p] + arcs[j, q]))
        for i, j in EDGES for p in levels for q in levels if p != q
    }
    best = None
    for assignment in permutations(levels, 6):
        debt = sum((F(e, 7 * (q * L - e)) for e, q in zip(H, assignment)), F(0))
        ranked = [
            (gain[i, j, assignment[i], assignment[j]], (i, j))
            for i, j in EDGES
        ]
        actual_gain, edge = max(ranked)
        row = (actual_gain - debt, assignment, edge, actual_gain, debt)
        if best is None or row < best:
            best = row
    return best


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-D", type=int, default=20)
    parser.add_argument("--m-multiples", type=int, default=1)
    parser.add_argument("--all-m", action="store_true")
    parser.add_argument("--brute-control", action="store_true")
    args = parser.parse_args()
    module = load()
    R = module.T.R
    ruler, safe = R.safe_cell_ranges(H)
    if ruler != L or not R.body_cell_is_safe(L, H, CELL):
        raise RuntimeError((ruler, safe, CELL))
    print("LRC H UPPER-MEDIAN TWO-STAR MATCHING -- FINITE-EXACT SCOUT", flush=True)
    print(f"body={H};ruler={L};cell={CELL};edges={EDGES}", flush=True)
    digest = hashlib.sha256()
    for D in range(6, args.max_D + 1):
        m_step = 1 if args.all_m else D
        for m in range(D, args.m_multiples * D + 1, m_step):
            certified, best, checks = scan(R, D, m)
            if args.brute_control and D <= 10:
                control = brute(R, D, m)
                if control[0] <= 0 or certified[0] > floor_scaled(control[0]):
                    raise RuntimeError(("brute certificate mismatch", D, m, best, control, certified))
                best = control
            margin_label = "minimum_margin" if args.brute_control and D <= 10 else "representative_margin"
            line = (f"D={D};m={m};threshold_checks={checks};"
                    f"certified_floor={certified[0]}/{SCALE};{margin_label}={best[0]};"
                    f"levels={best[1]};edge={best[2]};gain={best[3]};debt={best[4]}")
            print(line, flush=True)
            digest.update((line + "\n").encode())
            if certified[0] <= 0:
                print("FIRST_FAILURE")
                print(f"semantic_sha256={digest.hexdigest()}")
                return
    print(f"semantic_sha256={digest.hexdigest()}")


if __name__ == "__main__":
    main()
