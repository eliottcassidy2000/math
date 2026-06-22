#!/usr/bin/env python3
"""Tournament-induction switchboard for the LRC14 proof route.

This audit compares two induction grammars that kept recurring in the repo:

* tournament strong-ear induction, where a new vertex is controlled by the
  boundary state (start, end, Q);
* LRC scale-separated induction, where a new large speed is controlled by the
  seed safe-set topology (measure and interval components).

The point is not that tournaments prove LRC directly.  The point is that the
successful induction attempts on both sides preserve a boundary carrier before
scalarizing.
"""

from __future__ import annotations

import importlib.util
import sys
from collections import Counter
from fractions import Fraction
from itertools import permutations, product
from pathlib import Path


def tournament_from_bits(n: int, bits: tuple[int, ...]) -> list[list[bool]]:
    adj = [[False] * n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits[k]:
                adj[i][j] = True
            else:
                adj[j][i] = True
            k += 1
    return adj


def reachable(adj: list[list[bool]], start: int) -> set[int]:
    seen = {start}
    stack = [start]
    while stack:
        v = stack.pop()
        for w, has_edge in enumerate(adj[v]):
            if has_edge and w not in seen:
                seen.add(w)
                stack.append(w)
    return seen


def is_strong(adj: list[list[bool]]) -> bool:
    n = len(adj)
    return all(len(reachable(adj, v)) == n for v in range(n))


def delete_vertex(adj: list[list[bool]], skip: int) -> list[list[bool]]:
    keep = [v for v in range(len(adj)) if v != skip]
    return [[adj[i][j] for j in keep] for i in keep]


def valid_path(adj: list[list[bool]], perm: tuple[int, ...]) -> bool:
    return all(adj[perm[i]][perm[i + 1]] for i in range(len(perm) - 1))


def boundary_state(adj: list[list[bool]]) -> tuple[list[int], list[int], list[list[int]]]:
    n = len(adj)
    start = [0] * n
    end = [0] * n
    q = [[0] * n for _ in range(n)]
    for perm in permutations(range(n)):
        if valid_path(adj, perm):
            start[perm[0]] += 1
            end[perm[-1]] += 1
        for i in range(n - 1):
            a, b = perm[i], perm[i + 1]
            if all(j == i or adj[perm[j]][perm[j + 1]] for j in range(n - 1)):
                q[a][b] += 1
    return start, end, q


def hamiltonian_paths(adj: list[list[bool]]) -> int:
    return sum(1 for perm in permutations(range(len(adj))) if valid_path(adj, perm))


def add_vertex(adj: list[list[bool]], sig: tuple[int, ...]) -> list[list[bool]]:
    """Add x=n with sig[v]=1 iff x -> v."""
    n = len(adj)
    out = [row + [False] for row in adj]
    out.append([False] * (n + 1))
    for v, bit in enumerate(sig):
        if bit:
            out[n][v] = True
        else:
            out[v][n] = True
    return out


def ear_formula_value(adj: list[list[bool]], sig: tuple[int, ...]) -> int:
    start, end, q = boundary_state(adj)
    starts = sum(start[b] for b, bit in enumerate(sig) if bit)
    ends = sum(end[a] for a, bit in enumerate(sig) if not bit)
    slots = sum(
        q[a][b]
        for a, bit_a in enumerate(sig)
        for b, bit_b in enumerate(sig)
        if not bit_a and bit_b
    )
    return starts + ends + slots


def audit_strong_ears(max_n: int = 5) -> dict[str, object]:
    rows: list[tuple[int, int, int, int, int]] = []
    total_formula_failures = 0
    total_strong_failures = 0
    deletion_mins: dict[int, int] = {}
    for n in range(3, max_n + 1):
        strong_parents = 0
        ears = 0
        formula_failures = 0
        strong_failures = 0
        delete_counts: list[int] = []
        for bits in product((0, 1), repeat=n * (n - 1) // 2):
            adj = tournament_from_bits(n, bits)
            if not is_strong(adj):
                continue
            strong_parents += 1
            if n >= 4:
                delete_counts.append(sum(is_strong(delete_vertex(adj, v)) for v in range(n)))
            for sig in product((0, 1), repeat=n):
                if all(bit == sig[0] for bit in sig):
                    continue
                child = add_vertex(adj, sig)
                ears += 1
                if hamiltonian_paths(child) != ear_formula_value(adj, sig):
                    formula_failures += 1
                if not is_strong(child):
                    strong_failures += 1
        if delete_counts:
            deletion_mins[n] = min(delete_counts)
        rows.append((n, strong_parents, ears, formula_failures, strong_failures))
        total_formula_failures += formula_failures
        total_strong_failures += strong_failures
    return {
        "rows": rows,
        "deletion_mins": deletion_mins,
        "formula_failures": total_formula_failures,
        "strong_failures": total_strong_failures,
    }


def load_lrc_scale_module():
    path = Path(__file__).with_name("lrc_scale_separated_induction_codex.py")
    spec = importlib.util.spec_from_file_location("lrc_scale", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("could not load lrc scale-separated induction module")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def audit_lrc_peel() -> dict[str, object]:
    mod = load_lrc_scale_module()
    seed = tuple(list(range(1, 12)) + [13])
    n = 14
    safe = mod.safe_set(seed, n)
    return {
        "seed": seed,
        "n": n,
        "measure": mod.measure(safe),
        "components": len(safe),
        "least_v": mod.sufficient_v(safe, n),
        "checks": [
            (v, mod.measure(mod.intersect_intervals(safe, mod.safe_set((v,), n))), mod.lower_bound_after_one(safe, n, v))
            for v in (30030, 60060, 510510)
        ],
    }


def carrier_tournament() -> tuple[Counter[int], int, list[str]]:
    features = {
        "strong_ear_boundary_Q": {
            "boundary_state",
            "exact_insertion",
            "resonance_ledger",
            "strong_reducibility",
            "finite_basis",
        },
        "finite_comb_peel": {
            "boundary_state",
            "positive_measure",
            "component_budget",
            "size_descent",
            "exact_comb",
        },
        "multi_large_union_bound": {
            "positive_measure",
            "arc_budget",
            "size_descent",
            "multi_large",
        },
        "half_tiling_address": {
            "boundary_state",
            "mirror_quotient",
            "parity_address",
            "fixed_line",
        },
        "depth_parity_newton": {
            "boundary_state",
            "parity_address",
            "newton_packets",
            "cap_slack",
        },
        "prime_omission_witness": {
            "direct_witness",
            "covering_filter",
            "small_prime",
        },
        "dilation_normalization": {
            "scale_invariant",
            "primitive_core",
        },
        "bounded_covering_core": {
            "bounded_core",
            "finite_atlas",
            "three_gap",
            "ap_gw_tight",
        },
        "raw_runner_deletion": {"size_descent"},
        "raw_tournament_minor": {"size_descent"},
    }
    names = list(features)
    scores = Counter({name: 0 for name in names})
    adj = {name: set() for name in names}
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if j <= i:
                continue
            key_a = (
                "boundary_state" in features[a],
                len(features[a]),
                "size_descent" in features[a],
                -i,
            )
            key_b = (
                "boundary_state" in features[b],
                len(features[b]),
                "size_descent" in features[b],
                -j,
            )
            if key_a >= key_b:
                adj[a].add(b)
                scores[a] += 1
            else:
                adj[b].add(a)
                scores[b] += 1
    cycles3 = 0
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if j <= i:
                continue
            for k, c in enumerate(names):
                if k <= j:
                    continue
                if b in adj[a] and c in adj[b] and a in adj[c]:
                    cycles3 += 1
                if c in adj[a] and b in adj[c] and a in adj[b]:
                    cycles3 += 1
    path = sorted(names, key=lambda name: (scores[name], len(features[name]), name), reverse=True)
    return Counter(scores.values()), cycles3, path


def fmt(q: Fraction) -> str:
    if q.denominator == 1:
        return str(q.numerator)
    return f"{q.numerator}/{q.denominator}"


def main() -> None:
    print("Tournament-induction switchboard for LRC14")
    print("=" * 78)
    print("Principle: successful inductions preserve a boundary state before scalarizing.")

    print("\nTournament strong-ear boundary audit")
    print("-" * 78)
    ear = audit_strong_ears()
    print("  n  strong_parents  nonconstant_ears  formula_failures  strong_failures")
    for n, parents, ears, formula_failures, strong_failures in ear["rows"]:
        print(f"  {n:1d}  {parents:14d}  {ears:16d}  {formula_failures:16d}  {strong_failures:15d}")
    print(f"  deletion_min_by_n={ear['deletion_mins']}")
    print(f"  total_formula_failures={ear['formula_failures']}")
    print(f"  total_strong_failures={ear['strong_failures']}")
    print("  readout: the induction step is exact only after retaining (start,end,Q).")

    print("\nLRC finite-comb peel audit")
    print("-" * 78)
    peel = audit_lrc_peel()
    print(f"  seed={peel['seed']}, n={peel['n']}, threshold=1/{peel['n']}")
    print(f"  seed_measure={fmt(peel['measure'])} components={peel['components']} least_certified_v={peel['least_v']}")
    for v, exact, lower in peel["checks"]:
        print(f"  add v={v:6d}: exact_after={fmt(exact):>12s}  comb_lower={fmt(lower):>14s}")
    print("  readout: the induction step is exact only after retaining (measure,components).")

    print("\nShared reduction tree")
    print("-" * 78)
    rows = [
        ("omit-prime", "direct witness t=1/p", "small-prime cover profile", "done"),
        ("remove-large", "descend to smaller LRC seed", "safe-set measure + component/arc budget", "HYP-2904/S31v/S31w"),
        ("multi-large r<=6", "union-bound survival", "core floor + arc budget", "S31v"),
        ("multi-large r>=7", "second moment", "divisibility/resonant-pair graph", "open"),
        ("bounded covering core", "finite extremality", "three-gap/AP-GW + depth-parity packets", "open Node 2"),
    ]
    for name, action, state, status in rows:
        print(f"  {name:22s} -> {action:28s} | state={state}; status={status}")

    print("\nTournament Analysis over proof carriers")
    print("-" * 78)
    hist, cycles3, path = carrier_tournament()
    print(f"  score_hist={dict(sorted(hist.items()))}")
    print(f"  directed_3cycles={cycles3}")
    print("  Hamiltonian path=" + " > ".join(path))

    print("\nFlush-out target")
    print("-" * 78)
    print("  The LRC14 proof should now be stated as boundary-state induction:")
    print("  (1) descend all unbounded/non-covering rows by omit-prime and remove-large;")
    print("  (2) close r>=7 by a resonant-pair/divisibility second-moment ledger;")
    print("  (3) close the bounded covering core by AP/GW tight-locus rigidity plus")
    print("      missing-depth parity packets, not by raw runner deletion or scalar minors.")


if __name__ == "__main__":
    main()
