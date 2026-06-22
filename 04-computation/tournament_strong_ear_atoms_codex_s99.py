#!/usr/bin/env python3
"""
tournament_strong_ear_atoms_codex_s99.py

Codex S99.  Extend the strong-component H-atom lens by replacing raw
enumeration with one-vertex strong ears.

If T is a tournament and x is a new vertex, write sig[v]=1 when x -> v and
sig[v]=0 when v -> x.  For a Hamiltonian path P of T, x can be inserted in
a slot exactly when the two incident edges around that slot match sig.  The
resulting formula is

  H(T+x) =
      sum_{b: sig[b]=1} start_T[b]
    + sum_{a: sig[a]=0} end_T[a]
    + sum_{a: sig[a]=0, b: sig[b]=1} Q_T[a,b],

where Q_T[a,b] counts old-vertex permutations with a immediately followed by b
and every other adjacent step a valid edge of T.  For strong T, every
nonconstant sig gives a strong child, because x can both reach T and be reached
from T.

This script audits three claims through the validated n<=8 isomorphism-class
tower:

  1. Strong-ear reducibility: every strong tournament on n=4..8 has a vertex
     whose deletion is strong.
  2. Ear completeness: nonconstant ears from strong (n-1)-atoms generate the
     full strong H-spectrum at n, for n=4..8.
  3. Finite ear bases: a small set of parent atoms and a small set of cut
     weights cover the strong H-spectrum, mirroring the finite denominator
     basis in the LRC rational-witness route.

Tournament Analysis declaration.
  Vertices: proof carriers for strong atom growth.
  Pairwise observable:
    (preserves strongness, exact H formula, completeness evidence, finite-basis
     compression, LRC/residue analogy value, formalization readiness).
  Gauge: lexicographic comparison; ties follow the listed Hamiltonian path.

Assumption challenge.
  I considered tournament vertices, strong components, deletion vertices,
  insertion cuts, cut weights, bridge-exposure matrices Q, OCF packets,
  even-graph shadows, and LRC rational denominators as possible vertices.
  The chosen quotient is the ear boundary state (starts, ends, Q) plus the
  nonconstant cut.  It preserves exact H under atom growth and destroys raw
  vertex labels except where needed for the cut sum.
"""

from __future__ import annotations

import importlib.util
from collections import Counter, defaultdict
from itertools import combinations, permutations
from pathlib import Path

HERE = Path(__file__).resolve().parent
S5_PATH = HERE / "strong_H_spectrum_m8_isoclass_monad_s5.py"
SPEC = importlib.util.spec_from_file_location("strong_s5", S5_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot load {S5_PATH}")
s5 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(s5)

MAX_N = 8
A000568 = {1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880}


def has_arc(adj: list[int], i: int, j: int) -> bool:
    return ((adj[i] >> j) & 1) == 1


def induced_delete(adj: list[int], drop: int) -> list[int]:
    old = [v for v in range(len(adj)) if v != drop]
    idx = {v: i for i, v in enumerate(old)}
    out = [0] * len(old)
    for a in old:
        for b in old:
            if a != b and has_arc(adj, a, b):
                out[idx[a]] |= 1 << idx[b]
    return out


def add_vertex(adj: list[int], sig: int) -> list[int]:
    """sig bit v = 1 means new vertex x -> v; bit 0 means v -> x."""
    n = len(adj)
    out = list(adj) + [0]
    for v in range(n):
        if (sig >> v) & 1:
            out[n] |= 1 << v
        else:
            out[v] |= 1 << n
    return out


def boundary_state(adj: list[int]) -> tuple[list[int], list[int], list[list[int]]]:
    """Return start counts, end counts, and bridge-exposure Q."""
    n = len(adj)
    starts = [0] * n
    ends = [0] * n
    q = [[0] * n for _ in range(n)]

    for p in permutations(range(n)):
        valid = [has_arc(adj, p[i], p[i + 1]) for i in range(n - 1)]
        all_valid = all(valid)
        if all_valid:
            starts[p[0]] += 1
            ends[p[-1]] += 1

        # q[a,b] permits the a,b slot to be repaired by the new vertex; all
        # other old adjacencies must already be valid.
        prefix_bad = [0] * n
        for i in range(n - 1):
            prefix_bad[i + 1] = prefix_bad[i] + (0 if valid[i] else 1)
        total_bad = prefix_bad[n - 1]
        for k in range(n - 1):
            bad_except_k = total_bad - (0 if valid[k] else 1)
            if bad_except_k == 0:
                q[p[k]][p[k + 1]] += 1
    return starts, ends, q


def insertion_h(starts: list[int], ends: list[int], q: list[list[int]], sig: int) -> int:
    n = len(starts)
    total = 0
    for b in range(n):
        if (sig >> b) & 1:
            total += starts[b]
    for a in range(n):
        if ((sig >> a) & 1) == 0:
            total += ends[a]
            for b in range(n):
                if (sig >> b) & 1:
                    total += q[a][b]
    return total


def cut_weight(sig: int) -> int:
    return sig.bit_count()


def is_nonconstant(sig: int, n: int) -> bool:
    return sig != 0 and sig != (1 << n) - 1


def strong_deletion_counts(reps: dict[int, list[list[int]]]) -> dict[int, Counter[int]]:
    out: dict[int, Counter[int]] = {}
    for n in range(3, MAX_N + 1):
        dist: Counter[int] = Counter()
        for adj in reps[n]:
            if not s5.is_strong(n, adj):
                continue
            good = 0
            for v in range(n):
                child = induced_delete(adj, v)
                if len(child) >= 3 and s5.is_strong(n - 1, child):
                    good += 1
            dist[good] += 1
        out[n] = dist
    return out


def strong_spectra(reps: dict[int, list[list[int]]]) -> dict[int, set[int]]:
    spectra: dict[int, set[int]] = {}
    for n in range(3, MAX_N + 1):
        vals = set()
        for adj in reps[n]:
            if s5.is_strong(n, adj):
                vals.add(s5.Hcount(n, adj))
        spectra[n] = vals
    return spectra


def greedy_cover(universe: set[int], covers: dict[str, set[int]]) -> list[tuple[str, int, int]]:
    remaining = set(universe)
    chosen: list[tuple[str, int, int]] = []
    while remaining:
        best_key = None
        best_cover: set[int] = set()
        for key, vals in covers.items():
            gain = vals & remaining
            if len(gain) > len(best_cover):
                best_key = key
                best_cover = gain
        if best_key is None or not best_cover:
            raise RuntimeError(f"uncovered values remain: {sorted(remaining)[:20]}")
        chosen.append((best_key, len(best_cover), len(remaining) - len(best_cover)))
        remaining -= best_cover
    return chosen


def ear_transition_report(
    n: int,
    parent_reps: list[list[int]],
    target_spectrum: set[int],
) -> dict[str, object]:
    """Analyze strong n-parent ears generating strong (n+1)-children."""
    by_parent: dict[str, set[int]] = {}
    by_parent_h: dict[str, int] = {}
    by_weight: dict[str, set[int]] = defaultdict(set)
    ear_values: set[int] = set()
    formula_failures = []
    strong_failures = []
    total_cuts = 0

    for parent_id, adj in enumerate(parent_reps):
        if not s5.is_strong(n, adj):
            continue
        parent_h = s5.Hcount(n, adj)
        starts, ends, q = boundary_state(adj)
        key = f"parent{parent_id}:H{parent_h}"
        by_parent_h[key] = parent_h
        vals_for_parent: set[int] = set()
        for sig in range(1 << n):
            if not is_nonconstant(sig, n):
                continue
            total_cuts += 1
            h_formula = insertion_h(starts, ends, q, sig)
            child = add_vertex(adj, sig)
            if not s5.is_strong(n + 1, child):
                strong_failures.append((parent_id, sig))
            # Direct DP audit for small sizes; n=7 -> 8 is still cheap enough
            # because there are only 353 strong parents and 126 nonconstant cuts.
            h_direct = s5.Hcount(n + 1, child)
            if h_formula != h_direct:
                formula_failures.append((parent_id, sig, h_formula, h_direct))
            ear_values.add(h_formula)
            vals_for_parent.add(h_formula)
            by_weight[f"w={cut_weight(sig)}"].add(h_formula)
        by_parent[key] = vals_for_parent

    target_missing = sorted(target_spectrum - ear_values)
    extra_values = sorted(ear_values - target_spectrum)
    parent_basis = greedy_cover(target_spectrum, by_parent)
    weight_basis = greedy_cover(target_spectrum, by_weight)
    return {
        "total_cuts": total_cuts,
        "ear_values": ear_values,
        "target_missing": target_missing,
        "extra_values": extra_values,
        "formula_failures": formula_failures,
        "strong_failures": strong_failures,
        "by_parent": by_parent,
        "by_parent_h": by_parent_h,
        "by_weight": by_weight,
        "parent_basis": parent_basis,
        "weight_basis": weight_basis,
    }


def summarize_basis(basis: list[tuple[str, int, int]], limit: int = 12) -> str:
    head = basis[:limit]
    parts = [f"{key} (+{gain}, rem={rem})" for key, gain, rem in head]
    if len(basis) > limit:
        parts.append(f"... {len(basis)-limit} more")
    return "; ".join(parts)


def tournament_analysis() -> None:
    carriers = [
        "ear-boundary-Q-cut",
        "strong-component-H-atoms",
        "OCF-packets",
        "finite-denominator-residues",
        "even-graph-shadow",
        "raw-H-scalar",
    ]
    score = {
        "ear-boundary-Q-cut": (1, 1, 1, 1, 1, 1),
        "strong-component-H-atoms": (1, 1, 1, 1, 0, 1),
        "OCF-packets": (1, 1, 0, 1, 0, 0),
        "finite-denominator-residues": (0, 1, 0, 1, 1, 0),
        "even-graph-shadow": (0, 0, 0, 1, 1, 0),
        "raw-H-scalar": (0, 0, 0, 0, 0, 1),
    }
    outscore = Counter()
    wins: dict[tuple[str, str], str] = {}
    order = {v: i for i, v in enumerate(carriers)}
    for a, b in combinations(sorted(carriers), 2):
        if score[a] > score[b]:
            winner = a
        elif score[b] > score[a]:
            winner = b
        else:
            winner = a if order[a] < order[b] else b
        wins[(a, b)] = winner
        outscore[winner] += 1
    hist = Counter(outscore[v] for v in carriers)
    ranked = sorted(carriers, key=lambda v: (outscore[v], score[v]), reverse=True)

    def edge(a: str, b: str) -> bool:
        return wins[(min(a, b), max(a, b))] == a

    hp = 0
    for p in permutations(carriers):
        if all(edge(p[i], p[i + 1]) for i in range(len(p) - 1)):
            hp += 1
    print("\nTournament Analysis: strong-atom growth carriers")
    print(f"  score histogram: {dict(sorted(hist.items()))}")
    print(f"  Hamiltonian path count: {hp}")
    print("  leading path:", " -> ".join(ranked))


def main() -> None:
    print("=" * 78)
    print("TOURNAMENT STRONG EAR ATOMS")
    print("=" * 78)
    print("Generating non-isomorphic tournaments through n=8 via validated S5 engine...")
    reps, counts = s5.generate(MAX_N)
    print("  counts:", [counts[n] for n in range(1, MAX_N + 1)])
    print("  A000568:", [A000568[n] for n in range(1, MAX_N + 1)])
    print("  counts match:", all(counts[n] == A000568[n] for n in range(1, MAX_N + 1)))

    spectra = strong_spectra(reps)
    print("\nStrong H-spectra summary:")
    for n in range(3, MAX_N + 1):
        vals = sorted(spectra[n])
        print(f"  n={n}: count={len(vals):3d}, min={vals[0]:3d}, max={vals[-1]:3d}, "
              f"has7={7 in spectra[n]}, has21={21 in spectra[n]}")

    print("\nStrong-deletion reducibility:")
    deletion = strong_deletion_counts(reps)
    for n in range(3, MAX_N + 1):
        dist = deletion[n]
        total = sum(dist.values())
        print(f"  n={n}: strong classes={total}, deletion-count distribution={dict(sorted(dist.items()))}")

    print("\nEar-generation transitions:")
    transition_reports = {}
    for n in range(3, MAX_N):
        report = ear_transition_report(n, reps[n], spectra[n + 1])
        transition_reports[n] = report
        print(f"\n  {n} -> {n+1}:")
        print(f"    nonconstant strong-ear cuts checked: {report['total_cuts']}")
        print(f"    formula failures: {len(report['formula_failures'])}")
        print(f"    strongness failures: {len(report['strong_failures'])}")
        print(f"    ear H-values: {len(report['ear_values'])}; target strong spectrum: {len(spectra[n+1])}")
        print(f"    missing target values: {report['target_missing'][:20]}")
        print(f"    extra values outside target: {report['extra_values'][:20]}")
        print(f"    parent basis size: {len(report['parent_basis'])}")
        print(f"      {summarize_basis(report['parent_basis'])}")
        print(f"    cut-weight basis size: {len(report['weight_basis'])}")
        print(f"      {summarize_basis(report['weight_basis'])}")

    print("\nCut-weight coverage in the largest transition 7 -> 8:")
    big = transition_reports[7]
    for key in sorted(big["by_weight"], key=lambda s: int(s.split("=")[1])):
        vals = big["by_weight"][key]
        print(f"  {key}: covers {len(vals):3d} of {len(spectra[8])} strong H-values")
    for key in ("w=3", "w=4"):
        missing = sorted(spectra[8] - big["by_weight"][key])
        print(f"  {key}: missing values {missing}")

    print("\nInterpretation:")
    print("  1. Strong atom growth is an exact ear insertion problem, not a raw minor move.")
    print("  2. The boundary state (starts, ends, Q) is the atom analogue of the LRC")
    print("     character-sum ledger: linear boundary terms plus a cut/resonance matrix.")
    print("  3. Nonconstant ears preserve strongness automatically when the parent is strong.")
    print("  4. Through n=8, strong spectra are fully generated by ears from smaller")
    print("     strong atoms; small parent/cut bases cover all observed values.")
    print("  5. The forbidden {7,21} values are not minor-closed shadows; they are")
    print("     absent solutions of these recursive cut polynomials before the")
    print("     strong-min boundary lets the spectrum re-enter densely.")
    tournament_analysis()
    print("\nDONE.")


if __name__ == "__main__":
    main()
