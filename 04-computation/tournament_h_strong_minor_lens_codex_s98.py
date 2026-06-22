#!/usr/bin/env python3
"""
tournament_h_strong_minor_lens_codex_s98.py

Codex S98.  Positive companion to the S96 even-graph minor obstruction.

Kuratowski/Wagner says planarity is controlled by graph simplifications whose
predicate is closed under those simplifications.  S96 showed that the natural
even-graph simplifications do not preserve the tournament Hamiltonian-path
count H.  This script tests the closure that actually does preserve H in the
fixed-Hamiltonian-path tournament cube:

  * suppress singleton strong components (factor H=1);
  * contract each nontrivial strong component only as a labeled H-atom.

The result is a Wagner-style "minor lens" for H, but the minors are
strong-component/OCF packet minors, not arbitrary even-graph minors.

Tournament Analysis declaration.
  Vertices: proof carriers/simplification orders.
  Pairwise observable:
    (preserves H predicate, has exact factor law, sees {7,21}, finite audit
     support, compatibility with even-graph bridge, Beurling/BS analogy value).
  Gauge: lexicographic comparison; ties follow the declared Hamiltonian path.
  Fingerprints: score histogram, directed 3-cycles, SCC sizes, Hamiltonian path
  count.

Assumption challenge.
  I considered runners, arcs, even-graph edges, degree-2 smoothing cores,
  arbitrary edge contractions, strong components, OCF conflict components,
  H-values, Fourier/Beurling majorants, and proof obligations.  The chosen
  quotient is strong-component H-atoms: it preserves the predicate "which H
  values are multiplicatively attainable" and destroys the GF(2) cycle-space
  geometry.  S96 proves that geometry is too lossy for H; here we keep only the
  factor ledger that H actually respects.
"""

from __future__ import annotations

import importlib.util
from collections import Counter, defaultdict
from itertools import combinations, permutations
from math import comb, prod
from pathlib import Path

HERE = Path(__file__).resolve().parent
S96_PATH = HERE / "tournament_even_minor_obstruction_codex_s96.py"
SPEC = importlib.util.spec_from_file_location("s96_even_minor", S96_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot load {S96_PATH}")
s96 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(s96)

MAX_N = 7
FORBIDDEN = {7, 21}
_SUB_H_CACHE: dict[tuple[tuple[bool, ...], ...], int] = {}


def induced_subtournament(adj: list[list[bool]], comp: tuple[int, ...]) -> list[list[bool]]:
    idx = {v: i for i, v in enumerate(comp)}
    m = len(comp)
    out = [[False] * m for _ in range(m)]
    for a in comp:
        for b in comp:
            if a != b:
                out[idx[a]][idx[b]] = adj[a][b]
    return out


def adj_key(adj: list[list[bool]]) -> tuple[tuple[bool, ...], ...]:
    return tuple(tuple(row) for row in adj)


def h_cached(adj: list[list[bool]]) -> int:
    key = adj_key(adj)
    if key not in _SUB_H_CACHE:
        _SUB_H_CACHE[key] = s96.hamiltonian_paths(adj)
    return _SUB_H_CACHE[key]


def strong_factor_signature(adj: list[list[bool]]) -> tuple[tuple[int, int], ...]:
    comps = s96.strong_components(adj)
    atoms: list[tuple[int, int]] = []
    for comp in comps:
        if len(comp) == 1:
            continue
        h = h_cached(induced_subtournament(adj, comp))
        atoms.append((len(comp), h))
    return tuple(sorted(atoms))


def factor_h(adj: list[list[bool]]) -> int:
    values = []
    for comp in s96.strong_components(adj):
        if len(comp) == 1:
            values.append(1)
        else:
            values.append(h_cached(induced_subtournament(adj, comp)))
    return prod(values, start=1)


def semigroup_closure(gens: set[int], limit: int) -> set[int]:
    seen = {1}
    changed = True
    while changed:
        changed = False
        for x in list(seen):
            for g in gens:
                y = x * g
                if y <= limit and y not in seen:
                    seen.add(y)
                    changed = True
    return seen


def analyze() -> dict[str, object]:
    factor_failures = []
    singleton_suppression_failures = []
    strong_atoms_by_size: dict[int, Counter[int]] = defaultdict(Counter)
    h_spectrum_by_n: dict[int, Counter[int]] = {}
    signature_to_h: dict[tuple[tuple[int, int], ...], set[int]] = defaultdict(set)
    signature_to_count: Counter[tuple[tuple[int, int], ...]] = Counter()
    rows_checked = 0

    for n in range(1, MAX_N + 1):
        m = comb(n - 1, 2)
        total = 1 << m
        h_counter: Counter[int] = Counter()
        for mask in range(total):
            bits = s96.bits_of_mask(n, mask)
            adj = s96.tournament_from_bits(n, bits)
            h = h_cached(adj)
            fh = factor_h(adj)
            sig = strong_factor_signature(adj)
            rows_checked += 1
            h_counter[h] += 1
            signature_to_h[sig].add(h)
            signature_to_count[sig] += 1

            if h != fh:
                factor_failures.append((n, mask, h, fh, sig))

            # Suppressing singleton strong components should not change H; it
            # leaves exactly the product of nontrivial atom values.
            suppressed_h = prod((atom_h for _, atom_h in sig), start=1)
            if h != suppressed_h:
                singleton_suppression_failures.append((n, mask, h, suppressed_h, sig))

            comps = s96.strong_components(adj)
            if len(comps) == 1:
                strong_atoms_by_size[n][h] += 1
        h_spectrum_by_n[n] = h_counter

    atom_values = {h for counter in strong_atoms_by_size.values() for h in counter if h != 1}
    max_h = max(max(counter) for counter in h_spectrum_by_n.values())
    closure = semigroup_closure(atom_values, max_h)
    actual = {h for counter in h_spectrum_by_n.values() for h in counter}
    odd_values = set(range(1, max_h + 1, 2))
    actual_gaps = sorted(odd_values - actual)
    semigroup_gaps = sorted(odd_values - closure)

    return {
        "rows_checked": rows_checked,
        "factor_failures": factor_failures,
        "singleton_suppression_failures": singleton_suppression_failures,
        "strong_atoms_by_size": strong_atoms_by_size,
        "h_spectrum_by_n": h_spectrum_by_n,
        "signature_to_h": signature_to_h,
        "signature_to_count": signature_to_count,
        "atom_values": atom_values,
        "max_h": max_h,
        "closure": closure,
        "actual": actual,
        "actual_gaps": actual_gaps,
        "semigroup_gaps": semigroup_gaps,
    }


def strongly_connected_components(vertices: list[str], winner: dict[tuple[str, str], str]) -> list[list[str]]:
    adj = {v: [] for v in vertices}
    radj = {v: [] for v in vertices}
    for (a, b), w in winner.items():
        loser = b if w == a else a
        adj[w].append(loser)
        radj[loser].append(w)
    seen: set[str] = set()
    order: list[str] = []

    def dfs(v: str) -> None:
        seen.add(v)
        for u in adj[v]:
            if u not in seen:
                dfs(u)
        order.append(v)

    for v in vertices:
        if v not in seen:
            dfs(v)
    seen.clear()
    comps: list[list[str]] = []

    def rdfs(v: str, comp: list[str]) -> None:
        seen.add(v)
        comp.append(v)
        for u in radj[v]:
            if u not in seen:
                rdfs(u, comp)

    for v in reversed(order):
        if v not in seen:
            comp: list[str] = []
            rdfs(v, comp)
            comps.append(sorted(comp))
    return comps


def directed_triangles(vertices: list[str], winner: dict[tuple[str, str], str]) -> int:
    total = 0
    for a, b, c in combinations(vertices, 3):
        out = Counter(
            [
                winner[(min(a, b), max(a, b))],
                winner[(min(a, c), max(a, c))],
                winner[(min(b, c), max(b, c))],
            ]
        )
        if sorted(out.values()) == [1, 1, 1]:
            total += 1
    return total


def edge(winner: dict[tuple[str, str], str], a: str, b: str) -> bool:
    return winner[(min(a, b), max(a, b))] == a


def hamiltonian_path_count(vertices: list[str], winner: dict[tuple[str, str], str]) -> int:
    return sum(
        1
        for path in permutations(vertices)
        if all(edge(winner, path[i], path[i + 1]) for i in range(len(path) - 1))
    )


def tournament_analysis() -> None:
    carriers = [
        "strong-component-H-atoms",
        "OCF-conflict-packets",
        "even-graph-cycle-space",
        "degree2-GF2-smoothing",
        "arbitrary-edge-contraction",
        "Beurling-Selberg-majorants",
    ]
    score = {
        "strong-component-H-atoms": (1, 1, 1, 1, 1, 0),
        "OCF-conflict-packets": (1, 1, 1, 0, 0, 0),
        "even-graph-cycle-space": (0, 0, 0, 1, 1, 0),
        "degree2-GF2-smoothing": (0, 0, 0, 1, 0, 0),
        "arbitrary-edge-contraction": (0, 0, 0, 1, 0, 0),
        "Beurling-Selberg-majorants": (0, 0, 0, 0, 0, 1),
    }
    order = {name: i for i, name in enumerate(carriers)}
    wins: dict[tuple[str, str], str] = {}
    outscore = Counter()
    for a, b in combinations(sorted(carriers), 2):
        if score[a] > score[b]:
            w = a
        elif score[b] > score[a]:
            w = b
        else:
            w = a if order[a] < order[b] else b
        wins[(a, b)] = w
        outscore[w] += 1
    hist = Counter(outscore[v] for v in carriers)
    ranked = sorted(carriers, key=lambda v: (outscore[v], score[v]), reverse=True)
    print("\nTournament Analysis: minor-lens carriers")
    print("  score histogram:", dict(sorted(hist.items())))
    print("  directed 3-cycles:", directed_triangles(carriers, wins))
    print("  SCC sizes:", [len(c) for c in strongly_connected_components(carriers, wins)])
    print("  Hamiltonian path count:", hamiltonian_path_count(carriers, wins))
    print("  leading Hamiltonian path:", " -> ".join(ranked))


def main() -> None:
    data = analyze()
    print("=" * 78)
    print("TOURNAMENT H STRONG-MINOR LENS")
    print("=" * 78)
    print(f"fixed-path rooted tournaments checked through n={MAX_N}: {data['rows_checked']}")
    print(f"H factorization failures: {len(data['factor_failures'])}")
    print(f"singleton strong-component suppression failures: {len(data['singleton_suppression_failures'])}")

    print("\nH spectra by rooted fixed-path size:")
    for n, counter in data["h_spectrum_by_n"].items():
        values = sorted(counter)
        print(
            f"  n={n}: rows={sum(counter.values())}, H_count={len(values)}, "
            f"min={min(values)}, max={max(values)}, has7={7 in counter}, has21={21 in counter}"
        )

    print("\nStrong atom H-values by size (rooted fixed-path representatives):")
    for n, counter in data["strong_atoms_by_size"].items():
        values = sorted(counter)
        print(f"  size={n}: count={len(values)}, values={values[:80]}{' ...' if len(values) > 80 else ''}")

    print("\nSemigroup comparison up to max rooted H:")
    max_h = data["max_h"]
    actual_gaps = data["actual_gaps"]
    semigroup_gaps = data["semigroup_gaps"]
    print(f"  max_h={max_h}")
    print(f"  observed strong atom generators <=max_h: {sorted(data['atom_values'])}")
    print(f"  actual rooted odd gaps <=max_h (first 40): {actual_gaps[:40]}")
    print(f"  semigroup gaps from observed strong atoms <=max_h (first 40): {semigroup_gaps[:40]}")
    print(f"  forbidden {sorted(FORBIDDEN)} in semigroup closure? {sorted(h for h in FORBIDDEN if h in data['closure'])}")
    print(f"  forbidden {sorted(FORBIDDEN)} in actual rooted spectrum? {sorted(h for h in FORBIDDEN if h in data['actual'])}")

    print("\nCore signatures:")
    signature_to_h = data["signature_to_h"]
    signature_to_count = data["signature_to_count"]
    multi_h = [(sig, hs) for sig, hs in signature_to_h.items() if len(hs) > 1]
    print(f"  signatures={len(signature_to_h)}, signatures with multiple H values={len(multi_h)}")
    if multi_h:
        sig, hs = max(multi_h, key=lambda item: len(item[1]))
        print(f"  worst signature collision: sig={sig}, count={signature_to_count[sig]}, Hs={sorted(hs)[:30]}")
    else:
        print("  H is determined by the strong atom signature in this audit.")

    print("\nInterpretation:")
    print("  1. Strong-component factorization is the H-preserving contraction law.")
    print("  2. Singleton strong components are the safe suppressions: they carry H=1.")
    print("  3. Nontrivial strong components cannot be forgotten; they must be contracted")
    print("     as labeled H-atoms. This is the Wagner-style carrier that S96's even")
    print("     graph smoothing/contraction lacked.")
    print("  4. Beurling-Selberg/trigonometric majorants belong on the analytic LRC side:")
    print("     they are extremal carriers with explicit defect terms, analogous in role")
    print("     to H-atoms, not a replacement for the H-preserving minor order.")
    tournament_analysis()
    print("\nDONE.")


if __name__ == "__main__":
    main()
