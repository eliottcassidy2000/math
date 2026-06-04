#!/usr/bin/env python3
"""
S626: Moser unit-distance spines and flip-gauge tournaments.

The user question asks whether the Hamiltonian path forced by tournament
analysis should be read as a path through unit-distance pairs or through the
non-unit background.  This script separates the intrinsic graph question from
the tournament gauge question:

1. Build dense Moser-carrier witnesses using the S622 frontier beam.
2. Test whether the unit-distance graph has a spanning path, a "unit spine".
3. Compare flip gauges: unit edges flipped vs non-unit edges flipped relative
   to a fixed tie path.

Tournament Analysis is declared over gauges/proof routes, not only over points.
The point tournament is a useful diagnostic, but it forgets embedding
provenance, direction-shell channels, and traceability witnesses unless those
labels are reattached.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.append(str(SCRIPT_DIR))

import unit_distance_impairment_lab_s622 as lab  # noqa: E402


Point = tuple[int, int, int, int]
Cluster = tuple[Point, ...]


def bit_index(x: int) -> int:
    return x.bit_length() - 1


def unit_graph(cluster: Cluster) -> tuple[list[int], list[int]]:
    """Return adjacency bitsets and unit-edge bitsets for a Moser cluster."""
    n = len(cluster)
    unit_set = set(lab.UNITS)
    adj = [0] * n
    for i, j in combinations(range(n), 2):
        if lab.sub(cluster[j], cluster[i]) in unit_set or lab.sub(cluster[i], cluster[j]) in unit_set:
            adj[i] |= 1 << j
            adj[j] |= 1 << i
    return adj, adj[:]


def edge_count(adj: list[int]) -> int:
    return sum(x.bit_count() for x in adj) // 2


def degree_hist(adj: list[int]) -> dict[int, int]:
    hist: dict[int, int] = {}
    for bits in adj:
        d = bits.bit_count()
        hist[d] = hist.get(d, 0) + 1
    return dict(sorted(hist.items()))


def hamiltonian_unit_spine(adj: list[int]) -> tuple[bool, tuple[int, ...]]:
    """Find one Hamiltonian path in an undirected graph, if present."""
    n = len(adj)
    full = (1 << n) - 1
    ends = [0] * (1 << n)
    for v in range(n):
        ends[1 << v] = 1 << v

    for mask in range(1 << n):
        e = ends[mask]
        while e:
            b = e & -e
            v = bit_index(b)
            e ^= b
            nxts = adj[v] & ~mask
            while nxts:
                nb = nxts & -nxts
                nxts ^= nb
                ends[mask | nb] |= nb

    if ends[full] == 0:
        return False, ()

    end = bit_index(ends[full] & -ends[full])
    mask = full
    path = [end]
    while mask != (1 << end):
        prev_mask = mask ^ (1 << end)
        candidates = ends[prev_mask] & adj[end]
        prev = bit_index(candidates & -candidates)
        path.append(prev)
        mask = prev_mask
        end = prev
    path.reverse()
    return True, tuple(path)


def count_unit_hamiltonian_paths(adj: list[int]) -> int:
    """Count undirected Hamiltonian paths exactly; intended for n <= 14."""
    n = len(adj)
    full = (1 << n) - 1
    dp: list[dict[int, int]] = [dict() for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        row = dp[mask]
        for last, count in list(row.items()):
            nxts = adj[last] & ~mask
            while nxts:
                nb = nxts & -nxts
                nxts ^= nb
                nxt = bit_index(nb)
                nxt_row = dp[mask | nb]
                nxt_row[nxt] = nxt_row.get(nxt, 0) + count
    return sum(dp[full].values()) // 2


def orient_by_gauge(
    cluster: Cluster,
    order: tuple[int, ...],
    flip_units: bool,
) -> tuple[list[int], list[int]]:
    """
    Complete tournament from a tie order.

    The tie order is transitive by default.  If flip_units is true, exactly unit
    pairs are reversed.  If false, non-unit pairs are reversed, so unit pairs
    remain transitive.
    """
    n = len(cluster)
    unit_adj, unit_bits = unit_graph(cluster)
    rank = {v: i for i, v in enumerate(order)}
    out = [0] * n
    for i, j in combinations(range(n), 2):
        unit = bool(unit_adj[i] & (1 << j))
        forward_i_to_j = rank[i] < rank[j]
        flipped = unit if flip_units else not unit
        if flipped:
            forward_i_to_j = not forward_i_to_j
        if forward_i_to_j:
            out[i] |= 1 << j
        else:
            out[j] |= 1 << i
    return out, unit_bits


def max_unit_edges_in_directed_hp(out: list[int], unit_bits: list[int]) -> int:
    """Maximum number of unit-labelled arcs along any directed Hamiltonian path."""
    n = len(out)
    full = (1 << n) - 1
    neg = -10**9
    dp: list[dict[int, int]] = [dict() for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 0
    for mask in range(1 << n):
        for last, score in list(dp[mask].items()):
            nxts = out[last] & ~mask
            while nxts:
                nb = nxts & -nxts
                nxts ^= nb
                nxt = bit_index(nb)
                add = 1 if (unit_bits[last] & nb) else 0
                nxt_row = dp[mask | nb]
                old = nxt_row.get(nxt, neg)
                if score + add > old:
                    nxt_row[nxt] = score + add
    return max(dp[full].values())


def directed_triangles(out: list[int]) -> int:
    total = 0
    n = len(out)
    for i, j, k in combinations(range(n), 3):
        if (out[i] >> j) & 1 and (out[j] >> k) & 1 and (out[k] >> i) & 1:
            total += 1
        if (out[i] >> k) & 1 and (out[k] >> j) & 1 and (out[j] >> i) & 1:
            total += 1
    return total


def scc_sizes(out: list[int]) -> tuple[int, ...]:
    n = len(out)

    def reach(starts: list[int], graph: list[int]) -> set[int]:
        seen: set[int] = set(starts)
        stack = starts[:]
        while stack:
            v = stack.pop()
            bits = graph[v]
            while bits:
                b = bits & -bits
                bits ^= b
                w = bit_index(b)
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        return seen

    rev = [0] * n
    for v in range(n):
        bits = out[v]
        while bits:
            b = bits & -bits
            bits ^= b
            rev[bit_index(b)] |= 1 << v

    remaining = set(range(n))
    sizes: list[int] = []
    while remaining:
        v = next(iter(remaining))
        comp = reach([v], out) & reach([v], rev)
        sizes.append(len(comp))
        remaining -= comp
    return tuple(sorted(sizes, reverse=True))


def score_hist(out: list[int]) -> dict[int, int]:
    hist: dict[int, int] = {}
    for bits in out:
        d = bits.bit_count()
        hist[d] = hist.get(d, 0) + 1
    return dict(sorted(hist.items()))


@dataclass(frozen=True)
class Gauge:
    name: str
    preserves_unit_spine: int
    order_invariant: int
    geometry_retained: int
    proof_power: int
    computation: int
    risk: int
    note: str

    def score_tuple(self) -> tuple[int, int, int, int, int, int]:
        return (
            self.preserves_unit_spine,
            self.order_invariant,
            self.geometry_retained,
            self.proof_power,
            self.computation,
            -self.risk,
        )


GAUGES = (
    Gauge("unit graph traceability", 5, 5, 5, 4, 4, 1, "intrinsic yes/no unit-spine predicate"),
    Gauge("spine-order flip gauge", 5, 2, 4, 4, 4, 1, "chooses tie path from a found unit spine"),
    Gauge("direction-pair quotient", 3, 4, 5, 4, 5, 2, "retains Moser shell channels but forgets point order"),
    Gauge("frontier-gain recursion", 3, 4, 4, 5, 5, 1, "state-local construction/deletion ledger"),
    Gauge("lexicographic point flip", 2, 1, 3, 2, 4, 3, "useful diagnostic, strongly order-dependent"),
    Gauge("raw pair-distance tournament", 1, 1, 2, 1, 4, 4, "pretty but loses embedding side channels"),
    Gauge("proof-obligation tournament", 4, 4, 4, 5, 3, 2, "vertices are cores, ears, directions, obstructions"),
)


def gauge_tournament() -> list[list[int]]:
    n = len(GAUGES)
    adj = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        ai = GAUGES[i].score_tuple()
        aj = GAUGES[j].score_tuple()
        iv = sum(x > y for x, y in zip(ai, aj))
        jv = sum(y > x for x, y in zip(ai, aj))
        if iv > jv or (iv == jv and i < j):
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj


def hp_count_tournament(adj: list[list[int]]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            val = dp[mask][last]
            if not val:
                continue
            for nxt in range(n):
                if not ((mask >> nxt) & 1) and adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += val
    return sum(dp[-1])


def print_gauge_tournament() -> None:
    adj = gauge_tournament()
    out = [sum(row) for row in adj]
    hist: dict[int, int] = {}
    for d in out:
        hist[d] = hist.get(d, 0) + 1
    order = sorted(range(len(GAUGES)), key=lambda i: (-out[i], GAUGES[i].name))
    print("Tournament Analysis over unit-distance-to-tournament gauges")
    print("rank | gauge | outscore | note")
    print("--- | --- | --- | ---")
    for rank, idx in enumerate(order, 1):
        print(f"{rank} | {GAUGES[idx].name} | {out[idx]} | {GAUGES[idx].note}")
    print(f"score histogram: {dict(sorted(hist.items()))}")
    print(f"directed 3-cycles: {sum(1 for i, j, k in combinations(range(len(adj)), 3) if (adj[i][j] and adj[j][k] and adj[k][i]) or (adj[i][k] and adj[k][j] and adj[j][i]))}")
    print(f"SCC sizes: {scc_sizes([sum((adj[i][j] << j) for j in range(len(adj))) for i in range(len(adj))])}")
    print(f"Hamiltonian path count: {hp_count_tournament(adj)}")
    print(
        "Challenged assumption: tournament vertices need not be points.  "
        "For proof work the better vertices can be unit directions, dense cores, "
        "frontier ears, traceability obligations, and totally-unfaithful filters."
    )
    print()


def analyze_cluster(size: int, cluster: Cluster, exact: int | None, count_paths: bool) -> dict[str, object]:
    adj, _ = unit_graph(cluster)
    traceable, path = hamiltonian_unit_spine(adj)
    row: dict[str, object] = {
        "n": size,
        "edges": edge_count(adj),
        "exact": exact,
        "traceable": traceable,
        "degree_hist": degree_hist(adj),
        "path": path,
    }
    if count_paths and size <= 14:
        row["unit_hp"] = count_unit_hamiltonian_paths(adj)

    if size <= 14 and traceable:
        lex = tuple(range(size))
        spine = path
        for label, order in (("lex", lex), ("spine", spine)):
            for flip_units in (True, False):
                out, unit_bits = orient_by_gauge(cluster, order, flip_units)
                key = f"{label}_{'unitflip' if flip_units else 'nonunitflip'}_max_unit_arcs"
                row[key] = max_unit_edges_in_directed_hp(out, unit_bits)
                if label == "spine" and size == 14:
                    row[f"{key}_triangles"] = directed_triangles(out)
                    row[f"{key}_scc"] = scc_sizes(out)
                    row[f"{key}_score_hist"] = score_hist(out)
    return row


def print_rows(rows: list[dict[str, object]]) -> None:
    print("Unit-spine scan on Moser-carrier beam witnesses")
    print("n | edges | target | unit spine? | unit HPs | deg hist | lex unit-flip max | lex nonunit-flip max")
    print("--- | --- | --- | --- | --- | --- | --- | ---")
    for row in rows:
        n = int(row["n"])
        unit_hp = row.get("unit_hp", "not counted")
        lex_u = row.get("lex_unitflip_max_unit_arcs", "skip")
        lex_nu = row.get("lex_nonunitflip_max_unit_arcs", "skip")
        print(
            f"{n} | {row['edges']} | {row['exact']} | {row['traceable']} | "
            f"{unit_hp} | {row['degree_hist']} | {lex_u} | {lex_nu}"
        )
    print()


def print_spine_detail(row: dict[str, object], cluster: Cluster) -> None:
    n = int(row["n"])
    path = tuple(row["path"])  # type: ignore[arg-type]
    coords = [cluster[i] for i in path[: min(8, len(path))]]
    suffix = " ..." if len(path) > 8 else ""
    print(f"Representative unit spine for n={n}")
    print(f"path indices: {path[: min(16, len(path))]}{suffix}")
    print(f"first coordinates on path: {coords}{suffix}")
    if n == 14:
        print("n=14 spine-order point-tournament fingerprints")
        for key, value in sorted(row.items()):
            if key.startswith("spine_") and (
                key.endswith("_max_unit_arcs")
                or key.endswith("_triangles")
                or key.endswith("_scc")
                or key.endswith("_score_hist")
            ):
                print(f"{key}: {value}")
    print()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", type=int, default=22)
    parser.add_argument("--width", type=int, default=1200)
    parser.add_argument("--count-through", type=int, default=14)
    args = parser.parse_args()

    print("S626 unit-distance unit-spine / flip-gauge tournament lab")
    print(f"Moser shell: {len(lab.UNITS)} directed unit vectors, {len(lab.PAIRS)} antipodal pairs")
    print(f"beam target={args.target}, width={args.width}")
    result = lab.run_beam(args.target, args.width)
    final_target = lab.U_EXACT.get(args.target, 60 if args.target == 22 else result.target)
    print(f"beam final: n={args.target}, best_edges={result.best_edges}, target={final_target}")
    print()

    rows: list[dict[str, object]] = []
    for size in range(2, min(args.target, 14) + 1):
        cluster = result.cluster_at(size)
        rows.append(analyze_cluster(size, cluster, lab.U_EXACT.get(size), size <= args.count_through))

    if args.target > 14:
        for size in [21, 22]:
            if size <= args.target:
                cluster = result.cluster_at(size)
                target = lab.U_EXACT.get(size)
                if size == 22:
                    target = 60
                rows.append(analyze_cluster(size, cluster, target, False))

    print_rows(rows)
    if rows:
        detail_size = 14 if args.target >= 14 else int(rows[-1]["n"])
        detail_cluster = result.cluster_at(detail_size)
        detail_row = next(r for r in rows if r["n"] == detail_size)
        print_spine_detail(detail_row, detail_cluster)
    if args.target >= 22:
        row22 = next((r for r in rows if r["n"] == 22), None)
        if row22:
            print_spine_detail(row22, result.cluster_at(22))

    print_gauge_tournament()
    print("Interpretation")
    print("- The intrinsic question is graph traceability of the unit-distance graph.")
    print("- A traceable unit graph lets either flip convention carry an all-unit mandatory path after choosing the tie order from the spine.")
    print("- A fixed external order, such as lexicographic coordinates, can obscure that spine; any apparent flip/flop there is a gauge artifact.")
    print("- A genuine flop would require an extremal unit-distance graph with no unit Hamiltonian path, likely forced by three or more essential boundary ears or separator branches.")


if __name__ == "__main__":
    main()
