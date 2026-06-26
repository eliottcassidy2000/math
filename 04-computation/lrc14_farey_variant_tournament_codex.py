#!/usr/bin/env python3
"""Farey-variant tournament scout for LRC(14).

This is an exploratory proof-carrier ledger, not a proof.  It tests the user's
four Farey-sequence mutations at the LRC14 order:

  f1(a/b) = a*b
  f2(a/b) = a+b
  f3(a/b) = b**a
  f4(a/b) = a**b

The useful question is not whether these sequences themselves imply LRC.  The
question is which finite binary relations they create after projection to the
two LRC14 clocks:

  mod 7   : sector / Paley-Fano clock
  mod 27  : shell clock C=2n-1, the n=14 signed/worry modulus

For each transform we report local Farey-sequence behavior, Paley-mod-7
tournament statistics, folded mod-27 shell transition graph statistics, and
octahedral induced-subgraph counts.  A final hand-scored route tournament ranks
the current proof topics by preservation of the LRC predicate and actionability.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import gcd


N = 14
C = 2 * N - 1
QR7 = {1, 2, 4}


def farey_strict(order: int) -> list[tuple[int, int]]:
    out = []
    for b in range(2, order + 1):
        for a in range(1, b):
            if gcd(a, b) == 1:
                out.append((a, b))
    return sorted(out, key=lambda ab: Fraction(ab[0], ab[1]))


def transforms(a: int, b: int) -> dict[str, int]:
    return {
        "num_times_den": a * b,
        "num_plus_den": a + b,
        "den_pow_num": b**a,
        "num_pow_den": a**b,
    }


def folded_shell(x: int, modulus: int = C) -> int:
    r = x % modulus
    return min(r, modulus - r)


def unit_stratum(x: int, modulus: int = C) -> str:
    r = x % modulus
    if r == 0:
        return "zero"
    g = gcd(r, modulus)
    if g == 1:
        return "unit"
    return f"gcd{g}"


def directed_triangles(adj: list[list[int]]) -> int:
    total = 0
    for a, b, c in combinations(range(len(adj)), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        ):
            total += 1
    return total


def paley7_tournament(values: list[int]) -> tuple[list[list[int]], int]:
    """Blow up the Paley T7 relation along value residues.

    Equal residues are broken by Farey order.  The returned tie count is useful:
    high tie count means the transform lost sector information.
    """
    n = len(values)
    adj = [[0] * n for _ in range(n)]
    ties = 0
    residues = [v % 7 for v in values]
    for i, j in combinations(range(n), 2):
        d = (residues[j] - residues[i]) % 7
        if d == 0:
            adj[i][j] = 1
            ties += 1
        elif d in QR7:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj, ties


def graph_stats(edges: set[tuple[int, int]], vertices: set[int]) -> dict[str, object]:
    deg = Counter()
    for a, b in edges:
        deg[a] += 1
        deg[b] += 1
    for v in vertices:
        deg[v] += 0
    tri = 0
    for a, b, c in combinations(sorted(vertices), 3):
        if (
            tuple(sorted((a, b))) in edges
            and tuple(sorted((a, c))) in edges
            and tuple(sorted((b, c))) in edges
        ):
            tri += 1
    return {
        "vertices": len(vertices),
        "edges": len(edges),
        "degree_hist": dict(sorted(Counter(deg.values()).items())),
        "triangles": tri,
        "octahedral_induced": count_induced_octahedra(edges, vertices),
    }


def count_induced_octahedra(edges: set[tuple[int, int]], vertices: set[int]) -> int:
    """Count induced K6 minus a perfect matching subgraphs."""
    return len(induced_octahedra(edges, vertices))


def induced_octahedra(edges: set[tuple[int, int]], vertices: set[int]) -> list[tuple[int, ...]]:
    """Return induced K6 minus a perfect matching subgraphs."""
    out = []
    for six in combinations(sorted(vertices), 6):
        local_edges = 0
        local_deg = Counter()
        for a, b in combinations(six, 2):
            if tuple(sorted((a, b))) in edges:
                local_edges += 1
                local_deg[a] += 1
                local_deg[b] += 1
        if local_edges == 12 and all(local_deg[v] == 4 for v in six):
            out.append(six)
    return out


def shell_transition_graph(values: list[int]) -> tuple[set[int], set[tuple[int, int]], Counter]:
    shells = [folded_shell(v) for v in values]
    vertices = set(shells)
    edges: set[tuple[int, int]] = set()
    diff_strata: Counter = Counter()
    for x, y in zip(values, values[1:]):
        sx, sy = folded_shell(x), folded_shell(y)
        if sx != sy:
            edges.add(tuple(sorted((sx, sy))))
        diff_strata[unit_stratum(y - x)] += 1
    return vertices, edges, diff_strata


def longest_monotone_runs(values: list[int]) -> tuple[int, int, int]:
    up = down = flat = 1
    best_up = best_down = best_flat = 1
    for x, y in zip(values, values[1:]):
        if y > x:
            up += 1
            down = flat = 1
        elif y < x:
            down += 1
            up = flat = 1
        else:
            flat += 1
            up = down = 1
        best_up = max(best_up, up)
        best_down = max(best_down, down)
        best_flat = max(best_flat, flat)
    return best_up, best_down, best_flat


def model_graphs() -> None:
    print("\nMODEL GRAPH CHECKS")
    # Octahedron = K6 minus three opposite edges.
    oct_vertices = set(range(6))
    opposite = {tuple(sorted(e)) for e in [(0, 1), (2, 3), (4, 5)]}
    oct_edges = {
        tuple(sorted((a, b)))
        for a, b in combinations(oct_vertices, 2)
        if tuple(sorted((a, b))) not in opposite
    }
    print(f"  octahedron L(K4): {graph_stats(oct_edges, oct_vertices)}")

    # Clebsch = folded 5-cube.  Half-cube H_5 is its complement.
    verts = {min(x, x ^ 31) for x in range(32)}
    cleb_edges = set()
    half_edges = set()
    for a, b in combinations(sorted(verts), 2):
        d = (a ^ b).bit_count()
        fd = min(d, 5 - d)
        if fd == 1:
            cleb_edges.add((a, b))
        else:
            half_edges.add((a, b))
    print(f"  Clebsch folded 5-cube: {graph_stats(cleb_edges, verts)}")
    print(f"  half-cube complement:  {graph_stats(half_edges, verts)}")

    # Paley T7 reference.
    adj = [[0] * 7 for _ in range(7)]
    for i, j in combinations(range(7), 2):
        if (j - i) % 7 in QR7:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    print(
        "  Paley T7: "
        f"score_hist={dict(Counter(sum(row) for row in adj))}, "
        f"directed_triangles={directed_triangles(adj)}"
    )


@dataclass(frozen=True)
class Route:
    name: str
    preserves_predicate: int
    keeps_metric: int
    finite_exactness: int
    actionability: int
    risk: int

    @property
    def score(self) -> int:
        return (
            2 * self.preserves_predicate
            + 2 * self.keeps_metric
            + self.finite_exactness
            + self.actionability
            - 2 * self.risk
        )


def route_tournament() -> None:
    routes = [
        Route("metric winding tournament: order+gaps", 5, 5, 4, 5, 1),
        Route("covering-set split + equidistribution", 5, 4, 4, 5, 1),
        Route("k=9 Delsarte razor certificate", 4, 4, 5, 4, 2),
        Route("octahedral Hodge tail on L(K4)", 4, 3, 4, 4, 2),
        Route("Clebsch/half-cube cut covariance", 3, 3, 4, 3, 2),
        Route("Paley/Fano sector design", 3, 2, 4, 3, 2),
        Route("four Farey variants as classifiers", 2, 2, 3, 3, 2),
        Route("scalar H/Jensen/additive-energy only", 1, 1, 2, 1, 5),
    ]
    n = len(routes)
    adj = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        a, b = routes[i], routes[j]
        if (a.score, -a.risk, a.name) >= (b.score, -b.risk, b.name):
            adj[i][j] = 1
        else:
            adj[j][i] = 1

    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v in range(n):
            if not dp[mask][v]:
                continue
            for w in range(n):
                if (mask >> w) & 1:
                    continue
                if adj[v][w]:
                    dp[mask | (1 << w)][w] += dp[mask][v]
    hp = sum(dp[-1])

    print("\nROUTE TOURNAMENT")
    for idx, route in sorted(
        enumerate(routes), key=lambda item: (sum(adj[item[0]]), item[1].score), reverse=True
    ):
        print(
            f"  out={sum(adj[idx])} raw={route.score:2d} "
            f"{route.name}  features="
            f"({route.preserves_predicate},{route.keeps_metric},"
            f"{route.finite_exactness},{route.actionability},risk={route.risk})"
        )
    print(f"  directed_triangles={directed_triangles(adj)} hamiltonian_paths={hp}")


def main() -> None:
    farey = farey_strict(N)
    print(f"LRC14 Farey-variant tournament scout: strict Farey F_{N}, {len(farey)} fractions")
    print(f"First terms: {' '.join(f'{a}/{b}' for a, b in farey[:10])}")
    print(f"Last terms:  {' '.join(f'{a}/{b}' for a, b in farey[-10:])}")

    values_by_name: dict[str, list[int]] = defaultdict(list)
    for a, b in farey:
        for name, value in transforms(a, b).items():
            values_by_name[name].append(value)

    for name, vals in values_by_name.items():
        adjacent_desc = sum(1 for x, y in zip(vals, vals[1:]) if y < x)
        pair_inv = sum(1 for i, j in combinations(range(len(vals)), 2) if vals[j] < vals[i])
        best_up, best_down, best_flat = longest_monotone_runs(vals)
        paley_adj, sector_ties = paley7_tournament(vals)
        shells, shell_edges, diff_strata = shell_transition_graph(vals)
        print("\n" + "=" * 82)
        print(name)
        print("=" * 82)
        print(
            f"  raw sequence: adjacent_descents={adjacent_desc}, "
            f"pair_inversions={pair_inv}, longest_runs(up/down/flat)="
            f"{best_up}/{best_down}/{best_flat}"
        )
        print(
            f"  mod7 residues={dict(sorted(Counter(v % 7 for v in vals).items()))}, "
            f"sector_ties={sector_ties}, "
            f"Paley7_blowup_triangles={directed_triangles(paley_adj)}"
        )
        print(
            f"  mod27 shell_hist={dict(sorted(Counter(folded_shell(v) for v in vals).items()))}"
        )
        print(f"  adjacent diff strata mod27={dict(sorted(diff_strata.items()))}")
        print(f"  shell transition graph={graph_stats(shell_edges, shells)}")
        octs = induced_octahedra(shell_edges, shells)
        if octs:
            print(f"  induced octahedral shell witnesses={octs[:5]}")

    model_graphs()
    route_tournament()

    print("\nREADING")
    print("  Farey order supplies the perturbation grammar; mod 7 supplies the Paley/Fano")
    print("  sector clock; mod 27 supplies the signed shell clock.  The transforms are")
    print("  useful as packet classifiers only after these projections.  The proof route")
    print("  should keep metric gaps, then use Clebsch/half-cube cut labels and the")
    print("  octahedral L(K4) current tail as finite carriers for the remaining decorrelation.")


if __name__ == "__main__":
    main()
