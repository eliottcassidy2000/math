#!/usr/bin/env python3
"""
infinite_go_partition_recursion_s616.py

S616 scout: use Hamkins's infinite-Go ladder construction as a partition-
function lens on the repo's n -> n+2 recursions.

The MO update says that giving the main group n+2 liberties yields game value
omega*n: the two terminal liberties are the source/sink boundary, while each
extra liberty is one serial finite-announcement fuel cell.  A finite cutoff K
turns each fuel cell into the geometric packet 1 + q + ... + q^K; the ordinal
limit is the pole-order/serial-fuel limit rather than a scalar count.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations


def binom2(n: int) -> int:
    return n * (n - 1) // 2


def boundary_tiles_for_n_to_nplus2(inner_n: int) -> int:
    """Fixed-path tiling variables gained when an n-tournament becomes n+2."""
    return binom2(inner_n + 1) - binom2(inner_n - 1)


def geometric_packet_coeffs(k: int) -> list[int]:
    return [1] * (k + 1)


def convolve(a: list[int], b: list[int]) -> list[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def go_truncated_partition(phases: int, cutoff: int) -> list[int]:
    """Z_{phases,cutoff}(q)=(1+q+...+q^cutoff)^phases."""
    coeffs = [1]
    packet = geometric_packet_coeffs(cutoff)
    for _ in range(phases):
        coeffs = convolve(coeffs, packet)
    return coeffs


@dataclass(frozen=True)
class Lens:
    name: str
    evidence: int
    preserves_predicate: int
    recursion_signal: int
    partition_function_explicit: int
    cross_domain: int
    risk: int


LENSES = (
    Lens("Infinite Go serial ladder fuel", 5, 5, 5, 4, 4, 1),
    Lens("Tournament n->n+2 boundary OCF", 5, 5, 5, 5, 5, 1),
    Lens("LRC depth p0 / C=2n-1 resonance", 5, 5, 4, 5, 5, 1),
    Lens("Strong-component H semigroup gaps", 5, 4, 4, 4, 5, 1),
    Lens("Unit-distance arithmetic carrier", 4, 4, 3, 4, 4, 2),
    Lens("Collatz two-block residue partition", 3, 4, 3, 3, 4, 2),
    Lens("Raw scalar analogy without carrier", 1, 1, 1, 0, 2, 5),
)


def lens_votes(a: Lens, b: Lens) -> tuple[int, int]:
    criteria = [
        (a.evidence > b.evidence, b.evidence > a.evidence),
        (a.preserves_predicate > b.preserves_predicate, b.preserves_predicate > a.preserves_predicate),
        (a.recursion_signal > b.recursion_signal, b.recursion_signal > a.recursion_signal),
        (
            a.partition_function_explicit > b.partition_function_explicit,
            b.partition_function_explicit > a.partition_function_explicit,
        ),
        (a.cross_domain > b.cross_domain, b.cross_domain > a.cross_domain),
        (a.risk < b.risk, b.risk < a.risk),
    ]
    av = sum(1 for x, y in criteria if x and not y)
    bv = sum(1 for x, y in criteria if y and not x)
    return av, bv


def lens_tournament(lenses: tuple[Lens, ...]) -> list[list[int]]:
    n = len(lenses)
    adj = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        iv, jv = lens_votes(lenses[i], lenses[j])
        if iv > jv or (iv == jv and i < j):
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj


def score_hist(adj: list[list[int]]) -> dict[int, int]:
    hist: dict[int, int] = {}
    for row in adj:
        score = sum(row)
        hist[score] = hist.get(score, 0) + 1
    return dict(sorted(hist.items()))


def directed_triangles(adj: list[list[int]]) -> int:
    total = 0
    for i, j, k in combinations(range(len(adj)), 3):
        if adj[i][j] and adj[j][k] and adj[k][i]:
            total += 1
        if adj[i][k] and adj[k][j] and adj[j][i]:
            total += 1
    return total


def scc_sizes(adj: list[list[int]]) -> list[int]:
    n = len(adj)

    def reach(starts: list[int], forward: bool) -> set[int]:
        seen = set(starts)
        stack = list(starts)
        while stack:
            v = stack.pop()
            for w in range(n):
                edge = adj[v][w] if forward else adj[w][v]
                if edge and w not in seen:
                    seen.add(w)
                    stack.append(w)
        return seen

    left = set(range(n))
    sizes: list[int] = []
    while left:
        v = next(iter(left))
        comp = reach([v], True) & reach([v], False)
        sizes.append(len(comp))
        left -= comp
    return sorted(sizes, reverse=True)


def hamiltonian_paths(adj: list[list[int]]) -> int:
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
                if not (mask >> nxt) & 1 and adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += val
    return sum(dp[-1])


def print_go_partition_ledger() -> None:
    print("Infinite-Go finite-cutoff partition ledger")
    print("  Hamkins update: n+2 liberties -> game value omega*n.")
    print("  S616 reading: two terminal liberties are source/sink boundary;")
    print("    phases = liberties - 2 are serial finite-announcement fuel cells.")
    cutoff = 6
    print(f"  Finite cutoff model: Z_r,K(q)=(1+q+...+q^K)^r with K={cutoff}.")
    print(f"  {'liberties':>9} {'phases':>6} {'ordinal':>10} {'deg Z':>6} {'Z(1)':>8} {'top coeffs':>18}")
    for liberties in range(3, 9):
        phases = liberties - 2
        coeffs = go_truncated_partition(phases, cutoff)
        top = coeffs[-min(5, len(coeffs)) :]
        print(
            f"  {liberties:>9} {phases:>6} {'omega*' + str(phases):>10} "
            f"{len(coeffs)-1:>6} {sum(coeffs):>8} {str(top):>18}"
        )
    print("  Pole-order reading: as K grows, each phase contributes one geometric factor;")
    print("    the serial count of factors is the ordinal coefficient in omega*r.")
    print()


def print_nplus2_ledger() -> None:
    print("n -> n+2 boundary and LRC resonance ledger")
    print(f"  {'inner n':>7} {'Go liberties':>12} {'Go phases':>9} {'new tiles':>9} {'LRC C=2n-1':>12}")
    for inner_n in range(3, 15):
        new_tiles = boundary_tiles_for_n_to_nplus2(inner_n)
        c = 2 * inner_n - 1
        print(f"  {inner_n:>7} {inner_n+2:>12} {inner_n:>9} {new_tiles:>9} {c:>12}")
    print("  Exact identity: boundary_tiles(n -> n+2) = C(n+1,2)-C(n-1,2) = 2n-1.")
    print("  This is the same odd shell C used by the LRC floor resonance observer.")
    print()


def print_partition_dictionary() -> None:
    print("Partition-function dictionary")
    rows = [
        ("Infinite Go", "Z_r,K(q)=(1+...+q^K)^r", "life/death ordinal fuel", "pole order r"),
        ("Tournaments", "I(Omega(T),2)", "Hamiltonian path count H", "forbidden H=7,21"),
        ("LRC", "Z_delta(y)=sum p_k y^k", "ground cell p0", "2n-1 shell cancellation"),
        ("Unit distance", "carrier incidence / energy sum", "many unit pairs", "CM split-prime amplification"),
        ("Collatz", "two-block residue packet", "no live cycle state", "2^E-3^k determinant wall"),
    ]
    for name, z, target, obstruction in rows:
        print(f"  {name:14s} Z={z:36s} target={target:24s} obstruction={obstruction}")
    print()


def print_lens_tournament() -> None:
    print("Tournament Analysis on partition-function lenses")
    adj = lens_tournament(LENSES)
    scores = [sum(row) for row in adj]
    for lens, score in sorted(zip(LENSES, scores), key=lambda x: (-x[1], x[0].name)):
        print(
            f"  score={score} lens={lens.name}; "
            f"features=(evidence={lens.evidence}, preserve={lens.preserves_predicate}, "
            f"recursion={lens.recursion_signal}, explicitZ={lens.partition_function_explicit}, "
            f"cross={lens.cross_domain}, risk={lens.risk})"
        )
    print(f"  score histogram: {score_hist(adj)}")
    print(f"  directed 3-cycles: {directed_triangles(adj)}")
    print(f"  SCC sizes: {scc_sizes(adj)}")
    print(f"  Hamiltonian path count: {hamiltonian_paths(adj)}")
    print("  Reading: the best transfer keeps the partition function and its boundary state,")
    print("    not the visible scalar alone.")
    print()


def main() -> None:
    print("==== S616 infinite Go / n+2 recursion / partition-function scout ====")
    print_go_partition_ledger()
    print_nplus2_ledger()
    print_partition_dictionary()
    print_lens_tournament()
    print("Hypothesis HYP-2185:")
    print("  Infinite Go supplies the open-game analogue of the repo n->n+2 boundary recursion.")
    print("  The shared invariant is a partition function with a retained boundary/coimage state.")
    print("  The equality boundary_tiles = 2n-1 ties THM-291 to the LRC floor resonance C=2n-1.")
    print("  Extra liberties beyond the two terminal liberties are serial fuel cells, omega*n.")


if __name__ == "__main__":
    main()
