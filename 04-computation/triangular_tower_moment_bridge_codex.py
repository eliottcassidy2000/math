"""Triangular tower moment bridge.

The user supplied two adjacent-block towers:

  A_n: n^2 + ... + (n^2+n) = (n^2+n+1) + ... + (n^2+2n)
       as an equality of first moments.

  B_n: T_{2n}^2 + ... + (T_{2n}+n)^2
       = (T_{2n}+n+1)^2 + ... + (T_{2n}+2n)^2
       as an equality of second moments.

This script records exact formulas, overlap/crossover data, a higher-moment
extension, and a Tournament Analysis over proof routes.  It is exploratory:
the exact identities are elementary; the LRC transfer is a carrier hypothesis.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import isclose

import sympy as sp


def tri(n: int) -> int:
    return n * (n + 1) // 2


def interval_sum(lo: int, hi: int, power: int = 1) -> int:
    return sum(k**power for k in range(lo, hi + 1))


@dataclass(frozen=True)
class Side:
    tower: str
    n: int
    side: str
    lo: int
    hi: int

    @property
    def label(self) -> str:
        return f"{self.tower}{self.n}.{self.side}"

    @property
    def length(self) -> int:
        return self.hi - self.lo + 1

    @property
    def interval(self) -> tuple[int, int]:
        return (self.lo, self.hi)


def tower_a_sides(n: int) -> list[Side]:
    return [
        Side("A", n, "L", n * n, n * n + n),
        Side("A", n, "R", n * n + n + 1, n * n + 2 * n),
    ]


def tower_b_sides(n: int) -> list[Side]:
    start = 2 * n * n + n
    return [
        Side("B", n, "L", start, start + n),
        Side("B", n, "R", start + n + 1, start + 2 * n),
    ]


def verify_towers(limit: int = 8) -> list[str]:
    lines: list[str] = []
    lines.append("EXACT TOWER IDENTITIES")
    lines.append("======================")
    lines.append("A_n splits the square shell [n^2,(n+1)^2-1] into equal first moments.")
    lines.append("B_n splits the triangular shell [T_{2n},T_{2n+1}-1] into equal second moments.")
    lines.append("")
    for n in range(1, limit + 1):
        a_l, a_r = tower_a_sides(n)
        b_l, b_r = tower_b_sides(n)
        a_left = interval_sum(a_l.lo, a_l.hi)
        a_right = interval_sum(a_r.lo, a_r.hi)
        b_left_sq = interval_sum(b_l.lo, b_l.hi, 2)
        b_right_sq = interval_sum(b_r.lo, b_r.hi, 2)
        b_left_raw = interval_sum(b_l.lo, b_l.hi)
        b_right_raw = interval_sum(b_r.lo, b_r.hi)
        assert a_left == a_right
        assert b_left_sq == b_right_sq
        assert b_left_raw - b_right_raw == n * (n + 1)
        lines.append(
            f"n={n:2d} A: {a_l.interval} = {a_r.interval}, "
            f"sum={a_left}; B^2: {b_l.interval} = {b_r.interval}, "
            f"sumsq={b_left_sq}, raw_defect={b_left_raw-b_right_raw}"
        )
    lines.append("")

    n = sp.symbols("n", positive=True, integer=True)
    a_common = n * (n + 1) * (2 * n + 1) / 2
    b_common = n * (n + 1) * (2 * n + 1) * (12 * n**2 + 12 * n + 1) / 6
    b_raw_l = n * (n + 1) * (4 * n + 3) / 2
    b_raw_r = n * (4 * n**2 + 5 * n + 1) / 2
    lines.append("Closed forms")
    lines.append(f"  A_common(n) = {sp.factor(a_common)} = 3 * sum_{{j<=n}} j^2")
    lines.append(f"  B_common_square(n) = {sp.factor(b_common)}")
    lines.append(f"  B_raw_left(n) = {sp.factor(b_raw_l)}")
    lines.append(f"  B_raw_right(n) = {sp.factor(b_raw_r)}")
    lines.append(f"  B_raw_left - B_raw_right = {sp.factor(b_raw_l-b_raw_r)} = 2*T_n")
    lines.append("")
    return lines


def overlap_atlas(limit: int = 80) -> list[str]:
    sides: list[Side] = []
    for n in range(1, limit + 1):
        sides.extend(tower_a_sides(n))
        sides.extend(tower_b_sides(n))

    exact: list[tuple[Side, Side]] = []
    containments: list[tuple[Side, Side, tuple[int, int]]] = []
    partials: list[tuple[Side, Side, tuple[int, int]]] = []

    for x, y in combinations(sides, 2):
        if x.tower == y.tower:
            continue
        lo = max(x.lo, y.lo)
        hi = min(x.hi, y.hi)
        if lo > hi:
            continue
        if x.interval == y.interval:
            exact.append((x, y))
        elif x.lo <= y.lo and y.hi <= x.hi:
            containments.append((y, x, (y.lo - x.lo, x.hi - y.hi)))
        elif y.lo <= x.lo and x.hi <= y.hi:
            containments.append((x, y, (x.lo - y.lo, y.hi - x.hi)))
        else:
            partials.append((x, y, (lo, hi)))

    lines: list[str] = []
    lines.append("OVERLAP / CROSSOVER ATLAS")
    lines.append("==========================")
    lines.append(f"Scan range: n<= {limit} in both towers.")
    lines.append(f"Exact equal side intervals: {len(exact)}")
    for x, y in exact:
        lines.append(f"  {x.label} = {y.label} = {x.interval}")
    lines.append("")
    lines.append("First containments (inner in outer, with left/right padding):")
    for inner, outer, pad in containments[:18]:
        lines.append(
            f"  {inner.label} {inner.interval} in {outer.label} {outer.interval}; padding={pad}"
        )
    lines.append(f"  ... total containments: {len(containments)}")
    lines.append("")
    lines.append("First partial overlaps:")
    for x, y, overlap in partials[:16]:
        lines.append(f"  {x.label} {x.interval} overlaps {y.label} {y.interval} on {overlap}")
    lines.append(f"  ... total partial overlaps: {len(partials)}")
    lines.append("")
    lines.append("Diophantine exactness check for whole-side equality")
    lines.append("  B_L(n)=A_R(m) forces m=n+1 and n^2-2n-3=0, so n=3.")
    lines.append("  The other length/start equations have no positive integer solution.")
    lines.append("  Thus the exact hinge block 21..24 is unique among whole sides.")
    lines.append("")
    return lines


def positive_real_root(poly: sp.Expr, var: sp.Symbol) -> float | None:
    roots = sp.nroots(poly, n=30, maxsteps=200)
    real_roots = sorted(float(sp.re(r)) for r in roots if abs(float(sp.im(r))) < 1e-18)
    for r in real_roots:
        if r >= 0:
            return r
    return real_roots[-1] if real_roots else None


def moment_extension(max_power: int = 8, n_probe: int = 30) -> list[str]:
    A, n, i, c = sp.symbols("A n i c", positive=True)
    lines: list[str] = []
    lines.append("HIGHER-MOMENT EXTENSION")
    lines.append("=======================")
    lines.append(
        "For power p, solve sum_{i=0}^n (A+i)^p = "
        "sum_{i=0}^{n-1} (A+n+1+i)^p."
    )
    lines.append("p=1 and p=2 give the two exact integer towers; p>=3 exposes a fractional address.")
    lines.append("")

    for p in range(1, max_power + 1):
        expr = sp.summation((A + i) ** p, (i, 0, n)) - sp.summation(
            (A + n + 1 + i) ** p, (i, 0, n - 1)
        )
        if p <= 4:
            lines.append(f"p={p} defect factor: {sp.factor(expr)}")

    lines.append("")
    lines.append("Positive root compared to p*n^2+(p-1)*n at n=30:")
    for p in range(1, max_power + 1):
        Ap = sp.symbols("Ap")
        poly = sum((Ap + j) ** p for j in range(n_probe + 1)) - sum(
            (Ap + n_probe + 1 + j) ** p for j in range(n_probe)
        )
        root = positive_real_root(poly, Ap)
        baseline = p * n_probe * n_probe + (p - 1) * n_probe
        asym_const = Fraction((p - 1) * (p - 2), 12 * p)
        if root is None:
            lines.append(f"  p={p}: no real root found")
        else:
            lines.append(
                f"  p={p}: root={root:.12f}, root-baseline={root-baseline:.12f}, "
                f"asymptotic constant={(float(asym_const)):.12f} ({asym_const})"
            )

    lines.append("")
    lines.append("Asymptotic check:")
    lines.append("  Substitute A = p*n^2 + (p-1)*n + c.")
    lines.append("  The leading coefficient vanishes at c = (p-1)(p-2)/(12p).")
    for p in range(3, max_power + 1):
        expr = sp.summation((A + i) ** p, (i, 0, n)) - sp.summation(
            (A + n + 1 + i) ** p, (i, 0, n - 1)
        )
        substituted = sp.expand(expr.subs(A, p * n**2 + (p - 1) * n + c))
        poly_n = sp.Poly(substituted, n)
        top_degree = poly_n.degree()
        top_coeff = sp.factor(poly_n.coeff_monomial(n**top_degree))
        root_c = sp.solve(top_coeff, c)[0]
        lines.append(f"  p={p}: top degree {top_degree}, top coeff {top_coeff}, c={root_c}")
    lines.append("")
    return lines


def directed_3cycles(adj: list[list[bool]]) -> int:
    total = 0
    n = len(adj)
    for a, b, c in combinations(range(n), 3):
        if adj[a][b] and adj[b][c] and adj[c][a]:
            total += 1
        if adj[a][c] and adj[c][b] and adj[b][a]:
            total += 1
    return total


def scc_sizes(adj: list[list[bool]]) -> list[int]:
    n = len(adj)

    def reach(start: int, graph: list[list[bool]]) -> set[int]:
        seen = {start}
        stack = [start]
        while stack:
            v = stack.pop()
            for w, edge in enumerate(graph[v]):
                if edge and w not in seen:
                    seen.add(w)
                    stack.append(w)
        return seen

    radj = [[adj[j][i] for j in range(n)] for i in range(n)]
    unseen = set(range(n))
    sizes: list[int] = []
    while unseen:
        v = next(iter(unseen))
        comp = reach(v, adj) & reach(v, radj)
        sizes.append(len(comp))
        unseen -= comp
    return sorted(sizes, reverse=True)


def count_hamiltonian_paths(adj: list[list[bool]]) -> int:
    n = len(adj)
    dp: dict[tuple[int, int], int] = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp.get((mask, last), 0)
            if not count:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + count
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def tournament_analysis() -> list[str]:
    criteria = (
        "exactness",
        "computability",
        "hidden_lift",
        "LRC_transfer",
        "fraction_address",
        "sequence_anchor",
    )
    routes: dict[str, tuple[int, ...]] = {
        "moment1_square_shell": (5, 5, 3, 4, 2, 5),
        "moment2_triangular_shell": (5, 5, 4, 4, 3, 5),
        "bernoulli_fraction_address": (4, 4, 5, 4, 5, 3),
        "unique_hinge_overlap": (4, 5, 4, 3, 4, 4),
        "convolution_factor_lift_bridge": (5, 4, 5, 5, 4, 4),
        "lrc_resource_ledger": (3, 4, 5, 5, 5, 3),
        "raw_oeis_sequence_match": (2, 5, 2, 2, 1, 5),
    }
    names = list(routes)
    n = len(names)
    adj = [[False] * n for _ in range(n)]

    def beats(a: tuple[int, ...], b: tuple[int, ...], tie: bool) -> bool:
        wins = sum(x > y for x, y in zip(a, b))
        losses = sum(x < y for x, y in zip(a, b))
        if wins != losses:
            return wins > losses
        if sum(a) != sum(b):
            return sum(a) > sum(b)
        return tie

    for i, j in combinations(range(n), 2):
        if beats(routes[names[i]], routes[names[j]], True):
            adj[i][j] = True
        else:
            adj[j][i] = True

    scores = [sum(row) for row in adj]
    ranking = [name for _, name in sorted(zip(scores, names), reverse=True)]

    lines: list[str] = []
    lines.append("TOURNAMENT ANALYSIS")
    lines.append("===================")
    lines.append(f"Pairwise observable: majority win over criteria {criteria}.")
    lines.append("Switch/gauge: orient toward the route retaining more proof-bearing side-channel data.")
    lines.append("Tie Hamiltonian path: listed route order, then total score.")
    lines.append("Vertices are proof routes, not runners or arcs.")
    lines.append(f"routes = {routes}")
    lines.append(f"score_hist = {dict(sorted(Counter(scores).items()))}")
    lines.append(f"directed_3cycles = {directed_3cycles(adj)}")
    lines.append(f"scc_sizes = {scc_sizes(adj)}")
    lines.append(f"hamiltonian_paths = {count_hamiltonian_paths(adj)}")
    lines.append(f"ranking = {ranking}")
    lines.append("")
    return lines


def assumption_challenge() -> list[str]:
    return [
        "ASSUMPTION CHALLENGE",
        "====================",
        "Alternate vertex sets considered: tower equations, block sides, shell gaps,",
        "moment powers, fractional starts, OEIS sequences, LRC denominators, runner",
        "resources, hidden convolution grid cells, and proof-route obligations.",
        "Chosen quotient: adjacent unequal block splits of one shell, decorated by",
        "the first moment that equalizes.",
        "Preserved: exact shell interval geometry, first/second moment equality,",
        "overlap/hinge data, and the fractional carry needed at higher powers.",
        "Destroyed: runner positions, actual LRC circle arcs, Galois/root data, and",
        "large-scale asymptotics beyond the leading Bernoulli address. These must",
        "be reattached before any proof claim.",
        "Challenged assumption: the useful LRC clock need not be time t or runner",
        "speed.  Here the clock is the moment order p and the hidden address that",
        "makes two adjacent unequal blocks balance.",
        "",
    ]


def main() -> None:
    lines: list[str] = []
    lines.append("TRIANGULAR TOWER MOMENT BRIDGE")
    lines.append("==============================")
    lines.append("HYP-2453/T797/OPEN-Q-075 exploratory computation.")
    lines.append("")
    lines.extend(verify_towers())
    lines.extend(overlap_atlas())
    lines.extend(moment_extension())
    lines.extend(tournament_analysis())
    lines.extend(assumption_challenge())
    lines.append("NEXT DIRECTIONS")
    lines.append("===============")
    lines.append("1. Prove the general higher-moment asymptotic expansion to several terms and")
    lines.append("   identify which Bernoulli denominators can ever become integral after a")
    lines.append("   shell change or an added address coordinate.")
    lines.append("2. Turn containments between B-blocks and A-square shells into a Pell/Beatty")
    lines.append("   address system; this is where multiplication-by-2 meets square-root drift.")
    lines.append("3. Transfer to LRC14 by replacing scalar q-blocked counts with moment/resource")
    lines.append("   ledgers: blocked twists, owners, fiber, and the fractional address needed")
    lines.append("   to lift a scalar shell equality to a support proof.")
    lines.append("4. Compare with HYP-2452: coefficient convolution lifts and moment towers are")
    lines.append("   both boundary-total shadows of hidden 2D/filtered support.")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
