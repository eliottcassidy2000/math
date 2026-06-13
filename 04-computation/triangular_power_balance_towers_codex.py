#!/usr/bin/env python3
"""Triangular power-balance towers and additive/square bridges.

codex-2026-06-12

The user's two towers are treated as interval power-balance problems.
For row n and a center C, compare

    (C-n)^p + ... + C^p       with       (C+1)^p + ... + (C+n)^p.

The ordinary tower is the p=1 solution C=2*T_n.  The consecutive-square
tower is the p=2 solution C=4*T_n.  This script verifies those identities,
searches crossover structure between the two interval families, and records
Tournament Analysis over possible proof carriers.

External anchors:
  OEIS A059270: ordinary consecutive-integer equal-sum tower.
  OEIS A059255: consecutive-square equal-sum tower.

The output is a reproducible hypothesis generator, not a proof of the
speculative transfer to LRC14 or the self-dual [72,36,16] code.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations
from math import isqrt

try:
    import sympy as sp
except ImportError:  # pragma: no cover - repo environment has sympy, keep fallback honest.
    sp = None


def T(n: int) -> int:
    return n * (n + 1) // 2


def sum_power(values: tuple[int, ...], p: int) -> int:
    return sum(v**p for v in values)


def S1(n: int) -> int:
    """Common ordinary sum in the first tower."""
    return n * (n + 1) * (2 * n + 1) // 2


def S2(n: int) -> int:
    """Common square sum in the second tower."""
    return n * (n + 1) * (2 * n + 1) * (12 * n * n + 12 * n + 1) // 6


def L2_unsquared(n: int) -> int:
    return n * (n + 1) * (4 * n + 3) // 2


def R2_unsquared(n: int) -> int:
    return n * (n + 1) * (4 * n + 1) // 2


@dataclass(frozen=True)
class TowerRow:
    name: str
    n: int
    center: int
    left: tuple[int, ...]
    right: tuple[int, ...]
    power: int

    @property
    def left_value(self) -> int:
        return sum_power(self.left, self.power)

    @property
    def right_value(self) -> int:
        return sum_power(self.right, self.power)

    @property
    def ordinary_left(self) -> int:
        return sum(self.left)

    @property
    def ordinary_right(self) -> int:
        return sum(self.right)


def first_row(n: int) -> TowerRow:
    center = 2 * T(n)
    return TowerRow(
        "F",
        n,
        center,
        tuple(range(center - n, center + 1)),
        tuple(range(center + 1, center + n + 1)),
        1,
    )


def second_row(n: int) -> TowerRow:
    center = 4 * T(n)
    return TowerRow(
        "Q",
        n,
        center,
        tuple(range(center - n, center + 1)),
        tuple(range(center + 1, center + n + 1)),
        2,
    )


def D_power(p: int, center: int, n: int) -> int:
    return sum((center - i) ** p for i in range(n + 1)) - sum(
        (center + i) ** p for i in range(1, n + 1)
    )


def symbolic_D_lines(max_p: int = 4) -> list[str]:
    if sp is None:
        return ["sympy unavailable; symbolic D_p formulas skipped"]
    C, n, i = sp.symbols("C n i", integer=True, positive=True)
    lines: list[str] = []
    for p in range(1, max_p + 1):
        expr = sp.summation((C - i) ** p, (i, 0, n)) - sp.summation(
            (C + i) ** p, (i, 1, n)
        )
        lines.append(f"p={p}: {sp.factor(expr)}")
    return lines


def check_power_centers(max_n: int = 40, max_p: int = 8) -> dict[int, dict[str, object]]:
    out: dict[int, dict[str, object]] = {}
    for p in range(1, max_p + 1):
        integer_centers: list[tuple[int, int]] = []
        bracket_failures: list[tuple[int, int, int, int]] = []
        for n in range(1, max_n + 1):
            if p == 1:
                center = 2 * T(n)
                if D_power(p, center, n) == 0:
                    integer_centers.append((n, center))
                continue
            if p == 2:
                center = 4 * T(n)
                if D_power(p, center, n) == 0:
                    integer_centers.append((n, center))
                continue
            base = 2 * p * T(n)
            low = D_power(p, base, n)
            high = D_power(p, base + 1, n)
            if not (low < 0 < high):
                bracket_failures.append((n, base, low, high))
        out[p] = {
            "integer_centers": integer_centers,
            "bracket_failures": bracket_failures,
        }
    return out


def first_side(n: int, side: str) -> tuple[int, ...]:
    row = first_row(n)
    return row.left if side == "FL" else row.right


def second_side(n: int, side: str) -> tuple[int, ...]:
    row = second_row(n)
    return row.left if side == "QL" else row.right


def boundary_first(n: int) -> dict[str, int]:
    row = first_row(n)
    return {
        "FL_start": row.left[0],
        "FL_end": row.left[-1],
        "FR_start": row.right[0],
        "FR_end": row.right[-1],
    }


def boundary_second(n: int) -> dict[str, int]:
    row = second_row(n)
    return {
        "QL_start": row.left[0],
        "QL_end": row.left[-1],
        "QR_start": row.right[0],
        "QR_end": row.right[-1],
    }


def interval_relation(a: tuple[int, ...], b: tuple[int, ...]) -> str:
    sa = set(a)
    sb = set(b)
    if sa == sb:
        return "equal"
    if sa <= sb:
        missing = [x for x in b if x not in sa]
        return f"left contained in right; right extra={missing}"
    if sb <= sa:
        missing = [x for x in a if x not in sb]
        return f"right contained in left; left extra={missing}"
    inter = [x for x in a if x in sb]
    return f"overlap={inter}"


def exact_interval_equalities(max_n: int = 100, max_m: int = 150) -> list[tuple[str, int, str, int, tuple[int, ...]]]:
    hits: list[tuple[str, int, str, int, tuple[int, ...]]] = []
    for n in range(1, max_n + 1):
        for q_side in ("QL", "QR"):
            q = second_side(n, q_side)
            for m in range(1, max_m + 1):
                for f_side in ("FL", "FR"):
                    f = first_side(m, f_side)
                    if q == f:
                        hits.append((q_side, n, f_side, m, q))
    return hits


def best_overlaps(max_n: int = 30, max_m: int = 80, take: int = 16) -> list[tuple[float, int, str, int, str, int, str]]:
    rows: list[tuple[float, int, str, int, str, int, str]] = []
    for n in range(1, max_n + 1):
        for q_side in ("QL", "QR"):
            q = second_side(n, q_side)
            sq = set(q)
            for m in range(1, max_m + 1):
                for f_side in ("FL", "FR"):
                    f = first_side(m, f_side)
                    sf = set(f)
                    inter = len(sq & sf)
                    if not inter:
                        continue
                    union = len(sq | sf)
                    rel = interval_relation(q, f)
                    if sq <= sf or sf <= sq or inter / union >= 0.50:
                        rows.append((inter / union, inter, q_side, n, f_side, m, rel))
    rows.sort(key=lambda r: (-r[0], -r[1], r[3], r[5], r[2], r[4]))
    return rows[:take]


def endpoint_hits(max_n: int = 100) -> list[tuple[int, str, int, str, int]]:
    hits: list[tuple[int, str, int, str, int]] = []
    for n in range(1, max_n + 1):
        second = boundary_second(n)
        max_value = max(second.values())
        max_m = isqrt(max_value) + 4
        for m in range(1, max_m + 1):
            first = boundary_first(m)
            for q_name, q_value in second.items():
                for f_name, f_value in first.items():
                    if q_value == f_value:
                        hits.append((q_value, q_name, n, f_name, m))
    return hits


def sum_crossover_hits(limit: int = 500) -> tuple[list[tuple[int, int, int]], list[tuple[int, int, int]]]:
    l_hits: list[tuple[int, int, int]] = []
    r_hits: list[tuple[int, int, int]] = []
    s1_index = {S1(m): m for m in range(1, limit + 1)}
    for n in range(1, limit + 1):
        if L2_unsquared(n) in s1_index:
            l_hits.append((n, s1_index[L2_unsquared(n)], L2_unsquared(n)))
        if R2_unsquared(n) in s1_index:
            r_hits.append((n, s1_index[R2_unsquared(n)], R2_unsquared(n)))
    return l_hits, r_hits


def special_value_hits(values: tuple[int, ...], limit: int = 500) -> dict[int, list[tuple[str, int]]]:
    functions = (
        ("S1", S1),
        ("L2", L2_unsquared),
        ("R2", R2_unsquared),
        ("S2", S2),
    )
    hits: dict[int, list[tuple[str, int]]] = {v: [] for v in values}
    for n in range(1, limit + 1):
        for name, func in functions:
            val = func(n)
            if val in hits:
                hits[val].append((name, n))
    return hits


def strongly_connected_components(vertices: list[str], edges: dict[tuple[str, str], str]) -> list[list[str]]:
    graph = {v: [] for v in vertices}
    for (a, b), winner in edges.items():
        loser = b if winner == a else a
        graph[winner].append(loser)

    index = 0
    stack: list[str] = []
    on_stack: set[str] = set()
    indices: dict[str, int] = {}
    lowlinks: dict[str, int] = {}
    sccs: list[list[str]] = []

    def visit(v: str) -> None:
        nonlocal index
        indices[v] = index
        lowlinks[v] = index
        index += 1
        stack.append(v)
        on_stack.add(v)
        for w in graph[v]:
            if w not in indices:
                visit(w)
                lowlinks[v] = min(lowlinks[v], lowlinks[w])
            elif w in on_stack:
                lowlinks[v] = min(lowlinks[v], indices[w])
        if lowlinks[v] == indices[v]:
            comp: list[str] = []
            while True:
                w = stack.pop()
                on_stack.remove(w)
                comp.append(w)
                if w == v:
                    break
            sccs.append(comp)

    for v in vertices:
        if v not in indices:
            visit(v)
    return sccs


def count_hamiltonian_paths(vertices: list[str], edges: dict[tuple[str, str], str]) -> int:
    n = len(vertices)
    idx = {v: i for i, v in enumerate(vertices)}
    can = [[False] * n for _ in range(n)]
    for (a, b), winner in edges.items():
        loser = b if winner == a else a
        can[idx[winner]][idx[loser]] = True

    dp: dict[tuple[int, int], int] = {}
    for i in range(n):
        dp[(1 << i, i)] = 1
    for mask in range(1 << n):
        for last in range(n):
            val = dp.get((mask, last), 0)
            if not val:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if can[last][nxt]:
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + val
    full = (1 << n) - 1
    return sum(dp.get((full, last), 0) for last in range(n))


def tournament_analysis() -> dict[str, object]:
    tie_path = [
        "power_center_rigidity",
        "78_90_code_shadow",
        "square_tower_multiplier_lift",
        "unsquared_defect_channel",
        "first_square_row_partition",
        "pell_endpoint_alignments",
        "overlap_atlas",
        "lrc14_worry_modulus_bridge",
        "convolution_boundary_lift",
    ]
    scores = {
        # exactness, novelty, transfer_to_lrc, transfer_to_72, computability, proof_potential
        "power_center_rigidity": (5, 4, 3, 2, 5, 5),
        "78_90_code_shadow": (4, 5, 3, 5, 5, 4),
        "square_tower_multiplier_lift": (5, 4, 3, 3, 5, 4),
        "unsquared_defect_channel": (5, 4, 4, 4, 5, 4),
        "first_square_row_partition": (5, 3, 4, 2, 5, 4),
        "pell_endpoint_alignments": (3, 4, 3, 2, 4, 3),
        "overlap_atlas": (3, 4, 3, 3, 5, 3),
        "lrc14_worry_modulus_bridge": (4, 3, 5, 2, 4, 4),
        "convolution_boundary_lift": (3, 5, 4, 4, 3, 5),
    }
    vertices = list(tie_path)
    edges: dict[tuple[str, str], str] = {}
    flips_vs_exactness = 0
    for a, b in combinations(vertices, 2):
        a_score = scores[a]
        b_score = scores[b]
        votes_a = sum(1 for x, y in zip(a_score, b_score) if x > y)
        votes_b = sum(1 for x, y in zip(a_score, b_score) if y > x)
        if votes_a > votes_b:
            winner = a
        elif votes_b > votes_a:
            winner = b
        else:
            winner = a if tie_path.index(a) < tie_path.index(b) else b
        edges[(a, b)] = winner
        exact_winner = a if a_score[0] >= b_score[0] else b
        if winner != exact_winner:
            flips_vs_exactness += 1

    outdegrees = Counter({v: 0 for v in vertices})
    for (a, b), winner in edges.items():
        outdegrees[winner] += 1
    score_hist = Counter(outdegrees.values())
    cycles = []
    for a, b, c in combinations(vertices, 3):
        wins = {
            (a, b): edges[tuple(sorted((a, b), key=vertices.index))],
            (a, c): edges[tuple(sorted((a, c), key=vertices.index))],
            (b, c): edges[tuple(sorted((b, c), key=vertices.index))],
        }
        if (
            wins[(a, b)] == a
            and wins[(b, c)] == b
            and wins[(a, c)] == c
        ) or (
            wins[(a, b)] == b
            and wins[(b, c)] == c
            and wins[(a, c)] == a
        ):
            cycles.append((a, b, c))
    sccs = strongly_connected_components(vertices, edges)
    ranking = sorted(vertices, key=lambda v: (-outdegrees[v], tie_path.index(v)))
    return {
        "criteria": (
            "exactness",
            "novelty",
            "transfer_to_lrc",
            "transfer_to_72",
            "computability",
            "proof_potential",
        ),
        "tie_path": tie_path,
        "scores": scores,
        "ranking": ranking,
        "outdegrees": dict(outdegrees),
        "score_hist": dict(sorted(score_hist.items())),
        "directed_3cycles": cycles,
        "scc_sizes": sorted((len(c) for c in sccs), reverse=True),
        "hamiltonian_paths": count_hamiltonian_paths(vertices, edges),
        "edge_flips_vs_exactness": flips_vs_exactness,
    }


def print_rows() -> None:
    print("[1] The two towers as exact interval balances")
    print("    OEIS anchors: A059270 (ordinary sums), A059255 (square sums).")
    print("    T_n=n(n+1)/2.  First row center C1=2T_n=n(n+1).")
    print("    Second row center C2=4T_n=2n(n+1).")
    print()
    for n in range(1, 8):
        f = first_row(n)
        q = second_row(n)
        print(
            f"n={n}: F_L={f.left} F_R={f.right} sum={f.left_value} "
            f"formula_S1={S1(n)}"
        )
        print(
            f"     Q_L={q.left} Q_R={q.right} square_sum={q.left_value} "
            f"formula_S2={S2(n)} ordinary_shadow=({q.ordinary_left},{q.ordinary_right})"
        )
        assert f.left_value == f.right_value == S1(n)
        assert q.left_value == q.right_value == S2(n)
        assert q.ordinary_left == L2_unsquared(n)
        assert q.ordinary_right == R2_unsquared(n)
        assert q.ordinary_left - q.ordinary_right == 2 * T(n)
        assert q.ordinary_left + q.ordinary_right == 4 * S1(n)
    print()
    print("    First-tower row n partitions the square shell [n^2,(n+1)^2-1].")
    print("    Second-tower unsquared defect: L2(n)-R2(n)=2T_n; total L2+R2=4*S1(n).")
    print()


def print_power_center_section() -> None:
    print("[2] Power-center rigidity")
    print("    D_p(C,n)=sum_{i=0..n}(C-i)^p - sum_{i=1..n}(C+i)^p")
    for line in symbolic_D_lines(4):
        print(f"    {line}")
    print()
    centers = check_power_centers()
    for p in range(1, 9):
        data = centers[p]
        if p in (1, 2):
            sample = data["integer_centers"][:5]
            tail = data["integer_centers"][-1]
            print(f"    p={p}: integer centers verified n<=40, first={sample}, last={tail}")
        else:
            failures = data["bracket_failures"]
            print(
                f"    p={p}: no positive integer centers; "
                f"D_p(2pT_n)<0<D_p(2pT_n+1) for n<=40; failures={len(failures)}"
            )
    print()
    print("    Candidate theorem: p=1 and p=2 are the exact integer triangular centers;")
    print("    for p>=3 the positive root is trapped between consecutive integers 2pT_n and 2pT_n+1.")
    print()


def print_crossover_section() -> None:
    print("[3] Crossovers between the square-row partition F and the square-balance tower Q")
    scripted = [
        ("Q_L(2)", second_side(2, "QL"), "F_L(3)", first_side(3, "FL")),
        ("Q_R(2)", second_side(2, "QR"), "F_R(3)", first_side(3, "FR")),
        ("Q_L(3)", second_side(3, "QL"), "F_R(4)", first_side(4, "FR")),
        ("Q_R(3)", second_side(3, "QR"), "F_L(5)", first_side(5, "FL")),
        ("Q_L(4)", second_side(4, "QL"), "F_L(6)", first_side(6, "FL")),
    ]
    for left_name, left, right_name, right in scripted:
        print(f"    {left_name}={left} vs {right_name}={right}: {interval_relation(left, right)}")
    print()
    exact = exact_interval_equalities()
    print(f"    Exact side equalities with n<=100,m<=150: {len(exact)}")
    for q_side, n, f_side, m, values in exact:
        print(f"      {q_side}({n}) = {f_side}({m}) = {values}")
    print("    Derivation note: full side equality forces the side lengths to match;")
    print("    the nontrivial solution is Q_L(3)=F_R(4).")
    print()
    print("    Best nontrivial overlaps/containments:")
    for j, inter, q_side, n, f_side, m, rel in best_overlaps():
        print(f"      J={j:.3f}, inter={inter}: {q_side}({n}) vs {f_side}({m}) -> {rel}")
    print()
    hits = endpoint_hits()
    print(f"    Endpoint boundary coincidences through Q-row n<=100: {len(hits)}")
    for value, q_name, n, f_name, m in hits[:30]:
        print(f"      {value}: {q_name}({n}) = {f_name}({m})")
    print("    These are Pell-style quadratic boundary hits, e.g. 2n^2+n=m^2 or")
    print("    2n^2+2n+1=m^2, depending on the endpoint labels.")
    print()


def print_sum_bridge_section() -> None:
    print("[4] Sum-shadow bridges and code/LRC beacons")
    l_hits, r_hits = sum_crossover_hits()
    print(f"    L2(n)=S1(m) hits n,m<=500: {l_hits}")
    print(f"    R2(n)=S1(m) hits n,m<=500: {r_hits}")
    print("    The unique checked hit is L2(3)=90=S1(4).")
    print("    At the same row, R2(3)=78=C(13,2), the existing [72,36,16] lambda_5 beacon.")
    values = (25, 27, 42, 78, 90, 91, 165, 365, 2030, 7230, 249849)
    hits = special_value_hits(values)
    for value in values:
        print(f"      value {value}: {hits[value]}")
    print()
    print("    Interpretation: row Q(3) is not a code construction, but it creates a")
    print("    tight additive/square shadow: 21+22+23+24 has ordinary sum 90 while")
    print("    25+26+27 has ordinary sum 78, and the squared equality binds them.")
    print()


def print_tournament_section() -> None:
    print("[5] Tournament Analysis")
    ta = tournament_analysis()
    print(f"    Criteria: {ta['criteria']}")
    print(f"    Tie Hamiltonian path: {ta['tie_path']}")
    print("    Vertex scores:")
    for name in ta["tie_path"]:
        print(f"      {name}: {ta['scores'][name]}")
    print(f"    Ranking by majority tournament: {ta['ranking']}")
    print(f"    Outdegrees: {ta['outdegrees']}")
    print(f"    score_hist={ta['score_hist']}")
    print(f"    directed_3cycles={len(ta['directed_3cycles'])}: {ta['directed_3cycles'][:6]}")
    print(f"    scc_sizes={ta['scc_sizes']}")
    print(f"    hamiltonian_paths={ta['hamiltonian_paths']}")
    print(f"    edge_flips_vs_exactness={ta['edge_flips_vs_exactness']}")
    print()
    print("    Assumption challenge:")
    print("      Considered vertices: rows, intervals, endpoints, centers, sums, defects,")
    print("      square-equality shadows, Pell boundary events, code design parameters,")
    print("      LRC shell moduli, Fourier modes, and proof obligations.")
    print("      Chosen quotient preserves the power-balance center, side endpoints,")
    print("      ordinary-sum defect, and transfer beacons. It destroys individual square")
    print("      residue data, support/incidence fibers, and any actual code/LRC witness.")
    print("      Challenged assumption: the tower vertices need not be numbers in the")
    print("      intervals; the more useful vertices may be obligations such as 'center")
    print("      root is integral', 'endpoint aligns', or 'support lift exists'.")
    print()


def main() -> None:
    print("Triangular power-balance towers")
    print("HYP-2454 / T798 / OPEN-Q-076 candidate")
    print()
    print_rows()
    print_power_center_section()
    print_crossover_section()
    print_sum_bridge_section()
    print_tournament_section()
    print("[6] Next proof tasks")
    print("    1. Prove the p>=3 bracket D_p(2pT_n)<0<D_p(2pT_n+1), or find the first exception.")
    print("    2. Solve the endpoint Pell families and classify infinite overlap/containment types.")
    print("    3. Turn Q(3)'s 78/90 shadow into a concrete [72,36,16] support-ledger constraint.")
    print("    4. Attach the same defect ledger to LRC14: compare 27-shell, 78, 90, and 91 resources.")
    print("    5. Lift from interval balances to convolution/support balances, following HYP-2452.")


if __name__ == "__main__":
    main()
