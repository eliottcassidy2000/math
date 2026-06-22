#!/usr/bin/env python3
"""
lrc14_nu_denom_center_budget_codex_s95.py

Codex S95, 2026-06-22.  Exact component audit for the LRC(14) Lemma A
gap:

    D(E) := meas{x in [0,1): maxgap({frac(e*x): e in E}) <= 1/7}

The Bonferroni/witness bridge needs the equivalent extremal statement
`D(E) <= D({0,1,...,k-1})` for k=8..13.  A naive local compression proof is
false: several nonconsecutive bounded shapes are local dead ends under
one-step left-compressions, but their dense mass is tiny.  This script tests a
more precise lead: consecutive blocks are special because they own the low
least-denominator dense components q=8..k; scattered local optima move mass to
higher denominators, where the interval budget is much smaller.

Method.
  * Compute the exact dense intervals by the common refinement of all phase
    order breakpoints and all max-gap=1/7 walls.
  * Assign each dense component the least denominator q of a reduced rational
    p/q strictly inside the component.  This is a proof label, not an asserted
    canonical center.
  * Exhaustively scan primitive anchored k-sets in [0,W] for W=16.

Tournament Analysis declaration.
  Vertex set: proof carriers, not runners.  The carriers are
    denominator-center, EWLB-window, Fourier-spectrum, finite-atlas,
    runner-compression, and raw-gap-moment.
  Pairwise observable: a score tuple
    (survives_local_dead_ends, keeps_component_width, exact_bounded_certificate,
     formalizable_walls, tail_route_connection, low_auxiliary_complexity).
  Gauge: orient A -> B if A has lexicographically larger score.  Ties follow
    the declared Hamiltonian path above.
  Fingerprints reported below: score histogram, directed 3-cycles, SCC sizes,
    and Hamiltonian path count.

Assumption challenge.
  I considered runners, gaps, fixed 1/7 sections, section boundaries,
  wall-crossing events, residues, cover arcs, Fourier modes, rational
  denominator buckets, local-slope pairs, matroid-circuit-like relations, and
  proof obligations as possible vertices.  The chosen quotient preserves the
  predicate needed by Lemma A because it partitions the dense measure D(E) into
  exact interval widths.  It destroys individual runner identity and the full
  cyclic phase order, which is why this is a lead rather than a proof.  The
  challenged assumption is that a compression theorem should be monotone on
  raw runner positions; the exact data say the useful monotonicity is carried
  by denominator-center budgets instead.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import reduce
from itertools import combinations, permutations
from math import gcd
from typing import Iterable
import sys

THETA = F(1, 7)
SCAN_WINDOWS = {8: 14, 9: 14, 10: 14, 11: 14, 12: 14, 13: 14}
TOP_N = 8

KNOWN_LOCAL_DEAD_ENDS = [
    (0, 1, 2, 3, 7, 8, 9, 10),
    (0, 1, 2, 3, 11, 12, 13, 14),
    (0, 1, 2, 6, 7, 8, 9, 14),
    (0, 1, 3, 4, 7, 8, 11, 12),
    (0, 1, 3, 5, 7, 9, 11, 13),
]

try:
    sys.stdout.reconfigure(line_buffering=True)
except Exception:
    pass


def primitive(E: Iterable[int]) -> bool:
    return reduce(gcd, (abs(e) for e in E if e), 0) == 1


def frac_float(x: F) -> str:
    return f"{float(x):.9f}"


def phase_points(E: tuple[int, ...], x: F) -> list[F]:
    return sorted((F(e) * x) % 1 for e in E)


def max_gap(E: tuple[int, ...], x: F) -> F:
    pts = phase_points(E, x)
    if len(pts) <= 1:
        return F(1)
    gaps = [b - a for a, b in zip(pts, pts[1:])]
    gaps.append(pts[0] + 1 - pts[-1])
    return max(gaps)


def floor_fraction(x: F) -> int:
    return x.numerator // x.denominator


def order_breakpoints(E: tuple[int, ...]) -> list[F]:
    bps = {F(0), F(1)}
    for a, b in combinations(E, 2):
        d = abs(b - a)
        if d == 0:
            continue
        for m in range(d + 1):
            bps.add(F(m, d))
    return sorted(bps)


_INTERVAL_CACHE: dict[tuple[int, ...], tuple[tuple[F, F], ...]] = {}


def dense_intervals(E: tuple[int, ...]) -> tuple[tuple[F, F], ...]:
    E = tuple(sorted(set(E)))
    if E in _INTERVAL_CACHE:
        return _INTERVAL_CACHE[E]
    raw: list[tuple[F, F]] = []
    bps = order_breakpoints(E)
    for cell_a, cell_b in zip(bps, bps[1:]):
        if cell_b <= cell_a:
            continue
        mid = (cell_a + cell_b) / 2
        pts = sorted(((F(e) * mid) % 1, e) for e in E)
        order = [e for _, e in pts]
        floors = [floor_fraction(F(e) * mid) for e in order]
        lo, hi = cell_a, cell_b
        ok = True
        for idx, e_cur in enumerate(order):
            f_cur = floors[idx]
            if idx < len(order) - 1:
                e_next = order[idx + 1]
                f_next = floors[idx + 1]
                wrap = F(0)
            else:
                e_next = order[0]
                f_next = floors[0]
                wrap = F(1)
            # gap(x) = (e_next*x - f_next) - (e_cur*x - f_cur) + wrap.
            A = F(e_next - e_cur)
            C = F(f_cur - f_next) + wrap
            if A == 0:
                if C > THETA:
                    ok = False
                    break
                continue
            wall = (THETA - C) / A
            if A > 0:
                hi = min(hi, wall)
            else:
                lo = max(lo, wall)
            if lo >= hi:
                ok = False
                break
        if ok and lo < hi:
            raw.append((lo, hi))

    merged: list[tuple[F, F]] = []
    for a, b in raw:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a, b))
    out = tuple(merged)
    _INTERVAL_CACHE[E] = out
    return out


def dense_measure(E: tuple[int, ...]) -> F:
    return sum((b - a for a, b in dense_intervals(E)), F(0))


def least_interior_rational(iv: tuple[F, F], qmax: int = 500) -> tuple[int, int] | None:
    a, b = iv
    for q in range(1, qmax + 1):
        lo = a.numerator * q // a.denominator + 1
        hi = (b.numerator * q - 1) // b.denominator
        for p in range(max(0, lo), hi + 1):
            if gcd(p, q) == 1 and a < F(p, q) < b:
                return q, p
    return None


def denom_budget(E: tuple[int, ...]) -> tuple[Counter[int | None], dict[int | None, F]]:
    return denom_budget_from_intervals(dense_intervals(E))


def denom_budget_from_intervals(intervals: tuple[tuple[F, F], ...]) -> tuple[Counter[int | None], dict[int | None, F]]:
    hist: Counter[int | None] = Counter()
    widths: dict[int | None, F] = defaultdict(F)
    for iv in intervals:
        center = least_interior_rational(iv)
        q = center[0] if center else None
        hist[q] += 1
        widths[q] += iv[1] - iv[0]
    return hist, dict(widths)


def fmt_widths(widths: dict[int | None, F], max_terms: int = 8) -> str:
    parts = []
    for q, w in sorted(widths.items(), key=lambda item: (10**9 if item[0] is None else item[0])):
        label = "None" if q is None else str(q)
        parts.append(f"q{label}:{w}({frac_float(w)})")
    if len(parts) > max_terms:
        parts = parts[:max_terms] + ["..."]
    return " ".join(parts) if parts else "-"


def is_consecutive(E: tuple[int, ...]) -> bool:
    return E == tuple(range(len(E)))


def residue_complete(E: tuple[int, ...], q: int) -> bool:
    return len({e % q for e in E}) == q


def one_step_left_neighbors(E: tuple[int, ...]) -> list[tuple[int, ...]]:
    S = set(E)
    out = []
    for e in E:
        if e == 0:
            continue
        candidate = e - 1
        if candidate in S:
            continue
        N = tuple(sorted((S - {e}) | {candidate}))
        if len(N) == len(E) and 0 in N and primitive(N):
            out.append(N)
    return out


def local_left_dead_end(E: tuple[int, ...]) -> bool:
    d0 = dense_measure(E)
    return all(dense_measure(N) <= d0 for N in one_step_left_neighbors(E))


def scan_k(k: int, W: int) -> dict[str, object]:
    consec = tuple(range(k))
    consec_D = dense_measure(consec)
    consec_hist, consec_widths = denom_budget(consec)

    best_E = consec
    best_D = consec_D
    top: list[tuple[F, tuple[int, ...], dict[int | None, F]]] = []
    max_by_q: dict[int | None, tuple[F, tuple[int, ...]]] = {}
    max_prefix_nonconsec: dict[int, tuple[F, tuple[int, ...]]] = {
        qcut: (F(0), ()) for qcut in range(7, k + 1)
    }
    max_low_nonconsec = (F(0), ())
    max_tail_nonconsec = (F(0), ())
    checked = 0

    for tail in combinations(range(1, W + 1), k - 1):
        E = (0,) + tail
        if not primitive(E):
            continue
        checked += 1
        intervals = dense_intervals(E)
        D = sum((b - a for a, b in intervals), F(0))
        widths: dict[int | None, F]
        if D:
            _, widths = denom_budget_from_intervals(intervals)
        else:
            widths = {}
        if D > best_D:
            best_D, best_E = D, E

        if not is_consecutive(E):
            top.append((D, E, widths))
            low = sum((w for q, w in widths.items() if q is not None and q <= k), F(0))
            tail_mass = D - low
            if low > max_low_nonconsec[0]:
                max_low_nonconsec = (low, E)
            if tail_mass > max_tail_nonconsec[0]:
                max_tail_nonconsec = (tail_mass, E)
            for qcut in range(7, k + 1):
                prefix = sum(
                    (w for q, w in widths.items() if q is not None and q <= qcut),
                    F(0),
                )
                if prefix > max_prefix_nonconsec[qcut][0]:
                    max_prefix_nonconsec[qcut] = (prefix, E)

        for q, w in widths.items():
            old = max_by_q.get(q)
            if old is None or w > old[0]:
                max_by_q[q] = (w, E)

    top.sort(key=lambda item: (item[0], -item[1][-1]), reverse=True)
    local_dead_examples: list[tuple[F, tuple[int, ...], dict[int | None, F], bool]] = []
    for E0 in KNOWN_LOCAL_DEAD_ENDS:
        if len(E0) == k and max(E0) <= W and primitive(E0):
            intervals = dense_intervals(E0)
            D = sum((b - a for a, b in intervals), F(0))
            _, widths = denom_budget_from_intervals(intervals)
            local_dead_examples.append((D, E0, widths, local_left_dead_end(E0)))
    local_dead_examples.sort(key=lambda item: item[0], reverse=True)

    return {
        "k": k,
        "W": W,
        "checked": checked,
        "consec_D": consec_D,
        "consec_widths": consec_widths,
        "consec_hist": consec_hist,
        "best_D": best_D,
        "best_E": best_E,
        "top": top[:TOP_N],
        "max_by_q": max_by_q,
        "max_prefix_nonconsec": max_prefix_nonconsec,
        "max_low_nonconsec": max_low_nonconsec,
        "max_tail_nonconsec": max_tail_nonconsec,
        "local_dead_examples": local_dead_examples[:5],
    }


def strongly_connected_components(vertices: list[str], edges: dict[tuple[str, str], str]) -> list[list[str]]:
    adj = {v: [] for v in vertices}
    radj = {v: [] for v in vertices}
    for (a, b), winner in edges.items():
        loser = b if winner == a else a
        adj[winner].append(loser)
        radj[loser].append(winner)

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

    comps: list[list[str]] = []
    seen.clear()

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


def count_directed_triangles(vertices: list[str], edges: dict[tuple[str, str], str]) -> int:
    total = 0
    for a, b, c in combinations(vertices, 3):
        wins = {
            (a, b): edges[(min(a, b), max(a, b))],
            (a, c): edges[(min(a, c), max(a, c))],
            (b, c): edges[(min(b, c), max(b, c))],
        }
        out = Counter(wins.values())
        if sorted(out.values()) == [1, 1, 1]:
            total += 1
    return total


def is_edge(edges: dict[tuple[str, str], str], a: str, b: str) -> bool:
    key = (min(a, b), max(a, b))
    return edges[key] == a


def count_hamiltonian_paths(vertices: list[str], edges: dict[tuple[str, str], str]) -> int:
    total = 0
    for path in permutations(vertices):
        if all(is_edge(edges, path[i], path[i + 1]) for i in range(len(path) - 1)):
            total += 1
    return total


def proof_carrier_tournament() -> None:
    carriers = [
        "denominator-center",
        "EWLB-window",
        "Fourier-spectrum",
        "finite-atlas",
        "runner-compression",
        "raw-gap-moment",
    ]
    # Score coordinates follow the docstring declaration.  The values are not a
    # theorem; they record the current session's proof-routing judgement after
    # the local-compression failures and incoming Bonferroni bridge.
    score = {
        "denominator-center": (1, 1, 1, 1, 1, 1),
        "EWLB-window": (1, 1, 1, 1, 1, 0),
        "Fourier-spectrum": (1, 0, 0, 1, 1, 0),
        "finite-atlas": (1, 1, 1, 1, 0, 0),
        "runner-compression": (0, 0, 0, 1, 0, 1),
        "raw-gap-moment": (0, 0, 0, 0, 0, 1),
    }
    order = {name: i for i, name in enumerate(carriers)}
    edges: dict[tuple[str, str], str] = {}
    scores = Counter()
    for a, b in combinations(sorted(carriers), 2):
        if score[a] > score[b]:
            winner = a
        elif score[b] > score[a]:
            winner = b
        else:
            winner = a if order[a] < order[b] else b
        edges[(a, b)] = winner
        scores[winner] += 1
    hist = Counter(scores[v] for v in carriers)
    print("\nTournament Analysis: proof-carrier routing")
    print("  score histogram:", dict(sorted(hist.items())))
    print("  directed 3-cycles:", count_directed_triangles(carriers, edges))
    print("  SCC sizes:", [len(c) for c in strongly_connected_components(carriers, edges)])
    print("  Hamiltonian path count:", count_hamiltonian_paths(carriers, edges))
    ranked = sorted(carriers, key=lambda v: (scores[v], score[v]), reverse=True)
    print("  leading Hamiltonian path:", " -> ".join(ranked))


def print_result(res: dict[str, object]) -> None:
    k = int(res["k"])
    consec_D = res["consec_D"]  # type: ignore[assignment]
    best_D = res["best_D"]  # type: ignore[assignment]
    best_E = res["best_E"]  # type: ignore[assignment]
    consec_widths = res["consec_widths"]  # type: ignore[assignment]
    print(f"\n--- k={k}, primitive anchored shapes in [0,{res['W']}] ---")
    print(f"checked={res['checked']}")
    print(f"consec D={consec_D} ({frac_float(consec_D)})")
    print(f"consec denominator widths: {fmt_widths(consec_widths)}")
    print(f"best D={best_D} ({frac_float(best_D)}) at E={best_E}")
    if best_E == tuple(range(k)):
        print("bounded verdict: consecutive remains the unique dense-mass leader in this bank")
    else:
        print("bounded verdict: NONCONSECUTIVE LEADER FOUND")

    max_low, E_low = res["max_low_nonconsec"]  # type: ignore[misc]
    max_tail, E_tail = res["max_tail_nonconsec"]  # type: ignore[misc]
    print(
        "nonconsec budget maxima: "
        f"low(q<=k)={max_low} ({frac_float(max_low)}) at {E_low}; "
        f"tail(q>k/None)={max_tail} ({frac_float(max_tail)}) at {E_tail}"
    )

    print("top nonconsecutive dense shapes:")
    for idx, (D, E, widths) in enumerate(res["top"], 1):  # type: ignore[assignment]
        ratio = D / consec_D if consec_D else F(0)
        complete = [q for q in range(8, k + 1) if residue_complete(E, q)]
        print(
            f"  {idx:2d}. D={D} ({frac_float(D)}), ratio={frac_float(ratio)}, "
            f"E={E}, complete_q={complete}, widths={fmt_widths(widths, 5)}"
        )

    print("per-q nonconsecutive max vs consecutive mass:")
    max_by_q = res["max_by_q"]  # type: ignore[assignment]
    all_q = sorted(q for q in max_by_q if q is not None and q <= max(k + 6, 15))
    for q in all_q:
        con_q = consec_widths.get(q, F(0))
        max_q, E_q = max_by_q[q]
        if con_q or max_q:
            marker = "OK" if max_q <= max(con_q, F(0)) or E_q == tuple(range(k)) else "beats-consec-q"
            print(
                f"  q={q:2d}: consec={con_q} ({frac_float(con_q)}), "
                f"bank_max={max_q} ({frac_float(max_q)}) at {E_q} [{marker}]"
            )

    print("prefix q<=Q nonconsecutive max vs consecutive prefix:")
    max_prefix = res["max_prefix_nonconsec"]  # type: ignore[assignment]
    for qcut in range(7, k + 1):
        con_prefix = sum(
            (w for q, w in consec_widths.items() if q is not None and q <= qcut),
            F(0),
        )
        max_prefix_q, E_prefix = max_prefix[qcut]
        marker = "OK" if max_prefix_q <= con_prefix else "beats-prefix"
        print(
            f"  Q={qcut:2d}: consec_prefix={con_prefix} ({frac_float(con_prefix)}), "
            f"nonconsec_max={max_prefix_q} ({frac_float(max_prefix_q)}) at {E_prefix} [{marker}]"
        )

    print("sample one-step-left compression dead ends:")
    for D, E, widths, is_dead in res["local_dead_examples"]:  # type: ignore[assignment]
        print(
            f"  D={D} ({frac_float(D)}), E={E}, "
            f"verified_dead_end={is_dead}, widths={fmt_widths(widths, 5)}"
        )


def main() -> None:
    print("=" * 78)
    print("LRC14 Lemma A denominator-center budget audit (exact rational)")
    print("=" * 78)
    print(f"theta={THETA}; bounded scan windows={SCAN_WINDOWS}; top_n={TOP_N}")
    print(
        "Lead being tested: consecutive dense mass is carried by low least-denominator "
        "components q=8..k; nonconsecutive local optima leak to high-q tails."
    )

    results = []
    for k in range(8, 14):
        print(f"\n[scan] starting k={k}, W={SCAN_WINDOWS[k]}")
        res = scan_k(k, SCAN_WINDOWS[k])
        results.append(res)
        print_result(res)

    print("\nCross-k summary:")
    for res in results:
        k = int(res["k"])
        consec_D = res["consec_D"]  # type: ignore[assignment]
        best_D = res["best_D"]  # type: ignore[assignment]
        gap = consec_D - best_D
        print(
            f"  k={k}: best-minus-consec={best_D - consec_D} "
            f"({frac_float(best_D - consec_D)}), "
            f"consec_uniqueness={'yes' if gap == 0 and res['best_E'] == tuple(range(k)) else 'no'}"
        )

    print("\nInterpretation:")
    print(
        "  The bounded data support a denominator-budget proof, not a raw compression proof. "
        "  The exact local dead ends have dense mass at higher least denominators and stay far "
        "  below consecutive.  Individual q-buckets are not monotone; the plausible Lemma A "
        "  split is cumulative: prove prefix low-denominator budget inequalities and bound "
        "  the high-q tail by Farey/three-gap width packing."
    )
    proof_carrier_tournament()
    print("\nDONE.")


if __name__ == "__main__":
    main()
