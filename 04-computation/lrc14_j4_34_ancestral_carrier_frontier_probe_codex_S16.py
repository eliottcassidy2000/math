#!/usr/bin/env python3
r"""Bounded exact frontier probe for the pure `(3,4)` j=4 flood tail.

Fix

    E = {3,4,8,9,10,11,12,13,14}

and order the four external speeds as `15 <= a < b < c < d`.  This probe does
not evaluate any final `d`-sweep.  It proves a reusable ancestor-carrier
inequality, applies it globally at the root and first exact layer, and audits
only the bounded window `15 <= a <= 40` through the `c` frontier.

If A is the exact current survivor, C is any ancestor good-set carrier with
`A subset C`, and `s` ordered speeds remain, THM-732 gives

    |A minus union_i D_{w_i}|
      >= |A| - s|C|/7 - S2*r(C)/7 * sum_i 1/w_i.          (1)

For a whole tail whose first possible speed is v, the reciprocal sum is
largest at v,v+1,...,v+s-1.  One may also apply P2 to the first needle and
use C only for the remaining s-1 needles:

    >= 6|A|/7 - 8r(A)/(49v)
       - (s-1)|C|/7 - S2*r(C)/7 * sum_{i=1}^{s-1}1/(v+i). (2)

Taking the maximum of (1),(2) over every ancestor is a sound sufficient
certificate.  It is a one-sided proof quotient, not an equality quotient.

The bounded window is intentionally small enough to rerun routinely.  At an
`a` or `b` node it first computes only the exact sparse measure; it constructs
the more fragmented child interval carrier only if no ancestor closes the
whole remaining tail.  At surviving `(a,b)` nodes it counts the exact standard
and carrier-envelope `c` horizons, but never subtracts D_c and never enters d.
"""

from __future__ import annotations

import hashlib
import importlib.util
from fractions import Fraction as F
from itertools import permutations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CORE_PATH = ROOT / "04-computation/lrc14_thm741_2002_body_j4_tree_kps_S128c5.py"
CORE_SHA256 = "5aa81d9d78273c8f9e3e7a6574091a3bc3f64ab6086c7024c15f9420c99dac96"
E = tuple(sorted((*range(8, 15), 3, 4)))
A_MIN = 15
A_WINDOW_MAX = 40
SEARCH_GUARD = 10_000_000


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_core():
    require(sha256(CORE_PATH) == CORE_SHA256, "THM-741 core hash changed")
    spec = importlib.util.spec_from_file_location("thm741_34_carrier_dependency", CORE_PATH)
    require(spec is not None and spec.loader is not None, "cannot load THM-741 core")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def harmonic_packet(start: int, count: int) -> F:
    require(start >= 1 and count >= 0, "bad harmonic packet")
    return sum((F(1, start + index) for index in range(count)), F(0))


def certificate_rows(
    core,
    current_m: F,
    current_r: int,
    carriers: tuple[tuple[int, F, int], ...],
    start: int,
    remaining: int,
    include_p2: bool = True,
) -> tuple[tuple[F, int, str], ...]:
    """All ancestor union and P2-hybrid lower bounds at one tail start."""
    require(current_m > 0 and current_r > 0, "empty current carrier")
    require(remaining >= 1, "certificate requires a remaining needle")
    packet = harmonic_packet(start, remaining)
    rows: list[tuple[F, int, str]] = []
    for depth, carrier_m, carrier_r in carriers:
        rows.append(
            (
                current_m
                - F(remaining, 7) * carrier_m
                - core.S2 * carrier_r * packet / 7,
                depth,
                "union",
            )
        )
        if include_p2:
            rows.append(
                (
                    F(6, 7) * current_m
                    - F(8 * current_r, 49 * start)
                    - F(remaining - 1, 7) * carrier_m
                    - core.S2
                    * carrier_r
                    * harmonic_packet(start + 1, remaining - 1)
                    / 7,
                    depth,
                    "P2+union",
                )
            )
    return tuple(rows)


def best_certificate(*args, **kwargs) -> tuple[F, int, str]:
    return max(certificate_rows(*args, **kwargs), key=lambda row: (row[0], -row[1], row[2]))


def first_positive(core, current_m, current_r, carriers, start, remaining, include_p2=True):
    """First integer tail start with a strict, monotone positive certificate."""

    def value(v: int) -> tuple[F, int, str]:
        return best_certificate(
            core,
            current_m,
            current_r,
            carriers,
            v,
            remaining,
            include_p2=include_p2,
        )

    first = value(start)
    if first[0] > 0:
        return start, first
    low = start
    high = max(start + 1, 2 * start)
    while value(high)[0] <= 0:
        low, high = high, 2 * high
        require(high <= SEARCH_GUARD, "certificate search guard exceeded")
    while low + 1 < high:
        middle = (low + high) // 2
        if value(middle)[0] > 0:
            high = middle
        else:
            low = middle
    winner = value(high)
    require(winner[0] > 0, "nonpositive cutoff")
    require(value(high - 1)[0] <= 0, "cutoff is not minimal")
    return high, winner


def exact_union_closes(core, current_m, carriers, start, remaining) -> bool:
    """Post-sparse screen: (1) only, because the child component count is unknown."""
    rows = certificate_rows(
        core,
        current_m,
        1,  # unused when include_p2 is false
        carriers,
        start,
        remaining,
        include_p2=False,
    )
    return max(row[0] for row in rows) > 0


def tournament_fingerprint(pair_wins, residual_sums):
    vertices = (0, 1, 2)
    adjacency = {vertex: set() for vertex in vertices}
    observables = {}
    for left in vertices:
        for right in vertices:
            if left >= right:
                continue
            left_wins, right_wins = pair_wins[(left, right)]
            observable = (
                left_wins - right_wins,
                residual_sums[right] - residual_sums[left],
                right - left,
            )
            observables[(left, right)] = observable
            if observable > (0, 0, 0):
                adjacency[left].add(right)
            else:
                adjacency[right].add(left)
    scores = {vertex: len(adjacency[vertex]) for vertex in vertices}
    cycle = int(
        (1 in adjacency[0] and 2 in adjacency[1] and 0 in adjacency[2])
        or (2 in adjacency[0] and 1 in adjacency[2] and 0 in adjacency[1])
    )
    reach = {vertex: {vertex} | set(adjacency[vertex]) for vertex in vertices}
    for middle in vertices:
        for left in vertices:
            if middle in reach[left]:
                reach[left] |= reach[middle]
    scc_sizes = sorted(
        {
            tuple(v for v in vertices if v in reach[u] and u in reach[v])
            for u in vertices
        },
        key=lambda block: block[0],
    )
    paths = tuple(
        order
        for order in permutations(vertices)
        if all(order[index + 1] in adjacency[order[index]] for index in range(2))
    )
    return observables, scores, cycle, tuple(map(len, scc_sizes)), paths


def main() -> None:
    core = load_core()
    require(core.S2 == F(99, 70) and core.S2 * core.S2 > 2, "bad sqrt(2) majorant")
    good0, r0, m0 = core.good_norm(E)
    require(r0 == 28 and m0 == F(433607, 2522520), "root geometry changed")
    carriers0 = ((0, m0, r0),)
    V1 = core.minV(4, *(3 * m0 / (core.S2 * r0)).as_integer_ratio())
    ordered_union_cut, ordered_union_winner = first_positive(
        core, m0, r0, carriers0, A_MIN, 4, include_p2=False
    )
    root_cut, root_winner = first_positive(core, m0, r0, carriers0, A_MIN, 4)

    # The whole exact first layer is cheap and gives the honest pure-tail E2
    # baseline.  Deeper work below is restricted to A_MIN..A_WINDOW_MAX.
    a_rows = {}
    baseline_b_global = 0
    for a in range(A_MIN, V1):
        r1, m1, good1 = core.subtract(good0, a)
        require(m1 > 0, f"empty first layer a={a}")
        V2 = core.minV(3, *(4 * m1 / (core.S2 * r1)).as_integer_ratio())
        baseline_b_global += max(0, V2 - a - 1)
        a_rows[a] = (r1, m1, good1, V2)

    a_sparse = a_full_needed = a_post_sparse_closed = 0
    for a in range(A_MIN, root_cut):
        r1, exact_m1, _, _ = a_rows[a]
        sparse_m1 = core.subtract_sparse(good0, a)
        require(sparse_m1 == exact_m1, f"first-layer sparse mismatch a={a}")
        a_sparse += 1
        if exact_union_closes(core, sparse_m1, carriers0, a + 1, 3):
            a_post_sparse_closed += 1
        else:
            a_full_needed += 1
    require(
        a_sparse + (V1 - root_cut) == V1 - A_MIN,
        "root tail partition changed",
    )
    require(a_full_needed + a_post_sparse_closed == a_sparse, "first-layer ledger mismatch")

    window_a = window_baseline_b = window_b_sparse = 0
    window_b_algebra_closed = window_b_post_sparse_closed = window_full_g2 = 0
    window_standard_c = window_envelope_c = 0
    cutoff_winners = {(depth, kind): 0 for depth in range(3) for kind in ("union", "P2+union")}
    pair_wins = {(0, 1): [0, 0], (0, 2): [0, 0], (1, 2): [0, 0]}
    residual_sums = {0: 0, 1: 0, 2: 0}

    for a in range(A_MIN, A_WINDOW_MAX + 1):
        require(a < root_cut, "window extends into root-closed tail")
        r1, m1, good1, V2 = a_rows[a]
        require(
            not exact_union_closes(core, m1, carriers0, a + 1, 3),
            f"window a={a} unexpectedly closes after sparse measure",
        )
        window_a += 1
        baseline_here = max(0, V2 - a - 1)
        window_baseline_b += baseline_here
        carriers1 = carriers0 + ((1, m1, r1),)
        b_cut, _ = first_positive(core, m1, r1, carriers1, a + 1, 3)
        b_stop = min(V2, b_cut)
        window_b_algebra_closed += baseline_here - max(0, b_stop - a - 1)
        for b in range(a + 1, b_stop):
            m2_sparse = core.subtract_sparse(good1, b)
            require(m2_sparse > 0, f"empty sparse second layer a={a},b={b}")
            window_b_sparse += 1
            if exact_union_closes(core, m2_sparse, carriers1, b + 1, 2):
                window_b_post_sparse_closed += 1
                continue
            r2, m2, good2 = core.subtract(good1, b)
            require(m2 == m2_sparse, f"second-layer sparse mismatch a={a},b={b}")
            require(len(good2) == r2, f"second-layer component mismatch a={a},b={b}")
            window_full_g2 += 1
            V3 = core.minV(2, *(5 * m2 / (core.S2 * r2)).as_integer_ratio())
            standard_count = max(0, V3 - b - 1)
            window_standard_c += standard_count
            carriers2 = carriers1 + ((2, m2, r2),)
            c_cut, winner = first_positive(core, m2, r2, carriers2, b + 1, 2)
            c_count = max(0, min(V3, c_cut) - b - 1)
            window_envelope_c += c_count
            cutoff_winners[(winner[1], winner[2])] += 1

            individual_counts = {}
            for depth, carrier_m, carrier_r in carriers2:
                one = ((depth, carrier_m, carrier_r),)
                cut, _ = first_positive(core, m2, r2, one, b + 1, 2)
                individual_counts[depth] = max(0, min(V3, cut) - b - 1)
                residual_sums[depth] += individual_counts[depth]
            for left, right in pair_wins:
                if individual_counts[left] < individual_counts[right]:
                    pair_wins[(left, right)][0] += 1
                elif individual_counts[right] < individual_counts[left]:
                    pair_wins[(left, right)][1] += 1

    require(window_b_sparse + window_b_algebra_closed == window_baseline_b, "b ledger mismatch")
    require(
        window_full_g2 + window_b_post_sparse_closed == window_b_sparse,
        "post-sparse b ledger mismatch",
    )
    require(window_envelope_c <= window_standard_c, "carrier envelope enlarged c frontier")
    observables, scores, cycles, scc_sizes, paths = tournament_fingerprint(pair_wins, residual_sums)

    print("THM-741 PURE (3,4) TAIL: ANCESTRAL-CARRIER BOUNDED FRONTIER PROBE")
    print("=" * 96)
    print(f"dependency_sha256={CORE_SHA256}")
    print(f"body E={E}; external order 15<=a<b<c<d; Q(E)=empty")
    print(f"root r={r0} m={m0}")
    print("proved carrier law: |A\\U D_w| >= |A|-s|C|/7-(S2*r(C)/7)sum(1/w), A subset C")
    print("proved P2 hybrid: first needle uses 6|A|/7-8r(A)/(49v); later needles use ancestor C")
    print(
        f"root cutoffs: common-threshold V1={V1}; ordered-union={ordered_union_cut} "
        f"({ordered_union_winner[2]} C{ordered_union_winner[1]}); "
        f"ordered-P2-envelope={root_cut} ({root_winner[2]} C{root_winner[1]})"
    )
    print(
        f"global first layer: baseline a-nodes={V1-A_MIN}; baseline E2 nodes={baseline_b_global}; "
        f"root-algebra-closed a={V1-root_cut}; exact-m1 sparse={a_sparse}; "
        f"post-sparse-closed a={a_post_sparse_closed}; modeled full-G1 needed={a_full_needed}"
    )
    print(f"bounded deeper window: {A_MIN}<=a<={A_WINDOW_MAX} ({window_a} exact a branches)")
    print(
        f"  b frontier: standard nodes={window_baseline_b}; algebra-closed={window_b_algebra_closed}; "
        f"exact-m2 sparse={window_b_sparse}; post-sparse-closed={window_b_post_sparse_closed}; "
        f"full-G2 needed={window_full_g2}"
    )
    print(
        f"  c frontier at surviving G2: standard nodes={window_standard_c}; "
        f"ancestor-envelope nodes={window_envelope_c}; "
        f"algebra-closed={window_standard_c-window_envelope_c}"
    )
    print("  c-cutoff winning certificate counts:")
    for key in sorted(cutoff_winners):
        print(f"    C{key[0]} {key[1]}: {cutoff_winners[key]}")
    print("Tournament Analysis vertices: C0=root carrier,C1=after a,C2=after a,b (not runners)")
    print("pair observable: nodewise lower-c-horizon wins; tie by smaller residual sum, then older carrier")
    for edge in sorted(observables):
        print(
            f"  C{edge[0]} vs C{edge[1]}: wins={tuple(pair_wins[edge])}; "
            f"residual_sums=({residual_sums[edge[0]]},{residual_sums[edge[1]]}); "
            f"observable={observables[edge]}"
        )
    score_hist = {score: tuple(scores.values()).count(score) for score in sorted(set(scores.values()))}
    print(
        f"fingerprint: score_hist={score_hist}; directed_3cycles={cycles}; "
        f"SCC_sizes={scc_sizes}; Hamiltonian_paths={len(paths)}"
    )
    print("tie Hamiltonian paths=" + ",".join("->".join(f"C{x}" for x in path) for path in paths))
    print("kept: one-sided LRC closure horizons and ancestor containment")
    print("destroyed: exact later-needle overlaps, final margins, and all c/d interval geometry")
    print("challenged vertices: runners/root edges/Fano flags are unnecessary; carrier levels preserve the certificate")
    print("alternate vertices audited: gaps/sections/endpoints/wall-events retain geometry but not this scalar tail horizon")
    print(
        "VERDICT: exact root/first-layer reductions plus an exact a=15..40 frontier census; "
        "no c measures or d sweeps were evaluated, so the (3,4) body and sporadic branch remain open"
    )
    print(f"source_sha256={sha256(Path(__file__).resolve())}")
    print("ALL BOUNDED FRONTIER CHECKS PASSED")


if __name__ == "__main__":
    main()
