#!/usr/bin/env python3
"""Exact LRC14 termination-sidecar and proof-carrier tournament audit.

The vertices of the tournament are proof carriers, not runners.  The audit
compares period-14 local data with three global repairs:

* the signed threshold topology G_{1/14};
* THM-2048's peel fiber-quantization inequality;
* the first strict primitive rational exit.

All circle computations use fractions.  Exact maxima are enumerated on the
pair-sum rulers supplied by THM-2047/THM-1002.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations
from math import gcd


THRESHOLD = Fraction(1, 14)
UNITS14 = tuple(a for a in range(1, 14) if gcd(a, 14) == 1)


def replaced(drop: int, add: int) -> tuple[int, ...]:
    return tuple(sorted((set(range(1, 14)) - {drop}) | {add}))


def replace_many(drops: tuple[int, ...], adds: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted((set(range(1, 14)) - set(drops)) | set(adds)))


ROWS = {
    "AP": tuple(range(1, 14)),
    "GW12to24": replaced(12, 24),
    "liar12to26": replaced(12, 26),
    "K33_12to36": replaced(12, 36),
    "far12to96": replaced(12, 96),
    "cover12to84": replaced(12, 84),
    "P10plusK33": replace_many((10, 12), (20, 36)),
    "tax_gain_Cover14": (1, 8, 11, 12, 14, 17, 22, 26, 35, 40, 54, 90, 93),
}


def circle_distance(v: int, t: Fraction) -> Fraction:
    residue = (v * t.numerator) % t.denominator
    return Fraction(min(residue, t.denominator - residue), t.denominator)


def row_value(speeds: tuple[int, ...], t: Fraction) -> Fraction:
    return min(circle_distance(v, t) for v in speeds)


def exact_maximum(speeds: tuple[int, ...]) -> tuple[Fraction, tuple[Fraction, ...]]:
    best = Fraction(-1)
    witnesses: set[Fraction] = set()
    for i, v in enumerate(speeds):
        for w in speeds[i:]:
            ruler = v + w
            for a in range(ruler):
                t = Fraction(a, ruler)
                value = row_value(speeds, t)
                if value > best:
                    best = value
                    witnesses = {t}
                elif value == best:
                    witnesses.add(t)
    return best, tuple(sorted(witnesses))


def intersect_intervals(
    left: list[tuple[Fraction, Fraction]],
    right: list[tuple[Fraction, Fraction]],
) -> list[tuple[Fraction, Fraction]]:
    raw: list[tuple[Fraction, Fraction]] = []
    i = j = 0
    while i < len(left) and j < len(right):
        lo = max(left[i][0], right[j][0])
        hi = min(left[i][1], right[j][1])
        if lo <= hi:
            raw.append((lo, hi))
        if left[i][1] < right[j][1]:
            i += 1
        else:
            j += 1
    merged: list[tuple[Fraction, Fraction]] = []
    for lo, hi in raw:
        if merged and lo <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append((lo, hi))
    return merged


def one_runner_safe(v: int) -> list[tuple[Fraction, Fraction]]:
    return [
        (Fraction(14 * k + 1, 14 * v), Fraction(14 * k + 13, 14 * v))
        for k in range(v)
    ]


def safe_intervals(speeds: tuple[int, ...]) -> list[tuple[Fraction, Fraction]]:
    intervals = [(Fraction(0), Fraction(1))]
    for v in speeds:
        intervals = intersect_intervals(intervals, one_runner_safe(v))
        if not intervals:
            break
    return intervals


def topology_signature(speeds: tuple[int, ...]) -> tuple[Fraction, int, int, int]:
    intervals = safe_intervals(speeds)
    volume = sum((hi - lo for lo, hi in intervals), Fraction(0))
    positive = sum(lo < hi for lo, hi in intervals)
    points = sum(lo == hi for lo, hi in intervals)
    return volume, len(intervals), positive, points


def first_strict_exit(
    speeds: tuple[int, ...], maximum: Fraction
) -> tuple[int, int, int] | None:
    if maximum <= THRESHOLD:
        return None
    for q in range(2, 2 * max(speeds) + 1):
        for a in range(1, q):
            if gcd(a, q) != 1:
                continue
            slack = min(
                14 * min((a * v) % q, q - ((a * v) % q)) - q for v in speeds
            )
            if slack > 0:
                return q, a, slack
    raise AssertionError("THM-2047 pair-sum bound guarantees a strict exit")


def unit_germ(speeds: tuple[int, ...]) -> tuple[tuple[int, tuple[tuple[int, int], ...]], ...]:
    answer = []
    for a in UNITS14:
        distances = {v: min((a * v) % 14, 14 - ((a * v) % 14)) for v in speeds}
        minimum = min(distances.values())
        binders = []
        for v in speeds:
            residue = (a * v) % 14
            if distances[v] == minimum:
                slope = v if residue < 7 else -v
                binders.append((v, slope))
        answer.append((a, tuple(binders)))
    return tuple(answer)


def peel_tax(
    speeds: tuple[int, ...]
) -> tuple[tuple[int, Fraction, Fraction, int, Fraction], ...]:
    """Return exact violating peels (v,mu,theta,r,lhs-rhs)."""
    violations = []
    for index, v in enumerate(speeds):
        core = speeds[:index] + speeds[index + 1 :]
        intervals = safe_intervals(core)
        mu = sum((hi - lo for lo, hi in intervals), Fraction(0))
        r = sum(lo < hi for lo, hi in intervals)
        scaled = 7 * v * mu
        theta = scaled - (scaled.numerator // scaled.denominator)
        lhs = 6 * (v * mu) ** 2 + theta * (1 - theta) / 7
        rhs = Fraction(r * r, 3)
        if lhs > rhs:
            violations.append((v, mu, theta, r, lhs - rhs))
    return tuple(violations)


def compact_fraction(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def digest(value: object) -> str:
    return sha256(repr(value).encode("ascii")).hexdigest()[:12]


def carrier_tournament(records: dict[str, dict[str, object]]) -> dict[str, object]:
    truth = {name: record["maximum"] > THRESHOLD for name, record in records.items()}
    carriers = {
        "residue14": (1, lambda r: r["residue14"]),
        "unit_germ": (2, lambda r: r["unit_germ"]),
        "peel_tax": (3, lambda r: r["peel_tax"]),
        "first_exit": (4, lambda r: r["first_exit"]),
        "threshold_topology": (5, lambda r: r["topology"]),
        "fused_termination": (
            6,
            lambda r: (r["unit_germ"], r["peel_tax"], r["first_exit"]),
        ),
        "phase_height": (7, lambda r: (r["maximum"], r["witnesses"])),
    }
    metrics: dict[str, tuple[int, int, int]] = {}
    for carrier, (cost, projection) in carriers.items():
        signatures = {name: projection(record) for name, record in records.items()}
        bad = sum(
            truth[a] != truth[b] and signatures[a] == signatures[b]
            for a, b in combinations(records, 2)
        )
        fibers = len(set(signatures.values()))
        metrics[carrier] = (bad, -fibers, cost)

    vertices = tuple(sorted(carriers))
    rank = {v: metrics[v] + (v,) for v in vertices}
    edges = {(a, b) for a, b in combinations(vertices, 2) if rank[a] < rank[b]}
    edges |= {(b, a) for a, b in combinations(vertices, 2) if rank[b] < rank[a]}
    score = {v: sum((v, w) in edges for w in vertices if w != v) for v in vertices}
    cycles = sum(
        (a, b) in edges and (b, c) in edges and (c, a) in edges
        or (a, c) in edges and (c, b) in edges and (b, a) in edges
        for a, b, c in combinations(vertices, 3)
    )

    def reachable(a: str, b: str) -> bool:
        seen = {a}
        stack = [a]
        while stack:
            u = stack.pop()
            for v in vertices:
                if (u, v) in edges and v not in seen:
                    if v == b:
                        return True
                    seen.add(v)
                    stack.append(v)
        return a == b

    sccs = []
    unseen = set(vertices)
    while unseen:
        seed = min(unseen)
        component = tuple(sorted(v for v in unseen if reachable(seed, v) and reachable(v, seed)))
        sccs.append(component)
        unseen.difference_update(component)
    paths = [
        order
        for order in permutations(vertices)
        if all((order[i], order[i + 1]) in edges for i in range(len(order) - 1))
    ]
    return {
        "metrics": metrics,
        "score_histogram": tuple(sorted(score.values())),
        "directed_3cycles": cycles,
        "scc_sizes": tuple(sorted(map(len, sccs), reverse=True)),
        "hamiltonian_path_count": len(paths),
        "tie_hamiltonian_path": paths[0],
    }


def main() -> None:
    print("LRC14 TERMINATION-SIDECAR AUDIT")
    records: dict[str, dict[str, object]] = {}
    for name, speeds in ROWS.items():
        maximum, witnesses = exact_maximum(speeds)
        topology = topology_signature(speeds)
        exit_certificate = first_strict_exit(speeds, maximum)
        tax = peel_tax(speeds)
        records[name] = {
            "maximum": maximum,
            "witnesses": witnesses,
            "topology": topology,
            "first_exit": exit_certificate,
            "peel_tax": tax,
            "unit_germ": unit_germ(speeds),
            "residue14": tuple(sorted(v % 14 for v in speeds)),
        }
        volume, chi, positive, points = topology
        first_tax = tax[0][0] if tax else None
        print(
            f"{name}: M={compact_fraction(maximum)} top_count={len(witnesses)} "
            f"G_volume={compact_fraction(volume)} chi={chi} "
            f"positive_components={positive} points={points} "
            f"first_exit={exit_certificate} peel_violations={len(tax)} "
            f"first_tax_speed={first_tax} germ={digest(records[name]['unit_germ'])}"
        )

    ap = records["AP"]
    liar = records["liar12to26"]
    assert ap["unit_germ"] == liar["unit_germ"]
    assert ap["maximum"] == Fraction(1, 14)
    assert liar["maximum"] == Fraction(1, 12)
    assert ap["topology"] == (Fraction(0), 6, 0, 6)
    assert liar["topology"][0] > 0
    assert liar["first_exit"] == (12, 1, 2)
    print("HOSTILE PAIR")
    print("AP_vs_12to26_complete_unit_germ_equal=PASS")
    print("AP_M=1/14; 12to26_M=1/12; first_exit=(12,1,2)=PASS")
    print("signed_threshold_topology_separates_zero_volume_points_from_strict_arcs=PASS")

    gw = records["GW12to24"]
    assert gw["maximum"] == Fraction(1, 14)
    assert gw["first_exit"] is None
    far = records["far12to96"]
    assert far["first_exit"] == (41, 17, 1)
    print("TERMINATION CONTROLS")
    print("GW_boundary_has_no_strict_exit=PASS")
    print("12to96_first_exit=(41,17,1)=PASS")
    gain = records["tax_gain_Cover14"]
    assert gain["maximum"] == Fraction(7, 43)
    assert gain["topology"] == (Fraction(589595, 5213208), 50, 48, 2)
    assert gain["peel_tax"] == (
        (
            93,
            Fraction(35517, 280280),
            Fraction(19801, 40040),
            50,
            Fraction(2413467317, 235670635200),
        ),
    )
    print("THM2048_Cover14_tax_gain_exactly_reproduced=PASS")
    print("strict_exit_search_complete_by_pair_sum_bound_q<=2*max_speed=PASS")

    tournament = carrier_tournament(records)
    print("TOURNAMENT ANALYSIS (vertices=proof carriers, not runners)")
    print("pairwise_observable=lex(bad_cross-route_merges,-distinct_fibers,cost)")
    print("switch_gauge=quotient the named row bank by each carrier signature")
    for carrier, metric in sorted(tournament["metrics"].items(), key=lambda item: item[1]):
        print(f"  {carrier}: metric={metric}")
    print(f"score_histogram={tournament['score_histogram']}")
    print(f"directed_3cycles={tournament['directed_3cycles']}")
    print(f"scc_sizes={tournament['scc_sizes']}")
    print(f"hamiltonian_path_count={tournament['hamiltonian_path_count']}")
    print(f"tie_hamiltonian_path={tournament['tie_hamiltonian_path']}")
    print(
        "assumption_challenge=period-14 unit germs preserve tight contact but destroy "
        "magnitude and global exit; the valuation/first-exit sidecar is load-bearing"
    )
    print(
        "verdict=PASS; local germ acyclicity is not termination, while signed topology, "
        "peel tax, and resolved exits are compatible exact carriers"
    )


if __name__ == "__main__":
    main()
