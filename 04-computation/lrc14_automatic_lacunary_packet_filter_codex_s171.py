#!/usr/bin/env python3
"""Automatic/lacunary packet filters for LRC14.

This is a synthesis scout, not a proof.  It tests the user's finite-automaton
prompt against the current LRC14 packet discipline: automatic languages may
tag packet fibers, but they are only proof-safe when exact LRC coordinates
such as denominator, residue, endpoint state, and safe intervals remain
attached.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import gcd

N = 14
THRESHOLD = Fraction(1, N)


def fmt(fr: Fraction) -> str:
    return str(fr.numerator) if fr.denominator == 1 else f"{fr.numerator}/{fr.denominator}"


def bits_lsb(n: int) -> list[int]:
    bits: list[int] = []
    while n:
        bits.append(n & 1)
        n >>= 1
    return bits or [0]


def fibbinary_accept(n: int) -> tuple[bool, str]:
    """DFA read least-significant bit first: forbid consecutive 1s."""
    last_one = False
    for bit in bits_lsb(n):
        if bit and last_one:
            return False, "dead:11"
        last_one = bool(bit)
    return True, "ok:last1" if last_one else "ok:last0"


def moser_accept(n: int) -> tuple[bool, str]:
    """DFA read least-significant bit first: 1s only in even bit positions."""
    for pos, bit in enumerate(bits_lsb(n)):
        if bit and pos % 2 == 1:
            return False, f"dead:oddpos{pos}"
    return True, "ok:even-bit-support"


def first_fibbinary(count: int) -> list[int]:
    out: list[int] = []
    x = 1
    while len(out) < count:
        if fibbinary_accept(x)[0]:
            out.append(x)
        x += 1
    return out


def first_moser(count: int) -> list[int]:
    out: list[int] = []
    x = 1
    while len(out) < count:
        if moser_accept(x)[0]:
            out.append(x)
        x += 1
    return out


def union_intervals(intervals: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    if not intervals:
        return []
    intervals = sorted(intervals)
    out: list[tuple[Fraction, Fraction]] = []
    cur_l, cur_r = intervals[0]
    for left, right in intervals[1:]:
        if left <= cur_r:
            cur_r = max(cur_r, right)
        else:
            out.append((cur_l, cur_r))
            cur_l, cur_r = left, right
    out.append((cur_l, cur_r))
    return out


def danger_intervals(speeds: list[int]) -> list[tuple[Fraction, Fraction]]:
    intervals: list[tuple[Fraction, Fraction]] = []
    for speed in speeds:
        den = N * speed
        for k in range(speed):
            left = Fraction(k * N - 1, den)
            right = Fraction(k * N + 1, den)
            if left < 0:
                intervals.append((Fraction(0), right))
                intervals.append((Fraction(1) + left, Fraction(1)))
            elif right > 1:
                intervals.append((left, Fraction(1)))
                intervals.append((Fraction(0), right - 1))
            else:
                intervals.append((left, right))
    return union_intervals(intervals)


def safe_components(speeds: list[int]) -> list[tuple[Fraction, Fraction]]:
    """Positive open components left after removing closed threshold danger arcs."""
    danger = danger_intervals(speeds)
    components: list[tuple[Fraction, Fraction]] = []
    cursor = Fraction(0)
    for left, right in danger:
        if cursor < left:
            components.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < 1:
        components.append((cursor, Fraction(1)))
    return components


def component_measure(components: list[tuple[Fraction, Fraction]]) -> Fraction:
    return sum((right - left for left, right in components), Fraction(0))


def dist_to_integer(x: Fraction) -> Fraction:
    y = x % 1
    return min(y, Fraction(1) - y)


def min_distance(speeds: list[int], t: Fraction) -> Fraction:
    return min(dist_to_integer(speed * t) for speed in speeds)


def boundary_witnesses(speeds: list[int], denominator: int = N) -> list[tuple[int, Fraction]]:
    witnesses: list[tuple[int, Fraction]] = []
    for a in range(1, denominator):
        t = Fraction(a, denominator)
        if min_distance(speeds, t) >= THRESHOLD:
            witnesses.append((a, min_distance(speeds, t)))
    return witnesses


def residue_word(speeds: list[int], modulus: int = N) -> str:
    return ",".join(str(speed % modulus) for speed in speeds)


def max_adjacent_ratio(speeds: list[int]) -> Fraction:
    ordered = sorted(speeds)
    return max(Fraction(b, a) for a, b in zip(ordered, ordered[1:]))


def row_summary(label: str, speeds: list[int]) -> dict[str, object]:
    components = safe_components(speeds)
    largest = max(components, key=lambda x: x[1] - x[0]) if components else None
    midpoint = (largest[0] + largest[1]) / 2 if largest else None
    return {
        "label": label,
        "speeds": speeds,
        "gcd": gcd(*speeds),
        "residues_mod14": residue_word(speeds),
        "zero_mod14": sum(1 for speed in speeds if speed % N == 0),
        "fibbinary_count": sum(1 for speed in speeds if fibbinary_accept(speed)[0]),
        "moser_count": sum(1 for speed in speeds if moser_accept(speed)[0]),
        "positive_components": len(components),
        "safe_measure": component_measure(components),
        "largest_component": largest,
        "largest_midpoint": midpoint,
        "midpoint_min_distance": min_distance(speeds, midpoint) if midpoint else None,
        "boundary_units": boundary_witnesses(speeds),
        "max_adjacent_ratio": max_adjacent_ratio(speeds),
    }


@dataclass(frozen=True)
class Carrier:
    name: str
    role: str
    lrc_use: str
    forgets: str
    retention: tuple[int, ...]


FEATURES = (
    "predicate",
    "finite_state",
    "exact_period",
    "gap_block",
    "exception_ledger",
    "certificate",
    "anti_scalar",
)


CARRIERS = [
    Carrier(
        "labelled_packet_dfa_filter",
        "packet -> binary word -> DFA state -> labelled residual section",
        "Adds finite-state tags without dropping exact M/q/endpoint labels.",
        "Unsafe only if DFA state is used as the final quotient.",
        (3, 3, 3, 2, 2, 3, 3),
    ),
    Carrier(
        "hurwitz_2adic_window_automaton",
        "continued-fraction window and parity state for multiply-by-2",
        "Models the 2-adic/even-fold branch as a local state machine.",
        "Forgets real gap geometry unless paired with endpoint certificates.",
        (2, 3, 2, 1, 2, 2, 3),
    ),
    Carrier(
        "fibbinary_no11_carry_filter",
        "regular language with no adjacent 1 bits",
        "Tags Zeckendorf/no-adjacent carry debts in packet fields.",
        "Forgets exact period and safe intervals if used by itself.",
        (2, 3, 1, 1, 1, 1, 3),
    ),
    Carrier(
        "moser_even_bit_support_filter",
        "regular sublanguage with 1 bits only in even positions",
        "Tests sparse base-4 support and square-like lacunary block rows.",
        "Forgets neighboring odd-position defects and real witnesses.",
        (2, 3, 1, 2, 1, 1, 3),
    ),
    Carrier(
        "hadamard_gap_block_certificate",
        "block gaps / lacunary tail ratio",
        "Routes true large-gap tails to scale-separated interval peelers.",
        "Forgets intra-block collisions unless DFA/residue fields stay attached.",
        (2, 1, 1, 3, 1, 2, 2),
    ),
    Carrier(
        "fermat_catalan_exception_ledger",
        "finite exceptional packet list after exponent/arity constraints",
        "Names the desired proof shape: infinite families discharged, finite packets audited.",
        "Analogy only; has no LRC predicate without packet certificates.",
        (1, 1, 1, 1, 3, 1, 2),
    ),
    Carrier(
        "raw_growth_or_sequence_name",
        "row -> sequence label or growth rate",
        "Negative control for automatic/lacunary scalarization.",
        "Forgets the proof-bearing packet coordinates.",
        (1, 0, 0, 1, 0, 0, 0),
    ),
]


def orient(a: Carrier, b: Carrier, tie_order: list[str]) -> tuple[str, str]:
    wins_a = sum(x > y for x, y in zip(a.retention, b.retention))
    wins_b = sum(y > x for x, y in zip(a.retention, b.retention))
    if wins_a > wins_b:
        return a.name, b.name
    if wins_b > wins_a:
        return b.name, a.name
    return (a.name, b.name) if tie_order.index(a.name) < tie_order.index(b.name) else (b.name, a.name)


def tournament_edges(carriers: list[Carrier]) -> set[tuple[str, str]]:
    order = [carrier.name for carrier in carriers]
    edges: set[tuple[str, str]] = set()
    for a, b in combinations(carriers, 2):
        edges.add(orient(a, b, order))
    return edges


def score_hist(edges: set[tuple[str, str]], names: list[str]) -> dict[int, int]:
    scores = {name: 0 for name in names}
    for winner, _loser in edges:
        scores[winner] += 1
    hist: dict[int, int] = {}
    for score in scores.values():
        hist[score] = hist.get(score, 0) + 1
    return dict(sorted(hist.items()))


def directed_3cycles(edges: set[tuple[str, str]], names: list[str]) -> int:
    edge = set(edges)
    total = 0
    for a, b, c in combinations(names, 3):
        if ((a, b) in edge and (b, c) in edge and (c, a) in edge) or (
            (a, c) in edge and (c, b) in edge and (b, a) in edge
        ):
            total += 1
    return total


def scc_sizes(edges: set[tuple[str, str]], names: list[str]) -> list[int]:
    graph = {name: [] for name in names}
    rev = {name: [] for name in names}
    for a, b in edges:
        graph[a].append(b)
        rev[b].append(a)

    seen: set[str] = set()
    order: list[str] = []

    def dfs(v: str) -> None:
        seen.add(v)
        for w in graph[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for name in names:
        if name not in seen:
            dfs(name)

    seen.clear()
    sizes: list[int] = []

    def rdfs(v: str) -> int:
        seen.add(v)
        size = 1
        for w in rev[v]:
            if w not in seen:
                size += rdfs(w)
        return size

    for name in reversed(order):
        if name not in seen:
            sizes.append(rdfs(name))
    return sizes


def hamiltonian_path_count(edges: set[tuple[str, str]], names: list[str]) -> int:
    idx = {name: i for i, name in enumerate(names)}
    edge_idx = {(idx[a], idx[b]) for a, b in edges}
    n = len(names)
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp[mask][last]
            if not count:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if (last, nxt) in edge_idx:
                    dp[mask | (1 << nxt)][nxt] += count
    return sum(dp[-1])


def main() -> None:
    rows = [
        ("AP", list(range(1, 14))),
        ("GW_12_to_24", list(range(1, 12)) + [13, 24]),
        ("fibbinary_first13", first_fibbinary(13)),
        ("moser_de_bruijn_first13", first_moser(13)),
    ]

    print("LRC14 AUTOMATIC / LACUNARY PACKET FILTER")
    print("status: exact finite scout and proof-interface guardrail, not a proof")
    print()
    print("Source checks imported into this scout:")
    print("  - arXiv:2506.04110: multiplication-by-2 via finite continued-fraction windows.")
    print("  - Ostrowski-Hadamard: true Hadamard gaps force analytic natural-boundary behavior.")
    print("  - Fibbinary/Moser-de Bruijn: regular/2-automatic bit languages give finite-state packet tags.")
    print("  - Fermat-Catalan: use finite-exception ledger as proof shape, not as an LRC theorem.")
    print()

    print("Exact threshold-open safe components at threshold 1/14")
    for row in map(lambda pair: row_summary(*pair), rows):
        print(f"  {row['label']}:")
        print(f"    speeds={row['speeds']}")
        print(f"    residues mod14={row['residues_mod14']}  zero_mod14={row['zero_mod14']}  gcd={row['gcd']}")
        print(
            f"    fibbinary_count={row['fibbinary_count']}/13  "
            f"moser_count={row['moser_count']}/13  "
            f"max_adjacent_ratio={fmt(row['max_adjacent_ratio'])}"
        )
        print(
            f"    positive_safe_components={row['positive_components']}  "
            f"safe_measure={fmt(row['safe_measure'])}"
        )
        if row["largest_component"] is not None:
            left, right = row["largest_component"]  # type: ignore[misc]
            print(
                f"    largest_component=({fmt(left)}, {fmt(right)})  "
                f"midpoint={fmt(row['largest_midpoint'])}  "
                f"min_distance_at_midpoint={fmt(row['midpoint_min_distance'])}"
            )
        else:
            units = ",".join(str(a) for a, _ in row["boundary_units"])  # type: ignore[index]
            print(f"    no positive open component; boundary denominator-14 witness units={units}")
    print()

    print("Finite automaton readout on the two automatic rows")
    for label, speeds in rows[2:]:
        print(f"  {label}:")
        for speed in speeds:
            fib_ok, fib_state = fibbinary_accept(speed)
            mos_ok, mos_state = moser_accept(speed)
            binary = format(speed, "b")
            print(
                f"    {speed:>3} bin={binary:<8} "
                f"fib={str(fib_ok):<5} {fib_state:<10} "
                f"moser={str(mos_ok):<5} {mos_state}"
            )
    print()

    names = [carrier.name for carrier in CARRIERS]
    edges = tournament_edges(CARRIERS)
    top_path = sorted(names, key=lambda name: -sum(1 for a, _b in edges if a == name))

    print("Tournament Analysis")
    print("  vertices: proof carriers / automatic packet filters, not runners")
    print("  pairwise observable: coordinate-retention vector over")
    print(f"    {', '.join(FEATURES)}")
    print("  switch/gauge: majority of strictly larger retention coordinates")
    print("  tie Hamiltonian path: order listed in CARRIERS")
    print(f"  score_hist={score_hist(edges, names)}")
    print(f"  directed_3cycles={directed_3cycles(edges, names)}")
    print(f"  scc_sizes={scc_sizes(edges, names)}")
    print(f"  hamiltonian_path_count={hamiltonian_path_count(edges, names)}")
    print("  top path:")
    for i, name in enumerate(top_path, 1):
        print(f"    {i}. {name}")
    print()

    print("Proof-program takeaway:")
    print("  Automatic languages are good LRC packet filters only before scalarization.")
    print("  In this finite scout, the first 13 fibbinary and Moser-de Bruijn rows")
    print("  are not boundary threats: they have positive open safe mass at 1/14.")
    print("  AP and GW stay as zero-open denominator-14 boundary atoms.  The useful")
    print("  next theorem target is a DFA-labelled packet field for HYP-2963: automaton")
    print("  state, residue mod14, exact-period label, gap-block label, and first")
    print("  safe-component certificate or named F7/THM-572 residual.")


if __name__ == "__main__":
    main()
