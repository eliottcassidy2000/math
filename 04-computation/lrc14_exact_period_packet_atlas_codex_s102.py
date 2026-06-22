#!/usr/bin/env python3
"""
Exact-period packet atlas for the LRC(14) witness problem.

This is the LRC-side extension of the HYP-2882/HYP-2883 lesson:
keep primitive packets before scalarizing.

For a denominator D, an exact-period packet is a unit a mod D.  It certifies
the row S at level 1/14 when

    14 * min(sa mod D, D - sa mod D) >= D

for every speed s in S.  This script compares scalar counts N(S,D) with
mod-7, chi_7, affine-pair, and CRT/exact-period quotients.

Tournament Analysis:
  vertices: quotient lenses on exact-period packets, not runners.
  pairwise observable: aggregate weighted variance explained for packet safety.
  switch/gauge: orient lens A -> B if A explains more safety variance.
  tie Hamiltonian path: fewer quotient cells, then the printed lens order.
  fingerprints: score histogram, directed 3-cycles, SCCs, Hamiltonian path count.

Assumption challenge:
  Candidate vertices considered: runners, denominators, numerator residues,
  reduced exact denominators, mod-7 residues, chi_7 classes, affine zero-lane
  pairs, CRT factors, support-six packet classes, and proof obligations.
  The chosen quotient is exact-period unit packets.  It preserves the LRC
  predicate "a/D is a level-1/14 witness" and the Euler phi packet count, while
  destroying raw speed magnitude except through residues mod D.  That is the
  right destruction: HYP-2865 shows fixed denominators can be divisor-loaded,
  so the scaled/first-unblocked packet atlas is the only viable finite-residue
  target.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from math import gcd

N = 14
QR7 = {1, 2, 4}
BASIS_D = (21, 41, 53, 83, 89)
SCAN_D = (
    14,
    15,
    17,
    19,
    21,
    23,
    27,
    29,
    31,
    37,
    41,
    43,
    47,
    53,
    55,
    59,
    61,
    65,
    67,
    71,
    73,
    79,
    83,
    89,
    91,
    97,
    98,
    101,
    105,
    109,
    113,
    117,
    121,
    127,
    131,
)


@dataclass(frozen=True)
class Row:
    name: str
    speeds: tuple[int, ...]
    note: str


def lcm(a: int, b: int) -> int:
    return a // gcd(a, b) * b


def lcm_upto(n: int) -> int:
    out = 1
    for k in range(1, n + 1):
        out = lcm(out, k)
    return out


def base_cover(last: int) -> tuple[int, ...]:
    return tuple(range(1, 12)) + (13, last)


ROWS = [
    Row("AP13_boundary", tuple(range(1, 14)), "non-covering AP boundary"),
    Row("cover_84", base_cover(84), "loosest known small covering row"),
    Row("tower_m6", base_cover(84 * 6), "first HYP-2865 row needing D=53"),
    Row("tower_m53", base_cover(84 * 53), "later tower row needing D=55 in S93"),
    Row("AP12_182", tuple(range(1, 13)) + (182,), "single-far AP repair"),
    Row(
        "floor_star",
        (1, 2, 3, 5, 7, 8, 9, 10, 11, 12, 13, 38, 42),
        "S3 realized floor row from kps-S4",
    ),
    Row(
        "divload_B60",
        base_cover(84 * lcm_upto(60)),
        "fixed denominators <=60 deliberately killed",
    ),
    Row(
        "divload_B90",
        base_cover(84 * lcm_upto(90)),
        "fixed basis 83/89 deliberately killed",
    ),
]

CRT_ROWS = [row for row in ROWS if row.name in {"cover_84", "tower_m6", "tower_m53", "AP12_182", "floor_star"}]


def phi(n: int) -> int:
    out = n
    m = n
    p = 2
    while p * p <= m:
        if m % p == 0:
            out -= out // p
            while m % p == 0:
                m //= p
        p += 1 if p == 2 else 2
    if m > 1:
        out -= out // m
    return out


def units_mod(d: int) -> list[int]:
    return [a for a in range(1, d) if gcd(a, d) == 1]


def safe_residue(r: int, d: int) -> bool:
    r %= d
    return N * min(r, d - r) >= d


def is_witness_unit(speeds: tuple[int, ...], d: int, a: int) -> bool:
    return gcd(a, d) == 1 and all(safe_residue(s * a, d) for s in speeds)


def safe_units(speeds: tuple[int, ...], d: int) -> list[int]:
    return [a for a in units_mod(d) if is_witness_unit(speeds, d, a)]


def chi7(r: int) -> str:
    r %= 7
    if r == 0:
        return "0"
    return "QR" if r in QR7 else "NQR"


def affine_pair(r: int) -> str:
    r %= 7
    s = (2 - r) % 7
    return "/".join(map(str, sorted((r, s))))


def lens_key(name: str, a: int) -> str:
    r7 = a % 7
    if name == "parity":
        return str(a % 2)
    if name == "mod7":
        return str(r7)
    if name == "chi7":
        return chi7(r7)
    if name == "affine_pair":
        return affine_pair(r7)
    if name == "chi_x_affine":
        return f"{chi7(r7)}:{affine_pair(r7)}"
    if name == "mod14":
        return str(a % 14)
    raise KeyError(name)


LENSES = ("parity", "chi7", "affine_pair", "mod7", "chi_x_affine", "mod14")


def variance_explained(units: list[int], labels: dict[int, int], lens: str) -> tuple[Fraction, Fraction]:
    """Return (between-group SS, total SS) for 0/1 safety labels."""
    if not units:
        return Fraction(0), Fraction(0)
    ys = [labels[a] for a in units]
    n = len(ys)
    total_safe = sum(ys)
    # total sum of squares around the global mean, scaled by n.
    total = Fraction(total_safe * (n - total_safe), n)
    if total == 0:
        return Fraction(0), Fraction(0)

    groups: dict[str, list[int]] = defaultdict(list)
    for a in units:
        groups[lens_key(lens, a)].append(labels[a])
    between = Fraction(0)
    for vals in groups.values():
        m = len(vals)
        s = sum(vals)
        between += Fraction(s * s, m)
    between -= Fraction(total_safe * total_safe, n)
    return between, total


def first_unit_witness(speeds: tuple[int, ...], dmax: int = 160) -> tuple[int, int, int] | None:
    for d in range(2, dmax + 1):
        got = safe_units(speeds, d)
        if got:
            return d, got[0], len(got)
    return None


def short_int(n: int) -> str:
    s = str(n)
    if len(s) <= 18:
        return s
    return f"{s[:8]}...{s[-6:]} ({len(s)} digits)"


def packet_distribution(row: Row, d: int) -> str:
    units = units_mod(d)
    labels = {a: int(is_witness_unit(row.speeds, d, a)) for a in units}
    safe = sum(labels.values())
    by_mod7 = Counter(a % 7 for a in units if labels[a])
    by_chi = Counter(chi7(a % 7) for a in units if labels[a])
    by_aff = Counter(affine_pair(a % 7) for a in units if labels[a])
    return (
        f"N={safe}/{len(units)} "
        f"mod7={dict(sorted(by_mod7.items()))} "
        f"chi={dict(sorted(by_chi.items()))} "
        f"affine={dict(sorted(by_aff.items()))}"
    )


def aggregate_lens_scores() -> tuple[dict[str, Fraction], int, int]:
    numer = {lens: Fraction(0) for lens in LENSES}
    denom = Fraction(0)
    mixed_cases = 0
    all_cases = 0
    for row in ROWS:
        for d in SCAN_D:
            units = units_mod(d)
            labels = {a: int(is_witness_unit(row.speeds, d, a)) for a in units}
            safe = sum(labels.values())
            all_cases += 1
            if safe == 0 or safe == len(units):
                continue
            mixed_cases += 1
            _b0, total = variance_explained(units, labels, "mod7")
            denom += total
            for lens in LENSES:
                b, _t = variance_explained(units, labels, lens)
                numer[lens] += b
    scores = {lens: (numer[lens] / denom if denom else Fraction(0)) for lens in LENSES}
    return scores, mixed_cases, all_cases


def tournament_fingerprints(scores: dict[str, Fraction]) -> tuple[dict[int, int], int, list[list[str]], int]:
    order_index = {lens: i for i, lens in enumerate(LENSES)}

    def beats(a: str, b: str) -> bool:
        if scores[a] != scores[b]:
            return scores[a] > scores[b]
        # tie: simpler quotient wins, then declaration order.
        ca = len({lens_key(a, x) for x in range(14)})
        cb = len({lens_key(b, x) for x in range(14)})
        if ca != cb:
            return ca < cb
        return order_index[a] < order_index[b]

    adj = {v: set() for v in LENSES}
    score_hist = Counter()
    for a in LENSES:
        score = 0
        for b in LENSES:
            if a == b:
                continue
            if beats(a, b):
                adj[a].add(b)
                score += 1
        score_hist[score] += 1

    cycles = 0
    lenses = list(LENSES)
    for i, a in enumerate(lenses):
        for j, b in enumerate(lenses):
            for k, c in enumerate(lenses):
                if i < j < k:
                    triples = [(a, b, c), (a, c, b)]
                    for x, y, z in triples:
                        if y in adj[x] and z in adj[y] and x in adj[z]:
                            cycles += 1

    # SCCs by mutual reachability in the tournament orientation.
    reach: dict[str, set[str]] = {}
    for root in LENSES:
        seen = {root}
        q = deque([root])
        while q:
            v = q.popleft()
            for u in adj[v]:
                if u not in seen:
                    seen.add(u)
                    q.append(u)
        reach[root] = seen
    remaining = set(LENSES)
    comps: list[list[str]] = []
    while remaining:
        root = min(remaining, key=order_index.get)
        comp = sorted(
            [v for v in remaining if v in reach[root] and root in reach[v]],
            key=order_index.get,
        )
        comps.append(comp)
        remaining -= set(comp)

    # Count Hamiltonian paths by brute force; n=6.
    def path_count(prefix: tuple[str, ...], rest: tuple[str, ...]) -> int:
        if not rest:
            return 1
        total = 0
        last = prefix[-1]
        for i, nxt in enumerate(rest):
            if nxt in adj[last]:
                total += path_count(prefix + (nxt,), rest[:i] + rest[i + 1 :])
        return total

    hp = 0
    for i, start in enumerate(LENSES):
        rest = LENSES[:i] + LENSES[i + 1 :]
        hp += path_count((start,), rest)
    return dict(sorted(score_hist.items())), cycles, comps, hp


def raw_q_witnesses(row: Row, q: int) -> list[int]:
    return [
        a
        for a in range(q)
        if all(safe_residue(s * a, q) for s in row.speeds)
    ]


def denominator_histogram(row: Row) -> tuple[int, Counter[int]]:
    q = N * max(row.speeds)
    good = raw_q_witnesses(row, q)
    hist: Counter[int] = Counter()
    for a in good:
        if a == 0:
            continue
        hist[q // gcd(a, q)] += 1
    return len(good), hist


def rate(row: Row, d: int) -> Fraction:
    return Fraction(len(safe_units(row.speeds, d)), phi(d))


def crt_defect(row: Row, d1: int, d2: int) -> Fraction | None:
    if gcd(d1, d2) != 1:
        return None
    return rate(row, d1 * d2) - rate(row, d1) * rate(row, d2)


def section(title: str) -> None:
    print("\n" + "=" * 92)
    print(title)
    print("=" * 92)


def print_first_witness_table() -> None:
    section("FIRST EXACT-PERIOD UNIT WITNESSES")
    print(f"{'row':<14} {'max speed':>24} {'first D,a,N':>18} {'basis N(21,41,53,83,89)':>32}")
    for row in ROWS:
        fw = first_unit_witness(row.speeds, 160)
        fw_s = "none<=160" if fw is None else f"{fw[0]},{fw[1]},{fw[2]}"
        basis = tuple(len(safe_units(row.speeds, d)) for d in BASIS_D)
        print(f"{row.name:<14} {short_int(max(row.speeds)):>24} {fw_s:>18} {str(basis):>32}")
    print()
    print("Readout:")
    print("  Fixed finite bases are atlas coordinates, not closures: divload_B90 kills")
    print("  every denominator in {21,41,53,83,89}.  The usable object is the first")
    print("  unblocked exact-period packet, or a scaled CRT packet family.")


def print_selected_packet_tables() -> None:
    section("SELECTED PACKET DISTRIBUTIONS")
    wanted: set[tuple[str, int]] = set()
    for row in ROWS:
        fw = first_unit_witness(row.speeds, 160)
        if fw:
            wanted.add((row.name, fw[0]))
        for d in BASIS_D:
            if safe_units(row.speeds, d):
                wanted.add((row.name, d))
    for row in ROWS:
        for d in sorted(d for name, d in wanted if name == row.name):
            print(f"{row.name:<14} D={d:<3} {packet_distribution(row, d)}")
    print()
    print("The mod-7 residue is often lopsided before scalarization.  The scalar")
    print("count N(S,D) hides whether safety sits in the affine pairs 0/2, 3/6,")
    print("4/5, or in the fixed r=1 lane.")


def print_lens_tournament() -> None:
    section("TOURNAMENT ANALYSIS ON QUOTIENT LENSES")
    scores, mixed, total = aggregate_lens_scores()
    print(f"mixed safety cases used: {mixed}/{total}")
    print(f"{'lens':<14} {'cells mod14':>11} {'variance explained':>20}")
    for lens in sorted(LENSES, key=lambda x: (-scores[x], len({lens_key(x, a) for a in range(14)}), LENSES.index(x))):
        cells = len({lens_key(lens, a) for a in range(14)})
        print(f"{lens:<14} {cells:>11} {float(scores[lens]):>19.6f}")
    score_hist, cycles, comps, hp = tournament_fingerprints(scores)
    print()
    print(f"score_histogram: {score_hist}")
    print(f"directed_3_cycles: {cycles}")
    print(f"SCCs: {comps}")
    print(f"Hamiltonian_path_count: {hp}")
    print()
    print("Reading:")
    print("  mod14 wins because the LRC band has a hard 2*7 boundary.  The important")
    print("  proof fact is not the win itself; it is that chi_7 and affine_pair retain")
    print("  nontrivial signal before the mod14 scalar quotient.  That is the exact")
    print("  layer where the HYP-2632 signed-current kernel lives.")


def print_crt_decomposition() -> None:
    section("SCALED q=14*V EXACT-DENOMINATOR DECOMPOSITION")
    for row in CRT_ROWS:
        q = N * max(row.speeds)
        if q > 65000:
            print(f"{row.name:<14} q={q:<7} skipped for compact runtime")
            continue
        total, hist = denominator_histogram(row)
        top = hist.most_common(8)
        print(f"{row.name:<14} q={q:<7} good_a={total:<5} top exact denominators={top}")
    print()
    print("Scaled CRT witnesses are not bounded-D witnesses.  They decompose into")
    print("many exact periods; the proof obligation is to show at least one primitive")
    print("packet survives after the speed residues act as forbidden classes.")


def print_crt_defects() -> None:
    section("CRT MULTIPLICATIVITY DEFECTS")
    pairs = [(3, 7), (5, 7), (7, 13), (3, 41), (7, 83)]
    for row in [ROWS[0], ROWS[1], ROWS[4], ROWS[5]]:
        print(f"\n{row.name}:")
        for d1, d2 in pairs:
            defect = crt_defect(row, d1, d2)
            if defect is None:
                continue
            print(
                f"  D={d1}*{d2:<2} rate={rate(row,d1*d2)} "
                f"product={rate(row,d1)*rate(row,d2)} defect={defect}"
            )
    print()
    print("Strong-component H is multiplicative after condensation.  Exact-period")
    print("packet capacity phi is multiplicative under CRT, but the LRC safe band is")
    print("archimedean and creates a defect.  That defect is the LRC analogue of the")
    print("strong-ear insertion profile: keep it as data before forming a scalar.")


def main() -> None:
    section("LRC14 EXACT-PERIOD PACKET ATLAS - CODEX S102")
    print_first_witness_table()
    print_selected_packet_tables()
    print_lens_tournament()
    print_crt_decomposition()
    print_crt_defects()
    section("S102 CONCLUSION")
    print("1. HYP-2865's no-fixed-cap theorem is visible at packet level: divisor")
    print("   loading kills any prescribed finite denominator basis.")
    print("2. The productive replacement is a scaled exact-period atlas.  Unit packets")
    print("   keep Euler phi multiplicity, mod-7 residue, chi_7, and affine-pair data")
    print("   until after the LRC witness predicate is evaluated.")
    print("3. The HYP-2883 signed-current graph should be lifted to these packet")
    print("   atlases: first prove local conservation on exact-period residue fibers,")
    print("   then handle CRT multiplicativity defects as strong-component atoms.")


if __name__ == "__main__":
    main()
