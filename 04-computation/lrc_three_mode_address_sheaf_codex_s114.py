#!/usr/bin/env python3
"""Address-sheaf audit for the three LRC tournament recursion modes.

The same letters A..G are local addresses, not global variables.  Mobius
applies to every size as the full inclusion-exclusion skeleton; Legendre is the
odd half-tiling geometry; Eisenstein is the even half-tiling geometry.  This
script keeps the slot labels before projecting to scalar cell counts, so the
odd C/D size-N-2 cancellation does not erase the geometric distinction.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import comb, gcd, lcm


@dataclass(frozen=True)
class Slot:
    label: str
    sign: int
    offset: int


MODES: dict[str, list[Slot]] = {
    "mobius_full_all": [
        Slot("A", +1, 1),
        Slot("B", +1, 1),
        Slot("C", +1, 1),
        Slot("D", -1, 2),
        Slot("E", -1, 2),
        Slot("F", -1, 2),
        Slot("G", +1, 3),
    ],
    "eisenstein_even_half": [
        Slot("A", +1, 1),
        Slot("B", +1, 1),
        Slot("C", -1, 2),
    ],
    "legendre_odd_half": [
        Slot("A", +1, 1),
        Slot("B", +1, 1),
        Slot("C", -1, 2),
        Slot("D", +1, 2),
        Slot("E", -1, 3),
        Slot("F", -1, 3),
        Slot("G", +1, 4),
    ],
}

LEGENDRE_REGIONS: dict[str, Counter[str]] = {
    "corner_A": Counter({"A": 1}),
    "corner_D": Counter({"D": 1}),
    "corner_B": Counter({"B": 1}),
    "edge_AD": Counter({"A": 1, "D": 1, "E": -1}),
    "edge_BD": Counter({"B": 1, "D": 1, "F": -1}),
    "edge_AB": Counter({"A": 1, "B": 1, "C": -1}),
    "center": Counter({"A": 1, "B": 1, "D": 1, "C": -1, "E": -1, "F": -1, "G": 1}),
}


def full_cells(n: int) -> int:
    return comb(n - 1, 2) if n >= 2 else 0


def half_cells(n: int) -> int:
    return ((n - 1) * (n - 1)) // 4 if n >= 1 else 0


def signed_mode_value(mode: str, n: int, cell_func) -> int:
    return sum(slot.sign * cell_func(n - slot.offset) for slot in MODES[mode])


def scalar_projection(mode: str) -> dict[int, int]:
    projection: Counter[int] = Counter()
    for slot in MODES[mode]:
        projection[slot.offset] += slot.sign
    return {offset: projection[offset] for offset in sorted(projection)}


def slot_table(mode: str, n: int) -> str:
    terms = []
    for slot in MODES[mode]:
        sign = "+" if slot.sign > 0 else "-"
        terms.append(f"{sign}{slot.label}:N-{slot.offset}={n-slot.offset}")
    return " ".join(terms).lstrip("+")


def expr_text(expr: Counter[str]) -> str:
    out = []
    for label in ["A", "B", "C", "D", "E", "F", "G"]:
        coeff = expr[label]
        if coeff == 0:
            continue
        sign = "+" if coeff > 0 else "-"
        atom = label if abs(coeff) == 1 else f"{abs(coeff)}{label}"
        out.append(f"{sign}{atom}")
    return " ".join(out).lstrip("+")


def region_size_projection(expr: Counter[str], n: int) -> dict[int, int]:
    """Map a labelled Legendre expression to signed subtournament sizes."""
    offsets = {slot.label: slot.offset for slot in MODES["legendre_odd_half"]}
    projected: Counter[int] = Counter()
    for label, coeff in expr.items():
        projected[n - offsets[label]] += coeff
    return {size: projected[size] for size in sorted(projected, reverse=True)}


def print_mode_audit() -> None:
    print("Three-mode local slot atlas")
    print("=" * 78)
    print("Scalar projections by size offset:")
    for mode in MODES:
        print(f"  {mode:22s}: {scalar_projection(mode)}")
    print()

    print("Exact cell-count checks on the proper carriers")
    mobius_fail = [n for n in range(5, 40) if signed_mode_value("mobius_full_all", n, full_cells) != full_cells(n)]
    eisen_fail = [
        n
        for n in range(4, 40, 2)
        if signed_mode_value("eisenstein_even_half", n, half_cells) != half_cells(n)
    ]
    legendre_fail = [
        n
        for n in range(5, 40, 2)
        if signed_mode_value("legendre_odd_half", n, half_cells) != half_cells(n)
    ]
    print(f"  Mobius/full all sizes n=5..39 failures: {mobius_fail}")
    print(f"  Eisenstein/even half n even failures:   {eisen_fail}")
    print(f"  Legendre/odd half n odd failures:       {legendre_fail}")
    print()

    print("Parity carrier versus formula coincidence")
    print("  n  h(n)  even A+B-C  odd A+B-C+D-E-F+G  carrier")
    for n in range(4, 16):
        even = signed_mode_value("eisenstein_even_half", n, half_cells)
        odd = signed_mode_value("legendre_odd_half", n, half_cells)
        carrier = "Eisenstein geometry" if n % 2 == 0 else "Legendre geometry"
        print(f"  {n:2d} {half_cells(n):5d} {even:11d} {odd:22d}  {carrier}")
    print("  Note: the seven-term scalar formula also holds at even n, but the")
    print("  three-corner Legendre half-tiling geometry is only present at odd n.")
    print()


def print_legendre_geometry(n: int = 15) -> None:
    print(f"Corrected Legendre odd address geometry at N={n}")
    print("=" * 78)
    print("Slots:")
    print(f"  {slot_table('legendre_odd_half', n)}")
    print("Regions:")
    for name, expr in LEGENDRE_REGIONS.items():
        print(f"  {name:9s}: {expr_text(expr):22s} sizes={region_size_projection(expr, n)}")
    print()
    print("Key point:")
    print("  C and D are both size N-2, so they cancel in scalar cardinality,")
    print("  but they are different addresses: C sits on the AB edge, D is the")
    print("  third corner and feeds both AD and BD edges plus the center.")
    print()


def print_lrc14_composition() -> None:
    print("LRC14 composition: two size coordinates must be kept")
    print("=" * 78)
    n = 14
    print(f"At N={n}, h(N)={half_cells(n)} = 7*6 is pronic.")
    print(f"  Mobius skeleton slots:     {slot_table('mobius_full_all', n)}")
    print(f"  Eisenstein even slots:     {slot_table('eisenstein_even_half', n)}")
    print("  Shape/apex coordinate:     k=N/2=7, the odd Legendre character coordinate.")
    print()
    print("Do not conflate these coordinates:")
    print("  recurrence children of N=14 are sizes 13 and 12;")
    print("  the apex-7 coordinate is the pronic fold parameter k=7;")
    print("  the Legendre chart at apex size 7 has slots")
    print(f"    {slot_table('legendre_odd_half', 7)}")
    print("  while the odd half chart at N=15 has slots")
    print(f"    {slot_table('legendre_odd_half', 15)}")
    print()
    print("Thus composition means a parity-stratified address sheaf:")
    print("  every node carries the Mobius full atlas;")
    print("  even half nodes carry the Eisenstein pronic quotient;")
    print("  odd half nodes carry the Legendre square/three-corner quotient;")
    print("  scalarizing before matching these local addresses loses the C/D seam.")
    print()


def next_prime_after(n: int) -> int:
    def is_prime(x: int) -> bool:
        if x < 2:
            return False
        d = 2
        while d * d <= x:
            if x % d == 0:
                return False
            d += 1 if d == 2 else 2
        return True

    q = n + 1
    while not is_prime(q):
        q += 1
    return q


def safe_units(speeds: list[int], d: int) -> list[int]:
    good = []
    for a in range(1, d):
        if gcd(a, d) != 1:
            continue
        ok = True
        for s in speeds:
            r = (s * a) % d
            if 14 * min(r, d - r) < d:
                ok = False
                break
        if ok:
            good.append(a)
    return good


def first_witness(speeds: list[int], limit: int) -> tuple[int | None, list[int]]:
    for d in range(2, limit + 1):
        good = safe_units(speeds, d)
        if good:
            return d, good
    return None, []


def print_lcm_family_audit() -> None:
    print("Exact-period obstruction for S_X={1..11,13,lcm(2..X)}")
    print("=" * 78)
    print("For every D<=X, lcm(2..X) is 0 mod D, so N(S_X,D)=0 exactly.")
    print("The scan below corrects the stronger nextprime wording: q_min>X is")
    print("rigorous and sufficient to refute any finite basis; q_min need not be")
    print("nextprime(X) because the fixed AP-core has its own safe-denominator ladder.")
    print()
    print("  X   nextprime(X)  first witness D  first safe units")
    for x in [14, 15, 16, 20, 24, 27, 30, 41, 60, 90]:
        tail = 1
        for d in range(2, x + 1):
            tail = lcm(tail, d)
        speeds = list(range(1, 12)) + [13, tail]
        d0, units = first_witness(speeds, max(200, 3 * x))
        print(f"  {x:2d} {next_prime_after(x):14d} {str(d0):16s} {units[:6]}")
    print()


def print_tournament_analysis() -> None:
    print("Tournament Analysis over proof carriers")
    print("=" * 78)
    features = {
        "exact_period_packets": {"LRC_predicate", "CRT", "phi", "unit_residue", "denominator"},
        "address_sheaf": {"slot_label", "parity", "mode", "size", "LRC_predicate"},
        "legendre_CD_seam": {"slot_label", "parity", "chi7", "corner_edge_center"},
        "eisenstein_pronic_fold": {"parity", "mode", "apex_k", "complement_fold"},
        "mobius_scalar_skeleton": {"mode", "size", "inclusion_exclusion"},
        "three_gap_AP_hull": {"AP_extremal", "gap_lengths", "finite_Node2"},
        "raw_runner_vertices": set(),
    }
    names = list(features)
    adj = {name: set() for name in names}
    scores = Counter({name: 0 for name in names})
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if j <= i:
                continue
            key_a = (len(features[a]), len(features[a] & features[b]), -i)
            key_b = (len(features[b]), len(features[a] & features[b]), -j)
            if key_a >= key_b:
                adj[a].add(b)
                scores[a] += 1
            else:
                adj[b].add(a)
                scores[b] += 1
    cycles3 = 0
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if j <= i:
                continue
            for k, c in enumerate(names):
                if k <= j:
                    continue
                if b in adj[a] and c in adj[b] and a in adj[c]:
                    cycles3 += 1
                if c in adj[a] and b in adj[c] and a in adj[b]:
                    cycles3 += 1
    path = sorted(names, key=lambda name: (scores[name], len(features[name]), name), reverse=True)
    print("  vertices are proof carriers, not runners or arcs")
    print("  observable=(labels preserved, common labels, declaration order)")
    print(f"  score_hist={dict(sorted(Counter(scores.values()).items()))}")
    print(f"  directed_3cycles={cycles3}")
    print("  Hamiltonian path=" + " > ".join(path))
    print()
    print("Assumption challenged:")
    print("  A scalar half-tile recurrence is not enough to preserve the LRC")
    print("  predicate.  The quotient must retain exact-period packets, local")
    print("  slot labels, and parity mode until the final cap/floor scalarization.")
    print()


def main() -> None:
    print_mode_audit()
    print_legendre_geometry()
    print_lrc14_composition()
    print_lcm_family_audit()
    print_tournament_analysis()
    print("Proof-route synthesis")
    print("=" * 78)
    print("  Node 2 remains finite/algebraic: prove AP-hull or three-gap")
    print("  majorization for the bounded cap side without losing sector labels.")
    print("  Node 3 is analytic: lcm-loaded rows force unbounded witness")
    print("  denominators, so Weyl/equidistribution or exact-period packet")
    print("  density is irreducible.")
    print("  The three recursion modes compose as an address system connecting")
    print("  these branches, not as a closed scalar recurrence for witness counts.")


if __name__ == "__main__":
    main()
