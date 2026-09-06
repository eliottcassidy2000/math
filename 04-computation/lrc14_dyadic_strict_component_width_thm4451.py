#!/usr/bin/env python3
"""Computer-assisted all-height component-width audit for odd q=2 tails.

The infinite reduction is an interval-capacity lemma plus a geometric treatment
of the finitely many small base pairs whose capacity surplus alone permits an
unbounded third tail.  Every remaining triple is checked exactly.
"""

from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import importlib.util


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "lrc14_dyadic_strict_component_topology_thm4451.py"
SPEC = importlib.util.spec_from_file_location("component_topology", SOURCE)
C = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(C)
M = C.M


def need(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def allowed(n: int, ternary_unit: bool) -> bool:
    return n > 0 and n % 2 == 1 and (not ternary_unit or n % 3 != 0)


def comb_capacity(n: int, length: Q) -> Q:
    """Maximum C_n occupancy of an interval, C_n=D_n union tau(D_n)."""
    z = 2 * n * length
    whole = z.numerator // z.denominator
    remainder = z - whole
    return (Q(2 * whole, 7) + min(remainder, Q(2, 7))) / (2 * n)


def surplus(n: int, length: Q) -> Q:
    return comb_capacity(n, length) - Q(2, 7) * length


def pair_essential_geometry(a: int, b: int) -> tuple[Q, Q, int]:
    """Largest a.e.-filled component, least positive circular gap, count."""
    pair_a = M.merge(M.danger(a) + M.danger(b))
    pair_f = M.intersect(pair_a, M.shift_half(pair_a))
    if not pair_f:
        return Q(0), Q(1), 0
    widths = [right - left for left, right in pair_f]
    gaps = [pair_f[i + 1][0] - pair_f[i][1] for i in range(len(pair_f) - 1)]
    gaps.append(1 + pair_f[0][0] - pair_f[-1][1])
    need(all(gap > 0 for gap in gaps), (a, b, gaps))
    return max(widths), min(gaps), len(pair_f)


def first_geometric_safe_c(a: int, b: int, length: Q, ternary_unit: bool) -> int:
    """First allowed c>b after which the P_ab union C_c bound is monotone."""
    width, gap, count = pair_essential_geometry(a, b)
    c = b + 2
    while not allowed(c, ternary_unit):
        c += 2
    while True:
        tooth = Q(1, 7 * c)
        if count == 0:
            safe = tooth < length
        else:
            safe = tooth < gap and width + 2 * tooth < length
        if safe:
            return c
        c += 2
        while not allowed(c, ternary_unit):
            c += 2


def audit(label: str, length: Q, ternary_unit: bool, expected: tuple[int, int, int]) -> None:
    deficit = Q(8, 7) * length
    a_bound = Q(15, 49) / deficit
    print(f"CASE={label} target={length} deficit={deficit} a_bound={a_bound}")

    # The universal surplus bound is s_n <= 5/(49n).  The displayed finite
    # residue replay is a hostile implementation check, not its proof.
    for n in range(1, 20001, 2):
        need(surplus(n, length) <= Q(5, 49 * n), (label, n, "surplus"))
        need(surplus(n, length) >= 0, (label, n, "negative surplus"))

    possible_a = [
        a
        for a in range(1, a_bound.numerator // a_bound.denominator + 1, 2)
        if allowed(a, ternary_unit) and Q(a) <= a_bound
    ]
    base_pairs: list[tuple[int, int]] = []
    for a in possible_a:
        rem_after_a = deficit - surplus(a, length)
        need(rem_after_a > 0, (label, a, "single-tail surplus closes"))
        b_bound = Q(10, 49) / rem_after_a
        for b in range(a + 2, b_bound.numerator // b_bound.denominator + 1, 2):
            if allowed(b, ternary_unit) and Q(b) <= b_bound:
                base_pairs.append((a, b))

    unbounded = [
        (a, b) for a, b in base_pairs
        if surplus(a, length) + surplus(b, length) >= deficit
    ]
    bounded = [pair for pair in base_pairs if pair not in unbounded]
    print(f"possible_a={possible_a} base_pairs={base_pairs}")
    print(f"capacity_unbounded_base_pairs={unbounded}")

    candidates: set[tuple[int, int, int]] = set()
    print("BOUNDED_CAPACITY_FAMILIES")
    for a, b in bounded:
        remaining = deficit - surplus(a, length) - surplus(b, length)
        c_bound = Q(5, 49) / remaining
        count = 0
        for c in range(b + 2, c_bound.numerator // c_bound.denominator + 1, 2):
            if allowed(c, ternary_unit) and Q(c) <= c_bound:
                candidates.add((a, b, c))
                count += 1
        if count:
            print(f"pair={a}:{b} c_le={c_bound} candidates={count}")

    print("GEOMETRIC_UNBOUNDED_FAMILIES")
    for a, b in unbounded:
        width, gap, count = pair_essential_geometry(a, b)
        safe_from = first_geometric_safe_c(a, b, length, ternary_unit)
        small = []
        for c in range(b + 2, safe_from, 2):
            if allowed(c, ternary_unit):
                candidates.add((a, b, c))
                small.append(c)
        print(
            f"pair={a}:{b} pair_width={width} min_gap={gap} pair_components={count} "
            f"safe_for_allowed_c_ge={safe_from} finite_c={small}"
        )

    rows = []
    literal_checks = 0
    for tails in sorted(candidates):
        longest, components, _ = C.strict_stats(tails)
        literal = C.literal_strict_stats(tails)
        need((longest, components) == literal[:2], (label, tails, "literal", literal))
        literal_checks += 1
        rows.append((longest, tails, components))
    rows.sort(reverse=True)
    violations = [row for row in rows if row[0] > length]
    leaders = [tails for longest, tails, _ in rows if longest == length]
    need(not violations, (label, violations[:5]))
    need(leaders == [expected], (label, leaders))
    print(
        f"finite_distinct_triples={len(rows)} literal_checks={literal_checks} maximum={rows[0][0]} "
        f"leaders={leaders} runner_up={rows[1]}"
    )
    print("TOP10")
    for row in rows[:10]:
        print(row)


def main() -> None:
    print("Q2_COMPONENT_WIDTH_ALL_HEIGHT_AUDIT")
    print("STATUS=COMPUTER_ASSISTED_ALL_HEIGHT_PROOF")
    audit("ALL_ODD", Q(17, 693), False, (1, 9, 11))
    audit("ODD_3UNIT_STRICT", Q(19, 1001), True, (1, 11, 13))
    print("PASS")


if __name__ == "__main__":
    main()
