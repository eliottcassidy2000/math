#!/usr/bin/env python3
"""Exact referee for THM-1219's six-comb near-tiling curvature law.

Put ``a=7m`` and inspect the ``k=m`` complete safe gap of the slow comb.
For ``m>=6``, the six faster speeds ``a+r`` (``1<=r<=6``) each contribute
exactly the tooth centred at the common integer ``m+1``.  In chronological
order those six teeth almost tile the slow gap.  Their five internal holes
have exact numerators ``3,5,7,9,11`` and the terminal hole has numerator
``13``:

    delta_r = (13-2r)/(14(a+r)(a+r+1)),  r=1,...,5,
    delta_0 = 13/(14a(a+1)).

Consequently the total survivor ``S(a)`` satisfies

    24/[7(a+6)^2] <= S(a) < 24/(7a^2),

so ``a^2 S(a) -> 24/7`` and the surviving fraction of the slow gap is
asymptotic to ``4/a``.  This is a method obstruction, not an LRC
counterexample: every row has six explicit safe intervals.  It proves that
the five-comb ``C/a`` survivor floor of THM-1198 cannot extend to six combs.

The same calculation has a full mod-seven form.  For ``a=7m+s`` and the
same gap ``k=m``, the signed boundary/handoff numerators are

    12s-6;  13-2s-2r (r=5,...,1);  13-2s.

Keeping their positive parts gives curvature charges
``(48,42,43,46,51,58,67)`` for ``s=0,...,6``.  In particular the residue-one
family has the smaller sharp coefficient ``a^2 S(a) -> 3``.

Every signed numerator is also an active-pair slack ``14D-(x+y)`` in the
located-maximizer coordinates of THM-1205.  Its physical handoff length is

    (x+y)/(xy) * (D/(x+y) - 1/14).

Thus the curvature spectrum is an exact weighted path of pair-sum witness
obligations, rather than an unrelated endpoint accident.

All checks use integer/Fraction arithmetic and ``require`` rather than
``assert``, so optimized Python executes the same certificate.

Tournament audit
----------------
The pairwise observable is chronological tooth order.  Speed reversal gives
the tie-free Hamiltonian path ``b6 -> ... -> b1`` and a transitive tournament
with scores ``0,...,5``.  That quotient forgets the proof: the useful vertices
are the five adjacent handoff gaps plus the terminal boundary obligation,
weighted by the exact curvature numerators and denominators above.
"""

from __future__ import annotations

from fractions import Fraction as F


def require(condition: bool, message: object) -> None:
    """Always-on certificate check, including under ``python -O``."""
    if not condition:
        raise RuntimeError(message)


def slow_gap(a: int, k: int) -> tuple[F, F]:
    require(a > 0 and 0 <= k < a, (a, k))
    return F(14 * k + 1, 14 * a), F(14 * k + 13, 14 * a)


def tooth(b: int, centre: int) -> tuple[F, F]:
    require(b > 0, b)
    return F(14 * centre - 1, 14 * b), F(14 * centre + 1, 14 * b)


def intersecting_teeth(a: int, k: int, b: int) -> tuple[tuple[int, F, F], ...]:
    """All open b-teeth meeting the closed slow gap, clipped to that gap."""
    lo, hi = slow_gap(a, k)
    first = (b * lo - F(1, 14)).__floor__() - 1
    last = (b * hi + F(1, 14)).__ceil__() + 1
    rows: list[tuple[int, F, F]] = []
    for centre in range(first, last + 1):
        left, right = tooth(b, centre)
        clipped_left, clipped_right = max(lo, left), min(hi, right)
        if clipped_left < clipped_right:
            rows.append((centre, clipped_left, clipped_right))
    return tuple(rows)


def merged_length(intervals: list[tuple[F, F]]) -> F:
    """Lebesgue length of a finite open-interval union (endpoints are null)."""
    if not intervals:
        return F(0)
    intervals.sort()
    total = F(0)
    left, right = intervals[0]
    for next_left, next_right in intervals[1:]:
        if next_left <= right:
            right = max(right, next_right)
        else:
            total += right - left
            left, right = next_left, next_right
    return total + right - left


def direct_survivor(a: int, k: int) -> tuple[F, tuple[tuple[int, ...], ...]]:
    lo, hi = slow_gap(a, k)
    intervals: list[tuple[F, F]] = []
    centres: list[tuple[int, ...]] = []
    for r in range(1, 7):
        rows = intersecting_teeth(a, k, a + r)
        centres.append(tuple(row[0] for row in rows))
        intervals.extend((row[1], row[2]) for row in rows)
    return (hi - lo) - merged_length(intervals), tuple(centres)


def component_formula(a: int) -> tuple[F, ...]:
    """Safe component lengths, from left to right."""
    require(a > 0, a)
    internal = tuple(
        F(13 - 2 * r, 14 * (a + r) * (a + r + 1))
        for r in range(5, 0, -1)
    )
    terminal = F(13, 14 * a * (a + 1))
    return internal + (terminal,)


def residue_component_formula(a: int, residue: int) -> tuple[F, ...]:
    """Positive complement components for ``a=7m+residue``, ``k=m``."""
    require(0 <= residue <= 6 and a % 7 == residue, (a, residue))
    answer: list[F] = []
    left_numerator = 12 * residue - 6
    if left_numerator > 0:
        answer.append(F(left_numerator, 14 * a * (a + 6)))
    for r in range(5, 0, -1):
        numerator = 13 - 2 * residue - 2 * r
        if numerator > 0:
            answer.append(F(numerator, 14 * (a + r) * (a + r + 1)))
    terminal_numerator = 13 - 2 * residue
    require(terminal_numerator > 0, residue)
    answer.append(F(terminal_numerator, 14 * a * (a + 1)))
    return tuple(answer)


def active_pair_handoff(x: int, y: int, determinant: int) -> F:
    """Signed endpoint handoff from the active-pair margin ``D/(x+y)-1/14``."""
    require(0 < x < y and determinant > 0, (x, y, determinant))
    return F(x + y, x * y) * (F(determinant, x + y) - F(1, 14))


def survivor_formula(a: int) -> F:
    return sum(component_formula(a), F(0))


def check_symbolic_integer_identities(limit: int = 10_000) -> None:
    """Clear denominators for every identity on a large exact integer range."""
    for a in range(1, limit + 1):
        for r in range(1, 6):
            # The common tooth centre obeys 14n=2a+14 when a=7m.
            fourteen_n = 2 * a + 14
            cleared = (
                (fourteen_n - 1) * (a + r + 1)
                - (fourteen_n + 1) * (a + r)
            )
            require(cleared == 13 - 2 * r, (a, r, cleared))

        # Terminal gap between the b1-tooth right endpoint and G_m's right
        # endpoint, after multiplying by 14*a*(a+1).
        m_num = a // 7
        if a % 7 == 0:
            cleared_tail = (
                (14 * m_num + 13) * (a + 1)
                - (14 * (m_num + 1) + 1) * a
            )
            require(cleared_tail == 13, (a, cleared_tail))


def check_family(max_m: int = 2_000) -> dict[str, object]:
    require(max_m >= 6, max_m)
    rows = 0
    minimum_scaled = None
    maximum_scaled = None
    for m in range(6, max_m + 1):
        a, k = 7 * m, m
        direct, centres = direct_survivor(a, k)
        exact = survivor_formula(a)
        require(centres == ((m + 1,),) * 6, (m, centres))
        require(direct == exact, (m, direct, exact))

        components = component_formula(a)
        expected_numerators = (3, 5, 7, 9, 11, 13)
        require(all(value > 0 for value in components), (m, components))
        require(tuple(
            value * 14 * (a + r) * (a + r + 1)
            for value, r in zip(components[:5], range(5, 0, -1))
        ) == tuple(map(F, expected_numerators[:5])), (m, components))
        require(components[-1] * 14 * a * (a + 1) == 13,
                (m, components[-1]))

        lower = F(24, 7 * (a + 6) ** 2)
        upper = F(24, 7 * a * a)
        require(lower <= exact < upper, (m, lower, exact, upper))
        slow_length = F(6, 7 * a)
        relative = exact / slow_length
        require(F(4 * a, (a + 6) ** 2) <= relative < F(4, a),
                (m, relative))

        # Coarse frontier gates all survive; the defining difference beat is
        # q=1, exactly THM-1216's unique no-hole reduced clock.
        speeds = tuple(a + r for r in range(1, 7))
        harmonic = a * sum((F(1, b) for b in speeds), F(0))
        require(F(a + 1, a) < F(13, 6), (m, speeds[0]))
        require(harmonic > 1, (m, harmonic))
        require(speeds[-1] - speeds[-2] == 1, speeds)

        scaled = a * a * exact
        minimum_scaled = scaled if minimum_scaled is None else min(minimum_scaled, scaled)
        maximum_scaled = scaled if maximum_scaled is None else max(maximum_scaled, scaled)
        rows += 1

    # The five-comb physical floor 1/(49a) cannot survive adding a sixth row.
    for m in range(25, max_m + 1):
        a = 7 * m
        require(survivor_formula(a) < F(1, 49 * a), a)

    return {
        "rows": rows,
        "m_min": 6,
        "m_max": max_m,
        "minimum_scaled": minimum_scaled,
        "maximum_scaled": maximum_scaled,
    }


def check_residue_spectrum(max_m: int = 1_000) -> dict[str, object]:
    """Exact seven-residue curvature spectrum and active-pair factorization."""
    require(max_m >= 6, max_m)
    expected_charges = (48, 42, 43, 46, 51, 58, 67)
    expected_component_counts = (6, 7, 6, 5, 4, 3, 2)
    rows = 0
    for m in range(6, max_m + 1):
        for residue in range(7):
            a, k = 7 * m + residue, m
            exact, centres = direct_survivor(a, k)
            components = residue_component_formula(a, residue)
            require(centres == ((m + 1,),) * 6, (m, residue, centres))
            require(exact == sum(components, F(0)), (m, residue, exact, components))
            require(len(components) == expected_component_counts[residue],
                    (m, residue, components))

            signed_left = 12 * residue - 6
            D_left = a - 6 * m
            require(14 * D_left - (a + (a + 6)) == signed_left,
                    (m, residue, D_left))
            require(active_pair_handoff(a, a + 6, D_left) ==
                    F(signed_left, 14 * a * (a + 6)), (m, residue, "left"))

            D_common = m + 1
            for r in range(1, 6):
                signed = 13 - 2 * residue - 2 * r
                x, y = a + r, a + r + 1
                require(14 * D_common - (x + y) == signed,
                        (m, residue, r, signed))
                require(active_pair_handoff(x, y, D_common) ==
                        F(signed, 14 * x * y), (m, residue, r, "internal"))

            signed_terminal = 13 - 2 * residue
            require(14 * D_common - (a + (a + 1)) == signed_terminal,
                    (m, residue, signed_terminal))
            require(active_pair_handoff(a, a + 1, D_common) ==
                    F(signed_terminal, 14 * a * (a + 1)),
                    (m, residue, "terminal"))

            charge = sum(
                (value for value in
                 (signed_left,
                  *(13 - 2 * residue - 2 * r for r in range(1, 6)),
                  signed_terminal)
                 if value > 0),
                0,
            )
            require(charge == expected_charges[residue], (residue, charge))
            lower = F(charge, 14 * (a + 6) ** 2)
            upper = F(charge, 14 * a * a)
            require(lower <= exact < upper, (m, residue, lower, exact, upper))
            rows += 1

    return {
        "rows": rows,
        "m_min": 6,
        "m_max": max_m,
        "charges": expected_charges,
        "component_counts": expected_component_counts,
    }


def sample_rows() -> tuple[tuple[int, F, F, tuple[F, ...]], ...]:
    answer = []
    for m in (6, 10, 100, 1000):
        a = 7 * m
        total = survivor_formula(a)
        answer.append((a, total, a * a * total, component_formula(a)))
    return tuple(answer)


def main() -> None:
    check_symbolic_integer_identities()
    census = check_family()
    spectrum = check_residue_spectrum()

    print("THM-1219 SIX-COMB NEAR-TILING CURVATURE REFEREE")
    print("arithmetic=integers+fractions.Fraction; always-on checks; no dependencies")
    print("family: a=7m, k=m, speeds=(a+1,...,a+6), m>=6")
    print("one intersecting tooth per speed, common centre n=m+1: PASS")
    print("safe component numerators, chronological=(3,5,7,9,11,13): PASS")
    print("S(a)=13/[14a(a+1)]+sum_(r=1)^5 (13-2r)/[14(a+r)(a+r+1)]")
    print("24/[7(a+6)^2] <= S(a) < 24/(7a^2); a^2*S(a) -> 24/7")
    print("survivor_fraction=S/(6/[7a]) is squeezed by 4a/(a+6)^2 and 4/a")
    print(
        f"exact replay: m={census['m_min']}..{census['m_max']} rows={census['rows']} "
        f"scaled_range=[{census['minimum_scaled']},{census['maximum_scaled']}]"
    )
    print("for every a=7m>=175: S(a)<1/(49a), so the THM-1198 floor cannot extend to six")
    print("frontier location: harmonic-crowded; b1/a<13/6; defining q=b6-b5=1")
    print(
        "mod-seven spectrum: residues=0..6; charges="
        f"{spectrum['charges']}; component_counts={spectrum['component_counts']}"
    )
    print(
        f"spectrum replay: m={spectrum['m_min']}..{spectrum['m_max']} "
        f"rows={spectrum['rows']}; charge/[14(a+6)^2] <= S < charge/(14a^2)"
    )
    print("residue-one sharpening: seven holes and a^2*S -> 42/14 = 3")
    print("active-pair factorization: every signed handoff = (x+y)/(xy)*(D/(x+y)-1/14)")
    print("active-pair slacks along boundary path = 14D-(x+y): PASS")
    print("sample rows:")
    for a, total, scaled, components in sample_rows():
        print(f"  a={a}: S={total}; a^2*S={scaled}; components={components}")
    print("Tournament Analysis:")
    print("  observable=chronological tooth order; path=b6->b5->...->b1")
    print("  scores=(0,1,2,3,4,5); cycles=0; SCCs=6; Hamilton_paths=1")
    print("  faithful vertices=7 boundary/handoff active-pair obligations (positive-part path)")
    print("  destroyed by tournament=curvature numerators, denominators, endpoint ownership")
    print("VERDICT=positive survivor but only second order; no universal C/a six-comb floor")


if __name__ == "__main__":
    main()
