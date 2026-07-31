#!/usr/bin/env python3
"""Central-arrival clock-return classification for p=2q-1.

Fix q>=2, p=2q-1, the half-open p-digit intervals

    I_m=[m/p,(m+1)/p),                 0<=m<p,

and the half-open nearest q-clock cells

    C_l=[(2l-1)/(2q),(2l+1)/(2q)) mod 1.

Equivalently the clock label is

    c_q(x)=floor(q*{x}+1/2) mod q,

so every boundary belongs to the cell on its right.  Put D(x)={p*x}.
For a same-arrival return z,Dz in I_m, the positive-measure clock-pair
support is completely classified:

* odd m: one diagonal pair;
* the central m=q-1: diagonal only;
* even m!=q-1: one diagonal pair and one adjacent off-diagonal pair.

Indeed the q-clock boundary b_l=(2l+1)/(2q) lies in exactly digit I_(2l).
Odd digits therefore lie in one cell.  If m=2r<q-1, the return has

    mass(r,r)   =(p-m)/(2q p^2),
    mass(r,r+1) =(m+1)/(2q p^2).

If m=2r>q-1, it has

    mass(r+1,r)   =(p-m)/(2q p^2),
    mass(r+1,r+1) =(m+1)/(2q p^2),

with labels modulo q.  When q is odd, the central digit is even and its
unique internal clock boundary is the branch fixed point 1/2; its two
diagonal masses are 1/(2p^2).  When q is even, the central digit is odd and
lies wholly in C_(q/2).  Thus the central-arrival diagonal lemma needs no
parity assumption.  The stronger phrase "unique boundary-straddling digit
whose return is diagonal" requires q odd.

The script also distinguishes varying arrivals z in I_a, Dz in I_b.  In the
LRC metric |a-6|+|b-6|, the nearest 13/7 off-diagonal returns are (5,6) and
(7,6), each of off-diagonal mass 1/338.  Holding the arrival fixed at every
event, the nearest evasions are m=4 and m=8.

Finally it audits three explicitly typed orientations.  The inverse central branch
T(u)=(q-1+u)/p preserves the side of 1/2.  The reflected map A=rho o D,
rho(x)={-x}, appears to swap the raw phase cells 3 and 4 when q=7.  But a
typed A-handoff transports a current owner to the *reflected* following
shallow clock and sends h to p-1-h.  The current stored shallow label is
c_q(rho z), not c_q(z), so the correctly relabelled support is diagonal
again.  Reflection changes finitely many half-open boundary assignments;
all reflection statements here concern positive interval support, for which
those null endpoints are explicitly discarded.

This is a phase/typing lemma.  An off-diagonal phase interval does not build a
positive LRC rail, word, unit, endpoint, or scalar exclusion.
"""

from collections import Counter
from fractions import Fraction


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def floor_fraction(value):
    return value.numerator // value.denominator


def clock_label(x, q):
    """Half-open nearest-cell label floor(q*{x}+1/2) modulo q."""
    x %= 1
    return floor_fraction(q * x + Fraction(1, 2)) % q


def clock_boundaries(q):
    return tuple(Fraction(2 * ell + 1, 2 * q) for ell in range(q))


def forward_cross_arrival_support(q, first_digit, second_digit):
    """Positive-length clock pairs for z in I_a and D(z) in I_b."""
    p = 2 * q - 1
    require(0 <= first_digit < p and 0 <= second_digit < p,
            "arrival digit outside its p-chart")
    # On I_a, D(z)=p*z-a.  The b-preimage is one interval of length 1/p^2.
    left = Fraction(first_digit * p + second_digit, p * p)
    right = Fraction(first_digit * p + second_digit + 1, p * p)
    cuts = {left, right}
    for boundary in clock_boundaries(q):
        if left < boundary < right:
            cuts.add(boundary)
        pulled = (first_digit + boundary) / p
        if left < pulled < right:
            cuts.add(pulled)
    cuts = sorted(cuts)
    support = Counter()
    for start, stop in zip(cuts, cuts[1:]):
        midpoint = (start + stop) / 2
        image = p * midpoint - first_digit
        support[clock_label(midpoint, q), clock_label(image, q)] += (
            stop - start
        )
    require(sum(support.values(), Fraction(0)) == Fraction(1, p * p),
            "cross-arrival support lost its return interval")
    return dict(sorted(support.items()))


def predicted_same_arrival_support(q, m):
    """Closed formula for the theorem in the module docstring."""
    p = 2 * q - 1
    denominator = 2 * q * p * p
    if m % 2 == 1:
        label = (m + 1) // 2
        return {(label, label): Fraction(1, p * p)}
    r = m // 2
    if m == q - 1:
        # This case occurs only for odd q; for even q the central digit is odd.
        return {
            (r, r): Fraction(1, 2 * p * p),
            (r + 1, r + 1): Fraction(1, 2 * p * p),
        }
    if m < q - 1:
        return {
            (r, r): Fraction(p - m, denominator),
            (r, (r + 1) % q): Fraction(m + 1, denominator),
        }
    return {
        ((r + 1) % q, r): Fraction(p - m, denominator),
        ((r + 1) % q, (r + 1) % q): Fraction(m + 1, denominator),
    }


def reflected_central_support(q):
    """Raw and correctly typed support under A(x)={-p*x}."""
    p = 2 * q - 1
    central = q - 1
    # Away from the left endpoint of I_c, A(z)=c+1-p*z.  Pulling I_c back
    # gives an interval with reversed endpoint inclusion.  Endpoint choices
    # have zero length, so the open interior suffices for support masses.
    left = Fraction((central + 1) * (p - 1), p * p)
    right = Fraction(p * (central + 1) - central, p * p)
    require(right - left == Fraction(1, p * p),
            "reflected central return length changed")
    cuts = {left, right}
    for boundary in clock_boundaries(q):
        if left < boundary < right:
            cuts.add(boundary)
        pulled = Fraction(central + 1, p) - boundary / p
        if left < pulled < right:
            cuts.add(pulled)
    cuts = sorted(cuts)
    raw = Counter()
    typed = Counter()
    for start, stop in zip(cuts, cuts[1:]):
        z = (start + stop) / 2
        image = central + 1 - p * z
        raw_pair = (clock_label(z, q), clock_label(image, q))
        # If y_1=A(y_0), then D(y_0)=rho(y_1).  Hence the current stored
        # shallow label is c_q(rho z), not the naive c_q(z).
        typed_pair = (clock_label(-z, q), clock_label(image, q))
        raw[raw_pair] += stop - start
        typed[typed_pair] += stop - start
    return dict(sorted(raw.items())), dict(sorted(typed.items()))


def off_diagonal_mass(support):
    return sum(
        mass for (left, right), mass in support.items() if left != right
    )


def normalized_support(support, unit):
    return tuple(
        (pair, int(mass / unit)) for pair, mass in sorted(support.items())
    )


def main():
    # Finite exact controls for the all-q formulas.  The proof is the boundary
    # calculation in the docstring; this sweep is a hostile implementation
    # check over both parities and every endpoint digit.
    q_values = tuple(range(2, 65))
    same_arrival_cases = 0
    boundary_digit_checks = 0
    central_checks = 0
    reflected_orientation_checks = 0
    nearest_varying_checks = 0
    for q in q_values:
        p = 2 * q - 1
        central = q - 1
        boundary_digits = {
            floor_fraction(p * boundary)
            for boundary in clock_boundaries(q)
        }
        require(boundary_digits == set(range(0, p, 2)),
                "q-clock boundaries are not exactly in the even digits")
        boundary_digit_checks += q

        for m in range(p):
            actual = forward_cross_arrival_support(q, m, m)
            expected = predicted_same_arrival_support(q, m)
            require(actual == expected,
                    f"same-arrival formula failed at q={q}, m={m}")
            has_off_diagonal = off_diagonal_mass(actual) > 0
            require(
                has_off_diagonal == (m % 2 == 0 and m != central),
                "parity classification of an off-diagonal return failed",
            )
            same_arrival_cases += 1

        central_support = forward_cross_arrival_support(
            q, central, central
        )
        require(off_diagonal_mass(central_support) == 0,
                "central arrival acquired an off-diagonal forward return")
        if q % 2 == 1:
            require(central in boundary_digits,
                    "odd-q central digit stopped straddling a boundary")
            require(
                sum(m in boundary_digits
                    and off_diagonal_mass(
                        forward_cross_arrival_support(q, m, m)
                    ) == 0
                    for m in range(p)) == 1,
                "central digit is not the unique diagonal boundary digit",
            )
        else:
            require(central not in boundary_digits,
                    "even-q central digit unexpectedly straddles a boundary")
        central_checks += 1

        # Inverse central branch T(u)=(central+u)/p merely reverses the same
        # diagonal relation.  The reflected handoff is also diagonal after
        # its mandatory clock reflection is included.
        inverse_support = {
            (right, left): mass
            for (left, right), mass in central_support.items()
        }
        require(off_diagonal_mass(inverse_support) == 0,
                "inverse central branch changed the q-clock label")
        raw_reflected, typed_reflected = reflected_central_support(q)
        require(off_diagonal_mass(typed_reflected) == 0,
                "typed reflected handoff acquired an off-diagonal return")
        if q % 2 == 1:
            require(off_diagonal_mass(raw_reflected) == Fraction(1, p * p),
                    "naive odd-q reflection hostile stopped looking off-diagonal")
        else:
            require(off_diagonal_mass(raw_reflected) == 0,
                    "even-q reflected central digit left its single clock")
        reflected_orientation_checks += 1

        # Exact nearest varying-arrival classification in the L1 digit metric.
        neighbours = tuple(
            (a, b)
            for a, b in (
                (central - 1, central),
                (central + 1, central),
                (central, central - 1),
                (central, central + 1),
            )
            if 0 <= a < p and 0 <= b < p
        )
        actual_nearest = {
            pair: off_diagonal_mass(
                forward_cross_arrival_support(q, *pair)
            )
            for pair in neighbours
        }
        actual_nearest = {
            pair: mass for pair, mass in actual_nearest.items() if mass
        }
        if q % 2 == 1:
            expected_nearest = {
                (central - 1, central): Fraction(1, 2 * p * p),
                (central + 1, central): Fraction(1, 2 * p * p),
            }
        else:
            expected_nearest = {
                (central - 1, central): Fraction(1, p * p),
                (central + 1, central): Fraction(1, p * p),
                (central, central - 1): Fraction(q + 1, 2 * q * p * p),
                (central, central + 1): Fraction(q + 1, 2 * q * p * p),
            }
        require(actual_nearest == expected_nearest,
                "nearest varying-arrival classification changed")
        nearest_varying_checks += 1

    require(same_arrival_cases == 4_095,
            "all-q same-arrival control universe changed")
    require(boundary_digit_checks == 2_079,
            "all-q clock-boundary control universe changed")

    # Full p=13,q=7 atlas, normalized in the common unit 1/(2q p^2).
    q = 7
    p = 13
    central = 6
    unit = Fraction(1, 2 * q * p * p)
    same_arrival_atlas = tuple(
        (m, normalized_support(
            forward_cross_arrival_support(q, m, m), unit
        ))
        for m in range(p)
    )
    expected_same_arrival_atlas = (
        (0, (((0, 0), 13), ((0, 1), 1))),
        (1, (((1, 1), 14),)),
        (2, (((1, 1), 11), ((1, 2), 3))),
        (3, (((2, 2), 14),)),
        (4, (((2, 2), 9), ((2, 3), 5))),
        (5, (((3, 3), 14),)),
        (6, (((3, 3), 7), ((4, 4), 7))),
        (7, (((4, 4), 14),)),
        (8, (((5, 4), 5), ((5, 5), 9))),
        (9, (((5, 5), 14),)),
        (10, (((6, 5), 3), ((6, 6), 11))),
        (11, (((6, 6), 14),)),
        (12, (((0, 0), 13), ((0, 6), 1))),
    )
    require(same_arrival_atlas == expected_same_arrival_atlas,
            "13/7 same-arrival atlas changed")
    same_arrival_evasions = tuple(
        m for m in range(p)
        if off_diagonal_mass(
            forward_cross_arrival_support(q, m, m)
        )
    )
    require(same_arrival_evasions == (0, 2, 4, 8, 10, 12),
            "13/7 same-arrival evasion digits changed")
    fixed_arrival_distance = min(
        abs(m - central) for m in same_arrival_evasions
    )
    fixed_arrival_nearest = tuple(
        m for m in same_arrival_evasions
        if abs(m - central) == fixed_arrival_distance
    )
    require(fixed_arrival_distance == 2
            and fixed_arrival_nearest == (4, 8),
            "13/7 nearest fixed-arrival evasions changed")

    # Exhaust all 13^2 varying arrival pairs for the exact nearest metric.
    varying_evasions = []
    for first in range(p):
        for second in range(p):
            mass = off_diagonal_mass(
                forward_cross_arrival_support(q, first, second)
            )
            if mass:
                varying_evasions.append((
                    abs(first - central) + abs(second - central),
                    first, second, mass,
                ))
    varying_distance = min(item[0] for item in varying_evasions)
    varying_nearest = tuple(
        item[1:] for item in varying_evasions
        if item[0] == varying_distance
    )
    require(varying_distance == 1 and varying_nearest == (
        (5, 6, Fraction(1, 338)),
        (7, 6, Fraction(1, 338)),
    ), "13/7 nearest varying-arrival evasions changed")

    raw_reflected, typed_reflected = reflected_central_support(q)
    require(normalized_support(raw_reflected, Fraction(1, 338)) == (
        ((3, 4), 1), ((4, 3), 1),
    ), "13/7 naive reflection hostile changed")
    require(normalized_support(typed_reflected, Fraction(1, 338)) == (
        ((3, 3), 1), ((4, 4), 1),
    ), "13/7 typed reflection repair changed")

    print("Central-arrival q-clock return classification for p=2q-1")
    print("convention=I_m=[m/p,(m+1)/p); C_l=[(2l-1)/(2q),(2l+1)/(2q)) mod1")
    print("clock_label=floor(q*{x}+1/2) mod q; claims concern positive interval support")
    print(
        f"general_exact_controls=q=2..64 same_arrival_cases={same_arrival_cases} "
        f"clock_boundary_digit_checks={boundary_digit_checks}"
    )
    print("general_same_arrival_offdiagonal_iff=m even and m!=q-1")
    print("central_diagonal_all_q=yes; unique_diagonal_boundary_digit_requires=q odd")
    print("general_fixed_arrival_nearest=q odd: distance2; q even: distance1")
    print(
        "general_varying_arrival_nearest="
        "q odd: (q-2,q-1),(q,q-1); "
        "q even: also (q-1,q-2),(q-1,q)"
    )
    print(f"p13_q7_same_arrival_atlas_units_1_over_2366={same_arrival_atlas}")
    print(
        f"p13_q7_same_arrival_evasions={same_arrival_evasions} "
        f"nearest={fixed_arrival_nearest} distance={fixed_arrival_distance}"
    )
    print(
        f"p13_q7_varying_arrival_nearest={varying_nearest} "
        f"distance={varying_distance}"
    )
    print(
        "p13_q7_reflection_raw_vs_typed="
        f"{normalized_support(raw_reflected, Fraction(1,338))} -> "
        f"{normalized_support(typed_reflected, Fraction(1,338))}"
    )
    print("inverse_central_branch=diagonal; typed_reflected_handoff=diagonal")
    print(
        "minimal_actual_LRC_evasions="
        "change arrival carrier (fixed 4/8 or varying 5->6/7->6), "
        "allow a repeat clock, or prove a genuinely different typed transporter"
    )
    print(
        "verdict=PASS: arrival six is the odd-q fixed-point-aligned boundary "
        "digit; forward D, the central inverse branch, and typed rho_o_D "
        "do not bypass its diagonal trap"
    )
    print(
        "scope=phase support only; off-diagonal intervals do not supply a "
        "positive rail, unit, endpoint, scalar exclusion, or LRC(14)"
    )


if __name__ == "__main__":
    main()
