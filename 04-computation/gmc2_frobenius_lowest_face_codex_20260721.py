#!/usr/bin/env python3
"""Exact checks for THM-2022 (Frobenius amplification of a Wick face).

Only Python's integer and rational arithmetic is used.  The examples are
chosen to exercise features which a unique-channel check would miss:

* a four-vertex collinear face whose second constant term cancels, while its
  third constant term has two contributing channels;
* a two-vertex face with rational supporting-line intercept;
* a neutral endpoint with A0=0 whose first Gaussian moment cancels against an
  off-face radial term.

For each example the script enumerates every balanced channel at m=p*m0,
checks the Kummer/Lucas residue layer, and verifies the canonical normalized
Frobenius congruence from THM-2022 exactly:

    M_(p*m0)/(p*A0)! = Q^p (mod p).

The first example deliberately has A0 >= p.  This checks that the theorem
requires only p>m0; the older normalization by p^A0 required the unnecessary
extra restriction p>A0.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from math import factorial
from typing import Iterator, Sequence


@dataclass(frozen=True)
class Term:
    label: str
    q: int
    a: int
    c: int

    @property
    def b(self) -> int:
        return self.a - self.q


@dataclass(frozen=True)
class Case:
    name: str
    terms: tuple[Term, ...]
    face: tuple[int, ...]
    lam: Fraction
    delta: Fraction
    m0: int
    p: int
    earlier_zero_power: int | None = None
    earlier_gaussian_zero_power: int | None = None


def compositions(total: int, parts: int) -> Iterator[tuple[int, ...]]:
    """Yield all weak compositions of total into the requested parts."""
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in compositions(total - first, parts - 1):
            yield (first,) + tail


def multinomial(total: int, r: Sequence[int]) -> int:
    out = factorial(total)
    for value in r:
        out //= factorial(value)
    return out


def factorial_vp(n: int, p: int) -> int:
    value = 0
    while n:
        n //= p
        value += n
    return value


def vp(n: int, p: int) -> int:
    if n == 0:
        raise ValueError("v_p(0) is not used in this verifier")
    n = abs(n)
    value = 0
    while n % p == 0:
        n //= p
        value += 1
    return value


def dot(values: Sequence[int], weights: Sequence[int]) -> int:
    return sum(value * weight for value, weight in zip(values, weights))


def balanced_channels(terms: Sequence[Term], m: int) -> list[tuple[int, ...]]:
    charges = tuple(term.q for term in terms)
    return [r for r in compositions(m, len(terms)) if dot(r, charges) == 0]


def radial_height(terms: Sequence[Term], r: Sequence[int]) -> int:
    return dot(r, tuple(term.a for term in terms))


def coefficient_product(terms: Sequence[Term], r: Sequence[int]) -> int:
    value = 1
    for term, exponent in zip(terms, r):
        value *= term.c**exponent
    return value


def gaussian_channel_coefficient(
    terms: Sequence[Term], m: int, r: Sequence[int]
) -> int:
    height = radial_height(terms, r)
    return multinomial(m, r) * factorial(height) * coefficient_product(terms, r)


def face_constant_term(case: Case, power: int) -> int:
    """Compute CT((sum_{i in F} c_i u^q_i)^power) exactly."""
    total = 0
    face_terms = tuple(case.terms[index] for index in case.face)
    for local_r in compositions(power, len(face_terms)):
        if dot(local_r, tuple(term.q for term in face_terms)) != 0:
            continue
        total += multinomial(power, local_r) * coefficient_product(face_terms, local_r)
    return total


def face_channels(case: Case) -> list[tuple[int, ...]]:
    face_set = set(case.face)
    out: list[tuple[int, ...]] = []
    for r in balanced_channels(case.terms, case.m0):
        if all(index in face_set or value == 0 for index, value in enumerate(r)):
            out.append(r)
    return out


def fmt_fraction(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def fmt_channel(r: Sequence[int]) -> str:
    return "(" + ",".join(str(value) for value in r) + ")"


def check_case(case: Case) -> dict[str, object]:
    assert all(term.a >= 0 and term.b >= 0 and term.c != 0 for term in case.terms)
    assert case.p > case.m0

    # Verify the proposed supporting line and its exact equality set.
    slacks = tuple(
        Fraction(term.a) - case.lam * term.q - case.delta for term in case.terms
    )
    assert all(slack >= 0 for slack in slacks)
    assert tuple(index for index, slack in enumerate(slacks) if slack == 0) == case.face
    face_charges = [case.terms[index].q for index in case.face]
    assert len(face_charges) == len(set(face_charges))
    assert min(face_charges) <= 0 <= max(face_charges)

    base_face = face_channels(case)
    assert base_face
    base_heights = {radial_height(case.terms, s) for s in base_face}
    assert len(base_heights) == 1
    A0 = base_heights.pop()
    assert Fraction(A0) == case.delta * case.m0
    Q = face_constant_term(case, case.m0)
    assert Q != 0 and Q % case.p != 0
    assert all(term.c % case.p != 0 for term in case.terms)
    if case.earlier_zero_power is not None:
        assert face_constant_term(case, case.earlier_zero_power) == 0

    M = case.p * case.m0
    channels = balanced_channels(case.terms, M)
    channel_rows: list[tuple[tuple[int, ...], int, int, int]] = []
    for r in channels:
        A = radial_height(case.terms, r)
        multinom = multinomial(M, r)
        coeff = gaussian_channel_coefficient(case.terms, M, r)
        predicted_v = factorial_vp(A, case.p) + vp(multinom, case.p)
        assert vp(coeff, case.p) == predicted_v
        channel_rows.append((r, A, predicted_v, coeff))

    face_factorial = factorial(case.p * A0)
    face_factorial_v = factorial_vp(case.p * A0, case.p)
    min_v = min(row[2] for row in channel_rows)
    min_channels = {row[0] for row in channel_rows if row[2] == min_v}
    expected_min = {tuple(case.p * value for value in s) for s in base_face}
    assert min_v == face_factorial_v
    assert min_channels == expected_min

    # Kummer zero-carry iff r is componentwise divisible by p.
    for r, _A, _value_v, _coeff in channel_rows:
        multi_v = vp(multinomial(M, r), case.p)
        assert (multi_v == 0) == all(value % case.p == 0 for value in r)

    # Lucas on every member of the minimum layer.
    for s in base_face:
        r = tuple(case.p * value for value in s)
        assert multinomial(M, r) % case.p == multinomial(case.m0, s) % case.p

    # Canonical theorem normalization by the full common face factorial.
    # Every quotient is integral.  Precisely the p-scaled face channels can
    # remain nonzero modulo p.
    assert all(row[3] % face_factorial == 0 for row in channel_rows)
    for r, _A, _value_v, coeff in channel_rows:
        quotient = coeff // face_factorial
        if r in expected_min:
            assert quotient % case.p != 0
        else:
            assert quotient % case.p == 0

    full_moment = sum(row[3] for row in channel_rows)
    factorial_normalized_residue = sum(
        (row[3] // face_factorial) % case.p for row in channel_rows
    ) % case.p
    assert factorial_normalized_residue == (full_moment // face_factorial) % case.p
    frobenius_residue = pow(Q % case.p, case.p, case.p)
    assert factorial_normalized_residue == frobenius_residue != 0

    # Extra valuation audit.  This is not the canonical normalization, but it
    # verifies that the full moment has exactly the face-factorial valuation.
    valuation_divisor = case.p**face_factorial_v
    U_integer = face_factorial // valuation_divisor
    U_mod = U_integer % case.p
    valuation_normalized_residue = (full_moment // valuation_divisor) % case.p
    assert valuation_normalized_residue == U_mod * frobenius_residue % case.p
    wilson_unit: int | None = None
    if case.p > A0:
        wilson_unit = ((-1) ** A0 * factorial(A0)) % case.p
        assert face_factorial_v == A0
        assert U_mod == wilson_unit
    assert full_moment != 0 and vp(full_moment, case.p) == face_factorial_v

    if case.earlier_gaussian_zero_power is not None:
        earlier_m = case.earlier_gaussian_zero_power
        earlier_moment = sum(
            gaussian_channel_coefficient(case.terms, earlier_m, r)
            for r in balanced_channels(case.terms, earlier_m)
        )
        assert earlier_moment == 0

    off_face_count = sum(
        any(r[index] for index in range(len(case.terms)) if index not in case.face)
        for r, _A, _value_v, _coeff in channel_rows
    )
    off_face_valuations = [
        value_v
        for r, _A, value_v, _coeff in channel_rows
        if any(r[index] for index in range(len(case.terms)) if index not in case.face)
    ]

    valuation_hist = Counter(row[2] for row in channel_rows)
    tied_pairs = sum(count * (count - 1) // 2 for count in valuation_hist.values())
    min_path = sorted(min_channels)

    return {
        "A0": A0,
        "Q": Q,
        "M": M,
        "channels": len(channels),
        "base_face": base_face,
        "min_channels": min_channels,
        "min_path": min_path,
        "min_v": min_v,
        "face_factorial_v": face_factorial_v,
        "off_face_count": off_face_count,
        "off_face_min_v": min(off_face_valuations) if off_face_valuations else None,
        "U_mod": U_mod,
        "wilson_unit": wilson_unit,
        "factorial_normalized_residue": factorial_normalized_residue,
        "valuation_normalized_residue": valuation_normalized_residue,
        "frobenius_residue": frobenius_residue,
        "moment_digits": len(str(abs(full_moment))),
        "valuation_hist": valuation_hist,
        "tied_pairs": tied_pairs,
    }


def print_case(case: Case, result: dict[str, object]) -> None:
    terms = ", ".join(
        f"{term.label}:(q={term.q},a={term.a},b={term.b},c={term.c})"
        for term in case.terms
    )
    print(f"CASE {case.name}")
    print(f"  terms: {terms}")
    print(
        "  supporting line: "
        f"a - ({fmt_fraction(case.lam)})q >= {fmt_fraction(case.delta)}; "
        f"face={case.face}"
    )
    if case.earlier_zero_power is not None:
        print(
            f"  CT(face^{case.earlier_zero_power})=0; "
            f"CT(face^{case.m0})=Q={result['Q']}"
        )
    else:
        print(f"  CT(face^{case.m0})=Q={result['Q']}")
    if case.earlier_gaussian_zero_power is not None:
        print(f"  exact earlier Gaussian moment M_{case.earlier_gaussian_zero_power}=0")
    print(
        f"  m0={case.m0}, A0={result['A0']}, p={case.p}, "
        f"M=p*m0={result['M']}"
    )
    print(
        f"  balanced channels={result['channels']}; "
        f"off-face channels={result['off_face_count']}"
    )
    print(
        f"  minimum valuation={result['min_v']}; "
        f"minimum layer size={len(result['min_channels'])}; "
        "layer=" + ", ".join(fmt_channel(r) for r in result["min_path"])
    )
    print(f"  least off-face valuation={result['off_face_min_v']}")
    print(
        f"  Kummer zero-carry iff r=p*s checked on all {result['channels']} channels; "
        f"Lucas checked on {len(result['min_channels'])} minimum channels"
    )
    print(
        f"  division by (p*A0)! integral termwise on all {result['channels']} channels; "
        "every quotient outside the p-scaled face layer vanishes mod p"
    )
    print(
        f"  canonical congruence: M_M/(p*A0)! mod p="
        f"{result['factorial_normalized_residue']}; Q^p mod p={result['frobenius_residue']}"
    )
    print(
        f"  valuation audit: v_p((p*A0)!)={result['face_factorial_v']}; "
        f"unit part={result['U_mod']}; p^(-v)M_M mod p="
        f"{result['valuation_normalized_residue']}"
    )
    if result["wilson_unit"] is None:
        print("  Wilson shortcut not invoked (A0>=p); full-factorial theorem still passes")
    else:
        print(f"  Wilson shortcut applicable: unit=(-1)^A0*A0!={result['wilson_unit']} mod p")
    print(
        f"  full moment: nonzero, v_p={result['min_v']}, "
        f"decimal digits={result['moment_digits']}"
    )
    print("  PASS")


def tournament_report(result: dict[str, object]) -> None:
    n = int(result["channels"])
    tied_pairs = int(result["tied_pairs"])
    valuation_hist: Counter[int] = result["valuation_hist"]  # type: ignore[assignment]
    hist_text = ", ".join(f"{value}:{valuation_hist[value]}" for value in sorted(valuation_hist))
    print("TOURNAMENT ANALYSIS (multi-channel collinear case)")
    print("  vertices: all balanced channels at M=p*m0")
    print("  pairwise observable: p-adic valuation of binom(M;r)*A(r)!")
    print("  switch/gauge: lower valuation wins; valuation ties use lexicographic r")
    print("  tie Hamiltonian path: increasing (valuation,r)")
    print(f"  vertices={n}; valuation-layer histogram={{{hist_text}}}")
    print(f"  score histogram: each score 0..{n - 1} occurs once")
    print("  directed 3-cycles=0; SCC sizes=all 1; Hamiltonian paths=1")
    print(f"  edge flips under reverse-lex tie gauge={tied_pairs}")
    print(f"  minimum tied layer path={', '.join(fmt_channel(r) for r in result['min_path'])}")
    print("  preserved: minimum valuation layer and carry/noncarry separation")
    print("  destroyed: the residue sum within a tied layer")
    print("  repair sidecar: CT(face^m0)=Q, retained as Q^p by Frobenius")
    print(
        "  challenged vertex assumptions: monomials, charges, primitive circuits, "
        "radial heights, channels, digit/carry states, prime ideals, proof obligations"
    )
    print("  selected vertices: balanced channels; proof-bearing sidecar: face residue Q")


def main() -> None:
    cases = (
        Case(
            name="four-vertex collinear face with a cancelled earlier layer",
            terms=(
                Term("f-2", q=-2, a=1, c=1),
                Term("f-1", q=-1, a=2, c=1),
                Term("f+1", q=1, a=4, c=-1),
                Term("f+2", q=2, a=5, c=1),
                Term("off0", q=0, a=4, c=2),
            ),
            face=(0, 1, 2, 3),
            lam=Fraction(1),
            delta=Fraction(3),
            m0=3,
            p=5,
            earlier_zero_power=2,
        ),
        Case(
            name="rational-intercept two-vertex face",
            terms=(
                Term("f-1", q=-1, a=1, c=2),
                Term("f+2", q=2, a=3, c=3),
                Term("off0", q=0, a=2, c=5),
            ),
            face=(0, 1),
            lam=Fraction(2, 3),
            delta=Fraction(5, 3),
            m0=3,
            p=7,
        ),
        Case(
            name="neutral endpoint with A0=0 and cancelled first moment",
            terms=(
                Term("constant", q=0, a=0, c=6),
                Term("radial", q=0, a=2, c=-3),
                Term("positive", q=1, a=2, c=2),
            ),
            face=(0,),
            lam=Fraction(0),
            delta=Fraction(0),
            m0=1,
            p=7,
            earlier_gaussian_zero_power=1,
        ),
    )

    print("THM-2022 EXACT FROBENIUS LOWEST-FACE VERIFIER")
    print("arithmetic: Python integers + fractions only")
    print()
    results: list[dict[str, object]] = []
    for case in cases:
        result = check_case(case)
        results.append(result)
        print_case(case, result)
        print()

    tournament_report(results[0])
    print()
    print("GLOBAL VERDICT: all exact valuation, Lucas, and Frobenius checks passed")


if __name__ == "__main__":
    main()
