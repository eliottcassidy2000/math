#!/usr/bin/env python3
"""Exact local-density audit for Sun's 2-4-6-8 counterexample.

For B_k(t)=binomial(t,k), the sequence modulo p^a has period
P=p^(a+floor(log_p(k))).  We enumerate one such period for each summand,
cyclically convolve the four exact histograms, and report

    sigma_(p,a)(n) = p^a * Prob(F(x)=n mod p^a).

The lower bounds on the original integer variables do not affect these
natural/p-adic densities.  This file uses standard-library exact integers and
Fractions only.  It does not verify global non-representability.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from itertools import product
from math import comb, gcd


N = 896_315_812_331_399
DEGREES = (2, 4, 6, 8)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def floor_log_p(value: int, prime: int) -> int:
    exponent = 0
    while value >= prime:
        value //= prime
        exponent += 1
    return exponent


def period(prime: int, exponent: int, degree: int) -> int:
    return prime ** (exponent + floor_log_p(degree, prime))


def histogram(prime: int, exponent: int, degree: int) -> tuple[int, ...]:
    modulus = prime**exponent
    answer = [0] * modulus
    for value in range(period(prime, exponent, degree)):
        answer[comb(value, degree) % modulus] += 1
    return tuple(answer)


def cyclic_convolution(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    modulus = len(left)
    require(len(right) == modulus, "cyclic convolution modulus drift")
    result = [0] * modulus
    right_nonzero = tuple((j, value) for j, value in enumerate(right) if value)
    for i, a in enumerate(left):
        if not a:
            continue
        for j, b in right_nonzero:
            result[(i + j) % modulus] += a * b
    return tuple(result)


def full_counts(prime: int, exponent: int) -> tuple[tuple[int, ...], int]:
    hists = tuple(histogram(prime, exponent, degree) for degree in DEGREES)
    counts = hists[0]
    for hist in hists[1:]:
        counts = cyclic_convolution(counts, hist)
    universe = 1
    for degree in DEGREES:
        universe *= period(prime, exponent, degree)
    require(sum(counts) == universe, "representation histogram mass drift")
    return counts, universe


def factor(count: int, universe: int, modulus: int) -> Fraction:
    return Fraction(modulus * count, universe)


def state_slopes(prime: int, exponent: int, degree: int) -> tuple[tuple[int, int], ...]:
    """Return (value mod p^a, first lift slope mod p) over the level-a states."""

    modulus = prime**exponent
    step = period(prime, exponent, degree)
    states = []
    for x in range(step):
        base = comb(x, degree)
        delta = comb(x + step, degree) - base
        require(delta % modulus == 0, "claimed period failed")
        slope = (delta // modulus) % prime
        # The p digit lifts must have the claimed affine first-order law.
        for digit in range(prime):
            lhs = comb(x + digit * step, degree) % (prime * modulus)
            rhs = (base + digit * modulus * slope) % (prime * modulus)
            require(lhs == rhs, "one-level affine lift law failed")
        states.append((base % modulus, slope))
    return tuple(states)


def critical_tuple_count(prime: int, exponent: int, target: int) -> tuple[int, tuple[int, ...]]:
    """Count target states for which every one-coordinate lift slope vanishes."""

    modulus = prime**exponent
    critical_hists = []
    critical_sizes = []
    for degree in DEGREES:
        hist = [0] * modulus
        size = 0
        for value, slope in state_slopes(prime, exponent, degree):
            if slope == 0:
                hist[value] += 1
                size += 1
        critical_hists.append(tuple(hist))
        critical_sizes.append(size)
    counts = critical_hists[0]
    for hist in critical_hists[1:]:
        counts = cyclic_convolution(counts, hist)
    return counts[target % modulus], tuple(critical_sizes)


def regular_prime_base(prime: int, target: int) -> tuple[int, tuple[int, ...]]:
    """For p>8 this is the ordinary simultaneous-gradient singular count."""

    require(prime > 8, "base derivative certificate requires p>8")
    return critical_tuple_count(prime, 1, target)


def fmt_fraction(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def audit_row(prime: int, exponent: int) -> dict[str, object]:
    modulus = prime**exponent
    counts, universe = full_counts(prime, exponent)
    factors = tuple(factor(count, universe, modulus) for count in counts)
    target = N % modulus
    minimum = min(factors)
    maximum = max(factors)
    return {
        "p": prime,
        "a": exponent,
        "q": modulus,
        "target": target,
        "image": sum(count > 0 for count in counts),
        "target_count": counts[target],
        "universe": universe,
        "target_factor": factors[target],
        "minimum": minimum,
        "minimum_residues": tuple(i for i, value in enumerate(factors) if value == minimum),
        "maximum": maximum,
    }


def main() -> None:
    # Exact powers chosen to expose stabilization and cover every prime used by
    # the transcript's refined CRT family.
    plan = {
        2: range(1, 9),
        3: range(1, 7),
        5: range(1, 5),
        7: range(1, 4),
        11: range(1, 4),
        13: range(1, 3),
        17: range(1, 3),
        19: range(1, 3),
        23: range(1, 3),
        29: range(1, 2),
        31: range(1, 3),
        43: range(1, 2),
    }
    rows = []
    for prime, exponents in plan.items():
        for exponent in exponents:
            rows.append(audit_row(prime, exponent))

    # Exact period formula is independently tested on a hostile bank of starts.
    period_gates = 0
    for prime in (2, 3, 5, 7, 11, 13, 17, 19, 23):
        for exponent in range(1, 4):
            modulus = prime**exponent
            for degree in DEGREES:
                shift = period(prime, exponent, degree)
                for x in range(0, 4 * degree + 17):
                    require(
                        (comb(x + shift, degree) - comb(x, degree)) % modulus == 0,
                        "period hostile failed",
                    )
                    period_gates += 1

    # Hensel regularity certificates.  At p>8, zero singular solutions modulo
    # p rigorously implies all prescribed lifts have p^3 lifts per solution,
    # hence the target local factor is the base-p factor at every power.
    regular = {}
    for prime in (11, 13, 17, 19, 23, 29, 31, 43):
        singular, sizes = regular_prime_base(prime, N % prime)
        regular[prime] = (singular, sizes)

    # For small primes the degree denominators change the state periods.  The
    # same finite lift-slope test at the displayed level certifies one lift and
    # records whether a critical target state exists.  We do not silently turn
    # this finite certificate into an all-level theorem.
    small_regular = {}
    for prime, exponent in ((2, 1), (3, 2), (5, 2), (7, 2)):
        small_regular[(prime, exponent)] = critical_tuple_count(
            prime, exponent, N % (prime**exponent)
        )

    by_key = {(row["p"], row["a"]): row for row in rows}
    require(by_key[(3, 1)]["target_factor"] == Fraction(22, 27), "p=3 base factor")
    require(by_key[(3, 2)]["target_factor"] == Fraction(68, 81), "p=3 target factor")
    require(by_key[(11, 1)]["target_factor"] == Fraction(72, 121), "p=11 target factor")
    require(by_key[(5, 2)]["target_factor"] == Fraction(566, 625), "p=5^2 target factor")
    require(all(row["image"] == row["q"] for row in rows), "local image hole found")
    # The target fibre is smooth at every audited p>8 except p=31.  The four
    # singular p=31 states are deliberately retained as a hostile rather than
    # hidden behind a blanket stabilization claim.
    require(
        all(regular[p][0] == 0 for p in regular if p != 31),
        "unexpected large-prime singular target fibre",
    )
    require(regular[31][0] == 4, "p=31 singular hostile drift")

    # CRT family actually containing N.  Note that the transcript summary's
    # exemplar '(86,2,5,5,14,13)' is not N's row; N is the row below.
    crt_moduli = (9, 11, 25, 49, 13, 17, 19, 23)
    crt_residues = tuple(N % modulus for modulus in crt_moduli)
    crt_factors = (
        by_key[(3, 2)]["target_factor"],
        by_key[(11, 1)]["target_factor"],
        by_key[(5, 2)]["target_factor"],
        by_key[(7, 2)]["target_factor"],
        by_key[(13, 1)]["target_factor"],
        by_key[(17, 1)]["target_factor"],
        by_key[(19, 1)]["target_factor"],
        by_key[(23, 1)]["target_factor"],
    )
    crt_product = Fraction(1)
    for value in crt_factors:
        crt_product *= value
    crt_modulus = 1
    for i, value in enumerate(crt_moduli):
        require(
            all(gcd(value, earlier) == 1 for earlier in crt_moduli[:i]),
            "CRT moduli are not pairwise coprime",
        )
        crt_modulus *= value

    # The exact minimizers in this same independent CRT packet give six most
    # suppressed classes.  This is a finite congruence-density statement, not
    # a claim that any of the six classes contains a global counterexample.
    minimum_residue_sets = (
        by_key[(3, 2)]["minimum_residues"],
        by_key[(11, 1)]["minimum_residues"],
        by_key[(5, 2)]["minimum_residues"],
        by_key[(7, 2)]["minimum_residues"],
        by_key[(13, 1)]["minimum_residues"],
        by_key[(17, 1)]["minimum_residues"],
        by_key[(19, 1)]["minimum_residues"],
        by_key[(23, 1)]["minimum_residues"],
    )
    minimum_factors = (
        by_key[(3, 2)]["minimum"],
        by_key[(11, 1)]["minimum"],
        by_key[(5, 2)]["minimum"],
        by_key[(7, 2)]["minimum"],
        by_key[(13, 1)]["minimum"],
        by_key[(17, 1)]["minimum"],
        by_key[(19, 1)]["minimum"],
        by_key[(23, 1)]["minimum"],
    )
    minimum_product = Fraction(1)
    for value in minimum_factors:
        minimum_product *= value

    def crt(residues: tuple[int, ...]) -> int:
        answer = 0
        for residue, modulus in zip(residues, crt_moduli):
            cofactor = crt_modulus // modulus
            answer += residue * cofactor * pow(cofactor, -1, modulus)
        return answer % crt_modulus

    minimum_classes = tuple(
        (residues, crt(residues)) for residues in product(*minimum_residue_sets)
    )
    require(len(minimum_classes) == 6, "minimum CRT multiplicity drift")

    # A compact semantic digest covers every exact factor and critical count.
    semantic_payload = "\n".join(
        f"{r['p']},{r['a']},{r['q']},{r['target']},{r['image']},"
        f"{fmt_fraction(r['target_factor'])},{fmt_fraction(r['minimum'])},"
        f"{','.join(map(str, r['minimum_residues']))}"
        for r in rows
    )
    semantic_hash = sha256(semantic_payload.encode("ascii")).hexdigest()

    print("SUN_2468_LOCAL_MODULAR_EXACT_AUDIT")
    print(f"N={N}")
    print(f"period_formula_gates={period_gates}")
    for row in rows:
        residues = row["minimum_residues"]
        residue_text = ",".join(map(str, residues[:12])) + ("..." if len(residues) > 12 else "")
        print(
            f"p={row['p']} a={row['a']} q={row['q']} target={row['target']} "
            f"image={row['image']}/{row['q']} target_factor={fmt_fraction(row['target_factor'])} "
            f"min={fmt_fraction(row['minimum'])}@[{residue_text}] "
            f"max={fmt_fraction(row['maximum'])}"
        )
    for prime in sorted(regular):
        singular, sizes = regular[prime]
        stability = "PROVED" if singular == 0 else "NOT_CLAIMED_SINGULAR_FIBRE"
        print(
            f"hensel_p={prime} target={N % prime} singular_target_states={singular} "
            f"critical_coordinate_sizes={sizes} all_power_stability={stability}"
        )
    for (prime, exponent), (singular, sizes) in sorted(small_regular.items()):
        print(
            f"small_lift_p={prime} level={exponent} target={N % (prime**exponent)} "
            f"singular_target_states={singular} critical_coordinate_sizes={sizes} "
            f"scope=ONE_LEVEL_CERTIFICATE"
        )
    print(f"crt_moduli={crt_moduli}")
    print(f"crt_residues={crt_residues}")
    print("crt_factors=" + ",".join(fmt_fraction(value) for value in crt_factors))
    print(f"crt_product_factor={fmt_fraction(crt_product)} approx={float(crt_product):.12f}")
    print(f"crt_modulus={crt_modulus} N_mod_crt={N % crt_modulus}")
    print(
        f"crt_packet_minimum_factor={fmt_fraction(minimum_product)} "
        f"approx={float(minimum_product):.12f}"
    )
    for residues, combined_residue in minimum_classes:
        print(f"crt_packet_minimizer residues={residues} combined={combined_residue}")
    print(f"semantic_sha256={semantic_hash}")
    print("local_insolubility=NONE_IN_EXACT_UNIVERSES")
    print("status=FINITE_EXACT_PLUS_LARGE_PRIME_HENSEL_THEOREM;_GLOBAL_ZERO_NOT_IMPLIED")
    print("PASS")


if __name__ == "__main__":
    main()
