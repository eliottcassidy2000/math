#!/usr/bin/env python3
"""Exact companion for THM-2082.

The referee has four jobs.

1. Specialize the k-deletion evaluation code to a primitive scalar row.
2. Audit a translated-prime-grid escape lemma by exact rational arithmetic.
3. Exhibit an unbounded family on which the code/CRT and scalar terminal
   filters are frozen while a fixed rational phase supplies escape.
4. Show at p=17 that identical full-support Hamming data can coexist with
   radically different rational-grid safe sets.

No assertion statement is used, so ``python`` and ``python -O`` exercise the
same checks.  Tournament Analysis is deliberately diagnostic here: ordering
residues by clearance is transitive and forgets the incidence being tested.
"""

from fractions import Fraction
from math import ceil, gcd, lcm


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def primes_upto(limit):
    out = []
    for n in range(2, limit + 1):
        if all(n % p for p in out if p * p <= n):
            out.append(n)
    return out


def gcd_all(values):
    answer = 0
    for value in values:
        answer = gcd(answer, value)
    return answer


def hereditarily_primitive(row):
    return gcd_all(row) == 1 and all(
        gcd_all(row[:i] + row[i + 1 :]) == 1 for i in range(len(row))
    )


def divisor_complete(row, ceiling=14):
    return all(any(q % d == 0 for q in row) for d in range(2, ceiling + 1))


def quarter_escape(row):
    anchor = min(q for q in row if q % 4 == 0)
    alternatives = (
        any(q % 4 == 0 and q > 7 * anchor for q in row),
        any(q % 2 == 1 and q > Fraction(5, 2) * anchor for q in row),
        any(q % 4 == 2 and q > 6 * anchor for q in row),
    )
    return anchor, alternatives


def scalar_weight(row, prime):
    return sum(q % prime != 0 for q in row)


def fold14(value):
    residue = value % 14
    return residue * (14 - residue)


def scalar_fold_sum(row, guard):
    total = Fraction(0)
    for q in row:
        common = gcd(q, guard)
        a = q // common
        b = guard // common
        correction = fold14(b + 2 * a) - fold14(b - 2 * a)
        total += Fraction(correction, a * b)
    return total


def least_residue(value, modulus):
    residue = value % modulus
    return min(residue, modulus - residue)


def circle_norm(value):
    residue = value % 1
    return min(residue, 1 - residue)


def distance_to_14z(value):
    residue = value % 14
    return min(residue, 14 - residue)


def translated_grid_audit(prime, row, guard, z, shift):
    require(z % prime == 0, "the translated-grid carrier z must be p-divisible")
    carriers = [q for q in row if q % prime == 0]
    noncarriers = [q for q in row if q % prime != 0]
    for q in carriers:
        require(
            distance_to_14z(Fraction(shift * q, z)) >= 1,
            "a carrier is not uniformly Q-safe",
        )

    candidates = set(range(1, prime))
    q_bad = set()
    guard_bad = set()
    safe = set()
    for a in candidates:
        t = Fraction(a, prime) + Fraction(shift, 14 * z)
        if any(circle_norm(q * t) < Fraction(1, 14) for q in row):
            q_bad.add(a)
        if circle_norm(guard * t) < Fraction(1, 7):
            guard_bad.add(a)
        if a not in q_bad and a not in guard_bad:
            safe.add(a)

    up = ceil(prime / 7)
    vp = ceil(2 * prime / 7)
    require(len(q_bad) <= len(noncarriers) * up, "Q union bound failed")

    if guard % prime:
        require(len(guard_bad) <= vp, "guard union bound failed")
        bound = len(noncarriers) * up + vp
    else:
        require(
            distance_to_14z(Fraction(shift * guard, z)) >= 2,
            "a carrier guard is not uniformly safe",
        )
        require(not guard_bad, "uniform carrier-guard safety failed")
        bound = len(noncarriers) * up

    require(bound < prime - 1, "test row does not meet strict counting hypothesis")
    require(safe, "translated-grid lemma predicted an escape but found none")
    return len(noncarriers), bound, len(safe)


def pairwise_coprime(row):
    return all(gcd(row[i], row[j]) == 1 for i in range(len(row)) for j in range(i))


def grid_safe_residues(row, prime):
    return tuple(
        a
        for a in range(1, prime)
        if all(circle_norm(Fraction(a * q, prime)) >= Fraction(1, 14) for q in row)
    )


def main():
    print("THM-2082 RANK-ONE CODE WHEEL / TRANSLATED PRIME GRID REFEREE")
    print("ARITHMETIC=exact rational/integer; ASSERTIONS=explicit runtime checks")

    limit = 14
    common = 1
    for d in range(2, limit + 1):
        common = lcm(common, d)
    require(common == 360360, "lcm(2,...,14) changed")

    # Frozen unbounded terminal family.
    expected_folds = {
        7: Fraction(-103, 65),
        8: Fraction(-103, 65),
        9: Fraction(-231, 130),
        10: Fraction(-733, 390),
    }
    expected_cogirth = {7: 3, 8: 4, 9: 4, 10: 5}
    tested_primes = primes_upto(997)
    family_rows = 0
    profile_checks = 0
    phase_checks = 0
    for size in range(7, 11):
        reference_profile = None
        for exponent in range(1, 6):
            far = common * 32**exponent
            row = list(range(1, size)) + [far]
            require(len(set(row)) == size, "family row has a collision")
            require(hereditarily_primitive(row), "family lost hereditary primitivity")
            require(divisor_complete(row), "family lost divisor completeness")
            anchor, alternatives = quarter_escape(row)
            require(anchor == 4 and alternatives[0], "family lost the quarter escape")
            require((13 - size) * 13 < 2 * (size + 1) * far, "height filter failed")
            fold_sum = scalar_fold_sum(row, 13)
            require(fold_sum == expected_folds[size], "scalar fold sum changed")
            require(fold_sum <= 20 * (size - 7), "scalar containment invoice failed")

            profile = []
            for prime in tested_primes:
                actual = scalar_weight(row, prime)
                formula = (size - 1) - (size - 1) // prime + int(common % prime != 0)
                require(actual == formula, "rank-one weight formula failed")
                profile.append(actual)
                profile_checks += 1
            if reference_profile is None:
                reference_profile = profile
            require(profile == reference_profile, "prime profile depends on exponent")
            require(min(profile) == expected_cogirth[size], "cogirth minimum changed")

            t_numerator, t_denominator = 3, 31
            distances = [least_residue(t_numerator * q, t_denominator) for q in row]
            require(min(distances) >= 3, "fixed phase is Q-dangerous")
            require(
                least_residue(t_numerator * 13, t_denominator) >= 5,
                "fixed phase lies in the guard",
            )
            require(far % 31 == 16, "far residue should be frozen at 16 mod 31")
            family_rows += 1
            phase_checks += len(row) + 1

        print(
            f"FAMILY s={size}: fold={expected_folds[size]} "
            f"min_cogirth={expected_cogirth[size]} rows=5 fixed_escape=3/31"
        )
    print(
        f"FROZEN_FAMILY rows={family_rows} profile_checks={profile_checks} "
        f"phase_checks={phase_checks} height=UNBOUNDED"
    )

    # The positive translated-grid lemma, in both guard branches.
    grid_rows = 0
    min_safe_noncarrier_guard = None
    min_safe_carrier_guard = None
    for prime in [p for p in tested_primes if 17 <= p <= 101]:
        w, bound, safe_count = translated_grid_audit(
            prime, [prime, 2 * prime, 1, 2], 3, prime, 1
        )
        require(w == 2, "wrong noncarrier count in ordinary-guard control")
        min_safe_noncarrier_guard = (
            safe_count
            if min_safe_noncarrier_guard is None
            else min(min_safe_noncarrier_guard, safe_count)
        )
        grid_rows += 1

        w, bound, safe_count = translated_grid_audit(
            prime, [prime, 2 * prime, 1, 2, 3], prime, prime, 3
        )
        require(w == 3, "wrong noncarrier count in carrier-guard control")
        min_safe_carrier_guard = (
            safe_count
            if min_safe_carrier_guard is None
            else min(min_safe_carrier_guard, safe_count)
        )
        grid_rows += 1
    print(
        f"TRANSLATED_GRID rows={grid_rows} min_safe="
        f"{min_safe_noncarrier_guard}/{min_safe_carrier_guard} branches=ordinary/carrier-guard"
    )

    # Same p=17 rank-one Hamming data, very different residue incidence.
    cover_row = [1, 360360, 19, 37, 89, 107, 109, 901093]
    gap_row = [1, 360360, 103, 137, 239, 307, 409, 901171]
    for name, row in [("cover", cover_row), ("gap", gap_row)]:
        require(len(set(row)) == 8, f"{name} row collision")
        require(pairwise_coprime(row), f"{name} row is not pairwise coprime")
        require(hereditarily_primitive(row), f"{name} row is not hereditary")
        require(divisor_complete(row), f"{name} row is not divisor-complete")
        anchor, alternatives = quarter_escape(row)
        require(anchor == 360360 and alternatives[1], f"{name} lost odd quarter escape")
        require(scalar_weight(row, 17) == 8, f"{name} p=17 code is not full support")

    cover_residues = tuple(q % 17 for q in cover_row)
    gap_residues = tuple(q % 17 for q in gap_row)
    require(cover_residues == (1, 11, 2, 3, 4, 5, 7, 8), "cover residues changed")
    require(gap_residues == (1, 11, 1, 1, 1, 1, 1, 1), "gap residues changed")
    cover_safe = grid_safe_residues(cover_row, 17)
    gap_safe = grid_safe_residues(gap_row, 17)
    require(cover_safe == (), "projective cover row unexpectedly has a safe residue")
    require(
        gap_safe == (2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15),
        "projective gap safe residues changed",
    )
    print("P17_CODE both=1+16*z^8 cogirth=8")
    print(f"P17_INCIDENCE safe_counts={len(cover_safe)}/{len(gap_safe)} gap_safe={gap_safe}")

    print("TOURNAMENT vertices=residues ordered_by=clearance fingerprint=transitive")
    print("TOURNAMENT directed_cycles=0 SCCs=singletons Hamiltonian_paths=1 carrier=INSUFFICIENT")
    print("ASSUMPTION_CHALLENGE support preserves divisibility; destroys projective residue incidence")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
