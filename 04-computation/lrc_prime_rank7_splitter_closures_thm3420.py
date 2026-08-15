#!/usr/bin/env python3
"""Exact companion for THM-3420, prime rank-seven splitter closures.

The proof is analytic.  This standard-library companion freezes the periodic
capacity defects, both Faulhaber/Newton calculations, the Lucas certificate
for the only factor above 64 bits, the two positive atoms, and hostile
compatibility graphs.  Every truth gate remains active under ``python -O``.
"""

from fractions import Fraction
from hashlib import sha256
from math import gcd


EXPECTED_SEMANTIC_DIGEST = "4dffcebc010848fc022824c46a916208d8b46c70ac3306168463b613ce58a7ed"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def danger_mask(q, residue, epsilon):
    modulus = 2 * q
    mask = 0
    for sheet in range(q):
        word = residue * (2 * sheet + epsilon) % modulus
        if 14 * min(word, modulus - word) < modulus:
            mask |= 1 << sheet
    return mask


def zero_count(order):
    return 1 + 2 * ((order - 1) // 14)


def odd_count(order):
    return 2 * ((((order - 1) // 7) + 1) // 2)


def lagrange_value(values, point):
    """Evaluate the polynomial through y(0),...,y(n) at ``point``."""

    total = Fraction(0)
    for index, value in enumerate(values):
        basis = Fraction(1)
        for other in range(len(values)):
            if other != index:
                basis *= Fraction(point - other, index - other)
        total += value * basis
    return total


def power_sum_value(power, point):
    # The sum 1^power+...+n^power has degree power+1.  Interpolation avoids
    # importing a symbolic algebra system into the proof artifact.
    values = []
    running = 0
    for n in range(power + 2):
        if n:
            running += n**power
        values.append(running)
    return lagrange_value(tuple(values), point)


def odd_power_sum_value(power, point):
    # 1^power+3^power+...+(2k+1)^power.
    return power_sum_value(power, 2 * point + 1) - 2**power * power_sum_value(power, point)


def factor_product(factors):
    value = 1
    for prime, exponent in factors.items():
        value *= prime**exponent
    return value


def is_prime_u64(value):
    """Deterministic Miller--Rabin on the unsigned 64-bit range."""

    if value < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for prime in small:
        if value % prime == 0:
            return value == prime
    require(value < 2**64, (value, "outside deterministic u64 gate"))
    odd = value - 1
    shifts = 0
    while odd % 2 == 0:
        shifts += 1
        odd //= 2
    for base in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
        if base % value == 0:
            continue
        witness = pow(base, odd, value)
        if witness in (1, value - 1):
            continue
        for _ in range(shifts - 1):
            witness = witness * witness % value
            if witness == value - 1:
                break
        else:
            return False
    return True


LARGE_PRIME = 7_779_713_652_980_688_586_832_393
LARGE_PRIME_MINUS_ONE = {
    2: 3,
    7: 2,
    13: 1,
    191: 1,
    619: 1,
    643: 1,
    4649: 1,
    4_319_561_459: 1,
}


def lucas_certificate():
    """Verify a complete Lucas primality certificate for ``LARGE_PRIME``."""

    candidate = LARGE_PRIME
    require(factor_product(LARGE_PRIME_MINUS_ONE) == candidate - 1, "N-1 factorization")
    require(all(is_prime_u64(prime) for prime in LARGE_PRIME_MINUS_ONE), "N-1 prime factors")
    base = 3
    require(pow(base, candidate - 1, candidate) == 1, "Lucas Fermat gate")
    gcd_rows = tuple(
        (prime, gcd(pow(base, (candidate - 1) // prime, candidate) - 1, candidate))
        for prime in LARGE_PRIME_MINUS_ONE
    )
    require(all(value == 1 for _, value in gcd_rows), gcd_rows)
    return candidate, tuple(LARGE_PRIME_MINUS_ONE.items()), base, gcd_rows


def defect_audit():
    zero_delta = (-7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6)
    odd_delta = (0, -1, -2, -3, -4, -5, -6, -7, 6, 5, 4, 3, 2, 1)
    for order in range(2, 1001):
        require(7 * zero_count(order) - order == zero_delta[order % 14], order)
        require(7 * odd_count(order) - order == odd_delta[order % 14], order)

    # These are the two infinite equality profiles left by capacity,
    # reflection, and the primitive parity/fixed-point invoice alone.
    formal_zero = tuple(
        (order, 7 * zero_count(order) - order - 6)
        for order in (29, 43, 71, 113)
    )
    formal_half = tuple(
        (
            order,
            6 * odd_count(order) + zero_count(order) - order,
        )
        for order in (13, 41, 83, 97)
    )
    require(all(tax == 0 for _, tax in formal_zero), formal_zero)
    require(all(excess == 0 for _, excess in formal_half), formal_half)
    return zero_delta, odd_delta, formal_zero, formal_half


FIXED_MOMENTS = (
    (2, Fraction(-13, 1372), {13: 1}),
    (4, Fraction(611, 268_912), {13: 1, 47: 1}),
    (6, Fraction(-605_371, 368_947_264), {13: 1, 46_567: 1}),
    (8, Fraction(23_768_771, 10_330_523_392), {13: 1, 17: 1, 131: 1, 821: 1}),
    (
        10,
        Fraction(-10_591_414_093, 2_024_782_584_832),
        {13: 1, 19: 1, 41: 1, 1_045_859: 1},
    ),
    (
        12,
        Fraction(6_936_376_509_731, 396_857_386_627_072),
        {211: 1, 32_873_822_321: 1},
    ),
    (
        14,
        Fraction(-6_266_545_343_716_333, 77_784_047_778_906_112),
        {13: 1, 29: 1, 16_622_136_190_229: 1},
    ),
)


def fixed_moment_audit():
    point = Fraction(-1, 14)
    rows = []
    candidates = set()
    for power, expected, factors in FIXED_MOMENTS:
        observed = power_sum_value(power, point)
        require(observed == expected, (power, observed, expected))
        require(factor_product(factors) == abs(observed.numerator), (power, factors))
        require(all(is_prime_u64(prime) for prime in factors), (power, "factor primality"))
        candidates.update(prime for prime in factors if prime % 14 == 1)
        rows.append((power, observed, tuple(factors.items())))
    require(candidates == {29, 211}, candidates)
    return tuple(rows), tuple(sorted(candidates))


HALF_S = (
    Fraction(13, 1372),
    Fraction(-611, 268_912),
    Fraction(605_371, 368_947_264),
    Fraction(-23_768_771, 10_330_523_392),
    Fraction(10_591_414_093, 2_024_782_584_832),
    Fraction(-6_936_376_509_731, 396_857_386_627_072),
    Fraction(6_266_545_343_716_333, 77_784_047_778_906_112),
)

HALF_R = (
    Fraction(-8, 343),
    Fraction(544, 16_807),
    Fraction(-588_416, 5_764_801),
    Fraction(23_601_664, 40_353_607),
    Fraction(-10_572_781_568, 1_977_326_743),
    Fraction(6_933_325_225_984, 96_889_010_407),
    Fraction(-6_265_856_151_093_248, 4_747_561_509_943),
)

HALF_Q = (
    Fraction(13, 32),
    Fraction(611, 8704),
    Fraction(605_371, 37_658_624),
    Fraction(23_768_771, 6_042_025_984),
    Fraction(10_591_414_093, 10_826_528_325_632),
    Fraction(6_936_376_509_731, 28_398_900_125_630_464),
    Fraction(216_087_770_472_977, 3_539_992_661_362_475_008),
)

HALF_E = (
    Fraction(13, 32),
    Fraction(1651, 34_816),
    Fraction(11_646_661, 5_121_572_864),
    Fraction(24_615_544_210_543, 513_729_978_022_494_208),
    Fraction(38_463_412_097_547_188_869, 84_868_044_415_082_372_683_268_096),
    Fraction(
        25_169_871_946_198_759_841_867_924_664_541,
        21_132_471_965_361_157_685_281_632_380_249_058_574_336,
    ),
)

HALF_RESIDUAL = Fraction(
    279_717_535_826_786_771_570_267_752_942_596_708_036_779,
    4_458_951_692_213_221_631_351_994_734_945_496_909_892_106_026_221_568,
)

HALF_RESIDUAL_NUMERATOR = {
    7: 10,
    13: 2,
    277: 1,
    2719: 1,
    LARGE_PRIME: 1,
}

HALF_RESIDUAL_DENOMINATOR = {
    2: 39,
    17: 3,
    31: 1,
    61: 1,
    83: 1,
    311: 1,
    1487: 1,
    4597: 2,
    84_631: 1,
    255_443: 1,
    49_785_481: 1,
}


def half_moment_audit():
    point = Fraction(-13, 14)
    s_rows = tuple(power_sum_value(2 * n, point) for n in range(1, 8))
    r_rows = tuple(odd_power_sum_value(2 * n, point) for n in range(1, 8))
    require(s_rows == HALF_S, s_rows)
    require(r_rows == HALF_R, r_rows)
    q_rows = tuple(-left / right for left, right in zip(s_rows, r_rows))
    require(q_rows == HALF_Q, q_rows)

    elementary = [Fraction(1)]
    for n in range(1, 7):
        value = sum(
            (-1) ** (index - 1) * elementary[n - index] * q_rows[index - 1]
            for index in range(1, n + 1)
        ) / n
        elementary.append(value)
    require(tuple(elementary[1:]) == HALF_E, tuple(elementary[1:]))
    predicted_seventh = sum(
        (-1) ** (index - 1) * elementary[index] * q_rows[6 - index]
        for index in range(1, 7)
    )
    residual = q_rows[6] - predicted_seventh
    require(residual == HALF_RESIDUAL, residual)

    require(
        factor_product(HALF_RESIDUAL_NUMERATOR) == abs(residual.numerator),
        "half residual numerator",
    )
    require(
        factor_product(HALF_RESIDUAL_DENOMINATOR) == residual.denominator,
        "half residual denominator",
    )
    lucas = lucas_certificate()
    require(
        all(
            prime == LARGE_PRIME or is_prime_u64(prime)
            for prime in HALF_RESIDUAL_NUMERATOR
        ),
        "half numerator primes",
    )
    require(all(is_prime_u64(prime) for prime in HALF_RESIDUAL_DENOMINATOR), "half denominator primes")
    require(
        {prime for prime in HALF_RESIDUAL_NUMERATOR if prime % 14 == 13} == {13},
        "half numerator residue filter",
    )
    exceptional = tuple(
        sorted(prime for prime in HALF_RESIDUAL_DENOMINATOR if prime % 14 == 13)
    )
    require(exceptional == (83, 255_443), exceptional)
    # Both exceptional denominators come from R_14=0, while S_14 is nonzero.
    for prime in exceptional:
        require(HALF_R[6].numerator % prime == 0, (prime, "R14"))
        require(HALF_S[6].numerator % prime != 0, (prime, "S14"))
    return (
        tuple(zip(range(1, 8), s_rows, r_rows, q_rows)),
        tuple(elementary[1:]),
        predicted_seventh,
        residual,
        tuple(HALF_RESIDUAL_NUMERATOR.items()),
        tuple(HALF_RESIDUAL_DENOMINATOR.items()),
        exceptional,
        lucas,
    )


def multiplicity_record(q, epsilon, residues):
    masks = tuple(danger_mask(q, residue, epsilon) for residue in residues)
    joined = 0
    for mask in masks:
        joined |= mask
    multiplicities = tuple(
        sum((mask >> sheet) & 1 for mask in masks)
        for sheet in range(q)
    )
    require(joined == (1 << q) - 1, (q, epsilon, residues))
    return (
        q,
        epsilon,
        residues,
        tuple(mask.bit_count() for mask in masks),
        tuple((value, multiplicities.count(value)) for value in sorted(set(multiplicities))),
        sum(mask.bit_count() for mask in masks) - q,
    )


def positive_atoms():
    fixed = multiplicity_record(29, 0, (1, 4, 5, 6, 7, 9, 13))
    half = multiplicity_record(13, 1, (1, 2, 3, 5, 7, 9, 11))
    require(fixed[4] == ((1, 28), (7, 1)), fixed)
    require(half[4] == ((1, 13),), half)

    residues = fixed[2]
    splitter = {
        sign * pow(residue, -1, 29) % 29
        for residue in residues
        for sign in (1, -1)
    }
    quadratic_residues = {value * value % 29 for value in range(1, 29)}
    require(splitter == quadratic_residues and len(splitter) == 14, splitter)
    require(pow(2, 14, 29) == 28, "2 is a nonsquare modulo 29")
    products = tuple(sorted((multiplier * value) % 29 for multiplier in (1, 2) for value in splitter))
    require(products == tuple(range(1, 29)), products)
    return fixed, half, tuple(sorted(splitter)), products


def q211_hostile():
    prime = 211
    interval = {value % prime for value in range(-15, 16) if value}
    ratios = {
        left * pow(right, -1, prime) % prime
        for left in interval
        for right in interval
    }
    require(ratios == set(range(1, prime)), "I/I is all of F_211^*")
    gap_rows = []
    for multiplier in range(1, prime):
        points = sorted(index * multiplier % prime for index in range(16))
        gaps = tuple(
            (points[(index + 1) % 16] - points[index]) % prime
            for index in range(16)
        )
        gap_rows.append(min(gaps))
    require(max(gap_rows) <= 13, max(gap_rows))
    return len(interval), len(ratios), max(gap_rows), sha256(repr(tuple(gap_rows)).encode("ascii")).hexdigest()


def maximum_clique(vertices, compatible):
    best = ()
    calls = 0

    def search(chosen, candidates):
        nonlocal best, calls
        calls += 1
        if len(chosen) + len(candidates) <= len(best):
            return
        while candidates:
            if len(chosen) + len(candidates) <= len(best):
                return
            vertex = candidates[0]
            rest = candidates[1:]
            search(
                chosen + (vertex,),
                tuple(other for other in rest if compatible(vertex, other)),
            )
            candidates = rest
        if len(chosen) > len(best):
            best = chosen

    search((), tuple(vertices))
    return best, calls


def fixed_clique_record(prime):
    grouped = {}
    for residue in range(1, prime):
        mask = danger_mask(prime, residue, 0)
        grouped.setdefault(mask, residue)
    bank = tuple((residue, mask) for mask, residue in grouped.items())
    best, calls = maximum_clique(
        tuple(range(len(bank))),
        lambda left, right: bank[left][1] & bank[right][1] == 1,
    )
    return (
        prime,
        len(bank),
        tuple(sorted({mask.bit_count() for _, mask in bank})),
        len(best),
        tuple(bank[index][0] for index in best),
        calls,
    )


def half_clique_record(prime):
    grouped = ({}, {})
    for residue in range(1, 2 * prime):
        if residue == prime:
            continue
        mask = danger_mask(prime, residue, 1)
        grouped[residue % 2].setdefault(mask, residue)
    even = tuple((residue, mask) for mask, residue in grouped[0].items())
    odd = tuple((residue, mask) for mask, residue in grouped[1].items())
    overall = ()
    selected_even = None
    total_calls = 0
    for even_residue, even_mask in even:
        candidates = tuple(
            index for index, (_, odd_mask) in enumerate(odd)
            if not even_mask & odd_mask
        )
        best, calls = maximum_clique(
            candidates,
            lambda left, right: not (odd[left][1] & odd[right][1]),
        )
        total_calls += calls
        if len(best) > len(overall):
            overall = best
            selected_even = even_residue
    return (
        prime,
        len(odd),
        len(even),
        tuple(sorted({mask.bit_count() for _, mask in odd})),
        tuple(sorted({mask.bit_count() for _, mask in even})),
        len(overall),
        selected_even,
        tuple(odd[index][0] for index in overall),
        total_calls,
    )


def hostile_controls():
    fixed = tuple(fixed_clique_record(prime) for prime in (29, 43, 71, 211))
    half = tuple(half_clique_record(prime) for prime in (13, 41, 83))
    require(tuple(row[3] for row in fixed) == (7, 6, 4, 1), fixed)
    require(tuple(row[5] for row in half) == (6, 5, 3), half)
    return fixed, half, q211_hostile()


def source_overlap_controls():
    rows = (
        multiplicity_record(13, 1, (1, 2, 3, 5, 7, 9, 11)),
        multiplicity_record(14, 1, (1, 3, 4, 5, 9, 11, 13)),
        multiplicity_record(29, 0, (1, 4, 5, 6, 7, 9, 13)),
        multiplicity_record(29, 1, (1, 5, 7, 8, 12, 13, 22)),
        multiplicity_record(38, 1, (1, 9, 17, 20, 21, 29, 37)),
        multiplicity_record(51, 1, (1, 11, 12, 18, 23, 34, 35)),
        multiplicity_record(68, 1, (8, 11, 23, 24, 45, 56, 57)),
        multiplicity_record(148, 1, (8, 33, 41, 100, 107, 115, 140)),
    )
    require(tuple(row[5] for row in rows) == (0, 0, 6, 2, 4, 16, 8, 8), rows)
    return rows


def main():
    defects = defect_audit()
    fixed_moments = fixed_moment_audit()
    half_moments = half_moment_audit()
    atoms = positive_atoms()
    hostiles = hostile_controls()
    overlaps = source_overlap_controls()

    semantic_payload = (
        defects,
        fixed_moments,
        half_moments,
        atoms,
        hostiles,
        overlaps,
    )
    semantic = sha256(repr(semantic_payload).encode("ascii")).hexdigest()
    if EXPECTED_SEMANTIC_DIGEST is not None:
        require(semantic == EXPECTED_SEMANTIC_DIGEST, (semantic, EXPECTED_SEMANTIC_DIGEST))

    print("THM-3420 exact prime rank-seven splitter companion")
    print(f"defect_tables_and_formal_tails={defects}")
    print(f"fixed_Faulhaber_rows_and_candidates={fixed_moments}")
    print(
        "half_Newton="
        f"(moment_rows,elementary,predicted7,residual,num_factors,den_factors,exceptions)="
        f"{half_moments[:7]}"
    )
    print(f"large_prime_Lucas_certificate={half_moments[7]}")
    print(f"positive_atoms_and_Paley_splitter={atoms}")
    print(f"compatibility_hostiles={hostiles}")
    print(f"rank7_source_overlap_controls={overlaps}")
    print(
        "scope=prime fixed-zero cap7 iff p29; prime half-twist cap7 in p=13mod14 iff p13; "
        "other half prime classes and composite rank7 remain open; no LRC14 decrement"
    )
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
