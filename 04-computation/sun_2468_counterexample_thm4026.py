#!/usr/bin/env python3
"""Exact audit of the proposed counterexample to Sun's 2-4-6-8 conjecture.

The search is exhaustive in (x,y,z).  The remaining triangular summand is
tested by the equivalence

    r = binom(w, 2)  <=>  8*r + 1 is an odd square.

Odd-prime residue masks are necessary filters only; every survivor is checked
again with exact integer arithmetic.  Binomial coefficients are formed over
the integers before reduction, so no factorial denominator is inverted modulo
a prime (the MISTAKE-363 failure mode).
"""

from __future__ import annotations

from bisect import bisect_right
from fractions import Fraction
from hashlib import sha256
from math import comb, isqrt
import json


TARGET = 896_315_812_331_399
CONTROL = 4_655
SIEVE_PRIMES = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43)
RARE_WHEEL_PRIMES = (3, 5, 7, 11, 17, 19, 23)
MOD_CHECK = 1_000_000_007


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def canonical_values(k: int, target: int) -> tuple[int, list[int]]:
    """Return the distinct value support {C(t,k): t>=2, C(t,k)<=target}.

    All 2 <= t < k give zero.  The single canonical zero representative is
    t=k-1, after which C(t,k) is strictly increasing.
    """

    start = max(2, k - 1)
    values: list[int] = []
    t = start
    while True:
        value = comb(t, k)
        if value > target:
            break
        values.append(value)
        t += 1
    require(values and values[0] == 0, f"missing canonical zero for k={k}")
    require(all(a < b for a, b in zip(values, values[1:])),
            f"non-increasing support for k={k}")
    return start, values


def support_boundary(k: int, target: int) -> tuple[int, int, int]:
    """Largest t with C(t,k)<=target, plus the two bracketing values."""

    lo, hi = k - 1, max(k, 2)
    while comb(hi, k) <= target:
        hi *= 2
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if comb(mid, k) <= target:
            lo = mid
        else:
            hi = mid
    return lo, comb(lo, k), comb(lo + 1, k)


def triangular_residues(p: int) -> frozenset[int]:
    # For odd p, 2 is invertible and C(w+p,2)=C(w,2) (mod p).
    return frozenset(comb(w, 2) % p for w in range(p))


def build_x_masks(x_values: list[int]) -> dict[int, tuple[int, ...]]:
    """masks[p][t] has bit i iff t-C(x_i,4) is triangular modulo p."""

    masks: dict[int, tuple[int, ...]] = {}
    for p in SIEVE_PRIMES:
        triangular = triangular_residues(p)
        x_mod = [value % p for value in x_values]
        rows: list[int] = []
        for target_residue in range(p):
            mask = 0
            for i, residue in enumerate(x_mod):
                if (target_residue - residue) % p in triangular:
                    mask |= 1 << i
            rows.append(mask)
        masks[p] = tuple(rows)
    return masks


def audit_target(target: int) -> dict[str, object]:
    x_start, x_values = canonical_values(4, target)
    y_start, y_values = canonical_values(6, target)
    z_start, z_values = canonical_values(8, target)
    masks = build_x_masks(x_values)
    all_x = (1 << len(x_values)) - 1

    yz_pairs = 0
    admissible_triples = 0
    modular_survivors = 0
    exact_tests = 0
    x_checksum = 0
    yz_checksum = 0
    modular_digest = sha256()
    exact_digest = sha256()
    representations: list[tuple[int, int, int, int]] = []

    for z_offset, c8 in enumerate(z_values):
        z = z_start + z_offset
        for y_offset, c6 in enumerate(y_values):
            y = y_start + y_offset
            base = target - c8 - c6
            if base < 0:
                break
            yz_pairs += 1
            admissible_triples += bisect_right(x_values, base)
            candidates = all_x
            for p in SIEVE_PRIMES:
                candidates &= masks[p][base % p]
                if not candidates:
                    break

            while candidates:
                low_bit = candidates & -candidates
                x_offset = low_bit.bit_length() - 1
                candidates -= low_bit
                x = x_start + x_offset
                modular_survivors += 1
                modular_digest.update(f"{x},{y},{z}\n".encode("ascii"))

                c4 = x_values[x_offset]
                if c4 > base:
                    continue
                exact_tests += 1
                remainder = base - c4
                x_checksum = (x_checksum + x) % MOD_CHECK
                yz_checksum = (yz_checksum + 1000 * y + z) % MOD_CHECK
                exact_digest.update(
                    f"{x},{y},{z},{remainder}\n".encode("ascii")
                )

                discriminant = 8 * remainder + 1
                root = isqrt(discriminant)
                if root * root == discriminant and remainder > 0:
                    w = (root + 1) // 2
                    require(w >= 2 and comb(w, 2) == remainder,
                            "square test produced an invalid triangular witness")
                    representations.append((w, x, y, z))

    return {
        "target": target,
        "ranges": {
            "x": [x_start, x_start + len(x_values) - 1],
            "y": [y_start, y_start + len(y_values) - 1],
            "z": [z_start, z_start + len(z_values) - 1],
        },
        "support_sizes": [len(x_values), len(y_values), len(z_values)],
        "yz_pairs": yz_pairs,
        "admissible_xyz_triples": admissible_triples,
        "rectangular_xyz_triples": len(x_values) * len(y_values) * len(z_values),
        "modular_survivors": modular_survivors,
        "exact_tests": exact_tests,
        "x_checksum_mod_1000000007": x_checksum,
        "yz_checksum_mod_1000000007": yz_checksum,
        "modular_survivor_sha256": modular_digest.hexdigest(),
        "exact_test_sha256": exact_digest.hexdigest(),
        "representations": representations,
    }


def cyclic_convolution(left: list[int], right: list[int], p: int) -> list[int]:
    out = [0] * p
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[(i + j) % p] += a * b
    return out


def prime_period(k: int, p: int) -> int:
    """A proved period modulo p: the least p-power strictly greater than k."""

    period = 1
    while period <= k:
        period *= p
    return period


def local_distribution(p: int) -> tuple[list[int], tuple[int, ...]]:
    total = [1] + [0] * (p - 1)
    periods: list[int] = []
    for k in (2, 4, 6, 8):
        period = prime_period(k, p)
        periods.append(period)
        counts = [0] * p
        for t in range(period):
            counts[comb(t, k) % p] += 1
        # Frobenius proves this period.  The explicit replay catches indexing
        # and exact-integer/modular implementation errors.
        for t in range(2 * period):
            require(comb(t + period, k) % p == comb(t, k) % p,
                    f"period replay failed for p={p}, k={k}")
        total = cyclic_convolution(total, counts, p)
    return total, tuple(periods)


def crt(residues: tuple[int, ...], moduli: tuple[int, ...]) -> tuple[int, int]:
    modulus = 1
    for m in moduli:
        modulus *= m
    value = 0
    for residue, m in zip(residues, moduli):
        cofactor = modulus // m
        value += residue * cofactor * pow(cofactor, -1, m)
    return value % modulus, modulus


def local_anatomy() -> dict[str, object]:
    rows: list[dict[str, object]] = []
    rare_residue_options: list[tuple[int, ...]] = []
    combined_factor = Fraction(1, 1)

    for p in (3, 5, 7, 11, 13, 17, 19, 23):
        counts, periods = local_distribution(p)
        total_words = 1
        for period in periods:
            total_words *= period
        require(total_words % p == 0, "nonintegral uniform local average")
        average = total_words // p
        residue = TARGET % p
        minimum = min(counts)
        minimizers = tuple(i for i, count in enumerate(counts) if count == minimum)
        factor = Fraction(counts[residue], average)
        rows.append({
            "p": p,
            "target_residue": residue,
            "periods": list(periods),
            "count": counts[residue],
            "uniform_average": average,
            "factor": str(factor),
            "minimizers": list(minimizers),
            "is_minimum": residue in minimizers,
        })
        if p in RARE_WHEEL_PRIMES:
            require(residue in minimizers,
                    f"target is not locally minimal at wheel prime {p}")
            combined_factor *= factor
            rare_residue_options.append(minimizers)

    # There are two wheel classes because p=19 has two minimizing residues.
    wheel_classes: list[int] = []
    for r19 in rare_residue_options[5]:
        residues = (
            rare_residue_options[0][0],
            rare_residue_options[1][0],
            rare_residue_options[2][0],
            rare_residue_options[3][0],
            rare_residue_options[4][0],
            r19,
            rare_residue_options[6][0],
        )
        value, modulus = crt(residues, RARE_WHEEL_PRIMES)
        wheel_classes.append(value)
    wheel_classes.sort()
    require(TARGET % modulus in wheel_classes, "target missed the rare CRT wheel")

    return {
        "rows": rows,
        "rare_wheel_primes": list(RARE_WHEEL_PRIMES),
        "rare_wheel_modulus": modulus,
        "rare_wheel_classes": wheel_classes,
        "target_wheel_class": TARGET % modulus,
        "combined_relative_density": str(combined_factor),
    }


def legendre(a: int, p: int) -> int:
    value = pow(a % p, (p - 1) // 2, p)
    return -1 if value == p - 1 else value


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


def quartic_hole_controls() -> dict[str, object]:
    # Exact integral identity behind the low-pair local hole.
    identity_checks = 0
    for w in range(0, 101):
        for x in range(0, 101):
            left = 384 * (comb(w, 2) + comb(x, 4)) + 64
            right = ((2 * x - 3) ** 2 - 5) ** 2 + 48 * (2 * w - 1) ** 2
            require(left == right, "quartic-hole identity failed")
            identity_checks += 1

    checked_primes: list[int] = []
    for p in range(7, 501):
        if not is_prime(p) or p % 30 not in (17, 23):
            continue
        require(legendre(-48, p) == -1 and legendre(5, p) == -1,
                f"quadratic-character premise failed at p={p}")
        forbidden = (-pow(6, -1, p)) % p
        triangular = {comb(w, 2) % p for w in range(p)}
        quartic = {comb(x, 4) % p for x in range(p)}
        require(forbidden not in {(a + b) % p for a in triangular for b in quartic},
                f"claimed low-pair hole filled at p={p}")
        checked_primes.append(p)

    congruences = {p: TARGET % p for p in (17, 23)}
    require(all(congruences[p] == (-pow(6, -1, p)) % p for p in congruences),
            "target does not occupy both proved low-pair holes")
    return {
        "identity_checks": identity_checks,
        "prime_control_bound": 500,
        "checked_primes": checked_primes,
        "target_hole_congruences": congruences,
    }


def combinadic_normal_form(n: int, rank: int) -> list[tuple[int, int, int]]:
    """Greedy Macaulay/combinatorial-number-system expansion through `rank`."""

    remainder = n
    upper: int | None = None
    expansion: list[tuple[int, int, int]] = []
    for k in range(rank, 0, -1):
        lo = k - 1  # C(k-1,k)=0 is always admissible.
        if upper is None:
            hi = max(k + 1, 2)
            while comb(hi, k) <= remainder:
                hi *= 2
        else:
            hi = upper
        while lo + 1 < hi:
            mid = (lo + hi) // 2
            if comb(mid, k) <= remainder:
                lo = mid
            else:
                hi = mid
        value = comb(lo, k)
        expansion.append((lo, k, value))
        remainder -= value
        upper = lo
    require(remainder == 0, "combinadic normal form left a remainder")
    require(all(a > b for (a, _, _), (b, _, _) in zip(expansion, expansion[1:])),
            "combinadic upper indices are not strictly descending")
    return expansion


def main() -> None:
    target = audit_target(TARGET)
    require(target["ranges"] == {"x": [3, 12112], "y": [5, 932], "z": [7, 281]},
            "candidate range boundary changed")
    require(target["yz_pairs"] == 248_160, "candidate yz universe changed")
    require(target["admissible_xyz_triples"] == 2_755_643_831,
            "candidate xyz universe changed")
    require(target["modular_survivors"] == 287_120,
            "candidate modular survivor count changed")
    require(target["exact_tests"] == 263_434, "candidate exact universe changed")
    require(target["representations"] == [], "the candidate is representable")

    control = audit_target(CONTROL)
    expected_control = [(85, 14, 9, 7), (94, 7, 9, 11)]
    require(sorted(control["representations"]) == sorted(expected_control),
            "known two-representation control failed")

    anatomy = local_anatomy()
    holes = quartic_hole_controls()
    normal_form = combinadic_normal_form(TARGET, 8)
    require(normal_form == [
        (281, 8, 871_896_500_955_975),
        (279, 7, 24_202_109_279_205),
        (234, 6, 213_748_248_714),
        (212, 5, 3_403_031_632),
        (188, 4, 50_404_915),
        (136, 3, 410_040),
        (43, 2, 903),
        (15, 1, 15),
    ], "target combinadic normal form changed")
    boundaries = {str(k): support_boundary(k, TARGET) for k in (2, 4, 6, 8)}
    require(boundaries == {
        "2": (42_339_481, 896_315_804_504_940, 896_315_846_844_421),
        "4": (12_112, 896_266_258_399_820, 896_562_324_551_740),
        "6": (932, 895_693_597_430_352, 901_490_966_993_008),
        "8": (281, 871_896_500_955_975, 897_353_333_100_675),
    }, "support boundaries changed")
    neighbor_witnesses = {
        str(TARGET - 1): (33_663_667, 9_433, 16, 9),
        str(TARGET + 1): (40_920_205, 6_138, 22, 13),
    }
    for value_text, witness in neighbor_witnesses.items():
        value = int(value_text)
        require(sum(comb(t, k) for t, k in zip(witness, (2, 4, 6, 8))) == value,
                f"neighbor witness failed at {value}")
    report = {
        "status": "PASS",
        "scope": "single-target exact counterexample; no minimality claim",
        "sieve_primes": list(SIEVE_PRIMES),
        "target_audit": target,
        "positive_control": control,
        "local_anatomy": anatomy,
        "quartic_hole_controls": holes,
        "rank_8_combinadic_normal_form": normal_form,
        "support_boundaries": boundaries,
        "represented_immediate_neighbors": neighbor_witnesses,
    }
    semantic = sha256(
        json.dumps(report, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    print(json.dumps(report, sort_keys=True, indent=2))
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
