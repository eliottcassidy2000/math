#!/usr/bin/env python3
"""Exact audit of the three 60-clocks and a four-dimensional Kakeya control.

The script deliberately separates:

* THM-4029's owner-selector law on n mod lcm(1,...,6);
* the Fibonacci state modulo 10;
* triangular numbers modulo 30; and
* a finite-field direction spine in P^3(F_61).

It checks both positive maps and hostile quotient losses.  Every semantic gate
uses ``require`` so the normal and ``python -O`` executions are identical.
"""

from __future__ import annotations

from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
from math import comb
import sys

from lrc14_ap_cover_finite_owner_formula_thm4029 import prove_phase_rational_law


sys.stdout.reconfigure(newline="\n")


def require(condition: bool, label) -> None:
    if not condition:
        raise RuntimeError(label)


def digest(obj) -> str:
    return sha256(repr(obj).encode("utf-8")).hexdigest()


def divisors(n: int):
    return tuple(d for d in range(1, n + 1) if n % d == 0)


def cyclic_periods(values):
    n = len(values)
    return tuple(
        d
        for d in divisors(n)
        if all(values[r] == values[(r + d) % n] for r in range(n))
    )


def fibonacci_modulus(modulus: int, count: int):
    values = [0, 1]
    while len(values) < count:
        values.append((values[-1] + values[-2]) % modulus)
    return tuple(values[:count])


def matmul2(left, right, modulus: int):
    return tuple(
        tuple(
            sum(left[i][k] * right[k][j] for k in range(2)) % modulus
            for j in range(2)
        )
        for i in range(2)
    )


def matpow2(matrix, exponent: int, modulus: int):
    result = ((1, 0), (0, 1))
    base = matrix
    while exponent:
        if exponent & 1:
            result = matmul2(result, base, modulus)
        base = matmul2(base, base, modulus)
        exponent //= 2
    return result


def matrix_order2(matrix, modulus: int, search_bound: int = 10000):
    identity = ((1, 0), (0, 1))
    value = identity
    for exponent in range(1, search_bound + 1):
        value = matmul2(value, matrix, modulus)
        if value == identity:
            return exponent
    raise RuntimeError((matrix, modulus, "matrix order search exhausted"))


def projective_pair(vector, modulus: int):
    for coordinate in vector:
        if coordinate % modulus:
            inverse = pow(coordinate, -1, modulus)
            return tuple((entry * inverse) % modulus for entry in vector)
    raise RuntimeError("zero vector has no projective class")


def triangular(n: int) -> int:
    return n * (n + 1) // 2


def triangular_period_formula(modulus: int) -> int:
    return modulus if modulus % 2 else 2 * modulus


def first_period(values, max_period: int):
    for period in range(1, max_period + 1):
        if all(values[n] == values[n + period] for n in range(len(values) - period)):
            return period
    return None


def aggregate_ap_laws():
    laws, constants = prove_phase_rational_law(period=60, minimum_n=12)
    signatures = []
    moments = []
    for phase in range(60):
        aggregate = {}
        for coefficient, shift in laws[phase][1]:
            aggregate[shift] = aggregate.get(shift, Q(0)) + coefficient
        signature = tuple(
            (shift, coefficient)
            for shift, coefficient in sorted(aggregate.items())
            if coefficient
        )
        signatures.append(signature)
        moments.append(
            (
                sum((coefficient * shift for shift, coefficient in signature), Q(0)) / 7,
                sum((coefficient * shift**2 for shift, coefficient in signature), Q(0)) / 7,
            )
        )
    return laws, constants, tuple(signatures), tuple(moments)


def fibres(values):
    answer = {}
    for phase, value in enumerate(values):
        answer.setdefault(value, []).append(phase)
    return {value: tuple(phases) for value, phases in answer.items()}


def det_mod(matrix, prime: int) -> int:
    work = [list(row) for row in matrix]
    determinant = 1
    n = len(work)
    for column in range(n):
        pivot = next((row for row in range(column, n) if work[row][column] % prime), None)
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            determinant = -determinant
        pivot_value = work[column][column] % prime
        determinant = determinant * pivot_value % prime
        inverse = pow(pivot_value, -1, prime)
        for row in range(column + 1, n):
            multiplier = work[row][column] * inverse % prime
            for j in range(column, n):
                work[row][j] = (work[row][j] - multiplier * work[column][j]) % prime
    return determinant % prime


def vandermonde(parameters, prime: int) -> int:
    value = 1
    for i in range(len(parameters)):
        for j in range(i + 1, len(parameters)):
            value = value * (parameters[j] - parameters[i]) % prime
    return value


def multiplicative_order(value: int, prime: int):
    current = 1
    for exponent in range(1, prime):
        current = current * value % prime
        if current == 1:
            return exponent
    raise RuntimeError((value, prime, "order search exhausted"))


def main():
    print("THM-4035 SIXTY-CLOCK AND FINITE KAKEYA AUDIT")

    # AP owner-selector phase law.
    laws, constants, signatures, moments = aggregate_ap_laws()
    require(constants == {Q(127, 35)}, "THM-4029 asymptotic constant changed")
    require(len(laws) == len(set(signatures)) == 60, "AP phase laws are not distinct")
    require(cyclic_periods(signatures) == (60,), "AP law has a proper phase period")
    first_moments = tuple(value[0] for value in moments)
    second_moments = tuple(value[1] for value in moments)
    require(len(set(first_moments)) == 25, "first phase-moment count changed")
    require(len(set(second_moments)) == 42, "second phase-moment count changed")
    require(len(set(moments)) == 60, "two phase moments do not recover phase")
    print(
        "ap_tail=(owners:12,owner_denominators:1..6,phase_period:60,"
        f"distinct_laws:{len(set(signatures))},law_sha256:{digest(signatures)})"
    )
    print(
        "ap_phase_fingerprint=(A1_distinct:25,A2_distinct:42,"
        f"pair_distinct:{len(set(moments))},pair_sha256:{digest(moments)})"
    )

    # Fibonacci clock modulo 10, independently through recurrence and matrix order.
    transition = ((0, 1), (1, 1))
    fib10 = fibonacci_modulus(10, 121)
    fib_states10 = tuple((fib10[r], fib10[r + 1]) for r in range(60))
    require(matrix_order2(transition, 2) == 3, "pi(2)")
    require(matpow2(transition, 5, 5) == ((3, 0), (0, 3)), "A^5 mod 5")
    require(matrix_order2(transition, 5) == 20, "pi(5)")
    require(matrix_order2(transition, 10) == 60, "pi(10)")
    require(len(set(fib_states10)) == 60, "Fibonacci state is not an injective clock")
    require(cyclic_periods(tuple(fib10[:60])) == (60,), "Fibonacci scalar period")
    fib_scalar_fibres = fibres(fib10[:60])
    require(sorted(map(len, fib_scalar_fibres.values())) == [4] * 5 + [8] * 5, "Fibonacci fibres")
    print(
        "fibonacci_mod10=(pi2:3,pi5:20,pi10:60,state_images:60,"
        f"scalar_fibre_sizes:{tuple(sorted(map(len, fib_scalar_fibres.values())))})"
    )

    # Triangular clock.  The finite audit checks the all-modulus formula on a
    # hostile bank; the theorem file supplies the two-line congruence proof.
    for modulus in range(1, 121):
        predicted = triangular_period_formula(modulus)
        values = tuple(triangular(n) % modulus for n in range(3 * predicted + 1))
        require(first_period(values, 2 * modulus) == predicted, (modulus, "triangular period"))
    tri30 = tuple(triangular(r) % 30 for r in range(121))
    tri_states30 = tuple((tri30[r], tri30[r + 1]) for r in range(60))
    require(cyclic_periods(tri30[:60]) == (60,), "triangular mod-30 period")
    require(len(set(tri_states30)) == 60, "two-term triangular state is lossy")
    tri_scalar_fibres = fibres(tri30[:60])
    require(len(tri_scalar_fibres) == 12, "triangular residue image count")
    require(triangular_period_formula(10) == 20, "triangular last-digit hostile")
    require(triangular_period_formula(60) == 120, "triangular mod-60 hostile")
    print(
        "triangular_period_theorem=(odd_M:M,even_M:2M,audit_M:1..120,"
        "tau10:20,tau30:60,tau60:120)"
    )
    print(
        f"triangular_mod30=(state_images:{len(set(tri_states30))},"
        f"scalar_images:{len(tri_scalar_fibres)},fibre_sizes:{tuple(sorted(map(len, tri_scalar_fibres.values())))})"
    )

    # Fibre test: being periodic with period 60 is weaker than determining the
    # AP owner law.  The combined scalars miss exactly one bit.
    combined = tuple((fib10[r], tri30[r]) for r in range(60))
    combined_parity = tuple((fib10[r], tri30[r], r % 2) for r in range(60))
    combined_fibres = fibres(combined)
    collisions = tuple(
        sorted((value, phases) for value, phases in combined_fibres.items() if len(phases) > 1)
    )
    require(len(combined_fibres) == 48, "combined scalar image count")
    require(len(collisions) == 12 and all(len(phases) == 2 for _, phases in collisions), "combined collisions")
    require(all((phases[0] - phases[1]) % 2 for _, phases in collisions), "parity does not split collision")
    require(len(set(combined_parity)) == 60, "one-bit parity repair failed")
    for values in (fib10[:60], tri30[:60], combined):
        require(
            any(len({signatures[r] for r in phases}) > 1 for phases in fibres(values).values()),
            "a purported lossy observable unexpectedly determines the AP law",
        )
    require(combined[0] == combined[15] == (0, 0), "hostile scalar collision moved")
    require(signatures[0] != signatures[15], "hostile AP laws collapsed")
    require(moments[0][0] == Q(37, 14) and moments[15][0] == Q(31, 14), "hostile A1 values")
    print(f"combined_scalars=(images:48,doubletons:{collisions})")
    print("parity_repair=(images:60,one_bit_sufficient:True)")
    print(
        "hostile_phase_pair=(phases:(0,15),F_mod10:(0,0),T_mod30:(0,0),"
        f"A1:({moments[0][0]},{moments[15][0]}),laws_equal:False)"
    )

    # The p=61 finite-field bridge.  sqrt(5)=26 and phi=44 give a Fibonacci
    # eigenvalue and a primitive generator of F_61^*.  Three independent clocks
    # chart the nonzero affine torus in P^3; one clock is only a curve.
    prime = 61
    sqrt_five = 26
    phi = (1 + sqrt_five) * pow(2, -1, prime) % prime
    require(sqrt_five * sqrt_five % prime == 5 and phi == 44, "golden root mod 61")
    require((phi * phi - phi - 1) % prime == 0, "golden polynomial mod 61")
    require(multiplicative_order(phi, prime) == 60, "golden root is not primitive")
    require(matrix_order2(transition, prime) == 60, "pi(61)")
    require(matpow2(transition, 15, prime) == ((11, 0), (0, 11)), "projective Fibonacci collapse")
    fib61 = fibonacci_modulus(prime, 61)
    fib_projective = tuple(projective_pair((fib61[r], fib61[r + 1]), prime) for r in range(60))
    require(len(set(fib_projective)) == 15, "Fibonacci projective orbit should have 15 points")

    parameters = tuple(pow(phi, phase, prime) for phase in range(60))
    require(len(set(parameters)) == 60, "phase parameters collide")
    require(
        all(parameters[r] == (fib61[r + 1] + 43 * fib61[r]) % prime for r in range(60)),
        "Fibonacci state does not recover the primitive phase parameter",
    )
    triangular_shadow = tuple(pow(phi, 2 * triangular(r), prime) for r in range(60))
    require(cyclic_periods(triangular_shadow) == (60,), "triangular multiplicative shadow period")
    require(len(set(triangular_shadow)) == 12, "triangular multiplicative shadow image")

    # The same clock charts nonzero input residues for the Sun 2-4-6-8
    # binomial roles, but the output sumset has already saturated at roles
    # two and four.  This is a hostile against treating p=61 as an obstruction.
    sun_images = tuple(
        frozenset(comb(top, degree) % prime for top in range(prime))
        for degree in (2, 4, 6, 8)
    )
    sun_partial_sumsets = []
    running_sumset = {0}
    for image in sun_images:
        running_sumset = {(left + right) % prime for left in running_sumset for right in image}
        sun_partial_sumsets.append(frozenset(running_sumset))
    sun_target = 896_315_812_331_399
    require(tuple(map(len, sun_images)) == (31, 24, 26, 24), "Sun role image sizes mod 61")
    require(tuple(map(len, sun_partial_sumsets)) == (31, 61, 61, 61), "Sun sumset saturation mod 61")
    require(sun_target % prime == 21, "Sun target residue mod 61")
    broad = tuple((1, t, t * t % prime, t * t * t % prime) for t in parameters)
    narrow = tuple((1, t, 1, t) for t in parameters)
    determinant_digest = sha256()
    broad_nonzero = 0
    narrow_zero = 0
    for phases in combinations(range(60), 4):
        ts = tuple(parameters[r] for r in phases)
        direct = det_mod(tuple(broad[r] for r in phases), prime)
        formula = vandermonde(ts, prime)
        require(direct == formula, (phases, "direct/Vandermonde disagreement"))
        require(direct != 0, (phases, "broad quartet is not transverse"))
        broad_nonzero += 1
        narrow_det = det_mod(tuple(narrow[r] for r in phases), prime)
        require(narrow_det == 0, (phases, "narrow control left its plane"))
        narrow_zero += 1
        determinant_digest.update(f"{phases}:{direct}\n".encode("ascii"))
    quartet_count = comb(60, 4)
    require(broad_nonzero == narrow_zero == quartet_count == 487635, "quartet universe")

    total_projective_directions = prime**3 + prime**2 + prime + 1
    torus_directions = (prime - 1) ** 3
    require(total_projective_directions == 230764 and torus_directions == 216000, "direction counts")
    broad_control = det_mod(tuple(broad[r] for r in (0, 1, 2, 3)), prime)
    narrow_control = det_mod(tuple(narrow[r] for r in (0, 1, 2, 3)), prime)
    require(broad_control != 0 and narrow_control == 0, "same-clock hostile control")
    print(
        "finite_field_clock=(p:61,sqrt5:26,phi:44,ord_phi:60,pi61:60,"
        f"fibonacci_projective_images:{len(set(fib_projective))},"
        "phase_from_fibonacci:F_(r+1)+43F_r)"
    )
    print("triangular_multiplicative_shadow=(period:60,images:12,formula:phi^(2T_r))")
    print(
        "sun_2468_mod61_hostile=(target_residue:21,role_image_sizes:(31,24,26,24),"
        "partial_sumset_sizes:(31,61,61,61),local_obstruction:False)"
    )
    print(
        "projective_direction_atlas=(P3_directions:230764,nonzero_torus:216000,"
        "torus_coordinates:3_independent_C60_clocks)"
    )
    print(
        "twisted_cubic_broad_control=(directions:60,quartets:487635,"
        f"nonzero_determinants:{broad_nonzero},det_sha256:{determinant_digest.hexdigest()})"
    )
    print(
        "same_clock_narrow_control=(directions:60,quartets:487635,"
        f"zero_determinants:{narrow_zero},phase0123_dets:({broad_control},{narrow_control}))"
    )
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
