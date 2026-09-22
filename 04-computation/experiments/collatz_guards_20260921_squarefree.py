"""Exact certificates for squarefree Collatz guards and the 101 resonance.

The universal sieve and affine-cycle proofs are in the companion note.
Every computational universe is explicit; no floating-point tests are used.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def primes_up_to(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            for multiple in range(p * p, limit + 1, p):
                sieve[multiple] = 0
    return [p for p in range(2, limit + 1) if sieve[p]]


def squarefree_factor_test(n: int) -> bool:
    p = 2
    while p * p <= n:
        if n % p == 0:
            n //= p
            if n % p == 0:
                return False
        p += 1
    return True


def squarefree_trial_squares(n: int) -> bool:
    return all(n % (d * d) != 0 for d in range(2, isqrt(n) + 1))


def step(n: int) -> tuple[int, int]:
    require(n > 0 and n % 2 == 1, "positive odd source required")
    value = 3 * n + 1
    exponent = 0
    while value % 2 == 0:
        value //= 2
        exponent += 1
    return value, exponent


def compose(left: tuple[int, int], right: tuple[int, int], modulus: int) -> tuple[int, int]:
    a, b = left
    c, d = right
    return a * c % modulus, (a * d + b) % modulus


def affine_power(a: int, b: int, exponent: int, modulus: int) -> tuple[int, int]:
    result = (1, 0)
    current = (a, b)
    while exponent:
        if exponent % 2:
            result = compose(current, result, modulus)
        current = compose(current, current, modulus)
        exponent //= 2
    return result


def resonance_certificate() -> dict:
    modulus = 101**2
    inv128 = pow(128, -1, modulus)
    a, b = 27 * inv128 % modulus, 19 * inv128 % modulus
    orbit = []
    indices = {}
    value = 0
    while value not in indices:
        indices[value] = len(orbit)
        orbit.append(value)
        value = (a * value + b) % modulus
    require(value == 0 and len(orbit) == modulus, "101-square full cycle failure")
    roots = [0, -pow(3, -1, modulus) % modulus,
             -5 * pow(9, -1, modulus) % modulus]
    anchors = [indices[root] for root in roots]
    require(roots == [0, 3400, 9067] and anchors == [0, 3795, 5214], "anchor failure")
    binary_certificates = []
    for root, exponent in zip(roots, anchors):
        powered = affine_power(a, b, exponent, modulus)
        require(powered[1] == root, "independent binary affine power failed")
        binary_certificates.append({"exponent": exponent, "slope": powered[0],
                                    "intercept": powered[1], "root": root})

    inv27 = pow(27, -1, modulus)

    def inverse_reset(x: int) -> int:
        return (128 * x - 19) * inv27 % modulus

    bad = {0}
    phase0, phase1, phase2 = roots
    counts = [modulus - 1]
    boundary = {}
    for repeats in range(1, 5101):
        phase0 = inverse_reset(phase0)
        bad.update([phase0, phase1, phase2])
        phase1 = inverse_reset(phase1)
        phase2 = inverse_reset(phase2)
        count = modulus - len(bad)
        formula = sum(max(gap - repeats, 0) for gap in [3795, 1419, 4986])
        require(count == formula, "all-length survivor formula mismatch")
        counts.append(count)
        if repeats in [1, 1419, 3795, 4985, 4986, 5100]:
            survivors = sorted(set(range(modulus)) - bad)
            boundary[str(repeats)] = {"count": count,
                                      "survivors_if_at_most_ten": survivors if count <= 10 else None}
    require(boundary["4985"]["survivors_if_at_most_ten"] == [7332], "last survivor failure")
    require(boundary["4986"]["count"] == 0, "cutoff failure")

    # Direct phase-by-phase modular evolution is independent of root unions.
    value = 7332
    first_zero = None
    for j, exponent in enumerate([1, 1, 5] * 4986, start=1):
        value = (3 * value + 1) * pow(2**exponent, -1, modulus) % modulus
        if value == 0 and first_zero is None:
            first_zero = j
    require(first_zero == 3 * 4986, "final endpoint off-by-one failure")
    return {"prime": 101, "modulus": modulus, "orbit_from_zero": orbit,
            "cycle_length": len(orbit), "phase_roots": roots, "anchor_indices": anchors,
            "binary_affine_certificates": binary_certificates,
            "survivor_counts_r_0_through_5100": counts, "boundary_cases": boundary,
            "sole_last_survivor_first_zero_step": first_zero}


def other_prime_controls() -> dict:
    rows = []
    for p in primes_up_to(251):
        if p < 5 or p == 101:
            continue
        modulus = p * p
        fixed = 19 * pow(101, -1, modulus) % modulus
        source = (fixed + 13) % modulus if p == 13 else fixed
        value = source
        for _ in range(4986):
            phases = [value, (3 * value + 1) * pow(2, -1, modulus) % modulus,
                      (9 * value + 5) * pow(4, -1, modulus) % modulus]
            require(all(phase != 0 for phase in phases), "nonresonant avoiding seed failure")
            if p == 13:
                require(phases[0] % 13 != 0 and phases[1] % 13 != 0,
                        "13 phase units failure")
                require(phases[2] % 13 == 0, "13 exact first valuation failure")
            value = (27 * value + 19) * pow(128, -1, modulus) % modulus
        require(value != 0, "nonresonant final endpoint failure")
        rows.append({"p": p, "modulus": modulus, "fixed_point": fixed,
                     "avoiding_source": source, "reset_blocks_checked": 4986})
    value = 50
    orbit13 = []
    while value not in orbit13:
        orbit13.append(value)
        value = (27 * value + 19) * pow(128, -1, 169) % 169
    require(value == 50 and len(orbit13) == 12, "13 repair period failure")
    return {"prime_universe": "5<=p<=251, p prime, p!=101", "rows": rows,
            "prime_13_repaired_orbit": orbit13}


def word_source_progression(repeats: int) -> tuple[int, int]:
    denominator = 128**repeats
    tripling = 27**repeats
    carry = 19 * (denominator - tripling) // 101
    binary_modulus = 2 * denominator
    residue = (denominator - carry) * pow(tripling, -1, binary_modulus) % binary_modulus
    source = next(residue + binary_modulus * t for t in range(9)
                  if (residue + binary_modulus * t) % 9 == 5)
    return source, 9 * binary_modulus


def actual_reset_witnesses() -> list[dict]:
    result = []
    for repeats in [1, 2, 3]:
        base, modulus = word_source_progression(repeats)
        witness = None
        for parameter in range(50):
            source = base + modulus * parameter
            value = source
            nodes = [value]
            for expected in [1, 1, 5] * repeats:
                value, exponent = step(value)
                require(exponent == expected, "word source cylinder failure")
                nodes.append(value)
            if all(squarefree_factor_test(n) for n in nodes):
                witness = parameter, nodes
                break
        require(witness is not None, "bounded small reset witness search failed")
        parameter, nodes = witness
        require(all(squarefree_trial_squares(n) for n in nodes), "independent square test failed")
        require(nodes[0] % 9 == 5 and nodes[-1] % 9 == 5, "row five reset failure")
        result.append({"r": repeats, "base": base, "source_modulus": modulus,
                       "first_success_parameter": parameter, "nodes": nodes,
                       "search_parameter_range": [0, 49]})
    return result


def growing_row_five_witnesses() -> list[dict]:
    result = []
    for length in [1, 2, 3, 6, 10]:
        modulus = 2 ** (length + 1)
        m = max(16, 2 * length + 4)
        q_residue = 6 * pow(modulus, -1, 9) % 9
        first_q = (2**m + 1 + modulus - 1) // modulus
        last_q = ((m + 1) * 2**m + m) // (m * modulus)
        first_compatible = first_q + (q_residue - first_q) % 9
        witness = None
        checked = 0
        for q in range(first_compatible, last_q + 1, 9):
            checked += 1
            nodes = [3**j * 2 ** (length + 1 - j) * q - 1 for j in range(length + 1)]
            if all(squarefree_factor_test(n) for n in nodes):
                witness = q, nodes
                break
        require(witness is not None, "row five growing witness missing")
        q, nodes = witness
        require(all(squarefree_trial_squares(n) for n in nodes), "row five independent square test failed")
        require(nodes[0] % 9 == 5, "source residue failure")
        for n, target in zip(nodes, nodes[1:]):
            require(step(n) == (target, 1), "growing exponent word failure")
        source = nodes[0]
        error_bound = 512 * (m + 1) * Fraction(source - 2**m, 2**m) ** 4
        require(error_bound < Fraction(1, 3), "energy excess bound failure")
        result.append({"L": length, "m": m, "q": q, "q_mod9": q_residue,
                       "q_interval": [first_q, last_q], "compatible_q_checked": checked,
                       "nodes": nodes,
                       "rational_energy_excess_bound": f"{error_bound.numerator}/{error_bound.denominator}"})
    return result


def prime_sum_control() -> dict:
    primes = primes_up_to(37)[1:]
    cumulative = [1]
    skipped = []
    prime_set = set(primes)
    for index, p in enumerate(primes, start=1):
        cumulative.append(cumulative[-1] + p)
        composites = sum(n not in prime_set for n in range(3, p + 1, 2))
        require(p == 2 * index + 1 + 2 * composites, "prime gap count identity failed")
        skipped.append(composites)
        require(cumulative[-1] == (index + 1) ** 2 + 2 * sum(skipped),
                "prime cumulative correction failed")
    require(cumulative[-1] == 196 == 14**2, "196 identity failed")
    require(len(cumulative) == 12 and sum(skipped) == 26, "196 index count failed")
    return {"odd_primes": primes, "cumulative_including_initial_one": cumulative,
            "odd_composite_counts": skipped, "sum_composite_counts": sum(skipped),
            "extra_seven_would_give": cumulative[-1] + 7,
            "claimed_Collatz_map_from_identity": False}


def constant_exponent_controls() -> list[dict]:
    result = []
    for exponent, prime in [(3, 5), (4, 13)]:
        modulus = prime**2
        inverse_denominator = pow(2**exponent, -1, modulus)
        orbit = []
        value = 0
        while value not in orbit:
            orbit.append(value)
            value = (3 * value + 1) * inverse_denominator % modulus
        require(value == 0 and len(orbit) == modulus, "constant exponent cycle failure")
        # The last surviving residue is F(0): starting there first hits
        # zero after p^2-1 steps, so the first p^2-1 displayed nodes avoid it.
        source = inverse_denominator
        value = source
        first_zero = None
        for j in range(1, modulus):
            value = (3 * value + 1) * inverse_denominator % modulus
            if value == 0 and first_zero is None:
                first_zero = j
        require(first_zero == modulus - 1, "constant exponent endpoint threshold failed")
        result.append({"exponent": exponent, "slope_gap": 2**exponent - 3,
                       "least_gap_prime": prime, "modulus": modulus,
                       "full_orbit_from_zero": orbit,
                       "maximum_squarefree_steps": modulus - 2,
                       "first_impossible_squarefree_steps": modulus - 1,
                       "sole_last_surviving_residue": source})
    return result


def prefix_zero_roots(word: list[int], prime: int) -> list[int]:
    require(prime >= 5 and prime in primes_up_to(prime), "prime >=5 required")
    require(all(exponent >= 1 for exponent in word), "positive exponents required")
    modulus = prime**2
    carry = 0
    tripling = 1
    binary = 1
    roots = [0]
    for exponent in word:
        carry = (3 * carry + binary) % modulus
        tripling = 3 * tripling % modulus
        binary = binary * pow(2, exponent, modulus) % modulus
        roots.append(-carry * pow(tripling, -1, modulus) % modulus)
    return roots


def finite_squarefree_word_criterion(word: list[int]) -> dict:
    require(all(exponent >= 1 for exponent in word), "positive exponents required")
    rows = []
    for prime in primes_up_to(isqrt(len(word) + 1)):
        if prime < 5:
            continue
        roots = sorted(set(prefix_zero_roots(word, prime)))
        rows.append({"prime": prime, "modulus": prime**2, "roots": roots,
                     "covers": len(roots) == prime**2})
    return {"word": word, "length": len(word), "node_count": len(word) + 1,
            "finite_prime_checks": rows,
            "positive_density_squarefree_realizations_in_row5": not any(row["covers"] for row in rows)}


def general_word_controls() -> list[dict]:
    result = []
    for word in [[3] * 23, [3] * 24, [1] * 24]:
        roots = sorted(set(prefix_zero_roots(word, 5)))
        direct_bad = []
        for source in range(25):
            value = source
            hits_zero = value == 0
            for exponent in word:
                value = (3 * value + 1) * pow(2**exponent, -1, 25) % 25
                hits_zero = hits_zero or value == 0
            if hits_zero:
                direct_bad.append(source)
        require(roots == direct_bad, "general carry roots disagree with direct source enumeration")
        criterion = finite_squarefree_word_criterion(word)
        expected = not (word == [3] * 24)
        require(criterion["positive_density_squarefree_realizations_in_row5"] == expected,
                "shortest obstruction criterion failed")
        criterion["independent_mod25_root_count"] = len(roots)
        criterion["independent_mod25_roots"] = roots
        result.append(criterion)
    require(result[0]["independent_mod25_root_count"] == 24, "length23 root count failure")
    require(result[1]["independent_mod25_root_count"] == 25, "length24 root count failure")
    return result


def run() -> dict:
    return {"status": "FINITE-EXACT certificates; universal theorems proved in companion note",
            "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
            "resonance_101": resonance_certificate(),
            "other_prime_controls": other_prime_controls(),
            "actual_all_squarefree_reset_witnesses": actual_reset_witnesses(),
            "actual_squarefree_growth_energy_row5_witnesses": growing_row_five_witnesses(),
            "constant_exponent_controls": constant_exponent_controls(),
            "general_finite_word_controls": general_word_controls(),
            "prime_sum": prime_sum_control(),
            "checks": "PASS; explicit checks remain active under python -O",
            "limitations": ["No Collatz convergence proof", "No infinite squarefree orbit claimed",
                            "No squarefree fixed-target basin completion claimed",
                            "The threshold includes all intermediate nodes and the final endpoint",
                            "Sieve limits are for fixed word length only"]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_suffix(".json"))
    args = parser.parse_args()
    output = run()
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8", newline="\n")
    print(json.dumps({"checks": output["checks"], "output": str(args.output),
                      "reset_threshold": 4986, "modular_orbit_length": 10201,
                      "survivor_count_cases": 5101}, sort_keys=True))
