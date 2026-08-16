#!/usr/bin/env python3
"""Exact companion for the proposed THM-3483 nondivisor residue compiler.

For p not dividing d this checks the one-digit residual formula for every
coefficient of the two resonant factorial rows on a bounded hostile universe,
then uses the proved formula (without constructing the large coefficients) to
certify five pair-only closures left by the THM-3475 divisor compiler.
"""

from fractions import Fraction
from math import comb
import hashlib
import importlib.util
import json
from pathlib import Path


COEFFICIENT_NAME = "factorial_multiplace_newton_degree_barcode_thm3152.py"
COEFFICIENT_SHA256 = (
    "f804d3996abe4df981dbf7db877af4aeca9218df64b0ac382af876a3cdca15a0"
)
DIVISOR_NAME = (
    "factorial_all_divisor_digit_pair_compiler_independent_audit_thm3475.py"
)
DIVISOR_SHA256 = (
    "9330ca1b991b9d5875779b9975fc88701ab36855a6527e1865e821e6cd3ea665"
)

CONTROL_D = tuple(range(4, 61))
CONTROL_PRIMES = (2, 3, 5, 7, 11, 13, 17, 19)
EXPECTED_CONTROL_COEFFICIENTS = 23395
EXPECTED_CONTROL_PROFILES = 746
EXPECTED_ADMISSIBLE_PROFILES = 665

EXPECTED_PACKETS = {
    2516: (503, 1006, 1509, 2012),
    2564: (466, 699, 1165, 1631, 1864, 2097, 2330),
    2571: (2056,),
    2576: tuple(range(103, 2473, 103)),
    2586: (
        47, 141, 188, 235, 282, 329,
        2209, 2256, 2303, 2350, 2397, 2444, 2491, 2538,
    ),
}

EXPECTED_CLOSURES = {
    (2516, 7): (
        ((0, 0), (1, 0), (15, 4), (113, 36), (2514, 836)),
        ((0, 0), (2, 0), (16, 4), (114, 36), (2515, 836)),
        (("0", 1, 1), ("2/7", 14, 7), ("16/49", 98, 49),
         ("800/2401", 2401, 2401)),
        36,
    ),
    (2564, 13): (
        ((0, 0), (1, 0), (27, 4), (365, 60), (2562, 426)),
        ((0, 0), (2, 0), (28, 4), (366, 60), (2563, 426)),
        (("0", 1, 1), ("2/13", 26, 13), ("28/169", 338, 169),
         ("366/2197", 2197, 2197)),
        36,
    ),
    (2571, 7): (
        ((0, 0), (21, 6), (168, 54), (2569, 854)),
        ((0, 0), (1, 0), (22, 6), (169, 54), (2570, 854)),
        (("2/7", 21, 7), ("16/49", 147, 49),
         ("800/2401", 2401, 2401)),
        32,
    ),
    (2576, 13): (
        ((0, 0), (39, 6), (377, 62), (2574, 428)),
        ((0, 0), (1, 0), (40, 6), (378, 62), (2575, 428)),
        (("2/13", 39, 13), ("28/169", 338, 169),
         ("366/2197", 2197, 2197)),
        24,
    ),
    (2586, 7): (
        ((0, 0), (1, 0), (22, 6), (169, 54), (2570, 854),
         (2584, 860)),
        ((0, 0), (2, 0), (23, 6), (170, 54), (2571, 854),
         (2585, 860)),
        (("0", 1, 1), ("2/7", 21, 7), ("16/49", 147, 49),
         ("800/2401", 2401, 2401), ("3/7", 14, 7)),
        96,
    ),
}

# Filled only after the complete semantic record is generated once.  It pins
# serialization as JSON with sorted keys and compact separators.
EXPECTED_SEMANTIC_SHA256 = (
    "f80c046942d62a8a6b6f3802d224cc47944568a7e9f3ef245d43baf91a4031c4"
)


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def digest(value):
    encoded = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def load_pinned(filename, expected_hash, module_name):
    source = Path(__file__).resolve().with_name(filename)
    require(source.is_file(), ("missing dependency", source))
    normalized = source.read_bytes().replace(b"\r\n", b"\n")
    actual_hash = hashlib.sha256(normalized).hexdigest()
    require(actual_hash == expected_hash, (source, actual_hash, expected_hash))
    spec = importlib.util.spec_from_file_location(module_name, source)
    require(spec is not None and spec.loader is not None, ("bad import", source))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def factorial_valuation(n, prime):
    value = 0
    while n:
        n //= prime
        value += n
    return value


def binomial_valuation(n, j, prime):
    return (
        factorial_valuation(n, prime)
        - factorial_valuation(j, prime)
        - factorial_valuation(n - j, prime)
    )


def raw_weight(n, j, prime):
    return binomial_valuation(n, j, prime) + factorial_valuation(2 * j, prime)


def rho(n, j, d, prime):
    """Return d^(-(n-j)) Z_(n,j) modulo a nondivisor prime."""
    require(d % prime != 0, ("rho requires p not dividing d", d, prime))
    m0 = (n - j) % prime
    j0 = j % prime
    inverse = pow(d % prime, -1, prime)
    total = 0
    rising = 1
    for ell in range(m0 + 1):
        if ell:
            rising = rising * (2 * j0 + ell) % prime
        term = comb(m0, ell) * pow(-inverse, ell, prime) * rising
        total = (total + term) % prime
    return total


def lower_hull(points):
    hull = []
    for point in points:
        while len(hull) >= 2:
            x0, y0 = hull[-2]
            x1, y1 = hull[-1]
            x2, y2 = point
            left = (y1 - y0) * (x2 - x1)
            right = (y2 - y1) * (x1 - x0)
            if left < right:
                break
            hull.pop()
        hull.append(point)
    require(len(hull) >= 2, ("short hull", points))
    return tuple(hull)


def raw_hull(n, prime):
    return lower_hull(tuple((j, raw_weight(n, j, prime)) for j in range(n + 1)))


def actual_hull(poly, prime, coefficient_module):
    points = []
    for j in range(poly.degree() + 1):
        height = coefficient_module.vp(poly[j], prime)
        if height is not None:
            points.append((j, height))
    return lower_hull(tuple(points))


def ledger(hull):
    answer = {}
    for (x0, y0), (x1, y1) in zip(hull, hull[1:]):
        slope = Fraction(y1 - y0, x1 - x0)
        require(slope not in answer, ("repeated slope", slope, hull))
        answer[slope] = x1 - x0
    return answer


def pair_degrees(f_hull, g_hull):
    f_ledger = ledger(f_hull)
    g_ledger = ledger(g_hull)
    degrees = {0}
    blocks = []
    for slope in sorted(set(f_ledger) & set(g_ledger)):
        capacity = min(f_ledger[slope], g_ledger[slope])
        denominator = slope.denominator
        require(capacity % denominator == 0,
                ("denominator does not divide capacity", slope, capacity))
        blocks.append((str(slope), capacity, denominator))
        degrees = {
            old + use
            for old in degrees
            for use in range(0, capacity + 1, denominator)
        }
    return tuple(sorted(degrees)), tuple(blocks)


def control_census(coefficient_module):
    coefficient_count = 0
    profile_count = 0
    admissible_count = 0
    records = []
    hostile = None
    for d in CONTROL_D:
        f_poly, g_poly = coefficient_module.pair(d)
        for prime in CONTROL_PRIMES:
            if d % prime == 0:
                continue
            for label, n, poly in (
                    ("F", d - 2, f_poly), ("G", d - 1, g_poly)):
                profile_count += 1
                raw = raw_hull(n, prime)
                actual = actual_hull(poly, prime, coefficient_module)
                vertex_residues = tuple(
                    (j, rho(n, j, d, prime)) for j, _ in raw
                )
                admissible = all(value != 0 for _, value in vertex_residues)
                for j in range(n + 1):
                    coefficient_count += 1
                    raw_height = raw_weight(n, j, prime)
                    actual_height = coefficient_module.vp(poly[j], prime)
                    residue = rho(n, j, d, prime)
                    require(
                        actual_height is None or actual_height >= raw_height,
                        ("raw lower bound", d, prime, label, j,
                         raw_height, actual_height),
                    )
                    if residue:
                        require(actual_height == raw_height,
                                ("rho unit exactness", d, prime, label, j))
                    else:
                        require(
                            actual_height is None or actual_height > raw_height,
                            ("rho zero must lift", d, prime, label, j),
                        )
                if admissible:
                    admissible_count += 1
                    require(actual == raw,
                            ("admissible raw/actual hull", d, prime, label))
                record = (
                    d, prime, label, admissible, raw, vertex_residues, actual,
                )
                records.append(record)
                if (d, prime, label) == (4, 5, "F"):
                    hostile = record

    require(coefficient_count == EXPECTED_CONTROL_COEFFICIENTS,
            ("coefficient count", coefficient_count))
    require(profile_count == EXPECTED_CONTROL_PROFILES,
            ("profile count", profile_count))
    require(admissible_count == EXPECTED_ADMISSIBLE_PROFILES,
            ("admissible count", admissible_count))
    require(hostile is not None, "missing hostile")
    require(hostile[4] == ((0, 0), (2, 0)), ("hostile raw", hostile))
    require(hostile[5] == ((0, 0), (2, 1)), ("hostile rho", hostile))
    require(hostile[6] == ((0, 1), (1, 0), (2, 0)),
            ("hostile actual", hostile))
    return tuple(records), hostile, coefficient_count, profile_count, admissible_count


def audit_binary_rule():
    record = []
    for d in range(3, 20, 2):
        for n in range(0, 40):
            for j in range(n + 1):
                value = rho(n, j, d, 2)
                predicted = 1 if (n - j) % 2 == 0 else 0
                require(value == predicted, ("binary parity", d, n, j, value))
        record.append((d, 40 * 41 // 2))
    return tuple(record)


def closure_census(divisor_module):
    require(divisor_module.EXPECTED_SURVIVORS == EXPECTED_PACKETS,
            ("THM-3475 packet mismatch", divisor_module.EXPECTED_SURVIVORS))
    records = []
    for (d, prime), expected in EXPECTED_CLOSURES.items():
        n = d - 1
        f_hull = raw_hull(n - 1, prime)
        g_hull = raw_hull(n, prime)
        f_rhos = tuple((j, rho(n - 1, j, d, prime)) for j, _ in f_hull)
        g_rhos = tuple((j, rho(n, j, d, prime)) for j, _ in g_hull)
        require(all(value != 0 for _, value in f_rhos),
                ("F inadmissible", d, prime, f_rhos))
        require(all(value != 0 for _, value in g_rhos),
                ("G inadmissible", d, prime, g_rhos))
        local, blocks = pair_degrees(f_hull, g_hull)
        packet = EXPECTED_PACKETS[d]
        intersection = tuple(sorted(set(local) & set(packet)))
        expected_f, expected_g, expected_blocks, expected_count = expected
        require(f_hull == expected_f, ("F hull", d, prime, f_hull))
        require(g_hull == expected_g, ("G hull", d, prime, g_hull))
        require(blocks == expected_blocks, ("blocks", d, prime, blocks))
        require(len(local) == expected_count, ("local count", d, prime, len(local)))
        require(not intersection, ("surviving packet", d, prime, intersection))
        # Both constant raw vertices are rho-units, so both exact constant
        # coefficients are nonzero.  The separate coordinate-root capacity is 0.
        require(f_rhos[0][0] == 0 and f_rhos[0][1] != 0,
                ("F constant", d, prime, f_rhos[0]))
        require(g_rhos[0][0] == 0 and g_rhos[0][1] != 0,
                ("G constant", d, prime, g_rhos[0]))
        coordinate_root_capacity = 0
        records.append((
            d, prime, f_hull, g_hull, f_rhos, g_rhos, blocks,
            packet, intersection, len(local), coordinate_root_capacity,
        ))
    return tuple(records)


def main():
    coefficient_module = load_pinned(
        COEFFICIENT_NAME, COEFFICIENT_SHA256, "thm3483_coefficient",
    )
    divisor_module = load_pinned(
        DIVISOR_NAME, DIVISOR_SHA256, "thm3483_divisor",
    )
    controls, hostile, coefficient_count, profile_count, admissible_count = (
        control_census(coefficient_module)
    )
    binary = audit_binary_rule()
    closures = closure_census(divisor_module)
    semantic_record = (controls, hostile, binary, closures)
    semantic_sha256 = digest(semantic_record)
    require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256, semantic_sha256)

    print("THM-3483 NONDIVISOR RESIDUE DIGIT PAIR COMPILER")
    print("dependencies=((%s,%s),(%s,%s))" % (
        COEFFICIENT_NAME, COEFFICIENT_SHA256, DIVISOR_NAME, DIVISOR_SHA256,
    ))
    print("formula=rho=sum_(ell=0)^(m mod p) C(m mod p,ell)(-d^-1)^ell(2(j mod p)+1)_ell mod p")
    print("control_universe=d=4..60; p=(2,3,5,7,11,13,17,19), p does not divide d; rows F=A_(d-2)^d,G=A_(d-1)^d")
    print("control_counts=(coefficients=%d,profiles=%d,rho_admissible_profiles=%d)" % (
        coefficient_count, profile_count, admissible_count,
    ))
    print("hostile=(d=4,p=5,F; raw=((0,0),(2,0)); rho=((0,0),(2,1)); actual=((0,1),(1,0),(2,0)))")
    print("binary_rule=rho=1 iff m=n-j is even; checked d odd 3..19 and all 0<=j<=n<40")
    print("closure_schema=(d,p,F_raw_hull,G_raw_hull,F_vertex_rhos,G_vertex_rhos,common_blocks,THM3475_packet,intersection,local_count,coordinate_root_capacity)")
    print("closures=%s" % (closures,))
    print("qualification=common slope zero is a unit-root block and is retained; coordinate-root capacity is separate and equals zero in all five rows")
    print("semantic_schema=(all_746_control_profiles,hostile,binary_checks,five_closure_records); json sort_keys compact separators")
    print("semantic_sha256=%s" % semantic_sha256)
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
