#!/usr/bin/env python3
"""Pinned-engine exact factorial pair census for d=4001..6000."""

import ast
from collections import Counter
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path


START, END = 4001, 6000
PRIMES_47 = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)
PRIMES_101 = PRIMES_47 + (53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101)
DIVISOR_NAME = "factorial_all_divisor_digit_pair_compiler_thm3475.py"
DIVISOR_SHA256 = "834d0913eb5cd5b15684c7fb88af60e42d2a6ef36feb821e261c0498f55027ab"
RHO_NAME = "factorial_nondivisor_residue_digit_pair_compiler_thm3483.py"
RHO_SHA256 = "9e37ead620f141617a9c6d51c182e09c034945793092e56e39fb061254662723"
INDEPENDENT_NAME = "factorial_adaptive_rho_block_6000_independent_audit.py"
INDEPENDENT_SHA256 = "d416cb2955fd745394cf1043ac8c2eba28a6a97beb264dd9cbe9919ed8c96724"
EXPECTED_COUNTS = (2000, 1272, 728, 600, 128, 0)
EXPECTED_RHO_HISTOGRAM = (
    (3, 1), (5, 16), (7, 46), (11, 42), (13, 10), (17, 9), (19, 4),
)
EXPECTED_SEMANTIC_SHA256 = "7f8ab74ae9fae027f32fd7eabaf0338c217319e274594bd603859a1bbcca28bd"


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def digest(value):
    data = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return hashlib.sha256(data).hexdigest()


def load_pinned(name, expected, module_name):
    path = Path(__file__).resolve().with_name(name)
    actual = hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    require(actual == expected, (name, actual, expected))
    spec = importlib.util.spec_from_file_location(module_name, path)
    require(spec is not None and spec.loader is not None, ("bad import", path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_engines():
    divisor = load_pinned(DIVISOR_NAME, DIVISOR_SHA256, "rho6000_divisor")
    rho_module = load_pinned(RHO_NAME, RHO_SHA256, "rho6000_rho")
    independent = load_pinned(INDEPENDENT_NAME, INDEPENDENT_SHA256, "rho6000_independent")
    coefficient = divisor.load_pinned(divisor.PRIMARY_NAME, divisor.PRIMARY_SHA256,
                                      "rho6000_coefficient")
    digital = divisor.load_pinned(divisor.DIGITAL_NAME, divisor.DIGITAL_SHA256,
                                  "rho6000_digital")
    return divisor, rho_module, independent, coefficient, digital


def factor_primes(coefficient, n):
    return tuple(p for p, _ in coefficient.factor(n))


def is_prime(coefficient, n):
    factors = coefficient.factor(n) if n >= 2 else ()
    return len(factors) == 1 and factors[0][1] == 1


def exit_label(coefficient, d):
    tests = (
        ("d_prime", is_prime(coefficient, d)),
        ("d_minus_1_prime_power", len(coefficient.factor(d - 1)) == 1),
        ("d_minus_2_prime", is_prime(coefficient, d - 2)),
        ("d_minus_3_prime", is_prime(coefficient, d - 3)),
        ("d_minus_4_prime", is_prime(coefficient, d - 4)),
        ("d_minus_5_prime", is_prime(coefficient, d - 5)),
        ("d_minus_6_prime", is_prime(coefficient, d - 6)),
    )
    return next((label for label, true in tests if true), None)


def divisor_barcode(divisor, coefficient, digital, n, p):
    _, _, blocks, degrees = divisor.predicted_pair_compiler(coefficient, digital, n, p)
    normalized = tuple((s.numerator, s.denominator, capacity)
                       for s, capacity, _ in blocks)
    return tuple(sorted(degrees)), normalized


def rho_barcode(rho_module, f_hull, g_hull):
    degrees, blocks = rho_module.pair_degrees(f_hull, g_hull)
    normalized = []
    for slope_text, capacity, denominator in blocks:
        slope = Fraction(slope_text)
        require(slope.denominator == denominator, (slope, denominator))
        normalized.append((slope.numerator, denominator, capacity))
    return degrees, tuple(normalized)


def scan_row(divisor, rho_module, coefficient, digital, d):
    n, packet, divisor_trace = d - 1, tuple(range(1, d - 1)), []
    for p in factor_primes(coefficient, n):
        degrees, blocks = divisor_barcode(divisor, coefficient, digital, n, p)
        pre, local = packet, set(degrees)
        packet = tuple(k for k in packet if k in local)
        divisor_trace.append((p, blocks, pre, packet))
        if not packet:
            return (d, "divisor", p, tuple(divisor_trace), (), ())

    divisor_packet, rho_trace = packet, []
    for p in PRIMES_101:
        if d % p == 0:
            rho_trace.append((p, "divides_d", (), (), (), packet, packet))
            continue
        f_hull = rho_module.raw_hull(n - 1, p)
        g_hull = rho_module.raw_hull(n, p)
        f_vertices = tuple((j, rho_module.rho(n - 1, j, d, p)) for j, _ in f_hull)
        g_vertices = tuple((j, rho_module.rho(n, j, d, p)) for j, _ in g_hull)
        pre, blocks, status = packet, (), "inadmissible"
        if all(value for _, value in f_vertices + g_vertices):
            degrees, blocks = rho_barcode(rho_module, f_hull, g_hull)
            packet = tuple(k for k in packet if k in set(degrees))
            status = "admissible"
        rho_trace.append((p, status, f_vertices, g_vertices, blocks, pre, packet))
        if not packet:
            return (d, "rho", p, tuple(divisor_trace), tuple(rho_trace), divisor_packet)
    return (d, "survivor", None, tuple(divisor_trace), tuple(rho_trace),
            divisor_packet, packet)


def build_semantic(divisor, rho_module, coefficient, digital):
    exits, residuals = [], []
    for d in range(START, END + 1):
        label = exit_label(coefficient, d)
        (exits if label else residuals).append(
            (d, label) if label else scan_row(divisor, rho_module, coefficient, digital, d)
        )
    return ((START, END), tuple(exits), tuple(residuals), PRIMES_47, PRIMES_101)


def main():
    divisor, rho_module, independent, coefficient, digital = load_engines()
    semantic = build_semantic(divisor, rho_module, coefficient, digital)
    independent_semantic, _ = independent.build_semantic()
    require(semantic == independent_semantic,
            (digest(semantic), digest(independent_semantic)))
    exits, divisor_rows, rho_rows, survivors, packets, killers, inadmissible = independent.summarize(semantic)
    counts = (END - START + 1, len(exits), len(semantic[2]), len(divisor_rows),
              len(rho_rows), len(survivors))
    exit_histogram = tuple(sorted(Counter(label for _, label in exits).items()))
    killer_histogram = tuple(sorted(Counter(p for _, p in killers).items()))
    require(counts == EXPECTED_COUNTS, counts)
    require(exit_histogram == independent.EXPECTED_EXIT_HISTOGRAM, exit_histogram)
    require(killer_histogram == EXPECTED_RHO_HISTOGRAM, killer_histogram)
    require(len(inadmissible) == 80 and not survivors, (len(inadmissible), survivors))
    require(digest(semantic) == EXPECTED_SEMANTIC_SHA256, digest(semantic))
    first_rho = rho_rows[0]
    require(first_rho[:3] == (4150, "rho", 11) and first_rho[5] == (3227,),
            first_rho[:3])
    require(not any(isinstance(node, ast.Assert)
                    for node in ast.walk(ast.parse(Path(__file__).read_text()))),
            "assert node")

    print("FACTORIAL ADAPTIVE RHO BLOCK d=4001..6000 PINNED-ENGINE CENSUS")
    print("dependencies=%s:%s;%s:%s;%s:%s" %
          (DIVISOR_NAME, DIVISOR_SHA256, RHO_NAME, RHO_SHA256,
           INDEPENDENT_NAME, INDEPENDENT_SHA256))
    print("implementation=THM-3475 divisor compiler + THM-3483 rho compiler; exact cross-record equality with self-contained referee")
    print("universe=d=4001..6000; canonical seven exits; F=A_(d-2)^d,G=A_(d-1)^d; rho bank=%s" % (PRIMES_101,))
    print("counts=%s exit_histogram=%s" % (counts, exit_histogram))
    print("rho_needed_divisor_packets=%s" % (packets,))
    print("rho_killers=%s histogram=%s" % (killers, killer_histogram))
    print("inadmissible_count=%d survivors=%s" % (len(inadmissible), survivors))
    print("first_residual=(d=4034,divisor_killer=109); first_rho=(d=4150,packet=(3227,),killer=11)")
    print("semantic_sha256=%s cross_record_equal=True" % EXPECTED_SEMANTIC_SHA256)
    print("consequence=FINITE-EXACT closure through d=6000, equivalently r<=5998; no first survivor in universe")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
