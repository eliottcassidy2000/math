#!/usr/bin/env python3
"""Pinned-engine exact THM-3483 block census through d=4000."""

import ast
from collections import Counter
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path

START, END = 2606, 4000
PRIMES_47 = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)
PRIMES_101 = PRIMES_47 + (53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101)
DIVISOR_NAME = "factorial_all_divisor_digit_pair_compiler_thm3475.py"
DIVISOR_SHA256 = "834d0913eb5cd5b15684c7fb88af60e42d2a6ef36feb821e261c0498f55027ab"
RHO_NAME = "factorial_nondivisor_residue_digit_pair_compiler_thm3483.py"
RHO_SHA256 = "9e37ead620f141617a9c6d51c182e09c034945793092e56e39fb061254662723"
INDEPENDENT_NAME = "factorial_adaptive_rho_block_4000_independent_audit_thm3483.py"
INDEPENDENT_SHA256 = "0b858f7b1154a3ee2dec43bf5238f7e6f24b524e9dc2b9f70f5b497fcef58934"
EXPECTED_SEMANTIC_SHA256 = "95d1c233d59d00c38ce456fa7c5f5e248414e01b5ba9dc2ae9f61725d6c19dbd"


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
    divisor = load_pinned(DIVISOR_NAME, DIVISOR_SHA256, "rho4000_divisor")
    rho_module = load_pinned(RHO_NAME, RHO_SHA256, "rho4000_rho")
    independent = load_pinned(INDEPENDENT_NAME, INDEPENDENT_SHA256, "rho4000_independent")
    coefficient = divisor.load_pinned(divisor.PRIMARY_NAME, divisor.PRIMARY_SHA256,
                                      "rho4000_coefficient")
    digital = divisor.load_pinned(divisor.DIGITAL_NAME, divisor.DIGITAL_SHA256,
                                  "rho4000_digital")
    return divisor, rho_module, independent, coefficient, digital


def factor_primes(coefficient, n):
    return tuple(p for p, _ in coefficient.factor(n))


def is_prime(coefficient, n):
    factors = coefficient.factor(n) if n >= 2 else ()
    return len(factors) == 1 and factors[0][1] == 1


def exit_label(coefficient, d):
    tests = (("d_prime", is_prime(coefficient, d)),
             ("d_minus_1_prime_power", len(coefficient.factor(d - 1)) == 1),
             ("d_minus_2_prime", is_prime(coefficient, d - 2)),
             ("d_minus_3_prime", is_prime(coefficient, d - 3)),
             ("d_minus_4_prime", is_prime(coefficient, d - 4)),
             ("d_minus_5_prime", is_prime(coefficient, d - 5)),
             ("d_minus_6_prime", is_prime(coefficient, d - 6)))
    return next((label for label, yes in tests if yes), None)


def normalized_divisor_barcode(divisor, coefficient, digital, n, p):
    _, _, blocks, degrees = divisor.predicted_pair_compiler(coefficient, digital, n, p)
    normalized = tuple((s.numerator, s.denominator, capacity)
                       for s, capacity, _ in blocks)
    return tuple(sorted(degrees)), normalized


def normalized_rho_barcode(rho_module, f_hull, g_hull):
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
        degrees, blocks = normalized_divisor_barcode(divisor, coefficient, digital, n, p)
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
        f_hull, g_hull = rho_module.raw_hull(n - 1, p), rho_module.raw_hull(n, p)
        f_vertices = tuple((j, rho_module.rho(n - 1, j, d, p)) for j, _ in f_hull)
        g_vertices = tuple((j, rho_module.rho(n, j, d, p)) for j, _ in g_hull)
        f_ok = all(value for _, value in f_vertices)
        g_ok = all(value for _, value in g_vertices)
        pre, blocks, status = packet, (), "inadmissible"
        if f_ok and g_ok:
            degrees, blocks = normalized_rho_barcode(rho_module, f_hull, g_hull)
            packet, status = tuple(k for k in packet if k in set(degrees)), "admissible"
        rho_trace.append((p, status, f_vertices, g_vertices, blocks, pre, packet))
        if not packet:
            return (d, "rho", p, tuple(divisor_trace), tuple(rho_trace), divisor_packet)
    return (d, "survivor", None, tuple(divisor_trace), tuple(rho_trace), divisor_packet, packet)


def build_semantic(divisor, rho_module, coefficient, digital):
    exits, residuals = [], []
    for d in range(START, END + 1):
        label = exit_label(coefficient, d)
        if label:
            exits.append((d, label))
        else:
            residuals.append(scan_row(divisor, rho_module, coefficient, digital, d))
    return ((START, END), tuple(exits), tuple(residuals), PRIMES_47, PRIMES_101)


def main():
    divisor, rho_module, independent, coefficient, digital = load_engines()
    semantic = build_semantic(divisor, rho_module, coefficient, digital)
    independent_semantic, _ = independent.build_semantic()
    require(semantic == independent_semantic, (digest(semantic), digest(independent_semantic)))
    require(digest(semantic) == EXPECTED_SEMANTIC_SHA256, digest(semantic))
    exits, divisor_rows, rho_rows, survivors, packets, killers, hostiles = independent.summarize(semantic)
    counts = (END - START + 1, len(exits), len(semantic[2]), len(divisor_rows), len(rho_rows), len(survivors))
    exit_hist = tuple(sorted(Counter(label for _, label in exits).items()))
    rho_hist = tuple(sorted(Counter(p for _, p in killers).items()))
    require(counts == independent.EXPECTED_COUNTS, counts)
    require(exit_hist == independent.EXPECTED_EXIT_HISTOGRAM, exit_hist)
    require(rho_hist == independent.EXPECTED_RHO_HISTOGRAM, rho_hist)
    require(len(hostiles) == 20 and not survivors, (len(hostiles), survivors))
    d2606 = next(r for r in rho_rows if r[0] == 2606)
    require(d2606[5] == independent.EXPECTED_D2606_PACKET and d2606[4][-1][4] == independent.EXPECTED_D2606_BLOCKS, d2606[:3])
    source = Path(__file__).read_text(encoding="utf-8")
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))), "assert node")
    print("THM-3483 ADAPTIVE RHO BLOCK THROUGH d=4000 PINNED-ENGINE CENSUS")
    print("dependencies=%s:%s;%s:%s;%s:%s" % (DIVISOR_NAME, DIVISOR_SHA256, RHO_NAME, RHO_SHA256, INDEPENDENT_NAME, INDEPENDENT_SHA256))
    print("implementation=THM-3475 predicted divisor compiler + THM-3483 rho/raw hull compiler; exact cross-record equality with self-contained referee")
    print("universe=d=2606..4000; seven exits in canonical order; rows F=A_(d-2)^d,G=A_(d-1)^d")
    print("counts=%s exit_histogram=%s" % (counts, exit_hist))
    print("rho_needed_divisor_packets=%s" % (packets,))
    print("rho_killers=%s histogram=%s extension_above_47=() survivors=()" % (killers, rho_hist))
    print("d2606=(divisor_packet=%s,rho_killer=3,common_blocks=%s,post=())" % (independent.EXPECTED_D2606_PACKET, independent.EXPECTED_D2606_BLOCKS))
    print("hostiles=(inadmissible_count=20,recorded_and_skipped); positive_controls=(prime-power exit,divisor denominator DP,zero-slope unit block)")
    print("semantic_schema=((start,end),exit_records,residual_records,ordered_primes_47,ordered_primes_101); divisor_trace=(p,blocks,pre,post); rho_trace=(p,status,F_vertex_rhos,G_vertex_rhos,blocks,pre,post)")
    print("serialization=json.dumps(value,separators=(',',':'),sort_keys=True).encode('ascii')")
    print("cross_record_equal=True semantic_sha256=%s" % EXPECTED_SEMANTIC_SHA256)
    print("consequence=FINITE-EXACT closure through d=4000, equivalently r<=3998; no first survivor in universe")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
