#!/usr/bin/env python3
"""Self-contained exact audit of the factorial pair frontier, d=4001..6000.

The policy is inherited verbatim from THM-3483: apply the seven elementary
exits in their canonical order, intersect every THM-3475 divisor-prime pair
barcode, and only then try rho-admissible THM-3483 raw pair polygons at the
ordered prime bank through 101.  No large factorial coefficient is built.
"""

import ast
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path


START, END = 4001, 6000
PRIMES_47 = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)
PRIMES_101 = PRIMES_47 + (53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101)
REFERENCE_HASHES = (
    ("factorial_all_divisor_digit_pair_compiler_independent_audit_thm3475.py",
     "9330ca1b991b9d5875779b9975fc88701ab36855a6527e1865e821e6cd3ea665"),
    ("factorial_nondivisor_residue_digit_pair_compiler_thm3483.py",
     "9e37ead620f141617a9c6d51c182e09c034945793092e56e39fb061254662723"),
)
EXPECTED_COUNTS = (2000, 1272, 728, 600, 128, 0)
EXPECTED_EXIT_HISTOGRAM = (
    ("d_minus_1_prime_power", 238), ("d_minus_2_prime", 193),
    ("d_minus_3_prime", 192), ("d_minus_4_prime", 159),
    ("d_minus_5_prime", 159), ("d_minus_6_prime", 98),
    ("d_prime", 233),
)
EXPECTED_RHO_HISTOGRAM = (
    (3, 1), (5, 16), (7, 46), (11, 42), (13, 10), (17, 9), (19, 4),
)
EXPECTED_SEMANTIC_SHA256 = "7f8ab74ae9fae027f32fd7eabaf0338c217319e274594bd603859a1bbcca28bd"


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def digest(value):
    encoded = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def verify_references():
    here = Path(__file__).resolve().parent
    answer = []
    for name, expected in REFERENCE_HASHES:
        actual = hashlib.sha256((here / name).read_bytes().replace(b"\r\n", b"\n")).hexdigest()
        require(actual == expected, (name, actual, expected))
        answer.append((name, actual))
    return tuple(answer)


def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    q = 3
    while q * q <= n:
        if n % q == 0:
            return False
        q += 2
    return True


def prime_divisors(n):
    answer, q = [], 2
    while q * q <= n:
        if n % q == 0:
            answer.append(q)
            while n % q == 0:
                n //= q
        q += 1
    if n > 1:
        answer.append(n)
    return tuple(answer)


def exit_label(d):
    tests = (
        ("d_prime", is_prime(d)),
        ("d_minus_1_prime_power", len(prime_divisors(d - 1)) == 1),
        ("d_minus_2_prime", is_prime(d - 2)),
        ("d_minus_3_prime", is_prime(d - 3)),
        ("d_minus_4_prime", is_prime(d - 4)),
        ("d_minus_5_prime", is_prime(d - 5)),
        ("d_minus_6_prime", is_prime(d - 6)),
    )
    return next((label for label, true in tests if true), None)


class WeightBank:
    """Repeated-factor valuation tables, independent of Legendre digit code."""

    def __init__(self, primes):
        self.tables = {}
        for p in primes:
            table = [0] * (2 * END + 1)
            for k in range(1, len(table)):
                x, valuation = k, 0
                while x % p == 0:
                    valuation += 1
                    x //= p
                table[k] = table[k - 1] + valuation
            self.tables[p] = table

    def weight(self, n, j, p):
        table = self.tables[p]
        return table[n] - table[j] - table[n - j] + table[2 * j]


def lower_hull(points):
    hull = []
    for point in points:
        while len(hull) >= 2:
            a, b = hull[-2], hull[-1]
            determinant = ((b[0] - a[0]) * (point[1] - b[1])
                           - (b[1] - a[1]) * (point[0] - b[0]))
            if determinant > 0:
                break
            hull.pop()
        hull.append(point)
    require(len(hull) >= 2, ("short hull", points[:5]))
    return tuple(hull)


def raw_hull(bank, n, p, start=0, step=1):
    return lower_hull(tuple((j, bank.weight(n, j, p))
                            for j in range(start, n + 1, step)))


def common_barcode(f_hull, g_hull, positive_only):
    def ledger(hull):
        answer = {}
        for a, b in zip(hull, hull[1:]):
            slope = Fraction(b[1] - a[1], b[0] - a[0])
            require(slope not in answer, ("repeated slope", slope, hull))
            answer[slope] = b[0] - a[0]
        return answer

    f_ledger, g_ledger = ledger(f_hull), ledger(g_hull)
    bits, blocks = 1, []
    for slope in sorted(set(f_ledger) & set(g_ledger)):
        if positive_only:
            require(slope > 0, ("nonpositive divisor suffix", slope))
        capacity = min(f_ledger[slope], g_ledger[slope])
        denominator = slope.denominator
        require(capacity % denominator == 0, ("capacity", slope, capacity))
        blocks.append((slope.numerator, denominator, capacity))
        old = bits
        for use in range(denominator, capacity + 1, denominator):
            bits |= old << use
    degrees = tuple(k for k in range(bits.bit_length()) if (bits >> k) & 1)
    return degrees, tuple(blocks)


def rho(n, j, d, p):
    require(d % p != 0, ("rho requires a nondivisor", d, p))
    m, inverse = (n - j) % p, pow(d % p, -1, p)
    choose = rising = negative_power = 1
    total, twice_j = 0, 2 * (j % p)
    for ell in range(m + 1):
        if ell:
            choose = choose * (m - ell + 1) * pow(ell, -1, p) % p
            rising = rising * (twice_j + ell) % p
            negative_power = negative_power * (-inverse) % p
        total = (total + choose * negative_power * rising) % p
    return total


def rho_profile(bank, n, d, p):
    hull = raw_hull(bank, n, p)
    vertices = tuple((j, rho(n, j, d, p)) for j, _ in hull)
    return hull, vertices, all(value for _, value in vertices)


def scan_row(bank, d):
    n, packet, divisor_trace = d - 1, tuple(range(1, d - 1)), []
    for p in prime_divisors(n):
        g_hull = raw_hull(bank, n, p)
        f_hull = (raw_hull(bank, n - 1, p, 1, 2)
                  if p == 2 else raw_hull(bank, n - 1, p, (p - 1) // 2))
        degrees, blocks = common_barcode(f_hull, g_hull, True)
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
        f_hull, f_vertices, f_ok = rho_profile(bank, n - 1, d, p)
        g_hull, g_vertices, g_ok = rho_profile(bank, n, d, p)
        pre, blocks, status = packet, (), "inadmissible"
        if f_ok and g_ok:
            degrees, blocks = common_barcode(f_hull, g_hull, False)
            packet = tuple(k for k in packet if k in set(degrees))
            status = "admissible"
        rho_trace.append((p, status, f_vertices, g_vertices, blocks, pre, packet))
        if not packet:
            return (d, "rho", p, tuple(divisor_trace), tuple(rho_trace), divisor_packet)
    return (d, "survivor", None, tuple(divisor_trace), tuple(rho_trace),
            divisor_packet, packet)


def build_semantic():
    factors = {p for n in range(START - 1, END) for p in prime_divisors(n)}
    bank = WeightBank(tuple(sorted(factors | set(PRIMES_101))))
    exits, residuals = [], []
    for d in range(START, END + 1):
        label = exit_label(d)
        (exits if label else residuals).append(
            (d, label) if label else scan_row(bank, d)
        )
    return ((START, END), tuple(exits), tuple(residuals), PRIMES_47, PRIMES_101), bank


def summarize(semantic):
    _, exits, residuals, _, _ = semantic
    divisor_rows = tuple(r for r in residuals if r[1] == "divisor")
    rho_rows = tuple(r for r in residuals if r[1] == "rho")
    survivors = tuple(r for r in residuals if r[1] == "survivor")
    packets = tuple((r[0], r[5]) for r in rho_rows)
    killers = tuple((r[0], r[2]) for r in rho_rows)
    inadmissible = tuple(
        (r[0], entry[0], entry[2], entry[3])
        for r in rho_rows + survivors
        for entry in r[4]
        if entry[1] == "inadmissible"
    )
    return exits, divisor_rows, rho_rows, survivors, packets, killers, inadmissible


def compact_row(row):
    divisor = tuple((p, blocks, len(pre), len(post), post if len(post) <= 32 else ())
                    for p, blocks, pre, post in row[3])
    rho_trace = tuple((p, status, blocks, pre, post)
                      for p, status, _, _, blocks, pre, post in row[4])
    return row[0], row[1], row[2], divisor, rho_trace, row[5]


def degree_witness(blocks, target):
    paths = {0: ()}
    for _, denominator, capacity in blocks:
        next_paths = {}
        for old, path in paths.items():
            for use in range(0, capacity + 1, denominator):
                next_paths.setdefault(old + use, path + (use,))
        paths = next_paths
    return paths.get(target)


def main():
    references = verify_references()
    semantic, bank = build_semantic()
    exits, divisor_rows, rho_rows, survivors, packets, killers, inadmissible = summarize(semantic)
    counts = (END - START + 1, len(exits), len(semantic[2]), len(divisor_rows),
              len(rho_rows), len(survivors))
    exit_histogram = tuple(sorted(Counter(label for _, label in exits).items()))
    killer_histogram = tuple(sorted(Counter(p for _, p in killers).items()))
    first_residual = semantic[2][0] if semantic[2] else ()
    first_rho = rho_rows[0] if rho_rows else ()

    # Positive and hostile controls inherited from the compiler boundary.
    require(exit_label(4001) == "d_prime", exit_label(4001))
    pair_6, blocks_6 = common_barcode(raw_hull(bank, 5, 3, 1),
                                      raw_hull(bank, 6, 3), True)
    require(pair_6 == (0, 3) and blocks_6 == ((2, 3, 3),), (pair_6, blocks_6))
    require(first_residual, "missing residual control")
    require(counts == EXPECTED_COUNTS, counts)
    require(exit_histogram == EXPECTED_EXIT_HISTOGRAM, exit_histogram)
    require(killer_histogram == EXPECTED_RHO_HISTOGRAM, killer_histogram)
    require(not survivors, survivors)
    require(digest(semantic) == EXPECTED_SEMANTIC_SHA256, digest(semantic))
    require(first_residual[:3] == (4034, "divisor", 109), first_residual[:3])
    require(first_rho[:3] == (4150, "rho", 11) and first_rho[5] == (3227,),
            first_rho[:3])
    p3 = next(entry for entry in first_rho[4] if entry[0] == 3)
    p7 = next(entry for entry in first_rho[4] if entry[0] == 7)
    p11 = next(entry for entry in first_rho[4] if entry[0] == 11)
    require(degree_witness(p3[4], 3227) == (0, 243, 729, 2187, 68), p3[4])
    require(degree_witness(p7[4], 3227) == (0, 7, 686, 2401, 126, 0, 7), p7[4])
    require(degree_witness(p11[4], 3227) is None and not p11[6], p11)
    require(not any(isinstance(node, ast.Assert)
                    for node in ast.walk(ast.parse(Path(__file__).read_text()))),
            "assert node")

    print("FACTORIAL ADAPTIVE RHO BLOCK d=4001..6000 INDEPENDENT AUDIT")
    print("reference_hashes=%s" % (references,))
    print("implementation=direct repeated-factor tables; determinant hull; bitset denominator DP; finite-field rho recurrence; no compiler imports")
    print("universe=d=4001..6000; canonical seven exits; F=A_(d-2)^d,G=A_(d-1)^d; divisor primes of d-1; ordered rho bank=%s" % (PRIMES_101,))
    print("counts=%s exit_histogram=%s" % (counts, exit_histogram))
    print("rho_needed_divisor_packets=%s" % (packets,))
    print("rho_killers=%s histogram=%s" % (killers, killer_histogram))
    print("survivors=%s" % (tuple((r[0], r[5], r[6]) for r in survivors),))
    print("inadmissible_count=%d inadmissible_profiles=%s" % (len(inadmissible), inadmissible))
    print("first_residual=%s" % (compact_row(first_residual),))
    print("first_rho_row=%s" % (compact_row(first_rho),))
    print("d4150_structure=(4149=3^2*461; divisor p3 then p461 leaves 3227=7*461=2187+729+243+68; rho p3 witness=%s; rho p7 witness=%s; rho p11 blocks=%s have low maximum 2817 and next positive branch 3993, so 3227 dies)" % (degree_witness(p3[4], 3227), degree_witness(p7[4], 3227), p11[4]))
    print("positive_controls=(d4001 is prime; divisor (N,p)=(6,3) gives degrees (0,3))")
    print("hostile_control=every inadmissible profile is recorded and skipped with its packet unchanged")
    print("semantic_schema=((start,end),exit_records,residual_records,ordered_primes_47,ordered_primes_101)")
    print("semantic_sha256=%s" % EXPECTED_SEMANTIC_SHA256)
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
