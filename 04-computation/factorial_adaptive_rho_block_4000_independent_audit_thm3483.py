#!/usr/bin/env python3
"""Self-contained exact referee for the THM-3483 block through d=4000."""

import ast
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path

START, END = 2606, 4000
PRIMES_47 = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)
PRIMES_101 = PRIMES_47 + (53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101)
REFERENCE_HASHES = (
    ("factorial_all_divisor_digit_pair_compiler_independent_audit_thm3475.py",
     "9330ca1b991b9d5875779b9975fc88701ab36855a6527e1865e821e6cd3ea665"),
    ("factorial_nondivisor_residue_digit_pair_compiler_thm3483.py",
     "9e37ead620f141617a9c6d51c182e09c034945793092e56e39fb061254662723"),
)
EXPECTED_COUNTS = (1395, 934, 461, 420, 41, 0)
EXPECTED_EXIT_HISTOGRAM = (
    ("d_minus_1_prime_power", 176), ("d_minus_2_prime", 143),
    ("d_minus_3_prime", 142), ("d_minus_4_prime", 112),
    ("d_minus_5_prime", 111), ("d_minus_6_prime", 78), ("d_prime", 172),
)
EXPECTED_RHO_HISTOGRAM = ((3, 6), (5, 1), (7, 17), (11, 9), (13, 6), (17, 2))
EXPECTED_OVERLAP = (1200, 782, 388, 30, ((3, 4), (5, 1), (7, 13), (11, 7), (13, 3), (17, 2)))
EXPECTED_D2606_PACKET = (521, 1042, 1563, 2084)
EXPECTED_D2606_BLOCKS = (
    (2, 3, 3), (8, 9, 9), (80, 81, 81), (242, 243, 243),
    (2186, 2187, 2187), (1, 1, 27), (55, 54, 54),
)
EXPECTED_SEMANTIC_SHA256 = "95d1c233d59d00c38ce456fa7c5f5e248414e01b5ba9dc2ae9f61725d6c19dbd"


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def digest(value):
    data = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return hashlib.sha256(data).hexdigest()


def verify_references():
    result = []
    here = Path(__file__).resolve().parent
    for name, expected in REFERENCE_HASHES:
        actual = hashlib.sha256((here / name).read_bytes().replace(b"\r\n", b"\n")).hexdigest()
        require(actual == expected, (name, actual, expected))
        result.append((name, actual))
    return tuple(result)


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
    result, q = [], 2
    while q * q <= n:
        if n % q == 0:
            result.append(q)
            while n % q == 0:
                n //= q
        q += 1
    if n > 1:
        result.append(n)
    return tuple(result)


def exit_label(d):
    tests = (("d_prime", is_prime(d)),
             ("d_minus_1_prime_power", len(prime_divisors(d - 1)) == 1),
             ("d_minus_2_prime", is_prime(d - 2)), ("d_minus_3_prime", is_prime(d - 3)),
             ("d_minus_4_prime", is_prime(d - 4)), ("d_minus_5_prime", is_prime(d - 5)),
             ("d_minus_6_prime", is_prime(d - 6)))
    return next((label for label, yes in tests if yes), None)


class WeightBank:
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
        t = self.tables[p]
        return t[n] - t[j] - t[n - j] + t[2 * j]


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
    return lower_hull(tuple((j, bank.weight(n, j, p)) for j in range(start, n + 1, step)))


def common_barcode(f_hull, g_hull, positive_only):
    def ledger(hull):
        answer = {}
        for a, b in zip(hull, hull[1:]):
            slope = Fraction(b[1] - a[1], b[0] - a[0])
            require(slope not in answer, ("repeated slope", slope))
            answer[slope] = b[0] - a[0]
        return answer
    f_ledger, g_ledger = ledger(f_hull), ledger(g_hull)
    bits, blocks = 1, []
    for slope in sorted(set(f_ledger) & set(g_ledger)):
        if positive_only:
            require(slope > 0, ("nonpositive divisor suffix", slope))
        capacity, denominator = min(f_ledger[slope], g_ledger[slope]), slope.denominator
        require(capacity % denominator == 0, ("capacity", slope, capacity))
        blocks.append((slope.numerator, denominator, capacity))
        old = bits
        for use in range(denominator, capacity + 1, denominator):
            bits |= old << use
    return tuple(k for k in range(bits.bit_length()) if (bits >> k) & 1), tuple(blocks)


def rho(n, j, d, p):
    require(d % p != 0, ("rho nondivisor", d, p))
    m, inverse = (n - j) % p, pow(d % p, -1, p)
    choose = rising = neg_power = 1
    total, twice_j = 0, 2 * (j % p)
    for ell in range(m + 1):
        if ell:
            choose = choose * (m - ell + 1) * pow(ell, -1, p) % p
            rising = rising * (twice_j + ell) % p
            neg_power = neg_power * (-inverse) % p
        total = (total + choose * neg_power * rising) % p
    return total


def rho_profile(bank, n, d, p):
    hull = raw_hull(bank, n, p)
    vertices = tuple((j, rho(n, j, d, p)) for j, _ in hull)
    return hull, vertices, all(value for _, value in vertices)


def scan_row(bank, d):
    n, packet, divisor_trace = d - 1, tuple(range(1, d - 1)), []
    for p in prime_divisors(n):
        g_hull = raw_hull(bank, n, p)
        f_hull = raw_hull(bank, n - 1, p, 1, 2) if p == 2 else raw_hull(bank, n - 1, p, (p - 1) // 2)
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
        pre, blocks = packet, ()
        status = "inadmissible"
        if f_ok and g_ok:
            degrees, blocks = common_barcode(f_hull, g_hull, False)
            local, status = set(degrees), "admissible"
            packet = tuple(k for k in packet if k in local)
        rho_trace.append((p, status, f_vertices, g_vertices, blocks, pre, packet))
        if not packet:
            return (d, "rho", p, tuple(divisor_trace), tuple(rho_trace), divisor_packet)
    return (d, "survivor", None, tuple(divisor_trace), tuple(rho_trace), divisor_packet, packet)


def build_semantic():
    factors = {p for n in range(START - 1, END) for p in prime_divisors(n)}
    bank = WeightBank(tuple(sorted(factors | set(PRIMES_101))))
    exits, residuals = [], []
    for d in range(START, END + 1):
        label = exit_label(d)
        (exits if label else residuals).append((d, label) if label else scan_row(bank, d))
    return ((START, END), tuple(exits), tuple(residuals), PRIMES_47, PRIMES_101), bank


def summarize(semantic):
    _, exits, residuals, _, _ = semantic
    divisor = tuple(r for r in residuals if r[1] == "divisor")
    rho_rows = tuple(r for r in residuals if r[1] == "rho")
    survivors = tuple(r for r in residuals if r[1] == "survivor")
    packets = tuple((r[0], r[5]) for r in rho_rows)
    killers = tuple((r[0], r[2]) for r in rho_rows)
    hostiles = tuple((r[0], e[0], e[2], e[3]) for r in rho_rows for e in r[4] if e[1] == "inadmissible")
    return exits, divisor, rho_rows, survivors, packets, killers, hostiles


def main():
    references = verify_references()
    semantic, bank = build_semantic()
    exits, divisor, rho_rows, survivors, packets, killers, hostiles = summarize(semantic)
    counts = (END - START + 1, len(exits), len(divisor) + len(rho_rows) + len(survivors), len(divisor), len(rho_rows), len(survivors))
    exit_hist = tuple(sorted(Counter(label for _, label in exits).items()))
    rho_hist = tuple(sorted(Counter(p for _, p in killers).items()))
    overlap = tuple(r for r in semantic[2] if r[0] >= 2801)
    overlap_rho = tuple(r for r in overlap if r[1] == "rho")
    overlap_summary = (1200, sum(exit_label(d) is not None for d in range(2801, 4001)),
                       sum(r[1] == "divisor" for r in overlap), len(overlap_rho),
                       tuple(sorted(Counter(r[2] for r in overlap_rho).items())))
    require(counts == EXPECTED_COUNTS, counts)
    require(exit_hist == EXPECTED_EXIT_HISTOGRAM, exit_hist)
    require(rho_hist == EXPECTED_RHO_HISTOGRAM, rho_hist)
    require(overlap_summary == EXPECTED_OVERLAP, overlap_summary)
    require(not survivors and not tuple((d, p) for d, p in killers if p > 47), (survivors, killers))
    d2606 = next(r for r in rho_rows if r[0] == 2606)
    require(d2606[5] == EXPECTED_D2606_PACKET and d2606[2] == 3, d2606[:3])
    require(d2606[4][-1][4] == EXPECTED_D2606_BLOCKS and not d2606[4][-1][6], d2606[4][-1])
    require(exit_label(2810) == "d_minus_1_prime_power", exit_label(2810))
    pair_6, blocks_6 = common_barcode(raw_hull(bank, 5, 3, 1), raw_hull(bank, 6, 3), True)
    require(pair_6 == (0, 3) and blocks_6 == ((2, 3, 3),), (pair_6, blocks_6))
    hostile = next(e for e in next(r for r in rho_rows if r[0] == 2608)[4] if e[0] == 5)
    require(hostile[1] == "inadmissible" and hostile[3][0] == (0, 0) and hostile[5] == hostile[6], hostile)
    zero_block = next(e for e in next(r for r in rho_rows if r[0] == 2608)[4] if e[0] == 7)
    require(zero_block[1] == "admissible" and (0, 1, 2) in zero_block[4] and not zero_block[6], zero_block)
    require(len(hostiles) == 20, len(hostiles))
    require(digest(semantic) == EXPECTED_SEMANTIC_SHA256, digest(semantic))
    source = Path(__file__).read_text(encoding="utf-8")
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))), "assert node")
    print("THM-3483 ADAPTIVE RHO BLOCK THROUGH d=4000 INDEPENDENT AUDIT")
    print("reference_hashes=%s" % (references,))
    print("implementation=direct Legendre tables; exact determinant hull; bitset denominator DP; finite-field rho recurrence; no imports")
    print("universe=d=2606..4000; seven exits in canonical order; rows F=A_(d-2)^d,G=A_(d-1)^d")
    print("counts=%s exit_histogram=%s" % (counts, exit_hist))
    print("rho_needed_divisor_packets=%s" % (packets,))
    print("rho_killers=%s histogram=%s extension_above_47=() survivors=()" % (killers, rho_hist))
    print("overlap_2801_4000=%s" % (overlap_summary,))
    print("d2606=(divisor_packet=%s,rho_killer=3,common_blocks=%s,post=())" % (EXPECTED_D2606_PACKET, EXPECTED_D2606_BLOCKS))
    print("hostiles=(inadmissible_count=%d,first=(d=2608,p=5,G_vertex=(0,0),packet_unchanged)); skips=(d=2606,p=2 divides d)" % len(hostiles))
    print("positive_controls=(d2810 exits by 2809=53^2; divisor (N,p)=(6,3) gives degrees (0,3); d2608,p7 retains zero-slope unit block (0,1,2) and closes)")
    print("semantic_schema=((start,end),exit_records,residual_records,ordered_primes_47,ordered_primes_101); exit=(d,first_exit_label); divisor_trace=(p,blocks,pre,post); rho_trace=(p,status,F_vertex_rhos,G_vertex_rhos,blocks,pre,post)")
    print("serialization=json.dumps(value,separators=(',',':'),sort_keys=True).encode('ascii')")
    print("semantic_sha256=%s" % EXPECTED_SEMANTIC_SHA256)
    print("consequence=FINITE-EXACT closure through d=4000, equivalently r<=3998; no first survivor in universe")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
