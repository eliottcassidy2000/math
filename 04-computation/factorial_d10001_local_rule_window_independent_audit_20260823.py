#!/usr/bin/env python3
"""Independent hostile audit of the d=10001 factorial-barcode neighbourhood.

This file deliberately imports no repository computation module.  It rebuilds
the 9980..10030 scan twice:

* route A uses repeated-factor valuation tables, a streaming determinant hull,
  bitset degree propagation, and a finite-field rho recurrence;
* route B uses Kummer carry counts/digit sums, a gift-wrapping lower hull,
  explicit set propagation, and the literal finite rho sum.

Small-degree integer coefficient expansions independently test the theorem
interfaces used by both large formula-only routes.
"""

from __future__ import annotations

import ast
from collections import Counter
from fractions import Fraction
import hashlib
import json
from math import comb, factorial, gcd
from pathlib import Path
import sys


sys.stdout.reconfigure(newline="\n")
START, END, ANCHOR = 9980, 10030, 10001
PRIMES_101 = (
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41,
    43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101,
)
PRIMES_47 = PRIMES_101[:15]
EXPECTED_RHO_ROWS = (
    (9982, 11, (2218, 6654)),
    (9986, 3, (3994, 5991)),
    (9988, 7, (6658,)),
    (9994, 11, (3331, 6662)),
    (9996, 19, (3998,)),
    (10030, 7, (3343, 6686)),
)
EXPECTED_SOURCE_SEMANTIC = (
    "f2d02ba3f51d6dd5b47a0976c775b307a3b73b64e01affaecf198f349a635213"
)
GATES = 0


def require(condition: bool, label: object) -> None:
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def digest(value: object) -> str:
    payload = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return hashlib.sha256(payload).hexdigest()


def packet_digest(packet: tuple[int, ...]) -> str:
    return hashlib.sha256(",".join(map(str, packet)).encode("ascii")).hexdigest()


def is_prime_trial(n: int) -> bool:
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


def trial_prime_divisors(n: int) -> tuple[int, ...]:
    answer = []
    q = 2
    while q * q <= n:
        if n % q == 0:
            answer.append(q)
            while n % q == 0:
                n //= q
        q += 1
    if n > 1:
        answer.append(n)
    return tuple(answer)


def build_spf(limit: int) -> tuple[list[int], tuple[bool, ...]]:
    spf = list(range(limit + 1))
    if limit >= 1:
        spf[1] = 1
    for p in range(2, int(limit**0.5) + 1):
        if spf[p] == p:
            for multiple in range(p * p, limit + 1, p):
                if spf[multiple] == multiple:
                    spf[multiple] = p
    primality = tuple(n >= 2 and spf[n] == n for n in range(limit + 1))
    return spf, primality


SPF, SIEVE_PRIME = build_spf(2 * END + 10)


def spf_prime_divisors(n: int) -> tuple[int, ...]:
    answer = []
    while n > 1:
        p = SPF[n]
        answer.append(p)
        while n % p == 0:
            n //= p
    return tuple(answer)


def exit_label_a(d: int) -> str | None:
    tests = (
        ("d_prime", is_prime_trial(d)),
        ("d_minus_1_prime_power", len(trial_prime_divisors(d - 1)) == 1),
        ("d_minus_2_prime", is_prime_trial(d - 2)),
        ("d_minus_3_prime", is_prime_trial(d - 3)),
        ("d_minus_4_prime", is_prime_trial(d - 4)),
        ("d_minus_5_prime", is_prime_trial(d - 5)),
        ("d_minus_6_prime", is_prime_trial(d - 6)),
    )
    return next((name for name, true in tests if true), None)


def exit_label_b(d: int) -> str | None:
    tests = (
        ("d_prime", SIEVE_PRIME[d]),
        ("d_minus_1_prime_power", len(spf_prime_divisors(d - 1)) == 1),
        ("d_minus_2_prime", SIEVE_PRIME[d - 2]),
        ("d_minus_3_prime", SIEVE_PRIME[d - 3]),
        ("d_minus_4_prime", SIEVE_PRIME[d - 4]),
        ("d_minus_5_prime", SIEVE_PRIME[d - 5]),
        ("d_minus_6_prime", SIEVE_PRIME[d - 6]),
    )
    return next((name for name, true in tests if true), None)


class RepeatedFactorWeights:
    """Factorial valuations accumulated from v_p(k), with no digit formula."""

    def __init__(self) -> None:
        self.tables: dict[int, tuple[int, ...]] = {}

    def table(self, p: int) -> tuple[int, ...]:
        if p not in self.tables:
            values = [0] * (2 * END + 2)
            for k in range(1, len(values)):
                x, order = k, 0
                while x % p == 0:
                    x //= p
                    order += 1
                values[k] = values[k - 1] + order
            self.tables[p] = tuple(values)
        return self.tables[p]

    def weight(self, n: int, j: int, p: int) -> int:
        table = self.table(p)
        return table[n] - table[j] - table[n - j] + table[2 * j]


def digit_sum(n: int, p: int) -> int:
    total = 0
    while n:
        total += n % p
        n //= p
    return total


def factorial_valuation_digits(n: int, p: int) -> int:
    return (n - digit_sum(n, p)) // (p - 1)


def binomial_valuation_carries(n: int, j: int, p: int) -> int:
    """Kummer: number of carries when j and n-j are added in base p."""
    a, b, carry, count = j, n - j, 0, 0
    while a or b or carry:
        total = a % p + b % p + carry
        carry = int(total >= p)
        count += carry
        a //= p
        b //= p
    return count


def weight_b(n: int, j: int, p: int) -> int:
    return binomial_valuation_carries(n, j, p) + factorial_valuation_digits(2 * j, p)


def streaming_lower_hull(points: tuple[tuple[int, int], ...]) -> tuple[tuple[int, int], ...]:
    hull: list[tuple[int, int]] = []
    for point in points:
        while len(hull) >= 2:
            a, b = hull[-2], hull[-1]
            cross = ((b[0] - a[0]) * (point[1] - b[1])
                     - (b[1] - a[1]) * (point[0] - b[0]))
            if cross > 0:
                break
            hull.pop()
        hull.append(point)
    require(len(hull) >= 2, ("short streaming hull", points[:4]))
    return tuple(hull)


def giftwrap_lower_hull(points: tuple[tuple[int, int], ...]) -> tuple[tuple[int, int], ...]:
    """Choose the least outgoing slope, farthest endpoint on a tie."""
    require(len(points) >= 2, ("short gift-wrap input", points))
    answer = [points[0]]
    i = 0
    while i < len(points) - 1:
        best = i + 1
        for j in range(i + 2, len(points)):
            dy = points[j][1] - points[i][1]
            dx = points[j][0] - points[i][0]
            best_dy = points[best][1] - points[i][1]
            best_dx = points[best][0] - points[i][0]
            comparison = dy * best_dx - best_dy * dx
            if comparison < 0 or (comparison == 0 and j > best):
                best = j
        answer.append(points[best])
        i = best
    return tuple(answer)


def slope_ledger(hull: tuple[tuple[int, int], ...]) -> dict[tuple[int, int], int]:
    answer: dict[tuple[int, int], int] = {}
    for a, b in zip(hull, hull[1:]):
        dx, dy = b[0] - a[0], b[1] - a[1]
        common = gcd(abs(dy), dx)
        slope = (dy // common, dx // common)
        require(slope not in answer, ("repeated slope", slope, hull))
        answer[slope] = dx
    return answer


def common_blocks(
    f_hull: tuple[tuple[int, int], ...],
    g_hull: tuple[tuple[int, int], ...],
    positive_only: bool,
) -> tuple[tuple[int, int, int], ...]:
    f_ledger, g_ledger = slope_ledger(f_hull), slope_ledger(g_hull)
    slopes = set(f_ledger) & set(g_ledger)
    if positive_only:
        require(all(num > 0 for num, _ in slopes), ("nonpositive divisor slope", slopes))
    ordered = sorted(slopes, key=lambda pair: Fraction(*pair))
    blocks = []
    for numerator, denominator in ordered:
        capacity = min(f_ledger[(numerator, denominator)], g_ledger[(numerator, denominator)])
        require(capacity % denominator == 0,
                ("capacity/denominator", numerator, denominator, capacity))
        blocks.append((numerator, denominator, capacity))
    return tuple(blocks)


def degrees_bitset(blocks: tuple[tuple[int, int, int], ...]) -> tuple[int, ...]:
    bits = 1
    for _, denominator, capacity in blocks:
        old = bits
        for use in range(denominator, capacity + 1, denominator):
            bits |= old << use
    return tuple(i for i in range(bits.bit_length()) if (bits >> i) & 1)


def degrees_sets(blocks: tuple[tuple[int, int, int], ...]) -> tuple[int, ...]:
    values = {0}
    for _, denominator, capacity in blocks:
        old = values
        values = {
            base + use
            for base in old
            for use in range(0, capacity + 1, denominator)
        }
    return tuple(sorted(values))


def rho_recurrence(n: int, j: int, d: int, p: int) -> int:
    require(d % p != 0, ("rho divisor", d, p))
    m, inverse = (n - j) % p, pow(d % p, -1, p)
    choose = rising = negative_power = 1
    total, twice_j = 0, 2 * (j % p)
    for ell in range(m + 1):
        if ell:
            choose = choose * (m - ell + 1) * pow(ell, -1, p) % p
            rising = rising * (twice_j + ell) % p
            negative_power = negative_power * (-inverse) % p
        total = (total + choose * rising * negative_power) % p
    return total


def rho_literal(n: int, j: int, d: int, p: int) -> int:
    require(d % p != 0, ("rho divisor literal", d, p))
    m, inverse = (n - j) % p, pow(d % p, -1, p)
    total = 0
    for ell in range(m + 1):
        rising = 1
        for shift in range(ell):
            rising = rising * (2 * (j % p) + 1 + shift) % p
        total += comb(m, ell) * pow(-inverse, ell, p) * rising
    return total % p


class Route:
    def __init__(self, name: str, weight, hull, degree_engine, rho_engine, factors, exit_label):
        self.name = name
        self.weight = weight
        self.hull_engine = hull
        self.degree_engine = degree_engine
        self.rho_engine = rho_engine
        self.factors = factors
        self.exit_label = exit_label
        self.hull_cache: dict[tuple[int, int, int, int], tuple[tuple[int, int], ...]] = {}
        self.divisor_cache: dict[tuple[int, int], tuple[tuple[int, ...], tuple[tuple[int, int, int], ...]]] = {}

    def hull(self, n: int, p: int, start: int = 0, step: int = 1):
        key = (n, p, start, step)
        if key not in self.hull_cache:
            points = tuple((j, self.weight(n, j, p)) for j in range(start, n + 1, step))
            self.hull_cache[key] = self.hull_engine(points)
        return self.hull_cache[key]

    def divisor_barcode(self, n: int, p: int):
        key = (n, p)
        if key not in self.divisor_cache:
            g_hull = self.hull(n, p)
            f_hull = self.hull(n - 1, p, 1, 2) if p == 2 else self.hull(n - 1, p, (p - 1) // 2)
            blocks = common_blocks(f_hull, g_hull, True)
            self.divisor_cache[key] = (self.degree_engine(blocks), blocks)
        return self.divisor_cache[key]

    def rho_barcode(self, d: int, p: int):
        n = d - 1
        f_hull, g_hull = self.hull(n - 1, p), self.hull(n, p)
        f_vertices = tuple((j, self.rho_engine(n - 1, j, d, p)) for j, _ in f_hull)
        g_vertices = tuple((j, self.rho_engine(n, j, d, p)) for j, _ in g_hull)
        if not all(value for _, value in f_vertices + g_vertices):
            return "inadmissible", f_vertices, g_vertices, (), ()
        blocks = common_blocks(f_hull, g_hull, False)
        return "admissible", f_vertices, g_vertices, self.degree_engine(blocks), blocks

    def scan_row(self, d: int):
        n = d - 1
        packet = tuple(range(1, d - 1))
        divisor_trace = []
        for p in self.factors(n):
            degrees, blocks = self.divisor_barcode(n, p)
            pre = packet
            local = set(degrees)
            packet = tuple(value for value in packet if value in local)
            divisor_trace.append((p, blocks, pre, packet))
            if not packet:
                return d, "divisor", p, tuple(divisor_trace), (), ()

        divisor_packet = packet
        rho_trace = []
        for p in PRIMES_101:
            if d % p == 0:
                rho_trace.append((p, "divides_d", (), (), (), packet, packet))
                continue
            status, fv, gv, degrees, blocks = self.rho_barcode(d, p)
            pre = packet
            if status == "admissible":
                local = set(degrees)
                packet = tuple(value for value in packet if value in local)
            rho_trace.append((p, status, fv, gv, blocks, pre, packet))
            if not packet:
                return d, "rho", p, tuple(divisor_trace), tuple(rho_trace), divisor_packet
        return d, "survivor", None, tuple(divisor_trace), tuple(rho_trace), divisor_packet, packet

    def build(self):
        exits, residuals = [], []
        for d in range(START, END + 1):
            label = self.exit_label(d)
            if label:
                exits.append((d, label))
            else:
                residuals.append(self.scan_row(d))
        return (START, END), tuple(exits), tuple(residuals), PRIMES_101


TABLES = RepeatedFactorWeights()
ROUTE_A = Route(
    "repeated-factor/streaming/bitset/recurrence",
    TABLES.weight,
    streaming_lower_hull,
    degrees_bitset,
    rho_recurrence,
    trial_prime_divisors,
    exit_label_a,
)
ROUTE_B = Route(
    "Kummer-carry/gift-wrap/set/literal-sum",
    weight_b,
    giftwrap_lower_hull,
    degrees_sets,
    rho_literal,
    spf_prime_divisors,
    exit_label_b,
)


def integer_coefficient(n: int, d: int, j: int) -> int:
    return comb(n, j) * sum(
        comb(n - j, ell) * d ** (n - j - ell) * (-1) ** ell * factorial(2 * j + ell)
        for ell in range(n - j + 1)
    )


def integer_valuation(value: int, p: int) -> int | None:
    if value == 0:
        return None
    value, order = abs(value), 0
    while value % p == 0:
        value //= p
        order += 1
    return order


def actual_hull(n: int, d: int, p: int):
    points = []
    for j in range(n + 1):
        valuation = integer_valuation(integer_coefficient(n, d, j), p)
        if valuation is not None:
            points.append((j, valuation))
    return streaming_lower_hull(tuple(points))


def direct_coefficient_controls() -> tuple[int, int, int, int]:
    coefficient_checks = profile_checks = admissible_checks = divisor_checks = 0
    control_primes = (2, 3, 5, 7, 11, 13, 17, 19)
    for d in range(4, 25):
        n = d - 1
        for row_n in (n - 1, n):
            coefficients = tuple(integer_coefficient(row_n, d, j) for j in range(row_n + 1))
            for p in control_primes:
                if d % p == 0:
                    continue
                raw = ROUTE_A.hull(row_n, p)
                vertices = tuple((j, rho_recurrence(row_n, j, d, p)) for j, _ in raw)
                actual = streaming_lower_hull(tuple(
                    (j, valuation)
                    for j, value in enumerate(coefficients)
                    for valuation in (integer_valuation(value, p),)
                    if valuation is not None
                ))
                for j, value in enumerate(coefficients):
                    valuation = integer_valuation(value, p)
                    raw_value = TABLES.weight(row_n, j, p)
                    residue = rho_literal(row_n, j, d, p)
                    require((valuation == raw_value) == (residue != 0),
                            ("coefficient iff rho", d, row_n, p, j, valuation, raw_value, residue))
                    coefficient_checks += 1
                if all(value for _, value in vertices):
                    require(actual == raw, ("admissible actual/raw", d, row_n, p, actual, raw))
                    admissible_checks += 1
                else:
                    require(actual != raw, ("rho-zero hostile unchanged", d, row_n, p))
                profile_checks += 1

        for p in trial_prime_divisors(n):
            f_actual, g_actual = actual_hull(n - 1, d, p), actual_hull(n, d, p)
            actual_blocks = tuple(
                block for block in common_blocks(f_actual, g_actual, False) if block[0] > 0
            )
            formula_degrees, formula_blocks = ROUTE_A.divisor_barcode(n, p)
            require(actual_blocks == formula_blocks,
                    ("divisor actual/formula", d, p, actual_blocks, formula_blocks))
            require(integer_valuation(integer_coefficient(n, d, 0), p) == 0,
                    ("G constant unit", d, p))
            require(formula_degrees == ROUTE_B.divisor_barcode(n, p)[0],
                    ("small divisor route degrees", d, p))
            divisor_checks += 1
    return coefficient_checks, profile_checks, admissible_checks, divisor_checks


def degree_witness(blocks: tuple[tuple[int, int, int], ...], target: int):
    paths: dict[int, tuple[int, ...]] = {0: ()}
    for _, denominator, capacity in blocks:
        next_paths: dict[int, tuple[int, ...]] = {}
        for old, path in paths.items():
            for use in range(0, min(capacity, target - old) + 1, denominator):
                next_paths.setdefault(old + use, path + (use,))
        paths = next_paths
    return paths.get(target)


def compact_row(row):
    divisors = tuple(
        (p, blocks, len(pre), len(post), post if len(post) <= 16 else ())
        for p, blocks, pre, post in row[3]
    )
    rhos = tuple(
        (
            p,
            status,
            blocks,
            len(pre),
            len(post),
            next(((name, j, value)
                  for name, vertices in (("F", fv), ("G", gv))
                  for j, value in vertices if value == 0), None),
        )
        for p, status, fv, gv, blocks, pre, post in row[4]
    )
    return row[:3] + (divisors, rhos) + row[5:]


def full_prime_hostile_profile(route: Route, d: int, target: int):
    result = []
    for p in PRIMES_101:
        if d % p == 0:
            result.append((p, "divides_d", None))
            continue
        status, _, _, degrees, blocks = route.rho_barcode(d, p)
        result.append((p, status, target in set(degrees) if status == "admissible" else None, blocks))
    return tuple(result)


def main() -> None:
    for d in range(START, END + 1):
        require(exit_label_a(d) == exit_label_b(d), ("exit routes", d))
        require(trial_prime_divisors(d - 1) == spf_prime_divisors(d - 1), ("factor routes", d))

    coefficient_controls = direct_coefficient_controls()
    semantic_a = ROUTE_A.build()
    semantic_b = ROUTE_B.build()
    require(semantic_a == semantic_b, "two independent complete records")

    _, exits, residuals, _ = semantic_a
    source_shaped_semantic = (semantic_a[0], exits, residuals, PRIMES_47, PRIMES_101)
    source_semantic_sha = digest(source_shaped_semantic)
    require(source_semantic_sha == EXPECTED_SOURCE_SEMANTIC,
            ("candidate full semantic digest", source_semantic_sha))
    divisor_rows = tuple(row for row in residuals if row[1] == "divisor")
    rho_rows = tuple(row for row in residuals if row[1] == "rho")
    survivors = tuple(row for row in residuals if row[1] == "survivor")
    counts = (END - START + 1, len(exits), len(residuals), len(divisor_rows), len(rho_rows), len(survivors))
    require(counts == (51, 9, 42, 36, 6, 0), ("census", counts))
    rho_summary = tuple((row[0], row[2], row[5]) for row in rho_rows)
    require(rho_summary == EXPECTED_RHO_ROWS, ("rho summary", rho_summary))

    rows = {row[0]: row for row in residuals}
    anchor = rows[ANCHOR]
    require(anchor[:3] == (10001, "divisor", 5), ("anchor", anchor[:3]))
    require(tuple(trace[0] for trace in anchor[3]) == (2, 5), "anchor divisor order")
    p2, p5 = anchor[3]
    expected_anchor_packet = (
        256, 512, 768, 1024, 1280, 1536, 1792,
        8192, 8448, 8704, 8960, 9216, 9472, 9728, 9984,
    )
    require(p2[3] == expected_anchor_packet, ("anchor p2 packet", p2[3]))
    require(p5[2] == expected_anchor_packet and p5[3] == (), "anchor p5 closes")
    require(ROUTE_A.divisor_barcode(10000, 5)[0] == (0, 3125, 6250), "anchor p5 barcode")

    hostile = rows[9996]
    require(hostile[:3] == (9996, "rho", 19), ("hostile", hostile[:3]))
    require(hostile[5] == (3998,), ("hostile divisor packet", hostile[5]))
    require(trial_prime_divisors(9995) == (5, 1999), "hostile factorization")
    p5_div, p1999_div = hostile[3]
    require(degree_witness(p5_div[1], 3998) == (0, 0, 125, 3125, 748), "p5 witness")
    require(degree_witness(p1999_div[1], 3998) == (3998,), "p1999 witness")
    rho_map = {entry[0]: entry for entry in hostile[4]}
    require(degree_witness(rho_map[11][4], 3998) == (5, 0, 0, 1331, 2662, 0), "p11 witness")
    require(rho_map[11][4][0] == (0, 1, 5), "p11 zero slope")
    require(next((j for j, value in rho_map[13][2] if value == 0), None) == 0, "p13 F zero")
    require(rho_map[13][5] == rho_map[13][6] == (3998,), "p13 skip")
    require(rho_map[19][4] == ((2, 19, 171), (40, 361, 2888), (762, 6859, 6859)), "p19 blocks")
    require(degree_witness(rho_map[19][4], 3998) is None, "p19 kills")
    require(171 + 2888 == 3059 < 3998 < 6859, "p19 gap")

    nearest = min(rho_rows, key=lambda row: (abs(row[0] - ANCHOR), row[0]))
    require((nearest[0], abs(nearest[0] - ANCHOR)) == (9996, 5), "nearest hostile")
    require(all(rows[d][1] == "divisor" for d in range(9997, 10006)), "radius-four divisor closure")

    for row in residuals:
        d = row[0]
        require(row[3][0][2] == tuple(range(1, d - 1)), ("initial exact-support packet", d))
        for _, blocks, _, post in row[3]:
            require(all(numerator > 0 for numerator, _, _ in blocks), ("divisor positive slopes", d))
            require(0 not in post, ("positive target packet", d))
        for _, status, _, _, _, pre, post in row[4]:
            if status == "inadmissible":
                require(pre == post, ("inadmissible changed packet", d))

    full_hostile = full_prime_hostile_profile(ROUTE_A, 9996, 3998)
    first_killer = next(p for p, status, retained, *_ in full_hostile
                        if status == "admissible" and retained is False)
    require(first_killer == 19, ("first ordered rho killer", first_killer))

    # The window plus the inherited canonical r<=9998 result gives a contiguous
    # finite-exact extension through d=10030, i.e. r=d-2<=10028.
    require(not survivors and semantic_a[0] == (9980, 10030), "finite boundary extension")

    source_checker = Path(__file__).with_name(
        "factorial_d10001_local_rule_window_20260823.py"
    )
    require(EXPECTED_SOURCE_SEMANTIC in source_checker.read_text(encoding="utf-8"), "source digest pin visible")
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(Path(__file__).read_text()))), "assert node")

    compact_semantic = (
        semantic_a[0],
        exits,
        tuple(compact_row(row) for row in residuals),
        PRIMES_101,
    )
    semantic_sha = digest(compact_semantic)
    exit_histogram = tuple(sorted(Counter(label for _, label in exits).items()))
    later_killers = tuple(
        p for p, status, retained, *_ in full_hostile
        if status == "admissible" and retained is False
    )

    print("FACTORIAL_D10001_LOCAL_RULE_INDEPENDENT_AUDIT_20260823")
    print("status=FINITE-EXACT;candidate_computation=PASS;scope_wording=REPAIR")
    print(f"universe=d={START}..{END};exact_support_degrees=1..d-2;orientation=F=A_(d-2)^d,G=A_(d-1)^d")
    print(f"routes=({ROUTE_A.name};{ROUTE_B.name});imports_repository_modules=False")
    print(f"direct_integer_controls={coefficient_controls}")
    print(f"counts={counts};exit_histogram={exit_histogram}")
    print("rho_rows=" + repr(rho_summary).replace(" ", ""))
    print("anchor=" + repr(compact_row(anchor)).replace(" ", ""))
    print("hostile=" + repr(compact_row(hostile)).replace(" ", ""))
    print("hostile_full_bank_first_killer=19;admissible_killers=" + repr(later_killers).replace(" ", ""))
    print("route_record_equal=True;candidate_full_semantic_sha256=" + source_semantic_sha
          + ";compact_audit_semantic_sha256=" + semantic_sha)
    print("stronger_finite_consequence=with_THM3483_r_le_9998,the_window_closes_all_exact-support_rows_through_d=10030,hence_r_le_10028")
    print("scope_repair=refutes_completeness_of_the_THM3475_all-divisor_Newton-barcode_empty-intersection_criterion;does_not_refute_every_possible_divisor-place_argument")
    print("p19_scope=first_admissible_killer_in_the_declared_ordered_prime_bank;not_a_uniqueness_claim")
    print(f"gates={GATES}")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
