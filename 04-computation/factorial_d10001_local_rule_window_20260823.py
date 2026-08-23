#!/usr/bin/env python3
"""Exact local-window hostile for the d=10001 factorial divisor barcode.

Two formula-only routes are compared record-for-record.  The scan preserves
the exact-support packet 1..d-2, the ordered F=A_(d-2)^d/G=A_(d-1)^d pair,
positive-slope typing at divisor places, and zero-slope typing at admissible
nondivisor places.
"""

from __future__ import annotations

import ast
from collections import Counter
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import sys


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")

START = 9980
END = 10030
ANCHOR = 10001
EXPECTED_SEMANTIC_SHA256 = "f2d02ba3f51d6dd5b47a0976c775b307a3b73b64e01affaecf198f349a635213"

ROOT = Path(__file__).resolve().parents[1]
COMPUTATION = ROOT / "04-computation"
DEPENDENCIES = (
    (
        "factorial_adaptive_rho_block_6000_independent_audit.py",
        "d416cb2955fd745394cf1043ac8c2eba28a6a97beb264dd9cbe9919ed8c96724",
        "window_engine",
    ),
    (
        "factorial_all_divisor_digit_pair_compiler_independent_audit_thm3475.py",
        "9330ca1b991b9d5875779b9975fc88701ab36855a6527e1865e821e6cd3ea665",
        "window_divisor_referee",
    ),
    (
        "factorial_nondivisor_residue_digit_pair_compiler_thm3483.py",
        "9e37ead620f141617a9c6d51c182e09c034945793092e56e39fb061254662723",
        "window_rho_referee",
    ),
)
GATES = 0


def require(condition: bool, label: object) -> None:
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def normalized_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def load_module(name: str, expected: str, module_name: str):
    path = COMPUTATION / name
    require(normalized_sha256(path) == expected, ("dependency hash", name))
    spec = importlib.util.spec_from_file_location(module_name, path)
    require(spec is not None and spec.loader is not None, ("loader", name))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def digest(value: object) -> str:
    payload = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return hashlib.sha256(payload).hexdigest()


def exit_label_ref(divisor, d: int) -> str | None:
    tests = (
        ("d_prime", divisor.is_prime(d)),
        ("d_minus_1_prime_power", len(divisor.prime_divisors(d - 1)) == 1),
        ("d_minus_2_prime", divisor.is_prime(d - 2)),
        ("d_minus_3_prime", divisor.is_prime(d - 3)),
        ("d_minus_4_prime", divisor.is_prime(d - 4)),
        ("d_minus_5_prime", divisor.is_prime(d - 5)),
        ("d_minus_6_prime", divisor.is_prime(d - 6)),
    )
    return next((label for label, value in tests if value), None)


def normalize_rho_blocks(blocks: tuple[tuple[str, int, int], ...]) -> tuple[
    tuple[int, int, int], ...
]:
    answer = []
    for slope_text, capacity, denominator in blocks:
        slope = Fraction(slope_text)
        require(slope.denominator == denominator,
                ("rho denominator", slope_text, denominator))
        answer.append((slope.numerator, denominator, capacity))
    return tuple(answer)


def scan_referee(divisor, rho_module, primes: tuple[int, ...], d: int):
    """Second route, built directly from the two theorem referees."""
    n = d - 1
    packet = tuple(range(1, d - 1))
    divisor_trace = []
    for prime in divisor.prime_divisors(n):
        degrees, blocks, _, _ = divisor.pair_barcode(n, prime)
        pre = packet
        packet = tuple(value for value in packet if value in set(degrees))
        divisor_trace.append((prime, blocks, pre, packet))
        if not packet:
            return (d, "divisor", prime, tuple(divisor_trace), (), ())

    divisor_packet = packet
    rho_trace = []
    for prime in primes:
        if d % prime == 0:
            rho_trace.append((prime, "divides_d", (), (), (), packet, packet))
            continue
        f_hull = rho_module.raw_hull(n - 1, prime)
        g_hull = rho_module.raw_hull(n, prime)
        f_vertices = tuple(
            (j, rho_module.rho(n - 1, j, d, prime)) for j, _ in f_hull
        )
        g_vertices = tuple(
            (j, rho_module.rho(n, j, d, prime)) for j, _ in g_hull
        )
        pre = packet
        blocks = ()
        status = "inadmissible"
        if all(value for _, value in f_vertices + g_vertices):
            degrees, raw_blocks = rho_module.pair_degrees(f_hull, g_hull)
            blocks = normalize_rho_blocks(raw_blocks)
            packet = tuple(value for value in packet if value in set(degrees))
            status = "admissible"
        rho_trace.append(
            (prime, status, f_vertices, g_vertices, blocks, pre, packet)
        )
        if not packet:
            return (
                d,
                "rho",
                prime,
                tuple(divisor_trace),
                tuple(rho_trace),
                divisor_packet,
            )
    return (
        d,
        "survivor",
        None,
        tuple(divisor_trace),
        tuple(rho_trace),
        divisor_packet,
        packet,
    )


def build_referee_semantic(divisor, rho_module, engine):
    exits = []
    residuals = []
    for d in range(START, END + 1):
        label = exit_label_ref(divisor, d)
        if label:
            exits.append((d, label))
        else:
            residuals.append(scan_referee(divisor, rho_module, engine.PRIMES_101, d))
    return (
        (START, END),
        tuple(exits),
        tuple(residuals),
        engine.PRIMES_47,
        engine.PRIMES_101,
    )


def factorization(n: int) -> tuple[tuple[int, int], ...]:
    factors = []
    prime = 2
    while prime * prime <= n:
        if n % prime == 0:
            exponent = 0
            while n % prime == 0:
                exponent += 1
                n //= prime
            factors.append((prime, exponent))
        prime += 1
    if n > 1:
        factors.append((n, 1))
    return tuple(factors)


def degree_witness(
    blocks: tuple[tuple[int, int, int], ...], target: int
) -> tuple[int, ...] | None:
    paths: dict[int, tuple[int, ...]] = {0: ()}
    for _, denominator, capacity in blocks:
        next_paths: dict[int, tuple[int, ...]] = {}
        for old, path in paths.items():
            maximum = min(capacity, target - old)
            for use in range(0, maximum + 1, denominator):
                next_paths.setdefault(old + use, path + (use,))
        paths = next_paths
    return paths.get(target)


def first_zero(entry: tuple[object, ...]) -> tuple[str, int, int] | None:
    for label, vertices in (("F", entry[2]), ("G", entry[3])):
        for j, value in vertices:
            if value == 0:
                return (label, j, value)
    return None


def compact_row(row: tuple[object, ...]) -> tuple[object, ...]:
    divisor = tuple(
        (
            prime,
            blocks,
            len(pre),
            len(post),
            post if len(post) <= 16 else (),
        )
        for prime, blocks, pre, post in row[3]
    )
    rho = tuple(
        (prime, status, blocks, pre, post, first_zero(entry))
        for entry in row[4]
        for prime, status, _, _, blocks, pre, post in (entry,)
    )
    suffix = row[5:] if row[1] == "survivor" else (row[5],)
    return row[:3] + (divisor, rho) + suffix


def main() -> None:
    loaded = [load_module(*dependency) for dependency in DEPENDENCIES]
    engine, divisor, rho_module = loaded
    # First route: repeated-factor tables and determinant hulls.
    engine.START = START
    engine.END = END
    semantic, bank = engine.build_semantic()

    # Second route: Legendre weights and the two formula-only theorem referees.
    referee_semantic = build_referee_semantic(divisor, rho_module, engine)
    require(semantic == referee_semantic, "two-route full-record equality")
    require(digest(semantic) == EXPECTED_SEMANTIC_SHA256,
            "frozen window semantic digest")

    exits, divisor_rows, rho_rows, survivors, packets, killers, inadmissible = (
        engine.summarize(semantic)
    )
    counts = (
        END - START + 1,
        len(exits),
        len(semantic[2]),
        len(divisor_rows),
        len(rho_rows),
        len(survivors),
    )
    exit_histogram = tuple(sorted(Counter(label for _, label in exits).items()))
    rho_summary = tuple((row[0], row[2], row[5]) for row in rho_rows)
    require(counts == (51, 9, 42, 36, 6, 0), ("counts", counts))
    require(rho_summary == (
        (9982, 11, (2218, 6654)),
        (9986, 3, (3994, 5991)),
        (9988, 7, (6658,)),
        (9994, 11, (3331, 6662)),
        (9996, 19, (3998,)),
        (10030, 7, (3343, 6686)),
    ), ("rho rows", rho_summary))
    require(not survivors, ("survivors", survivors))

    rows = {row[0]: row for row in semantic[2]}

    # Routed positive control d=10001, independently rederived here.
    anchor = rows[ANCHOR]
    require(anchor[:3] == (10001, "divisor", 5), ("anchor status", anchor[:3]))
    require(tuple(entry[0] for entry in anchor[3]) == (2, 5), "anchor divisor order")
    p2_anchor, p5_anchor = anchor[3]
    expected_p2_packet = (
        256, 512, 768, 1024, 1280, 1536, 1792,
        8192, 8448, 8704, 8960, 9216, 9472, 9728, 9984,
    )
    require(p2_anchor[3] == expected_p2_packet, "anchor p=2 packet")
    require(p5_anchor[2] == expected_p2_packet and not p5_anchor[3],
            "anchor p=5 closure")
    anchor_p5_degrees = divisor.pair_barcode(10000, 5)[0]
    require(anchor_p5_degrees == (0, 3125, 6250), "anchor p=5 local barcode")

    # Minimal nearby hostile to completeness of the THM-3475 barcode test.
    hostile = rows[9996]
    require(hostile[:3] == (9996, "rho", 19), ("hostile status", hostile[:3]))
    require(hostile[5] == (3998,), ("hostile divisor packet", hostile[5]))
    require(tuple(entry[0] for entry in hostile[3]) == (5, 1999),
            "hostile divisor order")
    p5_divisor, p1999_divisor = hostile[3]
    require(degree_witness(p5_divisor[1], 3998) == (0, 0, 125, 3125, 748),
            "p=5 divisor witness")
    require(degree_witness(p1999_divisor[1], 3998) == (3998,),
            "p=1999 divisor witness")

    rho_entries = {entry[0]: entry for entry in hostile[4]}
    require(tuple((entry[0], entry[1]) for entry in hostile[4]) == (
        (2, "divides_d"),
        (3, "divides_d"),
        (5, "admissible"),
        (7, "divides_d"),
        (11, "admissible"),
        (13, "inadmissible"),
        (17, "divides_d"),
        (19, "admissible"),
    ), "hostile rho status trace")
    p11 = rho_entries[11]
    p13 = rho_entries[13]
    p19 = rho_entries[19]
    require(degree_witness(p11[4], 3998) == (5, 0, 0, 1331, 2662, 0),
            "zero-slope p=11 preservation witness")
    require(p11[4][0] == (0, 1, 5), "zero slope retained and typed")
    require(first_zero(p13) == ("F", 0, 0), "p=13 inadmissible skip")
    require(p13[5] == p13[6] == (3998,), "inadmissible packet unchanged")
    require(p19[4] == (
        (2, 19, 171), (40, 361, 2888), (762, 6859, 6859),
    ), "p=19 hostile blocks")
    require(degree_witness(p19[4], 3998) is None, "p=19 kills degree 3998")
    require(171 + 2888 == 3059 < 3998 < 6859,
            "p=19 exact barcode gap")

    nearest = min(rho_rows, key=lambda row: (abs(row[0] - ANCHOR), row[0]))
    require(nearest[0] == 9996 and abs(nearest[0] - ANCHOR) == 5,
            ("nearest divisor-rule hostile", nearest[0]))
    require(all(
        row[1] == "divisor"
        for d, row in rows.items()
        if 0 < abs(d - ANCHOR) < 5
    ), "all closer residual rows divisor-close")

    # Typing gates over the full residual census.
    for row in semantic[2]:
        d = row[0]
        require(row[3][0][2] == tuple(range(1, d - 1)),
                ("exact-support initial packet", d))
        for _, blocks, _, post in row[3]:
            require(all(numerator > 0 for numerator, _, _ in blocks),
                    ("positive divisor slopes", d))
            require(0 not in post, ("positive exact-support packet", d))
        for entry in row[4]:
            if entry[1] == "inadmissible":
                require(entry[5] == entry[6], ("inadmissible unchanged", d, entry[0]))

    positive_control = divisor.pair_barcode(6, 3)
    require(positive_control[:2] == ((0, 3), ((2, 3, 3),)),
            "positive divisor control N=6,p=3")
    require(not any(isinstance(node, ast.Assert)
                    for node in ast.walk(ast.parse(Path(__file__).read_text()))),
            "assert node")

    row_statuses = []
    exit_map = dict(exits)
    for d in range(START, END + 1):
        if d in exit_map:
            row_statuses.append((d, "exit", exit_map[d]))
        else:
            row = rows[d]
            row_statuses.append(
                (d, row[1], row[2], row[5] if row[1] == "rho" else ())
            )

    print("FACTORIAL_D10001_LOCAL_DIVISOR_RULE_HOSTILE_20260823")
    print("status=FINITE-EXACT;THM3475_all_divisor_barcode_completeness=REFUTED;no_general_d_force")
    print(
        f"universe=d={START}..{END};orientation=F=A_(d-2)^d,G=A_(d-1)^d;"
        "exact_support=degrees_1..d-2;divisor_slopes=positive_only;"
        "rho_zero_slope=retained_distinct_from_coordinate_root"
    )
    print("dependencies=" + repr(tuple((name, sha) for name, sha, _ in DEPENDENCIES)).replace(" ", ""))
    print(f"counts={counts};exit_histogram={exit_histogram}")
    print("row_statuses=" + repr(tuple(row_statuses)).replace(" ", ""))
    print("rho_rows=" + repr(rho_summary).replace(" ", ""))
    print(
        "anchor_d10001="
        + repr((factorization(10000), compact_row(anchor), anchor_p5_degrees)).replace(" ", "")
    )
    print(
        "minimal_hostile_d9996="
        + repr((factorization(9995), compact_row(hostile))).replace(" ", "")
    )
    print(
        "hostile_witnesses="
        + repr((
            ("p5_divisor", (0, 0, 125, 3125, 748)),
            ("p1999_divisor", (3998,)),
            ("p11_rho_with_zero_slope", (5, 0, 0, 1331, 2662, 0)),
            ("p13_skip", ("F", 0, 0)),
            ("p19_gap", 3059, 3998, 6859),
        )).replace(" ", "")
    )
    print("positive_control_N6_p3=" + repr(positive_control[:2]).replace(" ", ""))
    print(f"semantic_sha256={EXPECTED_SEMANTIC_SHA256};cross_record_equal=True")
    print(f"gates={GATES}")
    print(
        "interpretation=d10001_closure_is_digit_sensitive_to_N=10000;"
        "five_steps_away_N=9995_leaves_degree3998_at_all_divisor_places;"
        "p19_is_the_first_admissible_killer_in_the_declared_bank"
    )
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
