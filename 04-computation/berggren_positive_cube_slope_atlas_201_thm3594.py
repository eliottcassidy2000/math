#!/usr/bin/env python3
"""Exact Berggren positive-two-cube slope atlas through denominator 201.

The verifier combines the finite prime screen of THM-3547 with the complete
generalized-Pell class and unit-orbit method of THM-3580.  Passing the prime
screen is not treated as solubility: every survivor is decided by its full
LMM class list modulo the exact compiler address 2*n.  Every admissible slope
is then moved into the positive invariant cone and compiled to a positive
pair of distinct cubes.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
import importlib.util
import json
from math import gcd
from pathlib import Path
import sys

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


ROOT = Path(__file__).resolve().parents[1]
SCREEN_PATH = ROOT / "04-computation/berggren_positive_cube_slope_atlas_101_kps_s183.py"
PELL_PATH = ROOT / "04-computation/berggren_positive_cube_slope_completion_101_thm3580.py"
EXPECTED_PARENT_HASHES = (
    "bdb8cd4fbd14235ee144d80c4766aed5117097321cae8ec87b1aad1ffff1212c",
    "09cd67a5bf0c2c0368e3e9ca4350d2cbb4c856ccc3bca4afa21254aaef26e184",
)

LIMIT = 201
PRIME_LIMIT = 499
EXPECTED_CANDIDATE_COUNT = 2072
EXPECTED_SURVIVORS = (
    (14, 23), (26, 29), (26, 47), (38, 47), (38, 53),
    (50, 53), (50, 71), (62, 95), (74, 95), (74, 101),
    (98, 101), (86, 125), (98, 125), (86, 149), (110, 149),
    (122, 149), (86, 167), (98, 167), (110, 167), (146, 167),
    (98, 173), (110, 173), (122, 173), (170, 173), (98, 191),
    (110, 191), (170, 191), (110, 197), (146, 197), (182, 197),
)
EXPECTED_ADMISSIBLE = (
    (14, 23), (26, 29), (26, 47), (98, 101),
    (86, 125), (86, 149), (110, 149), (122, 149),
    (86, 167), (98, 167), (110, 167), (146, 167),
    (170, 191), (110, 197), (182, 197),
)

# (number of LMM classes, common orbit-period tuple, good class indices)
EXPECTED_CLASS_ORBITS = {
    (14, 23): (3, (11,) * 3, (0, 1)),
    (26, 29): (6, (10,) * 6, (0, 3)),
    (26, 47): (9, (46,) * 9, (0, 1, 2, 3, 6, 7)),
    (38, 47): (0, (), ()),
    (38, 53): (9, (2,) * 9, ()),
    (50, 53): (0, (), ()),
    (50, 71): (0, (), ()),
    (62, 95): (1, (60,), ()),
    (74, 95): (1, (10,), ()),
    (74, 101): (6, (34,) * 6, ()),
    (98, 101): (6, (34,) * 6, (0, 1, 2, 3)),
    (86, 125): (3, (50,) * 3, (0, 1)),
    (98, 125): (0, (), ()),
    (86, 149): (6, (50,) * 6, (1, 2)),
    (110, 149): (6, (25,) * 6, (0, 3)),
    (122, 149): (5, (150,) * 5, (0, 1, 2, 3)),
    (86, 167): (6, (83,) * 6, (0, 1, 2, 3)),
    (98, 167): (2, (83,) * 2, (0, 1)),
    (110, 167): (12, (166,) * 12, tuple(range(8))),
    (146, 167): (6, (83,) * 6, (0, 1, 2, 3)),
    (98, 173): (0, (), ()),
    (110, 173): (0, (), ()),
    (122, 173): (0, (), ()),
    (170, 173): (6, (58,) * 6, ()),
    (98, 191): (0, (), ()),
    (110, 191): (0, (), ()),
    (170, 191): (6, (95,) * 6, (0, 1, 2, 3)),
    (110, 197): (4, (198,) * 4, tuple(range(4))),
    (146, 197): (3, (66,) * 3, ()),
    (182, 197): (6, (11,) * 6, (1, 2)),
}

# Pinned after the first exact construction; both are unconditional gates in
# the maintained version.
EXPECTED_SEED_LEDGER_SHA256 = "8ff45843ce97b3ea6bdeeeab691bf1c06b53f7b58ea82abb1fa28edc14952a40"
EXPECTED_SEMANTIC_SHA256 = "2389ac8015a1f1394f938fe5ebf87933175dcdc7b4480901bdeba2c3a6b57de0"


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def load_module(path: Path, name: str, expected_hash: str):
    observed = lf_sha256(path)
    require(observed == expected_hash, (name, "parent hash", observed, expected_hash))
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, (name, "loader"))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


SCREEN = load_module(SCREEN_PATH, "thm3594_screen_parent", EXPECTED_PARENT_HASHES[0])
PELL = load_module(PELL_PATH, "thm3594_pell_parent", EXPECTED_PARENT_HASHES[1])


def candidates() -> tuple[tuple[int, int], ...]:
    return tuple(
        (m, n)
        for n in range(3, LIMIT + 1, 2)
        for m in range(n // 2 + 1, n)
        if m % 2 == 0 and gcd(m, n) == 1
    )


def pell_multiply(left: tuple[int, int], right: tuple[int, int], k: int):
    a, b = left
    c, d = right
    return a * c + k * b * d, a * d + b * c


def pell_power(unit: tuple[int, int], exponent: int, k: int):
    require(exponent >= 0, "negative Pell exponent")
    answer = (1, 0)
    base = unit
    while exponent:
        if exponent & 1:
            answer = pell_multiply(answer, base, k)
        base = pell_multiply(base, base, k)
        exponent >>= 1
    return answer


def compile_seed(m: int, n: int, w: int, h: int) -> tuple[int, ...]:
    require(h % n == 0, ((m, n), "h address", h % n))
    u = h // n
    require(w > 0 and u > 0 and w % 2 == u % 2 == 1, ((m, n), "odd positive"))
    d = n * n * u * u + 2
    e = u * w
    a = m * n * n * u**3 + (2 * m + n) * u
    q = m * m * n * n * u**4 + 2 * m * (m + n) * u * u + 1
    require(0 < e < d, ((m, n), "positive chamber"))
    require((d + e) % 2 == (d - e) % 2 == 0, ((m, n), "cube parity"))
    x, y = (d + e) // 2, (d - e) // 2
    r = (a - 1) // 2
    require(a % 2 == 1 and r >= 1, ((m, n), "Berggren depth"))
    require(a * a + 2 == d * q, ((m, n), "norm factor"))
    require(4 * q - d * d == 3 * e * e, ((m, n), "cube discriminant"))
    require(x > y > 0, ((m, n), "positive distinct cubes"))
    require(x**3 + y**3 == a * a + 2 == (2 * r + 1) ** 2 + 2,
            ((m, n), "cube identity"))
    return w, h, u, d, e, a, q, x, y, r


def positive_seed(
    slope: tuple[int, int],
    k: int,
    constant: int,
    representatives: tuple[tuple[int, int], ...],
    unit: tuple[int, int],
    orbit_rows: tuple[tuple[int, tuple[int, ...]], ...],
    good: tuple[int, ...],
):
    m, n = slope
    class_index = good[0]
    period, hits = orbit_rows[class_index]
    exponent = hits[0]
    phase_unit = pell_power(unit, exponent, k)
    w, h = pell_multiply(representatives[class_index], phase_unit, k)
    # Independent sign and conjugation are part of the complete LMM orbit.
    w, h = abs(w), abs(h)
    period_unit = pell_power(unit, period, k)
    cycles = 0
    while not (0 < w < n * h):
        w, h = pell_multiply((w, h), period_unit, k)
        cycles += 1
        require(cycles <= 4, (slope, "positive cone took too many periods"))
    require(w * w - k * h * h == constant, (slope, "seed norm"))
    require(w % 2 == 1 and h % (2 * n) == n, (slope, "compiler residue"))
    lower = n * h - w
    upper = n * w - k * h
    require(lower > 0 and upper > 0, (slope, "invariant cone"))

    # One period step preserves the address and the positive cone exactly.
    next_w, next_h = pell_multiply((w, h), period_unit, k)
    next_lower = n * next_h - next_w
    next_upper = n * next_w - k * next_h
    require(next_lower == period_unit[0] * lower + period_unit[1] * upper,
            (slope, "lower cone recurrence"))
    require(next_upper == period_unit[0] * upper + k * period_unit[1] * lower,
            (slope, "upper cone recurrence"))
    require(next_lower > 0 and next_upper > 0 and next_h > h,
            (slope, "strict positive ray"))
    require(next_w % 2 == 1 and next_h % (2 * n) == n,
            (slope, "period address"))
    compiled = compile_seed(m, n, w, h)
    compile_seed(m, n, next_w, next_h)
    return class_index, exponent, period, cycles, period_unit, compiled


def digest_json(value: object) -> str:
    return sha256(json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")).hexdigest()


def main() -> None:
    primes = SCREEN.primes_up_to(PRIME_LIMIT)
    require(primes[0] == 2 and primes[-1] == 499 and len(primes) == 95,
            "prime table")
    slopes = candidates()
    require(len(slopes) == EXPECTED_CANDIDATE_COUNT, "candidate count")
    obstruction = {
        slope: SCREEN.first_obstruction(*slope, primes) for slope in slopes
    }
    survivors = tuple(slope for slope in slopes if obstruction[slope] is None)
    require(survivors == EXPECTED_SURVIVORS, "screen survivor list")
    obstruction_counts = tuple(sorted(Counter(
        value for value in obstruction.values() if value is not None
    ).items()))
    require(sum(count for _prime, count in obstruction_counts)
            == EXPECTED_CANDIDATE_COUNT - len(EXPECTED_SURVIVORS),
            "screen exclusion count")

    orbit_ledger = []
    seed_ledger = []
    admissible = []
    for slope in survivors:
        m, n = slope
        k, constant = PELL.conic_data(m, n)
        require(k > 0 and constant > 0 and PELL.isqrt(k) ** 2 != k,
                (slope, "positive nonsquare conic"))
        representatives, audit = PELL.lmm_classes(k, constant)
        unit_rows = PELL.pell_pm_one(k, 1)
        require(len(unit_rows) == 1, (slope, "fundamental unit"))
        unit = unit_rows[0]
        require(unit[0] > 1 and unit[1] > 0
                and unit[0] ** 2 - k * unit[1] ** 2 == 1,
                (slope, "unit norm"))
        orbit_rows = tuple(
            PELL.unit_orbit(rep, k, unit, 2 * n) for rep in representatives
        )
        periods = tuple(period for period, _hits in orbit_rows)
        good = tuple(i for i, (_period, hits) in enumerate(orbit_rows) if hits)
        expected_classes, expected_periods, expected_good = EXPECTED_CLASS_ORBITS[slope]
        require(len(representatives) == expected_classes, (slope, "class count"))
        require(periods == expected_periods, (slope, "orbit periods"))
        require(good == expected_good, (slope, "good classes"))
        orbit_ledger.append((slope, k, constant, len(representatives), periods, good,
                             audit["root_branches"], audit["negative_pell_fallbacks"]))
        if good:
            admissible.append(slope)
            seed_ledger.append(positive_seed(
                slope, k, constant, representatives, unit, orbit_rows, good
            ))

    require(tuple(admissible) == EXPECTED_ADMISSIBLE, "admissible list")
    require(len(seed_ledger) == 15, "positive ray count")
    seed_digest = digest_json(seed_ledger)
    semantic = digest_json((
        LIMIT, PRIME_LIMIT, EXPECTED_PARENT_HASHES, slopes, obstruction_counts,
        survivors, tuple(orbit_ledger), tuple(admissible), seed_digest,
    ))
    require(seed_digest == EXPECTED_SEED_LEDGER_SHA256,
            ("seed ledger digest", seed_digest, EXPECTED_SEED_LEDGER_SHA256))
    require(semantic == EXPECTED_SEMANTIC_SHA256,
            ("semantic digest", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== THM-3594 Berggren positive two-cube slope atlas through 201 ==")
    print(f"universe=(limit={LIMIT},candidate_count={len(slopes)})")
    print(f"prime_screen=(limit={PRIME_LIMIT},prime_count={len(primes)},"
          f"excluded={len(slopes)-len(survivors)},survivors={len(survivors)})")
    print("screen_survivors=" + repr(survivors).replace(" ", ""))
    print("admissible_slopes=" + repr(tuple(admissible)).replace(" ", ""))
    print("new_beyond_101=" + repr(tuple(s for s in admissible if s[1] > 101)).replace(" ", ""))
    for row, seed in zip((row for row in orbit_ledger if row[0] in admissible), seed_ledger):
        slope, k, constant, classes, periods, good, branches, fallbacks = row
        class_index, exponent, period, cycles, period_unit, compiled = seed
        print(
            f"row={slope[0]}/{slope[1]}|K={k}|C={constant}|classes={classes}|"
            f"periods={periods}|good={good}|branches={branches}|fallbacks={fallbacks}|"
            f"seed_class={class_index}|phase={exponent}|period={period}|cycles={cycles}|"
            f"period_unit_bits=({period_unit[0].bit_length()},{period_unit[1].bit_length()})|"
            f"seed_bits=({compiled[0].bit_length()},{compiled[2].bit_length()},"
            f"{compiled[7].bit_length()},{compiled[8].bit_length()})|"
            f"seed_sha256={digest_json(compiled)}"
        )
    print(f"obstruction_counts={obstruction_counts}")
    print(f"seed_ledger_sha256={seed_digest}")
    print(f"semantic_sha256={semantic}")
    print(f"script_sha256_lf={lf_sha256(Path(__file__).resolve())}")
    print("status=FINITE-EXACT+VERIFIED-EXACT;15 positive infinite rays")
    print("scope=n<=201 only;no unbounded parametrization;no density/asymptotic")
    print("PASS")


if __name__ == "__main__":
    main()
