#!/usr/bin/env python3
"""Exact THM-3640 companion: positive two-cube slopes through n=401.

The script extends the independently audited THM-3594 decision pipeline.  A
prime obstruction is only an exclusion gate; every local survivor is decided
by complete generalized-Pell class enumeration and the full unit orbit modulo
the compiler address 2n.  One positive invariant-cone seed is then compiled
for every admissible slope.
"""

from __future__ import annotations

import ast
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
CANON = ROOT / "01-canon/theorems"
COMPUTATION = ROOT / "04-computation"
RESULTS = ROOT / "05-knowledge/results"
PARENT_FILES = (
    (
        "theorem",
        CANON / "THM-3594-berggren-positive-cube-slope-atlas-through-201.md",
        "03f4bff0642cd0f84553f21e46ffc600c102c41d8ef0a9aeaaf974fdfe28881c",
    ),
    (
        "script",
        COMPUTATION / "berggren_positive_cube_slope_atlas_201_thm3594.py",
        "d3466965aa2ccaca81c3f4f4369dbca8f54f03c8472c7d5094459e38ea8977bd",
    ),
    (
        "output",
        RESULTS / "berggren_positive_cube_slope_atlas_201_thm3594.out",
        "dd86e7548773cf59113969bb7f0f7d99d345f98adbdf1de8b722533d15a02856",
    ),
)

LIMIT = 401
PRIME_LIMIT = 997
EXPECTED_CANDIDATES = 8195
EXPECTED_SURVIVORS = 104
EXPECTED_ADMISSIBLE = (
    (14, 23), (26, 29), (26, 47), (98, 101),
    (86, 125), (86, 149), (110, 149), (122, 149),
    (86, 167), (98, 167), (110, 167), (146, 167),
    (170, 191), (110, 197), (182, 197),
    (122, 215), (182, 215), (146, 239), (182, 239),
    (158, 263), (218, 263), (230, 263),
    (146, 269), (230, 269), (242, 269),
    (170, 293), (182, 293), (194, 293), (254, 293), (278, 293),
    (158, 311), (278, 317), (218, 335), (254, 335),
    (182, 359), (194, 359), (206, 359), (278, 365),
    (194, 383), (338, 383), (350, 383), (206, 389),
)
EXPECTED_DIGESTS = {
    "obstruction_histogram": "c430037d87fe8b94adb4a1b73ec1e3adf813215d50a7dcd41579581c9acbc990",
    "survivors": "2b1e24ba10e59deac190e1b4d20a0afe3792c138ecc84e06fa22e1671448ff8e",
    "orbit_ledger": "5390c4b137d3fdd90d90212a7b3c26cac2b819ff1fd1b7e0950149cf6e593265",
    "seed_ledger": "30af440ae6817b71045a3351d8480aef55c5ce62ac46580659738266f74fd6e5",
    "semantic": "879a912b16752dc8b97393affeff61387eedc45d23c490a7c9ce5f54c1836374",
}

CHECKS = 0


def require(condition: bool, payload: object) -> None:
    global CHECKS
    if condition is not True:
        raise RuntimeError(payload)
    CHECKS += 1


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def load_parent():
    observed = tuple((label, lf_sha256(path)) for label, path, _expected in PARENT_FILES)
    expected = tuple((label, expected) for label, _path, expected in PARENT_FILES)
    require(observed == expected, ("parent drift", observed, expected))
    path = PARENT_FILES[1][1]
    spec = importlib.util.spec_from_file_location("thm3640_parent", path)
    require(spec is not None and spec.loader is not None, "parent loader")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module, observed


def candidates() -> tuple[tuple[int, int], ...]:
    return tuple(
        (m, n)
        for n in range(3, LIMIT + 1, 2)
        for m in range(n // 2 + 1, n)
        if m % 2 == 0 and gcd(m, n) == 1
    )


def main() -> None:
    parent, parent_hashes = load_parent()
    primes = parent.SCREEN.primes_up_to(PRIME_LIMIT)
    require(len(primes) == 168 and primes[0] == 2 and primes[-1] == 997,
            ("prime table", len(primes), primes[:1], primes[-1:]))
    slopes = candidates()
    require(len(slopes) == EXPECTED_CANDIDATES, ("candidate count", len(slopes)))

    obstruction = {
        slope: parent.SCREEN.first_obstruction(*slope, primes) for slope in slopes
    }
    survivors = tuple(slope for slope in slopes if obstruction[slope] is None)
    obstruction_histogram = tuple(sorted(Counter(
        value for value in obstruction.values() if value is not None
    ).items()))
    require(len(survivors) == EXPECTED_SURVIVORS,
            ("survivor count", len(survivors)))
    require(sum(count for _prime, count in obstruction_histogram)
            == EXPECTED_CANDIDATES - EXPECTED_SURVIVORS,
            "obstruction partition")
    require(all((n - m) % 3 == 0 and (2 * m - n) % 3 != 0
                for m, n in survivors), "survivor u/v residue law")

    orbit_ledger = []
    seed_ledger = []
    admissible = []
    for slope in survivors:
        m, n = slope
        k, constant = parent.PELL.conic_data(m, n)
        require(k > 0 and constant > 0 and parent.PELL.isqrt(k) ** 2 != k,
                (slope, "positive nonsquare conic"))
        representatives, audit = parent.PELL.lmm_classes(k, constant)
        units = parent.PELL.pell_pm_one(k, 1)
        require(len(units) == 1, (slope, "fundamental unit", units))
        unit = units[0]
        require(unit[0] > 1 and unit[1] > 0
                and unit[0] ** 2 - k * unit[1] ** 2 == 1,
                (slope, "unit norm"))
        orbits = tuple(
            parent.PELL.unit_orbit(representative, k, unit, 2 * n)
            for representative in representatives
        )
        good = tuple(
            index for index, (_period, hits) in enumerate(orbits) if hits
        )
        orbit_ledger.append((
            slope, k, constant, representatives, unit, orbits, good,
            audit["root_branches"], audit["negative_pell_fallbacks"],
        ))
        if good:
            admissible.append(slope)
            seed_ledger.append(parent.positive_seed(
                slope, k, constant, representatives, unit, orbits, good
            ))

    admissible = tuple(admissible)
    require(admissible == EXPECTED_ADMISSIBLE,
            ("admissible list", admissible))
    require(len(seed_ledger) == len(EXPECTED_ADMISSIBLE) == 42,
            ("seed count", len(seed_ledger)))
    require(tuple(slope for slope in admissible if slope[1] <= 201)
            == parent.EXPECTED_ADMISSIBLE, "parent atlas prefix")
    require(len(tuple(slope for slope in admissible if slope[1] > 201)) == 27,
            "new-ray count")
    require(all(len(row[6]) % 2 == 0 for row in orbit_ledger),
            "parity-admissible classes occur in pairs")

    observed_digests = {
        "obstruction_histogram": digest(obstruction_histogram),
        "survivors": digest(survivors),
        "orbit_ledger": digest(orbit_ledger),
        "seed_ledger": digest(seed_ledger),
    }
    semantic_record = (
        LIMIT, PRIME_LIMIT, parent_hashes, len(primes), slopes,
        obstruction_histogram, survivors, tuple(orbit_ledger),
        admissible, tuple(seed_ledger), observed_digests,
    )
    observed_digests["semantic"] = digest(semantic_record)
    for label, observed in observed_digests.items():
        expected = EXPECTED_DIGESTS[label]
        if expected != "TO_BE_PINNED":
            require(observed == expected, (label, observed, expected))

    source = Path(__file__).resolve()
    source_bytes = source.read_bytes()
    require(b"\r\n" not in source_bytes, "source raw LF")
    require(not any(isinstance(node, ast.Assert)
                    for node in ast.walk(ast.parse(source_bytes.decode("utf-8")))),
            "Python assert node present")

    transformed = tuple((n - m, 2 * m - n) for m, n in admissible)
    class_counts = tuple(len(row[3]) for row in orbit_ledger)
    good_counts = tuple(len(row[6]) for row in orbit_ledger)
    period_sets = tuple(tuple(sorted({period for period, _hits in row[5]}))
                        for row in orbit_ledger)
    print("== THM-3640 Berggren positive two-cube slope atlas through 401 ==")
    print(f"parent_sha256_lf={parent_hashes}")
    print(f"universe=(limit={LIMIT},candidates={len(slopes)},primes={len(primes)},prime_limit={primes[-1]})")
    print(f"screen=(excluded={len(slopes)-len(survivors)},survivors={len(survivors)},histogram_digest={observed_digests['obstruction_histogram']})")
    print(f"survivors_sha256={observed_digests['survivors']}")
    print("admissible=" + repr(admissible).replace(" ", ""))
    print("new_beyond_201=" + repr(tuple(s for s in admissible if s[1] > 201)).replace(" ", ""))
    print("transformed_uv=" + repr(transformed).replace(" ", ""))
    print(f"class_histogram={tuple(sorted(Counter(class_counts).items()))}")
    print(f"good_class_histogram={tuple(sorted(Counter(good_counts).items()))}")
    print(f"period_set_histogram={tuple(sorted(Counter(period_sets).items()))}")
    print(f"orbit_ledger_sha256={observed_digests['orbit_ledger']}")
    print(f"seed_ledger_sha256={observed_digests['seed_ledger']}")
    print(f"semantic_sha256={observed_digests['semantic']}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print(f"CHECKS={CHECKS}")
    print("scope=3<=n<=401 only;no slope parametrization,density,or asymptotic")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
