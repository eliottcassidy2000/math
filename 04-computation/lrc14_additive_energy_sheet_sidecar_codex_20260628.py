#!/usr/bin/env python3
"""HYP-3425: additive-energy needs a sheet sidecar on the LRC14 covering floor.

This scout sits directly under HYP-3424's transfer rule

    additive / multiplicative energy
    -> low-SPEC penalty or named odd phase debt.

The question here is whether a raw additive-energy scalar can already see the
covering-floor bias, or whether the HYP-3140 sheet packet is genuinely
load-bearing.

The answer from the exact rows below is negative: the old energy scalars
collide across covering packets with different Rprime behavior.  The smallest
useful repair is an energy-plus-sheet packet, not another scalar.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import gcd
from pathlib import Path
from random import Random
import importlib.util
import sys


F = Fraction
MOD = 27


def load_fiber_pgf_module():
    path = Path(__file__).with_name("lrc14_fiber_pgf_certificate_codex_s273.py")
    spec = importlib.util.spec_from_file_location("hyp3140_fiber_pgf", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("could not load HYP-3140 fiber PGF module")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


HYP3140 = load_fiber_pgf_module()


def add_energy(speeds: tuple[int, ...], mod: int = MOD) -> tuple[int, int]:
    residues = [speed % mod for speed in speeds]
    counts = Counter((a + b) % mod for a, b in combinations(residues, 2))
    return sum(value * value for value in counts.values()), counts.get(0, 0)


def v2(n: int) -> int:
    depth = 0
    while n and n % 2 == 0:
        n //= 2
        depth += 1
    return depth


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for speed in speeds:
        g = gcd(g, speed)
    return g == 1


def is_covering(speeds: tuple[int, ...]) -> bool:
    return all(any(speed % q == 0 for speed in speeds) for q in range(2, 15))


def random_covering(rng: Random, max_speed: int = 90) -> tuple[int, ...]:
    for _ in range(5000):
        speeds: set[int] = set()
        for q in rng.sample(range(2, 15), 13):
            if not any(speed % q == 0 for speed in speeds):
                choices = [q * k for k in range(1, max_speed // q + 1)]
                speeds.add(rng.choice(choices))
        while len(speeds) < 13:
            speeds.add(rng.randint(1, max_speed))
        row = tuple(sorted(speeds))
        if len(row) == 13 and primitive(row) and is_covering(row):
            return row
    raise RuntimeError("failed to generate covering row")


def split_row(row: tuple[int, ...]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    r_speeds = tuple(speed for speed in row if speed % 14)
    q_speeds = tuple(sorted(speed // 14 for speed in row if speed % 14 == 0))
    return r_speeds, q_speeds


def corr(xs: list[float], ys: list[float]) -> float:
    x_mean = sum(xs) / len(xs)
    y_mean = sum(ys) / len(ys)
    num = sum((x - x_mean) * (y - y_mean) for x, y in zip(xs, ys))
    den_x = sum((x - x_mean) ** 2 for x in xs)
    den_y = sum((y - y_mean) ** 2 for y in ys)
    if den_x == 0 or den_y == 0:
        return 0.0
    return num / (den_x * den_y) ** 0.5


def fmt_frac(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


CURATED_ROWS: dict[str, tuple[int, ...]] = {
    "canonical_r1_drop12": tuple([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14]),
    "canonical_r2_drop13": tuple([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 28]),
    "canonical_r3": tuple([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 28, 42]),
    "canonical_r4": tuple([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 14, 28, 42, 56]),
    "canonical_r5": tuple([1, 2, 3, 4, 5, 6, 7, 8, 9, 14, 28, 42, 56, 70]),
    "canonical_r6": tuple([1, 2, 3, 4, 5, 6, 7, 8, 14, 28, 42, 56, 70, 84]),
    "spread_far": tuple([1, 2, 3, 4, 5, 6, 7, 8, 9, 14, 42, 70]),
    "gcd_even": tuple([2, 4, 6, 8, 10, 12, 14, 42, 70]),
    "gcd_triples": tuple([3, 6, 9, 12, 14, 28, 70]),
    "coprime_sparse": tuple([5, 9, 11, 13, 14, 28, 42]),
    "covering_AP_with_84": tuple(list(range(1, 12)) + [13, 84]),
    "covering_AP_with_12_and_84": tuple(list(range(1, 13)) + [84]),
}


def curated_metrics() -> dict[str, dict[str, object]]:
    out: dict[str, dict[str, object]] = {}
    for name, row in CURATED_ROWS.items():
        r_speeds, q_speeds = split_row(row)
        fiber = HYP3140.fiber_pgf_metrics(r_speeds, q_speeds)
        odd = tuple(speed for speed in row if speed % 2)
        even = tuple(speed for speed in row if speed % 2 == 0)
        full_e, full_shell = add_energy(row)
        r_e, _r_shell = add_energy(r_speeds)
        odd_e, _odd_shell = add_energy(odd) if len(odd) >= 2 else (0, 0)
        even_e, _even_shell = add_energy(even) if len(even) >= 2 else (0, 0)
        out[name] = {
            "row": row,
            "Rprime": fiber["rprime_fiber"],
            "delta": fiber["q_mean_sheets"] - fiber["mean_sheets"],
            "fullE": full_e,
            "fullShell": full_shell,
            "RE": r_e,
            "oddE": odd_e,
            "evenE": even_e,
            "Qsize": len(q_speeds),
            "q_zero_mass": fiber["q_zero_mass"],
            "q_range_lo": fiber["min_q_count"],
            "q_range_hi": fiber["max_q_count"],
            "sumv2_even": sum(v2(speed) for speed in even),
            "sign": 1 if fiber["q_mean_sheets"] > fiber["mean_sheets"] else -1,
        }
    return out


def print_collision_audit(metrics: dict[str, dict[str, object]]) -> None:
    collision_specs = (
        ("fullE",),
        ("RE",),
        ("oddE",),
        ("RE", "oddE"),
    )
    print("EXACT COLLISION AUDIT ON THE CURATED COVERING BANK")
    for spec in collision_specs:
        print(f"  key={spec}")
        classes: dict[tuple[object, ...], list[str]] = {}
        for name, data in metrics.items():
            key = tuple(data[field] for field in spec)
            classes.setdefault(key, []).append(name)
        found = False
        for key, names in classes.items():
            if len(names) < 2:
                continue
            found = True
            print(f"    collision {key}")
            for name in names:
                data = metrics[name]
                print(
                    "      "
                    f"{name}: Rprime={float(data['Rprime']):.6f} "
                    f"delta={float(data['delta']):+.6f} "
                    f"q_zero={fmt_frac(data['q_zero_mass'])} "
                    f"q_range={data['q_range_lo']}..{data['q_range_hi']} "
                    f"Q={data['Qsize']} "
                    f"sumv2_even={data['sumv2_even']}"
                )
        if not found:
            print("    none")


def print_separator_audit(metrics: dict[str, dict[str, object]]) -> None:
    candidate_fields = (
        "fullE",
        "RE",
        "oddE",
        "q_zero_mass",
        "q_range_hi",
        "Qsize",
        "sumv2_even",
    )

    sign_pairs: list[tuple[str, str]] = []
    full_pairs: list[tuple[str, str]] = []
    names = tuple(metrics)
    for a, b in combinations(candidate_fields, 2):
        seen_sign: dict[tuple[object, object], set[int]] = {}
        seen_full: dict[tuple[object, object], str] = {}
        sign_ok = True
        full_ok = True
        for name in names:
            key = (metrics[name][a], metrics[name][b])
            seen_sign.setdefault(key, set()).add(metrics[name]["sign"])
            if len(seen_sign[key]) > 1:
                sign_ok = False
            if key in seen_full and seen_full[key] != name:
                full_ok = False
            else:
                seen_full[key] = name
        if sign_ok:
            sign_pairs.append((a, b))
        if full_ok:
            full_pairs.append((a, b))

    print("MINIMAL CURATED-BANK REPAIRS")
    print("  sign-separating field pairs:")
    for pair in sign_pairs:
        print(f"    {pair}")
    print("  full row-separating field pairs:")
    for pair in full_pairs:
        print(f"    {pair}")


def random_sample_metrics(sample_size: int = 10, seed: int = 2) -> list[dict[str, object]]:
    rng = Random(seed)
    rows: list[dict[str, object]] = []
    while len(rows) < sample_size:
        row = random_covering(rng)
        r_speeds, q_speeds = split_row(row)
        if not (1 <= len(q_speeds) <= 3):
            continue
        fiber = HYP3140.fiber_pgf_metrics(r_speeds, q_speeds)
        odd = tuple(speed for speed in row if speed % 2)
        even = tuple(speed for speed in row if speed % 2 == 0)
        rows.append(
            {
                "row": row,
                "Rprime": float(fiber["rprime_fiber"]),
                "delta": float(fiber["q_mean_sheets"] - fiber["mean_sheets"]),
                "RE": add_energy(r_speeds)[0],
                "oddE": add_energy(odd)[0] if len(odd) >= 2 else 0,
                "evenE": add_energy(even)[0] if len(even) >= 2 else 0,
                "fullE": add_energy(row)[0],
                "q_zero_mass": float(fiber["q_zero_mass"]),
                "q_range_hi": fiber["max_q_count"],
                "Qsize": len(q_speeds),
            }
        )
    return rows


def print_random_sample(rows: list[dict[str, object]]) -> None:
    print("RANDOM SMALL COVERING SAMPLE (seed=2, max_speed=90, Qsize<=3)")
    print("  Rprime      delta      fullE   RE   oddE  evenE  q_zero  q_hi  Q  row")
    for row in rows:
        print(
            "  "
            f"{row['Rprime']:.6f}  {row['delta']:+.6f}  "
            f"{row['fullE']:>5}  {row['RE']:>4}  {row['oddE']:>4}  {row['evenE']:>5}  "
            f"{row['q_zero_mass']:.6f}  {row['q_range_hi']:>4}  {row['Qsize']:>1}  {row['row']}"
        )

    features = ("fullE", "RE", "oddE", "evenE", "q_zero_mass", "q_range_hi", "Qsize")
    print("  correlations against delta and Rprime")
    for feature in features:
        values = [row[feature] for row in rows]
        deltas = [row["delta"] for row in rows]
        rprimes = [row["Rprime"] for row in rows]
        print(
            "    "
            f"{feature}: corr(delta)={corr(values, deltas):+.3f} "
            f"corr(Rprime)={corr(values, rprimes):+.3f}"
        )


def print_carrier_ranking() -> None:
    carriers = (
        ("P00_energy_plus_sheet_packet", 9, 9, 8, 8, 9, 0),
        ("P01_R_energy_scalar", 6, 4, 5, 2, 3, 2),
        ("P02_sheet_range_sidecar", 4, 8, 7, 7, 5, 1),
        ("P03_zero_sheet_mass_sidecar", 3, 8, 7, 7, 5, 1),
        ("P04_odd_energy_scalar", 2, 1, 1, 0, 2, 4),
        ("P05_raw_additive_shadow", 1, 0, 0, 0, 1, 5),
    )
    scored = []
    for code, floor_feed, sheet_payload, exact_collision, random_signal, packet_legality, scalar_risk in carriers:
        score = (
            5 * floor_feed
            + 4 * sheet_payload
            + 3 * exact_collision
            + 2 * random_signal
            + 3 * packet_legality
            - 5 * scalar_risk
        )
        scored.append((score, code))
    scored.sort(reverse=True)
    print("PROOF-CARRIER RANKING")
    for score, code in scored:
        print(f"  {code}: score={score}")
    print(
        "  priority_path="
        "P00_energy_plus_sheet_packet -> P01_R_energy_scalar -> "
        "P02_sheet_range_sidecar -> P03_zero_sheet_mass_sidecar -> "
        "P04_odd_energy_scalar -> P05_raw_additive_shadow"
    )


def main() -> None:
    print("HYP-3425 / codex-2026-06-28")
    print("Additive energy needs a sheet sidecar on the LRC14 covering floor")
    print()
    metrics = curated_metrics()
    print_collision_audit(metrics)
    print()
    print_separator_audit(metrics)
    print()
    random_rows = random_sample_metrics()
    print_random_sample(random_rows)
    print()
    print_carrier_ranking()
    print()
    print("VERDICT")
    print("  Raw additive-energy scalars collide across covering packets.")
    print("  The first legal repair is an energy-plus-sheet packet, not another scalar.")
    print("  Read HYP-2272 through HYP-3140/HYP-3424, not around them.")


if __name__ == "__main__":
    main()
