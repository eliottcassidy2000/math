#!/usr/bin/env python3
"""Independent exact audit of the five j6 recursive flag survivors.

This is scratch-only evidence.  It imports the existing rank-one H4 probe,
reconstructs its 306 pair rows, selects the five rows missed by both the
top-three and B2+q1 certificates, applies THM-2893 with
``(k,s,ell)=(3,2,2)``, and closes the last singleton by the geometric
longest-interval criterion rather than the discrepancy horizon.
"""

from __future__ import annotations

import importlib.util
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
H4_PATH = (
    ROOT
    / ".scratch/lrc_sparse_j6_gate_audit_20260729/h4_pair_residual_probe.py"
)
FIRST_EXTERNAL = 15


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_h4():
    spec = importlib.util.spec_from_file_location("j6_h4_geometric_audit", H4_PATH)
    require(spec is not None and spec.loader is not None, "cannot load H4 probe")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


H = load_h4()


def floor_fraction(value: F) -> int:
    return value.numerator // value.denominator


def ceil_fraction(value: F) -> int:
    return -((-value.numerator) // value.denominator)


def mass(intervals: list[tuple[F, F]]) -> F:
    return sum((right - left for left, right in intervals), F(0))


def interval_in_one_tooth(left: F, right: F, speed: int) -> bool:
    """Exact containment up to measure-zero interval endpoints."""

    lower_center = ceil_fraction(speed * right - F(1, 14))
    upper_center = floor_fraction(speed * left + F(1, 14))
    return lower_center <= upper_center


def singleton_contains(
    residual: list[tuple[F, F]],
    speed: int,
) -> bool:
    return all(
        interval_in_one_tooth(left, right, speed)
        for left, right in residual
    )


def recursive_row(
    root: dict[str, object],
    pair: dict[str, object],
) -> dict[str, object]:
    residual = pair["residual"]
    residual_mass = pair["m"]
    q1 = pair["q1"]

    # THM-2893 with (k,s,ell)=(3,2,2): the complementary one-set
    # cap must satisfy q1 < (7-2)m/7=5m/7.
    require(q1 < F(5, 7) * residual_mass, "q1 reaches the H2 boundary")
    theta = residual_mass - q1
    level = theta / 2
    require(level > residual_mass / 7, "H2 core is not finite")

    threshold = (
        H.S2
        * pair["r"]
        / (7 * (level - residual_mass / 7))
    )
    tail_first = max(FIRST_EXTERNAL, H.ceiling(threshold))
    excluded = {root["apex"], *pair["pair"]}
    by_speed = {
        speed: value
        for value, speed in pair["pair_data"]["ranked"]
    }
    if tail_first > H.HORIZON + 1:
        by_speed.update(
            {
                speed: value
                for value, speed in H.T.coverages_many(
                    residual,
                    [
                        speed
                        for speed in range(H.HORIZON + 1, tail_first)
                        if speed not in excluded
                    ],
                )
            }
        )
    require(
        residual_mass / 7
        + H.S2 * pair["r"] / (7 * tail_first)
        <= level,
        "strict discrepancy did not seal H2",
    )
    core = tuple(
        sorted(
            speed
            for speed, value in by_speed.items()
            if speed not in excluded and value >= level
        )
    )

    heavy = []
    heavy_residuals = []
    for edge in combinations(core, 2):
        after = H.R.subtract_local_multi(residual, edge)
        if residual_mass - mass(after) >= theta:
            heavy.append(edge)
            heavy_residuals.append(after)

    horizons = []
    checks = 0
    covers = []
    for edge, after in zip(heavy, heavy_residuals):
        after_mass = mass(after)
        require(after and after_mass > 0, "heavy edge itself covers residual")

        family = tuple(
            sorted((*root["body"], root["apex"], *pair["pair"], *edge))
        )
        direct, direct_r, direct_m = H.T.CORE.good_norm(family)
        require(
            after == direct
            and len(after) == direct_r
            and after_mass == direct_m,
            "recursive literal/direct residual mismatch",
        )

        longest = max(right - left for left, right in after)
        horizon = floor_fraction(F(1, 7) / longest)
        horizons.append(horizon)
        edge_excluded = excluded | set(edge)
        require(
            longest > F(1, 7 * (horizon + 1)),
            "geometric tail equality mishandled",
        )
        for speed in range(FIRST_EXTERNAL, horizon + 1):
            if speed in edge_excluded:
                continue
            geometric_cover = singleton_contains(after, speed)
            scalar_cover = H.T.coverage(after, speed) == after_mass
            require(
                geometric_cover == scalar_cover,
                "interval/scalar singleton predicate mismatch",
            )
            checks += 1
            if geometric_cover:
                covers.append((edge, speed))

    return {
        "H": core,
        "heavy": tuple(heavy),
        "horizons": tuple(horizons),
        "checks": checks,
        "covers": tuple(covers),
    }


def main() -> None:
    roots = [H.h4_core(body, apex) for body, apex in H.CASES]
    rows = [
        (root, H.probe_pair(root, pair))
        for root in roots
        for pair in combinations(root["H"], 2)
    ]
    failures = [
        (root, pair)
        for root, pair in rows
        if pair["direct3_margin"] <= 0 and pair["cheap_margin"] <= 0
    ]
    require(len(rows) == 306, "rank-one H4 pair universe changed")
    require(len(failures) == 5, "adaptive failure count changed")

    recursive = [
        (root, pair, recursive_row(root, pair))
        for root, pair in failures
    ]
    require(
        tuple(len(row["H"]) for _, _, row in recursive) == (3, 3, 3, 2, 2)
        and all(len(row["heavy"]) == 1 for _, _, row in recursive)
        and tuple(row["horizons"][0] for _, _, row in recursive)
        == (41, 38, 32, 31, 33)
        and sum(row["checks"] for _, _, row in recursive) == 86
        and all(not row["covers"] for _, _, row in recursive),
        "recursive geometric census changed",
    )

    print("J6 RECURSIVE FLAG GEOMETRIC INDEPENDENT AUDIT")
    print(
        "first_flag=(k,s,ell)=(5,4,2);"
        f"H4_pairs={len(rows)};"
        f"adaptive_closed={len(rows)-len(failures)};"
        f"adaptive_failures={len(failures)}"
    )
    for root, pair, row in recursive:
        print(
            f"E={root['body']};apex={root['apex']};"
            f"H4pair={pair['pair']};H2={row['H']};"
            f"heavy={row['heavy']};"
            f"geometric_horizon={row['horizons'][0]};"
            f"checks={row['checks']};covers={row['covers']}"
        )
    print(
        "second_flag=(k,s,ell)=(3,2,2);"
        "H2_sizes=(3,3,3,2,2);heavy_edges=5;"
        "singleton_checks=86;covers=0"
    )
    print(
        "boundaries=q1<3h/7 then q1<5m/7;"
        "heavy/high use >=;geometric equality is scanned;"
        "empty residual is failure"
    )
    print("verdict=PASS")


if __name__ == "__main__":
    main()
