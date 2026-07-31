#!/usr/bin/env python3
r"""Exact global top-fifteen gate for every eight-speed body in ``{1,...,14}``.

For an eight-speed body ``E``, let ``G_E`` be its gap-``1/14`` good set,
``m_E=|G_E|``, and

    c_E(w)=|G_E intersect D_w|,              w>=15.

The strict THM-735(ii) discrepancy cap is

    c_E(w) < m_E/7 + (99/70) r_E/(7w),

where ``r_E`` is the number of interval components of ``G_E``.  This
companion reconstructs the global top fifteen coverages for all
``C(14,8)=3003`` bodies and proves

    q_11(E)+q_12(E)+q_13(E)+q_14(E)+q_15(E) < m_E.

Thus every five-speed set whose individual coverages can exhaust ``G_E``
meets the chosen global top ten.  This is a hitting reduction only; it does
not close the resulting apex carriers or prove the five-slot rung.
"""

from __future__ import annotations

import hashlib
import importlib.util
from fractions import Fraction as F
from itertools import combinations
from math import lcm
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CORE_PATH = ROOT / "04-computation/lrc14_thm741_2002_body_j4_tree_kps_S128c5.py"
CORE_SHA256 = "5aa81d9d78273c8f9e3e7a6574091a3bc3f64ab6086c7024c15f9420c99dac96"
FIRST_EXTERNAL = 15
BASE_HORIZON = 600
S2 = F(99, 70)
BODIES = tuple(combinations(range(1, 15), 8))

EXPECTED_RAW_TOP5 = (16, 2987)
EXPECTED_MINIMUM_HITTING = (
    F(61, 12936),
    (1, 3, 5, 7, 8, 9, 11, 13),
    (36, 42, 84, 25, 72),
)
EXPECTED_MAXIMUM_THRESHOLD = (
    F(1002456, 953),
    (2, 8, 9, 10, 11, 12, 13, 14),
    1052,
)
EXPECTED_MAXIMUM_TOP15_SPEED = 243
EXPECTED_LEDGER_DIGEST = (
    "f1c478fc386171e4bfbdab52100808d7a37853ae755a49c9e45c6412708bb20a"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_core():
    require(file_sha256(CORE_PATH) == CORE_SHA256, "THM-741 core hash changed")
    spec = importlib.util.spec_from_file_location("thm2885_eight_body_core", CORE_PATH)
    require(spec is not None and spec.loader is not None, "cannot load THM-741 core")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CORE = load_core()


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def ceiling(value: F) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def tooth_primitive(value: F) -> F:
    """Primitive of the unit-periodic radius-``1/14`` tooth indicator."""

    integer = value.numerator // value.denominator
    remainder = value - integer
    return (
        F(integer, 7)
        + min(remainder, F(1, 14))
        + max(F(0), remainder - F(13, 14))
    )


def coverage(good: list[tuple[F, F]], speed: int) -> F:
    """Evaluate one tooth coverage by the scalar exact primitive."""

    return sum(
        (
            tooth_primitive(speed * right)
            - tooth_primitive(speed * left)
        )
        / speed
        for left, right in good
    )


def coverages_many(
    good: list[tuple[F, F]],
    speeds: range | list[int],
    batch_size: int = 20_000,
) -> list[tuple[F, int]]:
    """Evaluate the same primitive exactly in guarded int64 batches."""

    speed_list = list(speeds)
    if not speed_list:
        return []
    require(good and batch_size >= 1, "bad vector primitive input")
    denominator = lcm(
        *(endpoint.denominator for interval in good for endpoint in interval)
    )
    maximum_speed = max(speed_list)
    # Endpoints lie in [0,1].  The first guard controls scaled endpoints and
    # primitive values; the second controls the nonnegative component sum
    # even under the crude coverage bound c_E(w)<=1.
    require(
        4 * len(good) * denominator * maximum_speed < 2**62,
        "vector primitive component arithmetic would overflow int64",
    )
    require(
        14 * denominator * maximum_speed < 2**63,
        "vector primitive axis reduction would overflow int64",
    )
    endpoint_ints = np.array(
        [
            endpoint.numerator * (denominator // endpoint.denominator)
            for interval in good
            for endpoint in interval
        ],
        dtype=np.int64,
    )
    out: list[tuple[F, int]] = []
    for start in range(0, len(speed_list), batch_size):
        speed_batch = speed_list[start : start + batch_size]
        speed_array = np.array(speed_batch, dtype=np.int64)
        scaled = endpoint_ints[:, None] * speed_array[None, :]
        quotient, remainder = np.divmod(scaled, denominator)
        primitive = (
            2 * denominator * quotient
            + np.minimum(14 * remainder, denominator)
            + np.maximum(0, 14 * remainder - 13 * denominator)
        )
        numerators = (primitive[1::2] - primitive[0::2]).sum(axis=0)
        out.extend(
            (F(int(numerator), 14 * denominator * speed), speed)
            for numerator, speed in zip(numerators, speed_batch)
        )
    return out


def profile_body(body: tuple[int, ...]) -> dict[str, object]:
    """Seal and return one global top-fifteen profile."""

    good, root_r, root_m = CORE.good_norm(body)
    require(root_r == len(good) and root_m > 0, f"bad root carrier {body}")

    rows = coverages_many(good, range(FIRST_EXTERNAL, BASE_HORIZON + 1))
    base_by_speed = {speed: value for value, speed in rows}
    require(
        base_by_speed[FIRST_EXTERNAL] == coverage(good, FIRST_EXTERNAL),
        f"base vector/scalar mismatch {body}",
    )
    ranked600 = sorted(rows, key=lambda item: (-item[0], item[1]))
    q15_600 = ranked600[14][0]
    require(q15_600 > root_m / 7, f"rank fifteen misses limit {body}")
    threshold = S2 * root_r / (7 * (q15_600 - root_m / 7))
    tail_first = ceiling(threshold)

    if tail_first > BASE_HORIZON + 1:
        rows.extend(
            coverages_many(
                good,
                range(BASE_HORIZON + 1, tail_first),
            )
        )
    control_speed = tail_first - 1 if tail_first > BASE_HORIZON + 1 else BASE_HORIZON
    by_speed = {speed: value for value, speed in rows}
    require(
        by_speed[control_speed] == coverage(good, control_speed),
        f"tail vector/scalar mismatch {body}",
    )

    # If threshold is integral, the majorant at tail_first can equal q15_600,
    # but the THM-735(ii) coverage inequality itself is strict.
    require(
        root_m / 7 + S2 * root_r / (7 * tail_first) <= q15_600,
        f"rank-fifteen tail seal failed {body}",
    )
    ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    top15 = tuple(ranked[:15])
    require(
        len(top15) == 15
        and len({speed for _, speed in top15}) == 15
        and all(speed >= FIRST_EXTERNAL for _, speed in top15),
        f"bad top-fifteen ledger {body}",
    )
    raw_top5_margin = root_m - sum(
        (value for value, _ in top15[:5]),
        F(0),
    )
    hitting_margin = root_m - sum(
        (value for value, _ in top15[10:15]),
        F(0),
    )
    require(hitting_margin > 0, f"top-ten hitting gate failed {body}")
    return {
        "body": body,
        "m": root_m,
        "r": root_r,
        "threshold": threshold,
        "tail_first": tail_first,
        "top15": top15,
        "raw_top5_margin": raw_top5_margin,
        "hitting_margin": hitting_margin,
        "scalar_controls": 2,
    }


def main() -> None:
    require(
        len(BODIES) == 3003 and len(set(BODIES)) == 3003,
        "eight-body universe changed",
    )
    require(S2**2 > 2, "sqrt(2) majorant changed")
    profiles = [profile_body(body) for body in BODIES]
    require(
        tuple(row["body"] for row in profiles) == BODIES,
        "profile order changed",
    )

    raw_positive = sum(row["raw_top5_margin"] > 0 for row in profiles)
    raw_nonpositive = len(profiles) - raw_positive
    require(
        (raw_positive, raw_nonpositive) == EXPECTED_RAW_TOP5,
        "raw top-five census changed",
    )
    minimum_hitting = min(
        (
            row["hitting_margin"],
            row["body"],
            tuple(speed for _, speed in row["top15"][10:15]),
        )
        for row in profiles
    )
    require(
        minimum_hitting == EXPECTED_MINIMUM_HITTING,
        "minimum hitting margin changed",
    )
    maximum_threshold = max(
        (row["threshold"], row["body"], row["tail_first"])
        for row in profiles
    )
    require(
        maximum_threshold == EXPECTED_MAXIMUM_THRESHOLD,
        "maximum rank-fifteen threshold changed",
    )
    maximum_top15_speed = max(
        speed
        for row in profiles
        for _, speed in row["top15"]
    )

    ledger_text = "THM2885/eight-body-global-top15-hitting/v1\n" + "".join(
        "E="
        + ",".join(map(str, row["body"]))
        + f";m={ftext(row['m'])};r={row['r']};"
        + f"T={ftext(row['threshold'])};tail={row['tail_first']};"
        + "top15="
        + ",".join(
            f"{speed}:{ftext(value)}"
            for value, speed in row["top15"]
        )
        + f";margin={ftext(row['hitting_margin'])}\n"
        for row in profiles
    )
    ledger_digest = hashlib.sha256(ledger_text.encode()).hexdigest()
    require(
        maximum_top15_speed == EXPECTED_MAXIMUM_TOP15_SPEED,
        "maximum top-fifteen speed changed",
    )
    require(ledger_digest == EXPECTED_LEDGER_DIGEST, "top-fifteen ledger changed")

    print("THM-2885 EIGHT-BODY GLOBAL TOP-FIFTEEN HITTING GATE")
    print("status=FINITE-EXACT+GLOBAL-TAIL-SEALED")
    print("universe=C(14,8)=3003;external_speeds=w>=15")
    print(
        f"raw_top5=positive:{raw_positive},nonpositive:{raw_nonpositive}"
    )
    print(
        "top10_hitting=3003/3003;"
        f"minimum_margin={ftext(minimum_hitting[0])};"
        f"minimum_body={minimum_hitting[1]};"
        "minimum_complement="
        + ",".join(map(str, minimum_hitting[2]))
    )
    print(
        "maximum_rank15_threshold="
        f"{ftext(maximum_threshold[0])};"
        f"body={maximum_threshold[1]};"
        f"tail_first={maximum_threshold[2]};"
        f"maximum_top15_speed={maximum_top15_speed}"
    )
    print(
        "vector_scalar_controls="
        + str(sum(int(row["scalar_controls"]) for row in profiles))
    )
    print(f"top15_ledger_sha256={ledger_digest}")
    print(
        "scope=global top-ten hitting reduction only;"
        "ranked-apex carrier closure and the j=5 rung remain open"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
