#!/usr/bin/env python3
"""Exact ranked-H1 closure of two omission-six LRC(14) roots.

THM-2899 leaves exactly one scalar-hard marked suffix on each target root.
For that suffix, let R4 be the sum of the four largest allowed singleton
coverages.  Since R4<6h/7, THM-2893 forces all five labels of any putative
cover into the finite core

    H1={w:c_C(w)>=h-R4}.

Order H1 by decreasing coverage on C, breaking ties by speed.  A five-set
has a unique earliest label x and four labels in its strict suffix.  For
each x this verifier proves that those four labels cannot cover C minus D_x,
using first their inherited parent coverages and otherwise their exact
literal-residual coverages.

The computation proves only the two named roots.  It does not claim a
uniform H1 closure or LRC(14).
"""

from __future__ import annotations

import hashlib
import importlib.util
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ENGINE_PATH = (
    ROOT
    / "04-computation"
    / "lrc14_j6_all_root_ranked_suffix_scalar_census_codex_20260729.py"
)
ENGINE_SHA256 = (
    "e0ff69252870f194549bba61289c1c5b15bef451e37a72836d71f9e71b1016e9"
)
TARGETS = (
    (1, 2, 3, 4, 5, 10, 12),
    (1, 2, 3, 4, 5, 12, 13),
)
PRIOR_ROOTS = (
    (2, 8, 9, 10, 11, 13, 14),
    (1, 3, 9, 10, 11, 12, 14),
    (2, 5, 9, 11, 12, 13, 14),
    (2, 3, 4, 5, 6, 7, 8),
    (1, 8, 10, 11, 12, 13, 14),
    (1, 2, 3, 4, 5, 6, 13),
    (1, 2, 3, 4, 6, 7, 14),
    (1, 2, 3, 4, 6, 11, 13),
    (1, 2, 3, 4, 6, 12, 13),
    (7, 8, 9, 11, 12, 13, 14),
    (1, 2, 3, 4, 6, 11, 12),
    (1, 2, 3, 5, 6, 10, 13),
    (1, 2, 4, 5, 6, 12, 13),
    (1, 3, 4, 5, 6, 7, 14),
    (5, 7, 8, 10, 11, 13, 14),
)

# Filled after discovery, then locked for ordinary and optimized replay.
EXPECTED_ROWS: tuple[tuple[object, ...], ...] | None = (
    (
        (1, 2, 3, 4, 5, 10, 12),
        10,
        (1,),
        1,
        22,
        (22,),
        "689/2310",
        28,
        (
            (26, "324/5005"),
            (39, "136/2145"),
            (27, "656/10395"),
            (28, "3/49"),
            (16, "1069/18480"),
        ),
        "18371/72765",
        "232/72765",
        "1333/29106",
        1774,
        1759,
        49,
        (
            26,
            39,
            27,
            28,
            16,
            78,
            18,
            65,
            32,
            54,
            21,
            52,
            56,
            40,
            69,
            63,
            130,
            91,
            46,
            120,
            108,
            42,
            45,
            41,
            68,
            25,
            286,
            220,
            81,
            92,
            156,
            97,
            107,
            82,
            125,
            182,
            264,
            275,
            121,
            260,
            250,
            195,
            104,
            242,
            162,
            175,
            50,
            23,
            237,
        ),
        43,
        2,
        4,
        95,
        "7949/3027024",
        "8e7695f12bd6fcf47a585a4fac96ddabc2ba4b61b4b9e5d616180be8e2337109",
        "e1b338a95d3e2a2344db0ec542d22441bd764835ba5398d7e0ac85229e16bc6b",
    ),
    (
        (1, 2, 3, 4, 5, 12, 13),
        6,
        (1,),
        1,
        22,
        (22,),
        "4127/15015",
        22,
        (
            (27, "1124/19305"),
            (28, "11/196"),
            (18, "38/693"),
            (16, "4339/80080"),
            (21, "158/3003"),
        ),
        "3380627/15135120",
        "26443/2162160",
        "59953/1164240",
        363,
        348,
        6,
        (27, 28, 18, 16, 21, 40),
        1,
        1,
        4,
        5,
        "22991/5045040",
        "cd0bd27d3259507850aac604d7752c1804be6e7db09271b92c301f5018861d82",
        "db3e2c332b80cbdeff5dc1716ae35ec2d0103a1e903136a71ba34c737deb3eb2",
    ),
)
EXPECTED_COMPLETE_DIGEST: str | None = (
    "7b173f9c5817ab6e73ddabba48bddf0f7d4714f7384ce034b23911e2056b5e96"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_engine():
    require(file_sha256(ENGINE_PATH) == ENGINE_SHA256, "THM-2899 engine changed")
    spec = importlib.util.spec_from_file_location("thm2902_engine", ENGINE_PATH)
    require(spec is not None and spec.loader is not None, "cannot load engine")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


S = load_engine()


def ftext(value: F | None) -> str | None:
    if value is None:
        return None
    return f"{value.numerator}/{value.denominator}"


def residual_mass(intervals: list[tuple[F, F]]) -> F:
    return sum((right - left for left, right in intervals), F(0))


def profile_target(body: tuple[int, ...]) -> tuple[tuple[object, ...], str]:
    """Reconstruct one target and close its unique scalar-hard branch."""

    root = S.profile_root(body)
    hard = tuple(
        row for row in root["branches"] if not row["closed"]
    )
    require(len(hard) == 1, f"target lost one-hard status: {body}")
    row = hard[0]
    require(
        row["rank"] == 1 and row["prefix"] == (row["apex"],),
        f"target hard branch is not rank one: {body}",
    )

    top4 = row["top5"][:4]
    r4 = sum((value for value, _ in top4), F(0))
    margin = F(6, 7) * row["m"] - r4
    require(margin > 0, f"target lost ranked r=4 eligibility: {body}")
    level = row["m"] - r4
    epsilon = level - row["m"] / 7
    require(epsilon == margin, f"H1 threshold identity changed: {body}")
    threshold = S.S2 * row["r"] / (7 * epsilon)
    cutoff = S.ceiling(threshold) - 1
    require(
        row["m"] / 7 + S.S2 * row["r"] / (7 * (cutoff + 1))
        <= level,
        f"H1 tail did not seal: {body}",
    )

    gate_root = S.G.profile_body(body)
    gate_root["good"] = S.G.T.CORE.good_norm(body)[0]
    carrier = S.R.subtract_local(gate_root["good"], row["apex"])
    direct, direct_r, direct_h = S.G.T.CORE.good_norm(
        tuple(sorted((*body, row["apex"])))
    )
    require(
        carrier == direct
        and len(carrier) == direct_r == row["r"]
        and residual_mass(carrier) == direct_h == row["m"],
        f"literal/direct carrier mismatch: {body}",
    )

    excluded = set(row["prefix"])
    speeds = tuple(
        speed
        for speed in range(S.FIRST_EXTERNAL, cutoff + 1)
        if speed not in excluded
    )
    coverages = S.G.T.coverages_many(carrier, speeds)
    for value, speed in coverages:
        require(
            value == S.G.T.coverage(carrier, speed),
            f"vector/scalar H1 mismatch: {body}, {speed}",
        )
    core_rows = tuple(
        sorted(
            (
                (value, speed)
                for value, speed in coverages
                if value >= level
            ),
            key=lambda item: (-item[0], item[1]),
        )
    )
    require(len(core_rows) >= 5, f"unexpected small H1 core: {body}")
    require(
        all(speed not in excluded for _, speed in core_rows),
        f"excluded prefix entered H1: {body}",
    )

    core_hash = hashlib.sha256(b"LRC14/THM2902/H1-core/v1\n")
    for value, speed in core_rows:
        core_hash.update(f"{speed}:{ftext(value)}\n".encode())

    parent_count = 0
    local_count = 0
    short_count = 0
    local_evaluations = 0
    local_certificates = []
    minimum_margin: F | None = None
    depth_hash = hashlib.sha256(b"LRC14/THM2902/depth-one/v1\n")
    for index, (_, speed) in enumerate(core_rows):
        residual = S.R.subtract_local(carrier, speed)
        mass = residual_mass(residual)
        suffix = core_rows[index + 1 :]
        if len(suffix) < 4:
            mode = "SHORT"
            used_margin = None
            upper = F(0)
            short_count += 1
        else:
            inherited_upper = sum(
                (value for value, _ in suffix[:4]),
                F(0),
            )
            inherited_margin = mass - inherited_upper
            if inherited_margin > 0:
                mode = "PARENT"
                used_margin = inherited_margin
                upper = inherited_upper
                parent_count += 1
            else:
                local_rows = []
                for _, other in suffix:
                    child = S.R.subtract_local(residual, other)
                    gain = mass - residual_mass(child)
                    local_rows.append((gain, other))
                    local_evaluations += 1
                local_rows.sort(key=lambda item: (-item[0], item[1]))
                local_upper = sum(
                    (value for value, _ in local_rows[:4]),
                    F(0),
                )
                local_margin = mass - local_upper
                require(
                    local_margin > 0,
                    f"depth-one branch remains open: {body}, {speed}",
                )
                local_certificates.append(
                    (
                        speed,
                        ftext(local_margin),
                        tuple(
                            (other, ftext(gain))
                            for gain, other in local_rows[:4]
                        ),
                    )
                )
                mode = "LOCAL"
                used_margin = local_margin
                upper = local_upper
                local_count += 1
        if used_margin is not None:
            minimum_margin = (
                used_margin
                if minimum_margin is None
                else min(minimum_margin, used_margin)
            )
        depth_hash.update(
            (
                f"{index};{speed};{mode};{ftext(mass)};"
                f"{ftext(upper)};{ftext(used_margin)}\n"
            ).encode()
        )

    require(
        parent_count + local_count + short_count == len(core_rows),
        f"depth-one partition incomplete: {body}",
    )
    require(minimum_margin is not None, f"no positive depth margin: {body}")

    summary = (
        body,
        root["K"],
        tuple(row["rank"] for row in hard),
        row["rank"],
        row["apex"],
        row["prefix"],
        ftext(row["m"]),
        row["r"],
        tuple((speed, ftext(value)) for value, speed in row["top5"]),
        ftext(r4),
        ftext(margin),
        ftext(level),
        cutoff,
        len(speeds),
        len(core_rows),
        tuple(speed for _, speed in core_rows),
        parent_count,
        local_count,
        short_count,
        local_evaluations,
        ftext(minimum_margin),
        core_hash.hexdigest(),
        depth_hash.hexdigest(),
    )
    detail = (
        f"E={body};K={root['K']};hard_ranks={(row['rank'],)};"
        f"rank={row['rank']};apex={row['apex']};prefix={row['prefix']};"
        f"h={ftext(row['m'])};r={row['r']};R4={ftext(r4)};"
        f"margin={ftext(margin)};level={ftext(level)};cutoff={cutoff};"
        f"scanned={len(speeds)};H1={len(core_rows)};"
        f"parent={parent_count};local={local_count};short={short_count};"
        f"local_evals={local_evaluations};"
        f"local_certificates={tuple(local_certificates)};"
        f"min_depth_margin={ftext(minimum_margin)};"
        f"core_digest={core_hash.hexdigest()};"
        f"depth_digest={depth_hash.hexdigest()};H1_labels="
        + ",".join(str(speed) for _, speed in core_rows)
    )
    return summary, detail


def main() -> None:
    require(
        set(TARGETS).isdisjoint(PRIOR_ROOTS),
        "target overlaps a previously proved root",
    )
    records = tuple(profile_target(body) for body in TARGETS)
    summaries = tuple(summary for summary, _ in records)
    digest = hashlib.sha256(b"LRC14/THM2902/omission-six/v1\n")
    for _, detail in records:
        digest.update((detail + "\n").encode())
    complete_digest = digest.hexdigest()

    if EXPECTED_ROWS is not None:
        require(summaries == EXPECTED_ROWS, "locked target summaries changed")
    if EXPECTED_COMPLETE_DIGEST is not None:
        require(
            complete_digest == EXPECTED_COMPLETE_DIGEST,
            "complete THM-2902 digest changed",
        )

    print("LRC14 j6 omission-six ranked H1 depth-one closure")
    for _, detail in records:
        print(detail)
    print(f"closed_roots={TARGETS}")
    print("disjoint_from_THM2895_THM2898_THM2899_THM2901=PASS")
    print(f"complete_digest={complete_digest}")
    if EXPECTED_ROWS is None:
        print("mode=DISCOVERY")
    else:
        print("mode=LOCKED")
    print("scope=two named seven-body roots;not the rung;not LRC14")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
