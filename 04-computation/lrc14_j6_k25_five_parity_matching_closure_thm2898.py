#!/usr/bin/env python3
r"""Locked exact closure of the unique K=25 seven-body j6 root.

THM-2896 gives a globally sealed six-cover hitting gate of least size 25
for

    E=(1,8,10,11,12,13,14).

For every gate apex this verifier builds the reusable THM-2893 scalar
activation profile.  It then performs an exact breadth-first search on the
scalar-closed masks

    X -> cl(X union {a}),

so the first full state has the exact minimum number of paid nonscalar
edges.  The closed-state logic is inlined here: no scratch search helper is
imported.

The verifier audits THM-2897's rank-selective pair cap, but the proof uses
no such activation.  Its intrinsic undirected two-edge matching repair
instead closes apex 17 from the smaller fixed-rank prefix, hence also from
the larger actual prefix by monotonicity.  The remaining five paid edges
are certified on their actual prefixes using THM-2895's
singleton-complement H4 flag.  Every H4-pair residual is checked by the two
locked child certificates (top-three singleton sum and B2+B1); none needs
the recursive H2 fallback on this root.

This is one finite-exact root closure.  It is not the all-root j6 rung and
not LRC(14).
"""

from __future__ import annotations

import hashlib
import importlib.util
from collections import deque
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
THM2895_PATH = (
    ROOT
    / "01-canon/theorems/"
    "THM-2895-singleton-complement-parity-descent-and-four-root-j6-closure.md"
)
THM2895_SHA256 = (
    "1bbec361c3bf0a570c2c87723ab11c7c51dbe9297974fc22fd37034bea6d5e56"
)
THM2896_PATH = (
    ROOT
    / "01-canon/theorems/"
    "THM-2896-seven-body-adaptive-six-cover-hitting-gate-atlas.md"
)
THM2896_SHA256 = (
    "44e1492e6fecaa51b901c01defaf00f940057324605cb40ece282e81d0e05615"
)
THM2897_PATH = (
    ROOT
    / "01-canon/theorems/"
    "THM-2897-partition-cap-tropical-convolution-and-alternating-pair-ladder.md"
)
THM2897_SHA256 = (
    "8eda5146ed10d62dddf793c613f52a26b2eb7198b37bff9c998b8daa4424b7fa"
)
ATLAS_PATH = (
    ROOT
    / "04-computation/lrc14_j6_all_root_adaptive_gate_atlas_codex_20260729.py"
)
ATLAS_SHA256 = (
    "fc36f26d4c8da5b005465696b954eec700c080376eef9ee5ba74a7111def99d7"
)
ATLAS_OUTPUT_PATH = (
    ROOT
    / "05-knowledge/results/lrc14_j6_all_root_adaptive_gate_atlas_codex_20260729.out"
)
ATLAS_OUTPUT_SHA256 = (
    "3081b93a870faacb31d205e43f7ca87872d7a9f196f4774d8740ced6a314d80b"
)
BOOTSTRAP_PATH = (
    ROOT
    / "04-computation/lrc14_j6_suffix_bootstrap_minseed_battery_codex_20260729.py"
)
BOOTSTRAP_SHA256 = (
    "a7a77dc433b21d94a54524064ccd62e553ed67ae8a3d3364bf79c41c36849d04"
)
BOOTSTRAP_OUTPUT_PATH = (
    ROOT
    / "05-knowledge/results/lrc14_j6_suffix_bootstrap_minseed_battery_codex_20260729.out"
)
BOOTSTRAP_OUTPUT_SHA256 = (
    "bdcf896d152d206b3ae77a3568609887190e1ba991909481769d7ad560a68835"
)
PARITY_PATH = (
    ROOT
    / "04-computation/lrc14_j6_suffix_parity_flag_closure_thm2895.py"
)
PARITY_SHA256 = (
    "970d77503f8d56d737e223dabb3c3562d7b19cd018ca75398e3deb054715e5f6"
)
PARITY_OUTPUT_PATH = (
    ROOT
    / "05-knowledge/results/lrc14_j6_suffix_parity_flag_closure_thm2895.out"
)
PARITY_OUTPUT_SHA256 = (
    "c11260f6544a319e1cc1862c9221b188a4314860422470e465b82e7ce492b1b4"
)

BODY = (1, 8, 10, 11, 12, 13, 14)
MATCHING_REFERENCE_PREFIX = (23, 27, 19, 46, 18, 17)
EXPECTED_GATE = (
    23,
    27,
    19,
    46,
    18,
    17,
    25,
    34,
    38,
    100,
    63,
    156,
    29,
    125,
    130,
    44,
    37,
    50,
    92,
    168,
    72,
    32,
    110,
    54,
    182,
)
EXPECTED_ROOT = (
    "66023/280280",
    38,
    25,
    "5703505/4933292364",
)
EXPECTED_SEARCH: tuple[int, int, int, int] | None = (
    6,
    68_407,
    323_776,
    68_407,
)
EXPECTED_PATH: tuple[object, ...] | None = (
    ((), 23, (), ()),
    ((23,), 27, (), ()),
    ((23, 27), 19, (), ()),
    ((23, 27, 19), 46, (), ()),
    (
        (23, 27, 19, 46),
        18,
        ((168, 182),),
        (("379903/214414200", "20510209/11363952600"),),
    ),
    (
        (23, 27, 19, 46, 18, 168, 182),
        17,
        (
            (25, 156, 125, 44),
            (130, 72, 54),
            (34, 100, 92, 32, 110),
            (38, 63, 29, 37, 50),
        ),
        (
            (
                "849923/6896289400",
                "276599/321621300",
                "8976649/1216870200",
                "808/1530765",
            ),
            (
                "22403/126266140",
                "888589/300232800",
                "17210524381/9826296173694",
            ),
            (
                "1115281/1554502950",
                "6990073943/5758026974700",
                "13781227/20443419360",
                "4232423633/2773490559840",
                "2216017/661110450",
            ),
            (
                "1720183/876922200",
                "899364751/98933234400",
                "1225034970487/111750034916520",
                "16613071/85711025400",
                "97703231/41377736400",
            ),
        ),
    ),
)
EXPECTED_SEED_BANK: tuple[int, ...] | None = (23, 27, 19, 46, 18, 17)
EXPECTED_PROOF_SCHEDULE: tuple[object, ...] | None = (
    ("PARITY", (), 23, "EXACT_FAIL", (), ()),
    ("PARITY", (23,), 27, "LOWER_REFUTED", (), ()),
    ("PARITY", (23, 27), 19, "LOWER_REFUTED", (), ()),
    ("PARITY", (23, 27, 19), 46, "LOWER_REFUTED", (), ()),
    (
        "PARITY",
        (23, 27, 19, 46),
        18,
        "EXACT_FAIL",
        ((168, 182),),
        (("379903/214414200", "20510209/11363952600"),),
    ),
    (
        "MATCHING",
        (23, 27, 19, 46, 18, 168, 182),
        17,
        "c703e55e854ddd83f09ccd919681a2542350041c502c2cafeaf1a904d921e219",
        (
            (25, 156, 125, 44),
            (130, 72, 54),
            (34, 100, 92, 32, 110),
            (38, 63, 29, 37, 50),
        ),
        (
            (
                "849923/6896289400",
                "276599/321621300",
                "8976649/1216870200",
                "808/1530765",
            ),
            (
                "22403/126266140",
                "888589/300232800",
                "17210524381/9826296173694",
            ),
            (
                "1115281/1554502950",
                "6990073943/5758026974700",
                "13781227/20443419360",
                "4232423633/2773490559840",
                "2216017/661110450",
            ),
            (
                "1720183/876922200",
                "899364751/98933234400",
                "1225034970487/111750034916520",
                "16613071/85711025400",
                "97703231/41377736400",
            ),
        ),
    ),
)
EXPECTED_STATE_SUMMARIES: tuple[object, ...] | None = (
    ((), 25, 18, 7, ()),
    ((23,), 24, 21, 3, ()),
    ((23, 27), 23, 19, 4, ()),
    ((23, 27, 19), 22, 20, 2, ()),
    ((23, 27, 19, 46), 21, 14, 7, ()),
)
EXPECTED_BRANCH_RESULTS: tuple[object, ...] | None = (
    (
        1,
        (),
        23,
        1,
        "247363/1381380",
        38,
        18,
        "62896/1450449",
        "968347/29008980",
        "171161/5045040",
        920,
        (17, 18, 19, 25, 27, 34, 36, 37, 38, 44, 47, 72, 125),
        78,
        78,
        78,
        78,
        0,
        765,
    ),
    (
        2,
        (23,),
        27,
        2,
        "170909/945945",
        44,
        19,
        "855419/17972955",
        "57752/1935549",
        "199321/5990985",
        1192,
        (17, 19, 25, 32, 34, 37, 38, 41, 46, 53, 57, 156, 182),
        78,
        77,
        78,
        78,
        0,
        569,
    ),
    (
        3,
        (23, 27),
        19,
        3,
        "59287/326040",
        34,
        18,
        "113389/2396394",
        "133391/4357080",
        "6447409/191711520",
        898,
        (17, 18, 25, 29, 32, 37, 46, 52, 81, 100),
        45,
        45,
        45,
        45,
        0,
        264,
    ),
    (
        4,
        (23, 27, 19),
        46,
        4,
        "3531007/19339320",
        42,
        34,
        "4628653/109589480",
        "6906637/191781590",
        "88733/2528988",
        943,
        (17, 18, 25, 32, 34, 38, 41, 44, 53, 72, 125),
        55,
        55,
        55,
        55,
        0,
        764,
    ),
    (
        5,
        (23, 27, 19, 46),
        18,
        5,
        "66481/360360",
        32,
        17,
        "1988057/42882840",
        "701237/21441420",
        "987197/28588560",
        791,
        (17, 25, 29, 34, 37, 38, 44),
        21,
        21,
        21,
        21,
        0,
        187,
    ),
)
EXPECTED_MATCHING_RESULT: tuple[object, ...] | None = (
    (23, 27, 19, 46, 18, 168, 182),
    (23, 27, 19, 46, 18, 17),
    17,
    "189841/1021020",
    34,
    (25, "1023641/23823800"),
    (125, "81891/2290750"),
    "33906683/440740300",
    (25, 37),
    7,
    "b296c6ae4e1407d8ba5e172b9091268141a1a34dbf788bed82bac216df1b22b7",
    "26834677/178678500",
    "121070701/1652776125",
    "1492091/400673000",
    1844,
    "523888817/17874997140",
    1824,
    "cdbc1748f16346880ba95d3912227784967471c4b31a0894237200d10ae589b4",
    42,
    "7b424e08965f373c0ff0850ea6ad23267c63e0f7e603a0da6c606b21c69a0570",
    (
        (25, 37, "33906683/440740300"),
        (25, 45, "3985831/53603550"),
        (37, 50, "65302143/881480600"),
    ),
    "babed5b1b371f0e046f0cb48157618b25998720acea098d33fa95b6333251e90",
    (),
    (
        (
            "600111641/3966662700",
            (25, 37),
            (25, 45),
            False,
        ),
        (
            "133115509/881480600",
            (25, 37),
            (37, 50),
            False,
        ),
        (
            "47104891/317333016",
            (25, 45),
            (37, 50),
            True,
        ),
    ),
    (
        (
            "600111641/3966662700",
            (25, 37),
            (25, 45),
            False,
        ),
        (
            "133115509/881480600",
            (25, 37),
            (37, 50),
            False,
        ),
    ),
    "47104891/317333016",
    (
        (
            "47104891/317333016",
            (25, 45),
            (37, 50),
            True,
        ),
    ),
    "69186919/39666627000",
)
EXPECTED_AGGREGATE: tuple[int, ...] | None = (
    5,
    0,
    0,
    1,
    277,
    276,
    277,
    277,
    0,
    2_549,
)
EXPECTED_RANKPAIR_AGGREGATE: tuple[int, ...] | None = (
    5,
    115,
    92,
    23,
    0,
    23,
    382,
)
EXPECTED_PROFILE_DIGEST: str | None = (
    "9e17ebe4b7be86ad4fce3dbadce2c1f3cc7957a232f3b336bed424be606b0e03"
)
EXPECTED_RANKPAIR_DIGEST: str | None = (
    "8d6acf68adf838a985548bcdcd774497540fceb409e7266c5be61167efd24e09"
)
EXPECTED_LEDGER_DIGEST: str | None = (
    "79c37deea9d831456e5ed913adadef7c465ff454a690fed5ff4aad7f0ed911a2"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path, expected: str, name: str):
    require(file_sha256(path) == expected, f"{path.name} changed")
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


require(file_sha256(ATLAS_OUTPUT_PATH) == ATLAS_OUTPUT_SHA256, "atlas output changed")
require(file_sha256(THM2895_PATH) == THM2895_SHA256, "THM-2895 changed")
require(file_sha256(THM2896_PATH) == THM2896_SHA256, "THM-2896 changed")
require(file_sha256(THM2897_PATH) == THM2897_SHA256, "THM-2897 changed")
require(
    file_sha256(BOOTSTRAP_OUTPUT_PATH) == BOOTSTRAP_OUTPUT_SHA256,
    "bootstrap output changed",
)
require(
    file_sha256(PARITY_OUTPUT_PATH) == PARITY_OUTPUT_SHA256,
    "parity output changed",
)
A = load(ATLAS_PATH, ATLAS_SHA256, "j6_k25_atlas")
M = load(BOOTSTRAP_PATH, BOOTSTRAP_SHA256, "j6_k25_bootstrap")
Q = load(PARITY_PATH, PARITY_SHA256, "j6_k25_parity")


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def ceiling(value: F) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def interval_mass(intervals: list[tuple[F, F]]) -> F:
    return sum((right - left for left, right in intervals), F(0))


def mask_labels(gate: tuple[int, ...], mask: int) -> tuple[int, ...]:
    return tuple(
        speed for bit, speed in enumerate(gate) if mask & (1 << bit)
    )


def scalar_closure(
    gate: tuple[int, ...],
    profiles: dict[int, dict[str, object]],
    start: int,
) -> tuple[
    int,
    tuple[tuple[int, ...], ...],
    tuple[tuple[F, ...], ...],
]:
    """Return the least scalar closure and deterministic simultaneous rounds."""

    index = {speed: bit for bit, speed in enumerate(gate)}
    active = start
    rounds: list[tuple[int, ...]] = []
    round_margins: list[tuple[F, ...]] = []
    while True:
        additions: list[int] = []
        margins: list[F] = []
        for bit, apex in enumerate(gate):
            if active & (1 << bit):
                continue
            values: list[F] = []
            for value, speed in profiles[apex]["ranked_candidates"]:
                speed_bit = index.get(speed)
                if speed_bit is not None and active & (1 << speed_bit):
                    continue
                values.append(value)
                if len(values) == 5:
                    break
            require(len(values) == 5, "activation candidate set too small")
            margin = profiles[apex]["m"] - sum(values, F(0))
            if margin > 0:
                additions.append(bit)
                margins.append(margin)
        if not additions:
            return active, tuple(rounds), tuple(round_margins)
        rounds.append(tuple(gate[bit] for bit in additions))
        round_margins.append(tuple(margins))
        for bit in additions:
            active |= 1 << bit


def minimum_closed_state_path(
    gate: tuple[int, ...],
    profiles: dict[int, dict[str, object]],
) -> dict[str, object]:
    """Exact unit-cost BFS on scalar-closed masks, stopping on first discovery."""

    full = (1 << len(gate)) - 1
    start, start_rounds, start_margins = scalar_closure(gate, profiles, 0)
    require(not start_rounds and not start_margins, "empty seed unexpectedly activates")
    queue = deque([start])
    distance = {start: 0}
    previous: dict[
        int,
        tuple[
            int,
            int,
            tuple[tuple[int, ...], ...],
            tuple[tuple[F, ...], ...],
        ],
    ] = {}
    cache = {
        start: (start, start_rounds, start_margins)
    }
    edges = 0
    found = start == full
    while queue and not found:
        state = queue.popleft()
        for bit, apex in enumerate(gate):
            if state & (1 << bit):
                continue
            raw = state | (1 << bit)
            if raw not in cache:
                cache[raw] = scalar_closure(gate, profiles, raw)
            target, rounds, margins = cache[raw]
            edges += 1
            if target in distance:
                continue
            distance[target] = distance[state] + 1
            previous[target] = (state, apex, rounds, margins)
            if target == full:
                found = True
                break
            queue.append(target)
    require(full in distance, "full gate is unreachable")

    steps: list[
        tuple[
            int,
            int,
            tuple[tuple[int, ...], ...],
            tuple[tuple[F, ...], ...],
        ]
    ] = []
    state = full
    while state != start:
        old, apex, rounds, margins = previous[state]
        steps.append((old, apex, rounds, margins))
        state = old
    steps.reverse()
    return {
        "minimum": distance[full],
        "states": len(distance),
        "edges": edges,
        "closure_calls": len(cache),
        "steps": tuple(steps),
    }


def available_q1(
    profile: dict[str, object],
    excluded: set[int],
) -> tuple[F, int]:
    for value, speed in profile["ranked_candidates"]:
        if speed not in excluded:
            return value, speed
    raise RuntimeError("empty apex candidate ledger")


def available_top(
    profile: dict[str, object],
    excluded: set[int],
    count: int,
) -> tuple[tuple[F, int], ...]:
    rows: list[tuple[F, int]] = []
    for value, speed in profile["ranked_candidates"]:
        if speed in excluded:
            continue
        rows.append((value, speed))
        if len(rows) == count:
            break
    require(len(rows) == count, "apex candidate ledger is too short")
    return tuple(rows)


def profile_ledger_line(
    profile: dict[str, object],
    gate: tuple[int, ...],
) -> str:
    return (
        f"E={','.join(map(str, profile['body']))};a={profile['apex']};"
        f"m={ftext(profile['m'])};r={profile['r']};"
        f"threshold={ftext(profile['threshold'])};"
        f"tail_first={profile['tail_first']};outside_top5="
        + ",".join(
            f"{speed}:{ftext(value)}"
            for value, speed in profile["outside_top5"]
        )
        + ";gate="
        + ",".join(
            f"{speed}:{ftext(profile['gate_coverages'][speed])}"
            for speed in gate
            if speed != profile["apex"]
        )
        + "\n"
    )


def rank_pair_state_audit(
    gate: tuple[int, ...],
    profiles: dict[int, dict[str, object]],
    carriers: dict[int, list[tuple[F, F]]],
    prefix: tuple[int, ...],
) -> tuple[list[dict[str, object]], tuple[int, ...]]:
    """Test THM-2897 on every inactive apex at one actual prefix."""

    excluded_prefix = set(prefix)
    rows: list[dict[str, object]] = []
    passing: list[int] = []
    for apex in gate:
        if apex in excluded_prefix:
            continue
        excluded = excluded_prefix | {apex}
        top5 = available_top(profiles[apex], excluded, 5)
        q5 = top5[4][0]
        witness_pair = (top5[0][1], top5[1][1])
        residual = Q.H.R.subtract_local_multi(
            carriers[apex],
            witness_pair,
        )
        residual_mass = interval_mass(residual)
        pair_union = profiles[apex]["m"] - residual_mass
        family = tuple(sorted((*BODY, apex, *witness_pair)))
        direct, direct_r, direct_m = Q.H.T.CORE.good_norm(family)
        require(
            residual == direct
            and len(residual) == direct_r
            and residual_mass == direct_m,
            f"rank-pair witness reconstruction changed: {apex}",
        )
        lower_margin = profiles[apex]["m"] - q5 - 2 * pair_union
        exact = None
        exact_margin = None
        if lower_margin > 0:
            exact = Q.H.pair_cap(carriers[apex], excluded)
            exact_margin = profiles[apex]["m"] - q5 - 2 * exact["cap"]
            status = "PASS" if exact_margin > 0 else "EXACT_FAIL"
            if exact_margin > 0:
                passing.append(apex)
        else:
            status = "LOWER_REFUTED"
        rows.append(
            {
                "prefix": prefix,
                "apex": apex,
                "q5_speed": top5[4][1],
                "q5": q5,
                "witness_pair": witness_pair,
                "witness_union": pair_union,
                "lower_margin": lower_margin,
                "status": status,
                "exact": exact,
                "exact_margin": exact_margin,
            }
        )
    return rows, tuple(passing)


def rank_pair_ledger_line(row: dict[str, object]) -> str:
    exact = row["exact"]
    exact_fields = (
        "none"
        if exact is None
        else (
            f"{ftext(exact['cap'])}:{ftext(row['exact_margin'])}:"
            f"{exact['maximizer'][0]},{exact['maximizer'][1]}:"
            f"{exact['paid']}:{exact['paid_digest']}"
        )
    )
    return (
        f"RP:P={','.join(map(str, row['prefix']))};a={row['apex']};"
        f"q5={row['q5_speed']}:{ftext(row['q5'])};"
        f"W={row['witness_pair'][0]},{row['witness_pair'][1]};"
        f"U={ftext(row['witness_union'])};"
        f"lower={ftext(row['lower_margin'])};"
        f"status={row['status']};exact={exact_fields}\n"
    )


def matching_certificate(
    profile: dict[str, object],
    carrier: list[tuple[F, F]],
    actual_prefix: tuple[int, ...],
    apex: int,
) -> dict[str, object]:
    """Lock the finite THM-2897 q5+M_(2,2) repair at apex 17."""

    require(apex == 17, "unexpected matching apex")
    reference = MATCHING_REFERENCE_PREFIX
    actual_excluded = set(actual_prefix) | {apex}
    excluded = set(reference)
    require(
        excluded <= actual_excluded,
        "matching reference prefix is not inherited by actual prefix",
    )
    direct, components, direct_mass = Q.H.T.CORE.good_norm(
        tuple(sorted((*BODY, apex)))
    )
    h = interval_mass(carrier)
    require(
        carrier == direct
        and len(carrier) == components == profile["r"]
        and h == direct_mass == profile["m"] > 0,
        "matching carrier reconstruction changed",
    )

    ranked, tail_single = Q.H.globally_ranked(carrier, excluded)
    q1, q1_speed = ranked[0]
    q5, q5_speed = ranked[4]
    pair = Q.H.pair_cap(carrier, excluded)
    require(
        pair["ranked"][:5] == ranked[:5]
        and pair["q1"] == q1
        and pair["tail_single"] == tail_single,
        "matching singleton and pair ledgers disagree",
    )
    b2 = pair["cap"]
    target = h - q5
    edge_floor = target - b2
    gamma = Q.H.S2 * components / 7
    delta = edge_floor - q1 - h / 7
    require(delta > 0, "matching finite-tail gap is not positive")
    cutoff = ceiling(gamma / delta) - 1
    require(
        cutoff >= Q.H.FIRST_EXTERNAL
        and q1 + h / 7 + gamma / (cutoff + 1) <= edge_floor,
        "matching endpoint cutoff did not seal",
    )

    by_speed = {speed: value for value, speed in ranked}
    require(
        cutoff <= Q.H.HORIZON,
        "matching cutoff unexpectedly exceeds imported exact horizon",
    )
    vertices = tuple(
        sorted(
            (
                (by_speed[speed], speed)
                for speed in range(Q.H.FIRST_EXTERNAL, cutoff + 1)
                if speed not in excluded
            ),
            key=lambda item: (-item[0], item[1]),
        )
    )
    vertex_hash = hashlib.sha256(b"LRC14/j6/K25/matching-vertices/v1\n")
    for value, speed in vertices:
        vertex_hash.update(f"V={speed};c={ftext(value)}\n".encode())

    candidates: list[tuple[F, int, int]] = []
    heavy_edges: list[tuple[F, int, int]] = []
    candidate_hash = hashlib.sha256(
        b"LRC14/j6/K25/matching-candidate-pairs/v1\n"
    )
    heavy_hash = hashlib.sha256(
        b"LRC14/j6/highK-matching-heavy-edges/v1\n"
    )
    for first_index, (first_value, first) in enumerate(vertices[:-1]):
        if first_value + vertices[first_index + 1][0] < edge_floor:
            break
        after_first = Q.H.R.subtract_local(carrier, first)
        for second_value, second in vertices[first_index + 1 :]:
            if first_value + second_value < edge_floor:
                break
            survivor = Q.H.R.subtract_local(after_first, second)
            weight = h - interval_mass(survivor)
            x, y = sorted((first, second))
            candidates.append((weight, x, y))
            candidate_hash.update(
                f"E={x},{y};W={ftext(weight)}\n".encode()
            )
            if weight >= edge_floor:
                reverse = Q.H.R.subtract_local(
                    Q.H.R.subtract_local(carrier, second),
                    first,
                )
                literal, literal_r, literal_m = Q.H.T.CORE.good_norm(
                    tuple(sorted((*BODY, apex, x, y)))
                )
                require(
                    survivor == reverse == literal
                    and len(survivor) == literal_r
                    and interval_mass(survivor) == literal_m,
                    f"matching edge reconstruction changed at {(x, y)}",
                )
                heavy_edges.append((weight, x, y))
                heavy_hash.update(
                    f"E={x},{y};W={ftext(weight)}\n".encode()
                )

    heavy_edges.sort(key=lambda edge: (-edge[0], edge[1], edge[2]))
    edge_pairs: list[
        tuple[F, tuple[int, int], tuple[int, int], bool]
    ] = []
    for first_index, first_edge in enumerate(heavy_edges[:-1]):
        for second_edge in heavy_edges[first_index + 1 :]:
            first_pair = (first_edge[1], first_edge[2])
            second_pair = (second_edge[1], second_edge[2])
            edge_pairs.append(
                (
                    first_edge[0] + second_edge[0],
                    first_pair,
                    second_pair,
                    set(first_pair).isdisjoint(second_pair),
                )
            )
    hostile_pairs = tuple(
        row for row in edge_pairs if row[0] >= target
    )
    hostile_matchings = tuple(row for row in hostile_pairs if row[3])
    require(
        not hostile_matchings,
        "matching threshold graph contains a hostile disjoint pair",
    )
    disjoint_pairs = tuple(row for row in edge_pairs if row[3])
    require(disjoint_pairs, "matching threshold graph has no control matching")
    maximum_disjoint = max(row[0] for row in disjoint_pairs)
    maximum_ties = tuple(
        row for row in disjoint_pairs if row[0] == maximum_disjoint
    )
    margin = target - maximum_disjoint
    require(margin > 0, "matching threshold-graph margin is not strict")

    weight_classes: dict[F, list[tuple[int, int]]] = {}
    for weight, first, second in heavy_edges:
        weight_classes.setdefault(weight, []).append((first, second))
    heavy_ties = tuple(
        (ftext(weight), tuple(edges))
        for weight, edges in sorted(
            weight_classes.items(),
            key=lambda item: (-item[0], item[1]),
        )
        if len(edges) > 1
    )
    return {
        "actual_prefix": actual_prefix,
        "reference_prefix": reference,
        "apex": apex,
        "h": h,
        "r": components,
        "q1": q1,
        "q1_speed": q1_speed,
        "q5": q5,
        "q5_speed": q5_speed,
        "B2": b2,
        "B2_maximizer": pair["maximizer"],
        "B2_paid": pair["paid"],
        "B2_paid_digest": pair["paid_digest"],
        "target": target,
        "edge_floor": edge_floor,
        "delta": delta,
        "cutoff": cutoff,
        "tail_single": tail_single,
        "vertex_count": len(vertices),
        "vertex_digest": vertex_hash.hexdigest(),
        "candidate_count": len(candidates),
        "candidate_digest": candidate_hash.hexdigest(),
        "heavy_edges": tuple(heavy_edges),
        "heavy_digest": heavy_hash.hexdigest(),
        "heavy_ties": heavy_ties,
        "edge_pairs": tuple(edge_pairs),
        "hostile_pairs": hostile_pairs,
        "maximum_disjoint": maximum_disjoint,
        "maximum_ties": maximum_ties,
        "margin": margin,
        "closed": True,
    }


def matching_result(row: dict[str, object]) -> tuple[object, ...]:
    return (
        row["actual_prefix"],
        row["reference_prefix"],
        row["apex"],
        ftext(row["h"]),
        row["r"],
        (row["q1_speed"], ftext(row["q1"])),
        (row["q5_speed"], ftext(row["q5"])),
        ftext(row["B2"]),
        row["B2_maximizer"],
        row["B2_paid"],
        row["B2_paid_digest"],
        ftext(row["target"]),
        ftext(row["edge_floor"]),
        ftext(row["delta"]),
        row["cutoff"],
        ftext(row["tail_single"]),
        row["vertex_count"],
        row["vertex_digest"],
        row["candidate_count"],
        row["candidate_digest"],
        tuple(
            (first, second, ftext(weight))
            for weight, first, second in row["heavy_edges"]
        ),
        row["heavy_digest"],
        row["heavy_ties"],
        tuple(
            (ftext(total), first, second, disjoint)
            for total, first, second, disjoint in row["edge_pairs"]
        ),
        tuple(
            (ftext(total), first, second, disjoint)
            for total, first, second, disjoint in row["hostile_pairs"]
        ),
        ftext(row["maximum_disjoint"]),
        tuple(
            (ftext(total), first, second, disjoint)
            for total, first, second, disjoint in row["maximum_ties"]
        ),
        ftext(row["margin"]),
    )


def matching_ledger_line(row: dict[str, object]) -> str:
    return f"MATCH:{matching_result(row)}\n"


def build_parity_branch(
    root: dict[str, object],
    gate: tuple[int, ...],
    profile: dict[str, object],
    prefix: tuple[int, ...],
    apex: int,
    parity_rank: int,
) -> tuple[
    dict[str, object],
    list[dict[str, object]],
    tuple[object, ...],
]:
    """Build and close one actual-prefix THM-2895 paid branch."""

    excluded = set(prefix) | {apex}
    carrier, direct_r, direct_m = M.M.T.CORE.good_norm(
        tuple(sorted((*BODY, apex)))
    )
    require(
        carrier == M.M.R.subtract_local(root["good"], apex)
        and len(carrier) == direct_r == profile["r"]
        and interval_mass(carrier) == direct_m == profile["m"] > 0,
        f"paid carrier reconstruction changed at apex {apex}",
    )
    ranked, tail_single = Q.H.globally_ranked(carrier, excluded)
    q1, q1_speed = ranked[0]
    require(
        (q1, q1_speed) == available_q1(profile, excluded),
        f"activation/parity q1 disagreement at apex {apex}",
    )
    singleton_margin = F(3, 7) * profile["m"] - q1
    require(singleton_margin > 0, f"H4 condition fails at apex {apex}")
    level = (profile["m"] - q1) / 4
    require(level > profile["m"] / 7, "H4 level is not finite")
    threshold = (
        Q.H.S2
        * profile["r"]
        / (7 * (level - profile["m"] / 7))
    )
    tail_first = max(Q.H.FIRST_EXTERNAL, ceiling(threshold))
    by_speed = {speed: value for value, speed in ranked}
    if tail_first > Q.H.HORIZON + 1:
        by_speed.update(
            {
                speed: value
                for value, speed in Q.H.T.coverages_many(
                    carrier,
                    [
                        speed
                        for speed in range(Q.H.HORIZON + 1, tail_first)
                        if speed not in excluded
                    ],
                )
            }
        )
    require(
        profile["m"] / 7
        + Q.H.S2 * profile["r"] / (7 * tail_first)
        <= level,
        "H4 tail did not seal",
    )
    core = tuple(
        sorted(
            speed
            for speed, value in by_speed.items()
            if speed not in excluded and value >= level
        )
    )
    branch = {
        "body": BODY,
        "rank": parity_rank,
        "root_rank": gate.index(apex) + 1,
        "apex": apex,
        "excluded_prefix": (*prefix, apex),
        "m": profile["m"],
        "r": profile["r"],
    }
    row = {
        "root": root,
        "branch": branch,
        "carrier": carrier,
        "excluded": excluded,
        "q1": q1,
        "tail_single": tail_single,
        "singleton_margin": singleton_margin,
        "level": level,
        "Htail": tail_first,
        "H": core,
    }
    pairs = [
        Q.pair_residual(row, hpair)
        for hpair in combinations(core, 2)
    ]
    require(
        all(pair["adaptive_closed"] for pair in pairs),
        f"hard H4-pair residual survives at apex {apex}",
    )
    result = (
        parity_rank,
        prefix,
        apex,
        gate.index(apex) + 1,
        ftext(profile["m"]),
        profile["r"],
        q1_speed,
        ftext(q1),
        ftext(singleton_margin),
        ftext(level),
        tail_first,
        core,
        len(pairs),
        sum(pair["direct_margin"] > 0 for pair in pairs),
        sum(pair["pair_margin"] > 0 for pair in pairs),
        sum(pair["adaptive_closed"] for pair in pairs),
        sum(not pair["adaptive_closed"] for pair in pairs),
        sum(pair["cap"]["paid"] for pair in pairs),
    )
    return row, pairs, result


def main() -> None:
    """Locked combined scalar/rank-pair/parity proof schedule."""

    base = A.profile_body(BODY)
    good, components, root_mass = A.T.CORE.good_norm(BODY)
    root = {**base, "good": good, "K": base["adaptive_k"]}
    require(
        components == root["r"] and root_mass == root["m"],
        "root reconstruction changed",
    )
    gate = tuple(speed for _, speed in root["top"][: root["K"]])
    require(gate == EXPECTED_GATE, "K25 gate changed")
    root_result = (
        ftext(root["m"]),
        root["r"],
        root["K"],
        ftext(root["adaptive_margin"]),
    )
    require(root_result == EXPECTED_ROOT, "K25 root profile changed")

    profiles = {
        apex: M.profile_apex(root, gate, apex) for apex in gate
    }
    carriers = {
        apex: M.M.T.CORE.good_norm(tuple(sorted((*BODY, apex))))[0]
        for apex in gate
    }
    profile_hash = hashlib.sha256(b"LRC14/j6/K25/apex-profiles/v1\n")
    for apex in gate:
        profile_hash.update(profile_ledger_line(profiles[apex], gate).encode())
    profile_digest = profile_hash.hexdigest()

    scalar_search = minimum_closed_state_path(gate, profiles)
    search_result = (
        scalar_search["minimum"],
        scalar_search["states"],
        scalar_search["edges"],
        scalar_search["closure_calls"],
    )
    scalar_path_result = tuple(
        (
            mask_labels(gate, state),
            apex,
            rounds,
            tuple(tuple(ftext(value) for value in row) for row in margins),
        )
        for state, apex, rounds, margins in scalar_search["steps"]
    )
    seed_bank = tuple(row[1] for row in scalar_search["steps"])

    index = {speed: bit for bit, speed in enumerate(gate)}
    full = (1 << len(gate)) - 1
    mask, initial_rounds, initial_margins = scalar_closure(gate, profiles, 0)
    require(
        not initial_rounds and not initial_margins,
        "empty prefix unexpectedly scalar-activates",
    )
    prefix: list[int] = []
    seed_cursor = 0
    parity_rows: list[dict[str, object]] = []
    all_pairs: list[tuple[dict[str, object], dict[str, object]]] = []
    branch_results: list[tuple[object, ...]] = []
    rank_audits: list[dict[str, object]] = []
    matching_rows: list[dict[str, object]] = []
    state_summaries: list[tuple[object, ...]] = []
    proof_schedule: list[tuple[object, ...]] = []

    while mask != full:
        require(
            set(prefix) == set(mask_labels(gate, mask)),
            "actual prefix and combined state disagree",
        )
        prior = tuple(prefix)
        if (
            17 not in prefix
            and set(MATCHING_REFERENCE_PREFIX) <= set(prefix) | {17}
        ):
            row = matching_certificate(
                profiles[17],
                carriers[17],
                prior,
                17,
            )
            require(row["closed"], "apex-17 matching certificate failed")
            matching_rows.append(row)
            prefix.append(17)
            mask |= 1 << index[17]
            target, rounds, margins = scalar_closure(gate, profiles, mask)
            for round_apices in rounds:
                for activated in round_apices:
                    require(
                        activated not in prefix,
                        "matching scalar round repeats an apex",
                    )
                    prefix.append(activated)
            proof_schedule.append(
                (
                    "MATCHING",
                    prior,
                    17,
                    hashlib.sha256(
                        matching_ledger_line(row).encode()
                    ).hexdigest(),
                    rounds,
                    tuple(
                        tuple(ftext(value) for value in margin_row)
                        for margin_row in margins
                    ),
                )
            )
            mask = target
            continue

        audit, passing = rank_pair_state_audit(
            gate,
            profiles,
            carriers,
            prior,
        )
        rank_audits.extend(audit)
        state_summaries.append(
            (
                prior,
                len(audit),
                sum(row["status"] == "LOWER_REFUTED" for row in audit),
                sum(row["status"] != "LOWER_REFUTED" for row in audit),
                passing,
            )
        )
        if passing:
            passing_data = tuple(
                (
                    row["apex"],
                    ftext(row["exact"]["cap"]),
                    ftext(row["exact_margin"]),
                    row["exact"]["maximizer"],
                    row["exact"]["paid"],
                )
                for row in audit
                if row["status"] == "PASS"
            )
            require(
                tuple(row[0] for row in passing_data) == passing,
                "rank-pair passing ledger changed order",
            )
            for apex in passing:
                require(apex not in prefix, "rank-pair repeats an apex")
                prefix.append(apex)
                mask |= 1 << index[apex]
            target, rounds, margins = scalar_closure(gate, profiles, mask)
            for round_apices in rounds:
                for apex in round_apices:
                    require(apex not in prefix, "scalar round repeats an apex")
                    prefix.append(apex)
            proof_schedule.append(
                (
                    "RANKPAIR",
                    prior,
                    passing_data,
                    rounds,
                    tuple(
                        tuple(ftext(value) for value in row)
                        for row in margins
                    ),
                )
            )
            mask = target
            continue

        while (
            seed_cursor < len(seed_bank)
            and mask & (1 << index[seed_bank[seed_cursor]])
        ):
            seed_cursor += 1
        require(
            seed_cursor < len(seed_bank),
            "scalar seed bank exhausted before combined closure",
        )
        apex = seed_bank[seed_cursor]
        seed_cursor += 1
        selected_audit = next(row for row in audit if row["apex"] == apex)
        require(
            selected_audit["status"] != "PASS",
            "paid parity apex still has a rank-pair certificate",
        )
        parity_rank = len(parity_rows) + 1
        row, pairs, result = build_parity_branch(
            root,
            gate,
            profiles[apex],
            prior,
            apex,
            parity_rank,
        )
        parity_rows.append(row)
        all_pairs.extend((row, pair) for pair in pairs)
        branch_results.append(result)

        prefix.append(apex)
        mask |= 1 << index[apex]
        target, rounds, margins = scalar_closure(gate, profiles, mask)
        for round_apices in rounds:
            for activated in round_apices:
                require(
                    activated not in prefix,
                    "scalar round repeats an apex",
                )
                prefix.append(activated)
        proof_schedule.append(
            (
                "PARITY",
                prior,
                apex,
                selected_audit["status"],
                rounds,
                tuple(
                    tuple(ftext(value) for value in row)
                    for row in margins
                ),
            )
        )
        mask = target

    require(
        len(prefix) == len(gate) and set(prefix) == set(gate),
        "combined proof schedule does not exhaust the gate",
    )
    hard = [
        (row, pair)
        for row, pair in all_pairs
        if not pair["adaptive_closed"]
    ]
    require(not hard, "combined proof path has a hard H4-pair residual")
    aggregate = (
        len(parity_rows),
        sum(step[0] == "RANKPAIR" for step in proof_schedule),
        sum(
            len(step[2])
            for step in proof_schedule
            if step[0] == "RANKPAIR"
        ),
        len(matching_rows),
        len(all_pairs),
        sum(pair["direct_margin"] > 0 for _, pair in all_pairs),
        sum(pair["pair_margin"] > 0 for _, pair in all_pairs),
        sum(pair["adaptive_closed"] for _, pair in all_pairs),
        len(hard),
        sum(pair["cap"]["paid"] for _, pair in all_pairs),
    )
    rankpair_aggregate = (
        len(state_summaries),
        len(rank_audits),
        sum(row["status"] == "LOWER_REFUTED" for row in rank_audits),
        sum(row["status"] != "LOWER_REFUTED" for row in rank_audits),
        sum(row["status"] == "PASS" for row in rank_audits),
        sum(row["status"] == "EXACT_FAIL" for row in rank_audits),
        sum(
            0 if row["exact"] is None else row["exact"]["paid"]
            for row in rank_audits
        ),
    )

    rank_hash = hashlib.sha256(b"LRC14/j6/K25/rank-pair-audit/v1\n")
    for row in rank_audits:
        rank_hash.update(rank_pair_ledger_line(row).encode())
    rankpair_digest = rank_hash.hexdigest()

    ledger = hashlib.sha256(b"LRC14/j6/K25/combined-closure/v3\n")
    ledger.update(
        (
            f"E={','.join(map(str, BODY))};root={root_result};gate={gate};"
            f"scalar_search={search_result};seed_bank={seed_bank};"
            f"scalar_path={scalar_path_result};"
            f"schedule={tuple(proof_schedule)};"
            f"ordered={tuple(prefix)};"
            f"states={tuple(state_summaries)}\n"
        ).encode()
    )
    for apex in gate:
        ledger.update(profile_ledger_line(profiles[apex], gate).encode())
    for row in rank_audits:
        ledger.update(rank_pair_ledger_line(row).encode())
    for row in matching_rows:
        ledger.update(matching_ledger_line(row).encode())
    for row in parity_rows:
        ledger.update(Q.branch_ledger_line(row).encode())
    for row, pair in all_pairs:
        ledger.update(Q.pair_ledger_line(row, pair).encode())
    ledger_digest = ledger.hexdigest()

    if EXPECTED_SEARCH is not None:
        require(search_result == EXPECTED_SEARCH, "scalar seed search changed")
    if EXPECTED_PATH is not None:
        require(scalar_path_result == EXPECTED_PATH, "scalar seed path changed")
    if EXPECTED_SEED_BANK is not None:
        require(seed_bank == EXPECTED_SEED_BANK, "hostile seed bank changed")
    if EXPECTED_PROOF_SCHEDULE is not None:
        require(
            tuple(proof_schedule) == EXPECTED_PROOF_SCHEDULE,
            "combined proof schedule changed",
        )
    if EXPECTED_STATE_SUMMARIES is not None:
        require(
            tuple(state_summaries) == EXPECTED_STATE_SUMMARIES,
            "combined state summaries changed",
        )
    if EXPECTED_BRANCH_RESULTS is not None:
        require(
            tuple(branch_results) == EXPECTED_BRANCH_RESULTS,
            "paid branch rows changed",
        )
    locked_matching_result = tuple(
        matching_result(row) for row in matching_rows
    )
    if EXPECTED_MATCHING_RESULT is not None:
        require(
            locked_matching_result == (EXPECTED_MATCHING_RESULT,),
            "matching certificate changed",
        )
    if EXPECTED_AGGREGATE is not None:
        require(aggregate == EXPECTED_AGGREGATE, "combined aggregate changed")
    if EXPECTED_RANKPAIR_AGGREGATE is not None:
        require(
            rankpair_aggregate == EXPECTED_RANKPAIR_AGGREGATE,
            "rank-pair aggregate changed",
        )
    if EXPECTED_PROFILE_DIGEST is not None:
        require(profile_digest == EXPECTED_PROFILE_DIGEST, "apex profiles changed")
    if EXPECTED_RANKPAIR_DIGEST is not None:
        require(
            rankpair_digest == EXPECTED_RANKPAIR_DIGEST,
            "rank-pair audit digest changed",
        )
    if EXPECTED_LEDGER_DIGEST is not None:
        require(ledger_digest == EXPECTED_LEDGER_DIGEST, "ledger digest changed")

    print("LRC14 j6 unique K25 five-parity matching closure proof")
    print(f"root={BODY};root_result={root_result}")
    print(f"gate={gate}")
    print(f"scalar_search={search_result};seed_bank={seed_bank}")
    proof_path = []
    for step in proof_schedule:
        if step[0] in {"PARITY", "MATCHING"}:
            activated = step[2]
        else:
            activated = tuple(row[0] for row in step[2])
        proof_path.append((step[0], step[1], activated))
    branch_summary = tuple(
        (
            row[0],
            row[2],
            len(row[11]),
            row[12],
            row[13],
            row[15],
            row[17],
        )
        for row in branch_results
    )
    print(f"proof_path={tuple(proof_path)}")
    print(f"state_summaries={tuple(state_summaries)}")
    print(f"ordered_prefix={tuple(prefix)}")
    print(f"branch_summary={branch_summary}")
    print(f"matching_summary={locked_matching_result}")
    print(f"aggregate={aggregate}")
    print(f"rankpair_aggregate={rankpair_aggregate}")
    print(f"apex_profile_sha256={profile_digest}")
    print(f"rankpair_audit_sha256={rankpair_digest}")
    print(f"canonical_ledger_sha256={ledger_digest}")
    if EXPECTED_PROOF_SCHEDULE is None:
        print(f"discovery_proof_schedule={tuple(proof_schedule)}")
    if EXPECTED_BRANCH_RESULTS is None:
        print(f"discovery_branch_results={tuple(branch_results)}")
    print(
        "mode="
        + (
            "DISCOVERY"
            if (
                EXPECTED_LEDGER_DIGEST is None
                or EXPECTED_MATCHING_RESULT is None
            )
            else "LOCKED"
        )
    )
    print("scope=one K25 root;not all-root j6;not LRC14")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
