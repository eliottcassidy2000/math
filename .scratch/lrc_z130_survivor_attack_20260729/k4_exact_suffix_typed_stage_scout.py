#!/usr/bin/env python3
"""Exact k=4 projected-suffix ordered typed-apex verifier.

It reconstructs the canonical 87,975 projected-suffix rows, partitions the
six-tail first-apex gate after ``z1`` into a bounded ``z2`` branch and a
complementary aligned-apex branch, and applies exact scalar/typed filters
before any three-drift projection scan.  The surviving rows are necessary
physical prefixes, not realized covers.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import inspect
import math
import multiprocessing as mp
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


def repo_root(start: Path) -> Path:
    for parent in (start, *start.parents):
        if (parent / "04-computation").is_dir() and (parent / "01-canon").is_dir():
            return parent
    raise RuntimeError("cannot locate repository root")


HERE = Path(__file__).resolve().parent
ROOT = repo_root(HERE)
GEOMETRIC_SOURCE = (
    ROOT / "04-computation/lrc14_j7_aligned_projected_arc_suffix_thm2941.py"
)
PROJECTED_SOURCE = (
    ROOT
    / "04-computation"
    / "lrc14_j7_five_aligned_two_drift_projected_closure_thm2941.py"
)
EXPECTED_GEOMETRIC_SHA256 = (
    "a003d287f618eb301edf6974d0b67dc128c4f380a169e7809ed5b5754e8b8303"
)
EXPECTED_PROJECTED_SHA256 = (
    "76f891edfcc029a08202481304a809e03e8bd81f247afaeabab685825c4d3662"
)
EXPECTED_SUPPORT_SHA256 = (
    "ae3377d67325f820e84a05d006c3c888f83472f1b99457f70237b966bad385b4"
)
DEFAULT_OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_four_aligned_three_drift_typed_apex_thm2941.out"
)
ETA = F(51, 1_183)
PROJECTED_ALPHA = F(2_366, 21_875)
EXCESS_IDENTITY_CONTROLS = (
    15,
    16,
    17,
    44,
    130,
    156,
    182,
    1_871,
    2_163,
    3_489,
    6_515,
    7_000,
)


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


G = load("k4_geometric_source", GEOMETRIC_SOURCE)
C = load("k4_projected_source", PROJECTED_SOURCE)
support_payload = (
    "\n".join(
        inspect.getsource(getattr(G, name))
        for name in ("suffix_upper", "top_insert")
    )
    + "\n"
    + "\n".join(
        inspect.getsource(getattr(C, name))
        for name in (
            "excess_numerator",
            "danger_fraction",
            "subtract_fraction",
            "interval_mass",
            "floor_fraction",
        )
    )
).encode()
def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def normalized_sha256(path: Path) -> str:
    payload = path.read_bytes()
    require(
        b"\r" not in payload.replace(b"\r\n", b""),
        f"{path.name}: lone carriage return",
    )
    return hashlib.sha256(payload.replace(b"\r\n", b"\n")).hexdigest()


require(
    normalized_sha256(GEOMETRIC_SOURCE) == EXPECTED_GEOMETRIC_SHA256,
    "pinned projected-suffix source changed",
)
require(
    normalized_sha256(PROJECTED_SOURCE) == EXPECTED_PROJECTED_SHA256,
    "pinned projected-closure source changed",
)
require(
    hashlib.sha256(support_payload).hexdigest() == EXPECTED_SUPPORT_SHA256,
    "pinned suffix/projected support changed",
)


def integer_apex(p: int, r: int, h: F) -> int:
    """Return the exact non-strict integer apex forced by a p-label cover.

    Suppose every label ``w`` obeys
    ``c_R(w) <= h/7 + gamma/w`` with ``gamma=6r/49``.  Put
    ``T=7*p*gamma/((7-p)*h)``.  If every one of the ``p`` integer labels
    exceeded ``floor(T)``, then each would exceed ``T`` (also when ``T`` is
    integral), so its contribution would be strictly below
    ``h/7+gamma/T``.  The p-sum would therefore be strictly below ``h``,
    contradicting that those labels cover R.  Thus some integer label is at
    most ``floor(T)=floor(6*p*r/[7*(7-p)*h])``.
    """

    require(1 <= p <= 6 and r > 0 and h > 0, "invalid apex parameters")
    return C.floor_fraction(F(6 * p * r, 7 * (7 - p)) / h)


def cached_excess_identity_audit(
    roots: tuple[tuple[int, ...], ...],
) -> tuple[int, str]:
    """Hostile-check the exact excess numerator reused by scalar screening."""

    digest = hashlib.sha256(b"LRC14/k4/cached-excess-identity/v1\n")
    comparisons = 0
    for body in roots:
        carrier = G.A.carrier_for(body)
        carrier_numerator = sum(right - left for left, right in carrier)
        h = F(carrier_numerator, G.A.RULER)
        canonical_l = 14 * math.lcm(*body)
        for speed in EXCESS_IDENTITY_CONTROLS:
            if speed % canonical_l == 0:
                continue
            delta = G.A.singleton_coverage(carrier, speed) - h / 7
            scaled = delta * 7 * G.A.RULER * speed
            require(
                scaled.denominator == 1,
                "hostile cached excess is not an integer numerator",
            )
            direct = C.excess_numerator(
                carrier,
                carrier_numerator,
                speed,
            )
            require(
                scaled.numerator == direct,
                "cached/direct exact excess identity failed",
            )
            digest.update(repr((body, speed, direct)).encode())
            digest.update(b"\n")
            comparisons += 1
    return comparisons, digest.hexdigest()


def suffix_rows(
    body: tuple[int, ...],
) -> tuple[tuple[dict[str, object], ...], dict[int, F]]:
    """Reconstruct the canonical k=4 rows and their already-paid excesses."""

    carrier = G.A.carrier_for(body)
    h = F(sum(right - left for left, right in carrier), G.A.RULER)
    components = len(carrier)
    canonical_l = 14 * math.lcm(*body)
    delta = {
        label: G.A.singleton_coverage(carrier, label) - h / 7
        for label in range(G.A.BASE_LABEL, G.HORIZON + 1)
        if label % canonical_l
    }
    cap_q = F(18 * components, 49) / (h * ETA)
    first_cap = cap_q.numerator // cap_q.denominator
    require(first_cap < G.HORIZON, "k4 first cap reached exact horizon")
    floor_q = PROJECTED_ALPHA * canonical_l
    require(
        floor_q > max(body),
        "projected wall does not force the third drift above the body",
    )
    high_floor = floor_q.numerator // floor_q.denominator + 1
    require(
        high_floor >= 15,
        "projected wall fell below the nonbody-label floor",
    )
    ordinary_tail = F(6 * components, 49 * (G.HORIZON + 1))
    high_tail = F(
        6 * components,
        49 * max(G.HORIZON + 1, high_floor),
    )

    arbitrary_top: list[tuple[F, int | None, str]] = []
    high_top: list[tuple[F, int | None, str]] = []
    rows: list[dict[str, object]] = []
    for first in range(G.HORIZON, G.A.BASE_LABEL - 1, -1):
        if first % canonical_l == 0:
            continue
        first_delta = delta[first]
        if first <= first_cap:
            suffix = G.suffix_upper(
                arbitrary_exact=arbitrary_top,
                high_exact=high_top,
                need=2,
                tail=ordinary_tail,
                high_tail=high_tail,
                constrained=first < high_floor,
            )
            lower = h * ETA
            upper = first_delta + sum(
                (value for value, _, _ in suffix),
                F(0),
            )
            if upper >= lower:
                rows.append(
                    {
                        "body": body,
                        "first": first,
                        "canonical_l": canonical_l,
                        "h": h,
                        "components": components,
                        "first_cap": first_cap,
                        "high_floor": high_floor,
                        "first_delta": first_delta,
                        "suffix": suffix,
                        "upper": upper,
                        "lower": lower,
                    }
                )
        item = (first_delta, first, "EXACT")
        G.top_insert(arbitrary_top, item, limit=4)
        if first >= high_floor:
            G.top_insert(high_top, item, limit=4)
    return tuple(sorted(rows, key=lambda row: row["first"])), delta


def profile(body: tuple[int, ...]) -> dict[str, object]:
    carrier_i = G.A.carrier_for(body)
    carrier_numerator = sum(right - left for left, right in carrier_i)
    carrier = tuple(
        (F(left, G.A.RULER), F(right, G.A.RULER))
        for left, right in carrier_i
    )
    rows, exact_delta = suffix_rows(body)
    if not rows:
        return {
            "body": body,
            "rows": (),
            "typed": (),
            "drift_candidates": 0,
            "scalar_killed": 0,
            "scalar_finite": 0,
            "no_positive_gap": 0,
            "aligned_states": (),
            "recursive_states": (),
            "z2_interval_rows": (),
            "merged_candidates": 0,
            "merged_scalar_killed": 0,
            "merged_scalar_finite": 0,
            "merged_no_positive_gap": 0,
        }

    canonical_l = rows[0]["canonical_l"]
    high_floor = rows[0]["high_floor"]
    typed: list[tuple[object, ...]] = []
    drift_candidates = scalar_killed = scalar_finite = no_positive_gap = 0
    aligned_states: list[tuple[object, ...]] = []
    recursive_states: list[tuple[object, ...]] = []
    z2_interval_rows: list[tuple[object, ...]] = []
    merged_candidates = 0
    merged_scalar_killed = 0
    merged_scalar_finite = 0
    merged_no_positive_gap = 0
    numerator_cache: dict[int, int] = {}
    danger_cache: dict[int, tuple[tuple[F, F], ...]] = {}

    def numerator(speed: int) -> int:
        if speed not in numerator_cache:
            if speed in exact_delta:
                scaled = (
                    exact_delta[speed]
                    * 7
                    * G.A.RULER
                    * speed
                )
                require(
                    scaled.denominator == 1,
                    "cached exact excess lost its integer numerator",
                )
                numerator_cache[speed] = scaled.numerator
            else:
                numerator_cache[speed] = C.excess_numerator(
                    carrier_i,
                    carrier_numerator,
                    speed,
                )
        return numerator_cache[speed]

    def danger(speed: int) -> tuple[tuple[F, F], ...]:
        if speed not in danger_cache:
            danger_cache[speed] = C.danger_fraction(speed)
        return danger_cache[speed]

    for row in rows:
        first = row["first"]
        residual1 = C.subtract_fraction(carrier, danger(first))
        h1 = C.interval_mass(residual1)
        r1 = len(residual1)
        apex1 = integer_apex(6, r1, h1)
        aligned_cap1 = apex1 // canonical_l
        typed.append(
            (
                body,
                first,
                canonical_l,
                high_floor,
                h1,
                r1,
                apex1,
                aligned_cap1,
            )
        )

        z2_intervals: list[tuple[int, int]] = []
        if first + 1 <= apex1:
            z2_intervals.append((first + 1, apex1))
        base_gap = row["lower"] - row["first_delta"]

        def scalar_type(second: int) -> int:
            """0=impossible, 1=finite z3 cap, 2=no positive scalar gap.

            Type 2 says only that the scalar excess inequality supplies no
            upper bound for ``z3``.  It does not assert that the eventual
            third-drift excess is nonnegative.
            """

            third_floor = max(second + 1, high_floor)
            second_numerator = numerator(second)
            if (
                second_numerator * base_gap.denominator
                >= base_gap.numerator * 7 * G.A.RULER * second
            ):
                return 2
            left = (
                7 * third_floor * second_numerator
                + 6 * row["components"] * G.A.RULER * second
            ) * base_gap.denominator
            right = (
                base_gap.numerator
                * 49
                * G.A.RULER
                * second
                * third_floor
            )
            return 0 if left < right else 1

        # Branch 2: z2>A1 forces a small aligned apex.  Choose the least
        # remaining actual aligned multiplier; subsequent choices are
        # strictly increasing, which retains THM-2893's forbidden-prefix
        # sidecar and makes the recursion disjoint.  Repeat until z2 is
        # bounded or all four aligned labels have been deleted.
        queue = [
            ((multiplier,), apex1 + 1)
            for multiplier in range(1, aligned_cap1 + 1)
        ]
        seen: set[tuple[tuple[int, ...], int]] = set()
        residual_cache: dict[
            tuple[int, ...],
            tuple[tuple[F, F], ...],
        ] = {(): residual1}
        while queue:
            selected, second_floor = queue.pop()
            state_key = (selected, second_floor)
            if state_key in seen:
                continue
            seen.add(state_key)
            if selected in residual_cache:
                residual2 = residual_cache[selected]
            else:
                residual2 = residual1
                for multiplier in selected:
                    residual2 = C.subtract_fraction(
                        residual2,
                        danger(canonical_l * multiplier),
                    )
                residual_cache[selected] = residual2
            h2 = C.interval_mass(residual2)
            r2 = len(residual2)
            require(h2 > 0 and r2 > 0, "aligned recursion erased residual")
            remaining = 6 - len(selected)
            require(2 <= remaining <= 5, "typed recursion depth escaped")
            apex2 = integer_apex(remaining, r2, h2)
            drift_possible = second_floor <= apex2
            if drift_possible:
                z2_intervals.append((second_floor, apex2))
            unused = tuple(
                multiplier
                for multiplier in range(1, apex2 // canonical_l + 1)
                if multiplier > selected[-1]
            )
            complement_closed = len(selected) == 4 or not unused
            recursive_states.append(
                (
                    body,
                    first,
                    canonical_l,
                    selected,
                    second_floor,
                    remaining,
                    h2,
                    r2,
                    apex2,
                    drift_possible,
                    unused,
                    complement_closed,
                )
            )
            if len(selected) == 1:
                aligned_states.append(
                    (
                        body,
                        first,
                        canonical_l,
                        high_floor,
                        apex1,
                        selected[0],
                        h2,
                        r2,
                        apex2,
                        second_floor,
                        apex2 < min(canonical_l, second_floor),
                        second_floor <= apex2 < canonical_l,
                        canonical_l <= apex2,
                        apex2 // canonical_l,
                    )
                )
            if complement_closed:
                continue
            next_floor = max(second_floor, apex2 + 1)
            for multiplier in unused:
                child = (*selected, multiplier)
                queue.append((child, next_floor))

        merged: list[list[int]] = []
        for left, right in sorted(z2_intervals):
            if left > right:
                continue
            if merged and left <= merged[-1][1] + 1:
                merged[-1][1] = max(merged[-1][1], right)
            else:
                merged.append([left, right])
        merged_tuple = tuple((left, right) for left, right in merged)
        pair_count = sum(
            right - left + 1 - (right // canonical_l - (left - 1) // canonical_l)
            for left, right in merged_tuple
        )
        row_killed = row_finite = row_no_gap = 0
        row_max_survivor: int | None = None
        for left, right in merged_tuple:
            for second in range(left, right + 1):
                if second % canonical_l == 0:
                    continue
                require(
                    first < second and second % canonical_l,
                    "merged physical drift ordering/alignment changed",
                )
                kind = scalar_type(second)
                if second <= apex1:
                    # The first bounded-apex branch is a prefix of the
                    # merged recursion bank, so account for it during this
                    # single scalar pass rather than integrating twice.
                    drift_candidates += 1
                    if kind == 2:
                        no_positive_gap += 1
                    elif kind == 0:
                        scalar_killed += 1
                    else:
                        scalar_finite += 1
                if kind == 0:
                    row_killed += 1
                elif kind == 1:
                    row_finite += 1
                    row_max_survivor = second
                else:
                    row_no_gap += 1
                    row_max_survivor = second
        require(
            row_killed + row_finite + row_no_gap == pair_count,
            "merged scalar partition changed",
        )

        def exact_z3_cap(second: int, kind: int) -> int | None:
            if kind != 1:
                return None
            second_delta = F(
                numerator(second),
                7 * G.A.RULER * second,
            )
            gap = base_gap - second_delta
            require(gap > 0, "finite scalar class lost its positive gap")
            cap = C.floor_fraction(F(6 * row["components"], 49) / gap)
            require(
                max(second + 1, high_floor) <= cap,
                "finite scalar cap fell below its required floor",
            )
            return cap

        raw_max = merged_tuple[-1][1]
        raw_kind = scalar_type(raw_max)
        row_max_kind = (
            None if row_max_survivor is None else scalar_type(row_max_survivor)
        )
        merged_candidates += pair_count
        merged_scalar_killed += row_killed
        merged_scalar_finite += row_finite
        merged_no_positive_gap += row_no_gap
        z2_interval_rows.append(
            (
                body,
                first,
                canonical_l,
                merged_tuple,
                pair_count,
                row_killed,
                row_finite,
                row_no_gap,
                row_max_survivor,
                row_max_kind,
                (
                    None
                    if row_max_survivor is None
                    else exact_z3_cap(row_max_survivor, row_max_kind)
                ),
                raw_max,
                raw_kind,
                exact_z3_cap(raw_max, raw_kind),
                high_floor,
            )
        )

    return {
        "body": body,
        "rows": rows,
        "typed": typed,
        "drift_candidates": drift_candidates,
        "scalar_killed": scalar_killed,
        "scalar_finite": scalar_finite,
        "no_positive_gap": no_positive_gap,
        "aligned_states": aligned_states,
        "recursive_states": recursive_states,
        "z2_interval_rows": z2_interval_rows,
        "merged_candidates": merged_candidates,
        "merged_scalar_killed": merged_scalar_killed,
        "merged_scalar_finite": merged_scalar_finite,
        "merged_no_positive_gap": merged_no_positive_gap,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(8, mp.cpu_count() or 1))
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    roots = tuple(combinations(range(1, 15), 6))
    if args.workers == 1:
        profiles = [profile(body) for body in roots]
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            profiles = list(pool.imap(profile, roots, chunksize=1))
    profiles.sort(key=lambda row: row["body"])
    excess_audit_count, excess_audit_digest = cached_excess_identity_audit(roots)
    rows = sorted(
        (row for item in profiles for row in item["rows"]),
        key=lambda row: (row["body"], row["first"]),
    )
    typed = sorted(row for item in profiles for row in item["typed"])
    aligned_states = sorted(
        row for item in profiles for row in item["aligned_states"]
    )
    recursive_states = sorted(
        row for item in profiles for row in item["recursive_states"]
    )
    z2_interval_rows = sorted(
        row for item in profiles for row in item["z2_interval_rows"]
    )
    aligned_types = Counter(
        "closed" if row[10] else "z2_forced" if row[11] else "aligned_possible"
        for row in aligned_states
    )
    digest = hashlib.sha256(
        b"LRC14/k4/exact-projected-suffix/typed-stage/v2\n"
        + repr(tuple(typed)).encode()
        + b"\n"
        + repr(tuple(aligned_states)).encode()
        + b"\n"
        + repr(tuple(recursive_states)).encode()
        + b"\n"
        + repr(tuple(z2_interval_rows)).encode()
    ).hexdigest()
    depth_counts = dict(
        sorted(Counter(len(row[3]) for row in recursive_states).items())
    )
    bounded_depth_counts = dict(
        sorted(
            Counter(len(row[3]) for row in recursive_states if row[9]).items()
        )
    )
    closed_depth_counts = dict(
        sorted(
            Counter(len(row[3]) for row in recursive_states if row[11]).items()
        )
    )
    suffix_roots = len({row["body"] for row in rows})
    maximum_z1 = max(row["first"] for row in rows)
    first_apex_max = max(row[6] for row in typed)
    aligned_first_rows = sum(row[7] > 0 for row in typed)
    aligned_first_states = sum(row[7] for row in typed)
    initial_candidates = sum(row["drift_candidates"] for row in profiles)
    initial_killed = sum(row["scalar_killed"] for row in profiles)
    initial_finite = sum(row["scalar_finite"] for row in profiles)
    initial_no_gap = sum(row["no_positive_gap"] for row in profiles)
    p5_max_apex = max(row[8] for row in aligned_states)
    p5_multiplier_max = max(row[13] for row in aligned_states)
    raw_cap = max(right for row in z2_interval_rows for _, right in row[3])
    raw_pairs = sum(row[4] for row in z2_interval_rows)
    final_finite = sum(row["merged_scalar_finite"] for row in profiles)
    final_no_gap = sum(row["merged_no_positive_gap"] for row in profiles)
    final_candidates = sum(row["merged_candidates"] for row in profiles)
    final_killed = sum(row["merged_scalar_killed"] for row in profiles)
    screened_cap = max(
        row[8] for row in z2_interval_rows if row[8] is not None
    )
    surviving_first_rows = sum(row[8] is not None for row in z2_interval_rows)
    killed_first_rows = sum(row[8] is None for row in z2_interval_rows)
    scalar_class = {
        0: "scalar_killed",
        1: "finite_z3_cap",
        2: "no_positive_gap_z3_cap",
    }
    raw_cap_witnesses = tuple(
        (
            row[0],
            row[1],
            row[11],
            row[2],
            max(row[11] + 1, row[14]),
            scalar_class[row[12]],
            row[13],
        )
        for row in z2_interval_rows
        if row[11] == raw_cap
    )
    screened_cap_witnesses = tuple(
        (
            row[0],
            row[1],
            row[8],
            row[2],
            max(row[8] + 1, row[14]),
            scalar_class[row[9]],
            row[10],
        )
        for row in z2_interval_rows
        if row[8] == screened_cap
    )
    cap_witness_digest = hashlib.sha256(
        b"LRC14/k4/raw-screened-cap-witness/v1\n"
        + repr((raw_cap_witnesses, screened_cap_witnesses)).encode()
    ).hexdigest()
    minimum_projected_wall = min(
        (PROJECTED_ALPHA * row[2], row[0], row[14])
        for row in z2_interval_rows
    )

    require(len(rows) == 87_975, "exact projected-suffix bank changed")
    require(
        excess_audit_count == 36_036
        and excess_audit_digest
        == "95c91236dfb3e99ac2a28b4e030e1397dfcf046e5debda329a801f1a6300f544",
        "cached/direct exact excess hostile audit changed",
    )
    require(suffix_roots == 3_003 and maximum_z1 == 182, "suffix support changed")
    require(
        (first_apex_max, aligned_first_rows, aligned_first_states)
        == (1_871, 706, 758),
        "first-apex split changed",
    )
    require(
        (initial_candidates, initial_killed, initial_finite, initial_no_gap)
        == (49_129_013, 40_630_772, 2_032_220, 6_466_021),
        "initial drift-apex scalar split changed",
    )
    require(
        len(aligned_states) == 758
        and aligned_types == Counter({"aligned_possible": 758})
        and p5_max_apex == 3_489
        and p5_multiplier_max == 9,
        "p5 aligned-apex state changed",
    )
    require(
        len(recursive_states) == 5_368
        and depth_counts == {1: 758, 2: 1_052, 3: 1_427, 4: 2_131}
        and bounded_depth_counts == {1: 758, 2: 990, 3: 1_293, 4: 543}
        and closed_depth_counts == {2: 54, 3: 684, 4: 2_131},
        "ordered aligned recursion changed",
    )
    require(
        all(
            all(left < right for left, right in zip(row[3], row[3][1:]))
            for row in recursive_states
        ),
        "selected aligned multipliers lost strict ordering",
    )
    require(
        len(
            {
                (row[0], row[1], row[3], row[4])
                for row in recursive_states
            }
        )
        == len(recursive_states),
        "ordered recursion state keys are not unique",
    )
    require(
        raw_cap == 6_515
        and raw_pairs == 50_285_016
        and len(z2_interval_rows) == 87_975
        and max(len(row[3]) for row in z2_interval_rows) == 1,
        "raw recursive z2 bank changed",
    )
    require(
        all(
            len(row[3]) == 1 and row[3][0][0] == row[1] + 1
            for row in z2_interval_rows
        ),
        "merged z2 bank lost its single prefix interval",
    )
    require(
        (final_candidates, final_killed, final_finite, final_no_gap)
        == (50_285_016, 41_770_842, 2_042_669, 6_471_505),
        "full recursive scalar split changed",
    )
    require(
        screened_cap == 2_163
        and final_finite + final_no_gap == 8_514_174
        and surviving_first_rows == 87_975
        and killed_first_rows == 0,
        "screened necessary pair bank changed",
    )
    require(
        minimum_projected_wall
        == (F(56_784, 3_125), (1, 2, 3, 4, 6, 12), 19),
        "minimum projected-wall witness changed",
    )
    require(
        raw_cap_witnesses
        == (
            (
                (2, 3, 4, 8, 9, 12),
                156,
                6_515,
                1_008,
                6_516,
                "scalar_killed",
                None,
            ),
        )
        and screened_cap_witnesses
        == (
            (
                (2, 3, 4, 6, 8, 12),
                44,
                2_163,
                336,
                2_164,
                "finite_z3_cap",
                23_477,
            ),
        )
        and cap_witness_digest
        == "9c3469a3600a928aebcdae07d51a1c3f16fc221740505502250e59ca716efefe",
        "raw/screened z2-cap witness ledger changed",
    )
    require(
        all(
            witness[1] < witness[2]
            and witness[2] % witness[3] != 0
            for witness in raw_cap_witnesses + screened_cap_witnesses
        ),
        "cap witness lost physical order or became aligned",
    )
    require(
        digest == "5303930eea27eb9713a81245efa0a975f7e2b21869edcd028e4b2f2bb56f8fb3",
        "typed-stage semantic digest changed",
    )

    lines = [
        "LRC14 k4 exact projected-suffix ordered typed-apex verifier",
        f"geometric_source_sha256={normalized_sha256(GEOMETRIC_SOURCE)}",
        f"projected_source_sha256={normalized_sha256(PROJECTED_SOURCE)}",
        f"support_sha256={hashlib.sha256(support_payload).hexdigest()}",
        f"cached_excess_identity_controls={EXCESS_IDENTITY_CONTROLS};"
        f"comparisons={excess_audit_count};digest={excess_audit_digest}",
        "physical_order=z1<z2<z3;"
        "denominator_sidecar=d(z)=L/gcd(z,L),not physically ordered;"
        "diagonal denominator triples are independently empty",
        "integer_apex_lemma=A_p(R)=floor(6*p*r/[7*(7-p)*h]);"
        "if_all_integer_w>A_then_w>T_and_the_p_sum<h",
        "projected_wall=alpha4*L>max(E),"
        "so_z3>=floor(alpha4*L)+1;"
        f"minimum_wall={minimum_projected_wall}",
        f"suffix_rows={len(rows)};suffix_roots={suffix_roots};"
        f"maximum_z1={maximum_z1}",
        f"first_apex_max={first_apex_max};"
        f"aligned_first_rows={aligned_first_rows};"
        f"aligned_first_states={aligned_first_states}",
        f"initial_drift_apex_candidates={initial_candidates};"
        f"initial_scalar_killed={initial_killed};"
        f"initial_scalar_finite_z3={initial_finite};"
        f"initial_no_positive_gap_z3_cap={initial_no_gap}",
        f"complementary_aligned_apex_states={len(aligned_states)};"
        f"p5_types={dict(sorted(aligned_types.items()))};"
        f"p5_max_apex={p5_max_apex};"
        f"p5_aligned_multiplier_max={p5_multiplier_max}",
        f"recursive_aligned_candidate_states={len(recursive_states)};"
        f"states_by_selected_depth={depth_counts};"
        f"bounded_z2_branches_by_depth={bounded_depth_counts};"
        f"closed_complements_by_depth={closed_depth_counts}",
        f"raw_uniform_z2_cap={raw_cap};raw_recursive_pair_bank={raw_pairs};"
        f"interval_rows={len(z2_interval_rows)};"
        f"maximum_intervals_per_row={max(len(row[3]) for row in z2_interval_rows)}",
        f"full_recursive_scalar_candidates={final_candidates};"
        f"full_recursive_scalar_killed={final_killed};"
        f"full_recursive_finite_z3={final_finite};"
        f"full_recursive_no_positive_gap_z3_cap={final_no_gap}",
        f"screened_uniform_z2_cap={screened_cap};"
        f"necessary_non_aligned_(E,z1,z2)_bank={final_finite+final_no_gap};"
        f"surviving_(E,z1)_rows={surviving_first_rows};"
        f"killed_(E,z1)_rows={killed_first_rows}",
        "cap_witness_schema=(E,z1,z2,L,z3_floor,scalar_class,z3_cap)",
        f"raw_cap_witnesses={raw_cap_witnesses}",
        f"screened_cap_witnesses={screened_cap_witnesses}",
        f"cap_witness_digest={cap_witness_digest}",
        f"typed_stage_digest={digest}",
        "all_exact_controls=PASS",
    ]
    output = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
