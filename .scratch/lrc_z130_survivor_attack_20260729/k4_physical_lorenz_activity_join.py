#!/usr/bin/env python3
"""Exploratory physical join for the k=4 scalar-surviving pair bank.

This is intentionally a scratch audit.  It joins the exact ordered physical
``(E,z1,z2)`` bank from ``k4_exact_suffix_typed_stage_scout.py`` to the
incoming Lorenz+activity residual of sorted denominator multisets.  A sorted
multiset is not treated as physical order: all ordered choices of two
positions are admitted as ``(d(z1),d(z2))``, with the third position retained
as a possible ``d(z3)``.

For a positive remaining scalar gap the join additionally requires an actual
integer ``z3`` in the exact interval

    [max(z2+1, floor(alpha4*L)+1), floor(6*r_E/(49*gap))]

having one of the permitted remaining denominators.  When the scalar gap is
nonpositive, only denominator existence is asserted.  The result is a
necessary screen, never a realized-cover census.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import math
import multiprocessing as mp
from collections import defaultdict
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCOUT_SOURCE = HERE / "k4_exact_suffix_typed_stage_scout.py"
DEFAULT_RESIDUAL = Path(
    "/private/tmp/math-wt-hunter-star-resume/.scratch/"
    "lrc14_three_drift_lorenz_activity_combined_residual.tsv"
)
EXPECTED_RESIDUAL_SHA256 = (
    "fceb560b80bbfd8f9995ffbfbbc8eede42c84503b89744fe32035a8adc576f33"
)
EXPECTED_CORRECTED_RESIDUAL_SHA256 = (
    "0c7bad48d32e3b8bce6fa0dd3b3f1d82d7f8ff31272c972c4f0035a0dfd335ef"
)
EXPECTED_SCOUT_SHA256 = ""
CORRECTED_DELETIONS = {
    (
        "1,4,5,7,9,11\t194040\t194040\t55392\t2\t2\t194040"
        "\t97020\t2\t3\t1\t2\t5"
    ),
    (
        "1,5,7,8,9,11\t388080\t388080\t109044\t2\t2\t388080"
        "\t194040\t2\t3\t1\t2\t5"
    ),
}
HEADER = (
    "F\tL\tD\tsupport\td1\td2\td3\tbest_q\tfiber_height"
    "\tmin_slack\tbest_s\tsupport_top_s\tneedle_top_s"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


K = load("k4_exact_typed_stage_join_source", SCOUT_SOURCE)


def residual_ledger(
    path: Path,
) -> tuple[
    dict[tuple[int, ...], dict[tuple[int, int], tuple[int, ...]]],
    tuple[int, int, int, int, int, int],
]:
    """Delete two audited bad rows and build multiset-deletion links."""

    require(sha256(path) == EXPECTED_RESIDUAL_SHA256, "residual TSV changed")
    expanded: dict[
        tuple[int, ...],
        dict[tuple[int, int], set[int]],
    ] = defaultdict(lambda: defaultdict(set))
    shapes: set[tuple[int, int, int, int]] = set()
    bodies: set[tuple[int, ...]] = set()
    body_divisors: set[tuple[tuple[int, ...], int]] = set()
    raw_lines = path.read_text(encoding="utf-8").splitlines()
    require(raw_lines and raw_lines[0] == HEADER, "residual TSV header changed")
    require(
        all(raw_lines.count(line) == 1 for line in CORRECTED_DELETIONS),
        "audited residual deletion multiplicity changed",
    )
    lines = [
        line
        for line in raw_lines
        if line not in CORRECTED_DELETIONS
    ]
    corrected_payload = ("\n".join(lines) + "\n").encode()
    require(
        hashlib.sha256(corrected_payload).hexdigest()
        == EXPECTED_CORRECTED_RESIDUAL_SHA256,
        "corrected residual byte stream changed",
    )
    for line in lines[1:]:
        fields = line.split("\t")
        require(len(fields) == 13, "residual TSV row width changed")
        body = tuple(map(int, fields[0].split(",")))
        canonical_l = int(fields[1])
        divisor_lcm = int(fields[2])
        denominators = tuple(map(int, fields[4:7]))
        require(
            len(body) == 6
            and tuple(sorted(body)) == body
            and canonical_l == 14 * math.lcm(*body),
            "residual physical body/L typing changed",
        )
        require(
            denominators == tuple(sorted(denominators))
            and all(
                d > 1 and canonical_l % d == 0
                for d in denominators
            )
            and math.lcm(*denominators) == divisor_lcm,
            "residual denominator-multiset typing changed",
        )
        bodies.add(body)
        body_divisors.add((body, divisor_lcm))
        shapes.add((divisor_lcm, *denominators))
        # Delete each possible third position from the sorted multiset.
        # Repeated values are deduplicated.  The surviving two-position
        # multiset is compared with sorted(d(z1),d(z2)); ledger slots carry
        # no physical ownership.
        for third_index in range(3):
            pair = tuple(
                sorted(
                    denominators[:third_index]
                    + denominators[third_index + 1 :]
                )
            )
            expanded[body][pair].add(denominators[third_index])
    frozen = {
        body: {
            pair: tuple(sorted(thirds))
            for pair, thirds in sorted(pairs.items())
        }
        for body, pairs in sorted(expanded.items())
    }
    require(
        len(lines) - 1 == 29_219
        and len(shapes) == 4_970
        and len(bodies) == 2_878
        and len(body_divisors) == 2_974,
        "residual TSV census changed",
    )
    link_count = sum(
        len(thirds)
        for pairs in frozen.values()
        for thirds in pairs.values()
    )
    key_count = sum(len(pairs) for pairs in frozen.values())
    require(
        key_count == 44_933 and link_count == 79_329,
        "corrected denominator deletion index changed",
    )
    return frozen, (
        len(lines) - 1,
        len(shapes),
        len(bodies),
        len(body_divisors),
        key_count,
        link_count,
    )


RESIDUAL_PATH = DEFAULT_RESIDUAL
ALLOWED, RESIDUAL_CENSUS = residual_ledger(RESIDUAL_PATH)


def first_label_with_denominator(
    canonical_l: int,
    denominator: int,
    lower: int,
    upper: int,
) -> int | None:
    """Return the first z in [lower,upper] with L/gcd(z,L)=denominator."""

    require(
        canonical_l % denominator == 0 and denominator > 1,
        "invalid requested physical denominator",
    )
    scale = canonical_l // denominator
    first_u = (lower + scale - 1) // scale
    last_u = upper // scale
    for unit in range(first_u, last_u + 1):
        if math.gcd(unit, denominator) == 1:
            label = scale * unit
            require(
                lower <= label <= upper
                and canonical_l // math.gcd(label, canonical_l)
                == denominator,
                "physical denominator constructor failed",
            )
            return label
    return None


def profile_join(body: tuple[int, ...]) -> tuple[object, ...]:
    """Join one body and return a compact independently hashable ledger."""

    result = K.profile(body)
    suffix_by_first = {
        row["first"]: row
        for row in result["rows"]
    }
    interval_rows = sorted(
        result["z2_interval_rows"],
        key=lambda row: row[1],
    )
    require(
        len(suffix_by_first) == len(interval_rows),
        "suffix/interval row pairing changed",
    )
    allowed = ALLOWED.get(body, {})
    carrier_i = K.G.A.carrier_for(body)
    carrier_numerator = sum(right - left for left, right in carrier_i)
    numerator_cache: dict[int, int] = {}

    def numerator(speed: int) -> int:
        if speed not in numerator_cache:
            numerator_cache[speed] = K.C.excess_numerator(
                carrier_i,
                carrier_numerator,
                speed,
            )
        return numerator_cache[speed]

    scalar_finite = scalar_no_gap = 0
    prefix_finite = prefix_no_gap = 0
    retained_finite = retained_no_gap = 0
    prefix_cap = retained_cap = finite_cap = no_gap_cap = 0
    prefix_rows = retained_rows = 0
    prefix_digest = hashlib.sha256(
        b"LRC14/k4/lorenz-activity/physical-prefix/body/v1\n"
    )
    retained_digest = hashlib.sha256(
        b"LRC14/k4/lorenz-activity/physical-realizability/body/v1\n"
    )

    for interval_row in interval_rows:
        first = interval_row[1]
        canonical_l = interval_row[2]
        high_floor = interval_row[14]
        row = suffix_by_first[first]
        base_gap = row["lower"] - row["first_delta"]
        row_prefix = row_retained = False
        require(
            len(interval_row[3]) == 1
            and interval_row[3][0][0] == first + 1,
            "joined z2 row is not a single prefix interval",
        )
        for left, right in interval_row[3]:
            for second in range(left, right + 1):
                if second % canonical_l == 0:
                    continue
                second_numerator = numerator(second)
                second_delta = Q(
                    second_numerator,
                    7 * K.G.A.RULER * second,
                )
                gap = base_gap - second_delta
                third_floor = max(second + 1, high_floor)

                # Independent exact cross-check against the source's
                # cross-multiplied three-way scalar classification.
                no_gap_cross = (
                    second_numerator * base_gap.denominator
                    >= (
                        base_gap.numerator
                        * 7
                        * K.G.A.RULER
                        * second
                    )
                )
                require(no_gap_cross == (gap <= 0), "gap-sign cross-check failed")
                if gap <= 0:
                    scalar_class = "no_positive_gap_z3_cap"
                    scalar_no_gap += 1
                    upper: int | None = None
                else:
                    upper = K.C.floor_fraction(
                        Q(6 * row["components"], 49) / gap
                    )
                    left_cross = (
                        7 * third_floor * second_numerator
                        + (
                            6
                            * row["components"]
                            * K.G.A.RULER
                            * second
                        )
                    ) * base_gap.denominator
                    right_cross = (
                        base_gap.numerator
                        * 49
                        * K.G.A.RULER
                        * second
                        * third_floor
                    )
                    require(
                        (third_floor <= upper) == (left_cross >= right_cross),
                        "finite scalar cap cross-check failed",
                    )
                    if upper < third_floor:
                        continue
                    scalar_class = "finite_z3_cap"
                    scalar_finite += 1

                first_denominator = canonical_l // math.gcd(
                    first,
                    canonical_l,
                )
                second_denominator = canonical_l // math.gcd(
                    second,
                    canonical_l,
                )
                possible_thirds = allowed.get(
                    tuple(
                        sorted(
                            (first_denominator, second_denominator)
                        )
                    ),
                    (),
                )
                if not possible_thirds:
                    continue

                row_prefix = True
                prefix_cap = max(prefix_cap, second)
                if upper is None:
                    prefix_no_gap += 1
                else:
                    prefix_finite += 1
                prefix_record = (
                    body,
                    first,
                    second,
                    scalar_class,
                    upper,
                    first_denominator,
                    second_denominator,
                    possible_thirds,
                )
                prefix_digest.update(repr(prefix_record).encode())
                prefix_digest.update(b"\n")

                chosen: tuple[int, int] | None = None
                if upper is None:
                    chosen = (possible_thirds[0], -1)
                else:
                    for third_denominator in possible_thirds:
                        third = first_label_with_denominator(
                            canonical_l,
                            third_denominator,
                            third_floor,
                            upper,
                        )
                        if third is not None:
                            candidate = (third_denominator, third)
                            if chosen is None or (third, third_denominator) < (
                                chosen[1],
                                chosen[0],
                            ):
                                chosen = candidate
                if chosen is None:
                    continue

                row_retained = True
                retained_cap = max(retained_cap, second)
                if upper is None:
                    retained_no_gap += 1
                    no_gap_cap = max(no_gap_cap, second)
                else:
                    retained_finite += 1
                    finite_cap = max(finite_cap, second)
                retained_record = (*prefix_record, chosen)
                retained_digest.update(repr(retained_record).encode())
                retained_digest.update(b"\n")
        prefix_rows += row_prefix
        retained_rows += row_retained

    require(
        scalar_finite == result["merged_scalar_finite"]
        and scalar_no_gap == result["merged_no_positive_gap"],
        "joined scalar bank disagrees with source profile",
    )
    return (
        body,
        scalar_finite,
        scalar_no_gap,
        prefix_finite,
        prefix_no_gap,
        retained_finite,
        retained_no_gap,
        prefix_cap,
        retained_cap,
        finite_cap,
        no_gap_cap,
        prefix_rows,
        retained_rows,
        prefix_digest.hexdigest(),
        retained_digest.hexdigest(),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(8, mp.cpu_count() or 1))
    parser.add_argument("--residual", type=Path, default=DEFAULT_RESIDUAL)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    require(
        args.residual.resolve() == DEFAULT_RESIDUAL.resolve(),
        "scratch verifier currently pins the default residual path",
    )
    if EXPECTED_SCOUT_SHA256:
        require(
            sha256(SCOUT_SOURCE) == EXPECTED_SCOUT_SHA256,
            "typed-apex scout source changed",
        )

    roots = tuple(combinations(range(1, 15), 6))
    if args.workers == 1:
        profiles = [profile_join(body) for body in roots]
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            profiles = list(pool.imap(profile_join, roots, chunksize=1))
    profiles.sort(key=lambda row: row[0])

    def total(index: int) -> int:
        return sum(row[index] for row in profiles)

    scalar_finite, scalar_no_gap = total(1), total(2)
    prefix_finite, prefix_no_gap = total(3), total(4)
    retained_finite, retained_no_gap = total(5), total(6)
    prefix_cap = max(row[7] for row in profiles)
    retained_cap = max(row[8] for row in profiles)
    finite_cap = max(row[9] for row in profiles)
    no_gap_cap = max(row[10] for row in profiles)
    prefix_rows, retained_rows = total(11), total(12)
    prefix_bodies = sum(row[11] > 0 for row in profiles)
    retained_bodies = sum(row[12] > 0 for row in profiles)
    prefix_digest = hashlib.sha256(
        b"LRC14/k4/lorenz-activity/physical-prefix/global/v1\n"
        + repr(tuple((row[0], row[13]) for row in profiles)).encode()
    ).hexdigest()
    retained_digest = hashlib.sha256(
        b"LRC14/k4/lorenz-activity/physical-realizability/global/v1\n"
        + repr(tuple((row[0], row[14]) for row in profiles)).encode()
    ).hexdigest()

    require(
        scalar_finite == 2_042_669
        and scalar_no_gap == 6_471_505,
        "global scalar pair bank changed",
    )
    lines = [
        "LRC14 k4 physical Lorenz+activity exploratory join",
        f"typed_apex_source_sha256={sha256(SCOUT_SOURCE)}",
        f"residual_sha256={sha256(args.residual)}",
        f"residual_rows={RESIDUAL_CENSUS[0]};"
        f"residual_shapes={RESIDUAL_CENSUS[1]};"
        f"residual_bodies={RESIDUAL_CENSUS[2]};"
        f"residual_body_D_rows={RESIDUAL_CENSUS[3]};"
        f"pair_multiset_keys={RESIDUAL_CENSUS[4]};"
        f"pair_to_third_links={RESIDUAL_CENSUS[5]}",
        f"corrected_residual_sha256={EXPECTED_CORRECTED_RESIDUAL_SHA256};"
        "audited_deleted_rows=2",
        "denominator_semantics=delete_each_possible_third_position_from_the_"
        "sorted_multiset;match_sorted(d(z1),d(z2));ledger_slots_have_no_"
        "physical_ownership",
        f"scalar_finite_z3_pairs={scalar_finite};"
        f"scalar_no_positive_gap_pairs={scalar_no_gap};"
        f"scalar_total={scalar_finite+scalar_no_gap}",
        f"denominator_prefix_finite_pairs={prefix_finite};"
        f"denominator_prefix_no_positive_gap_pairs={prefix_no_gap};"
        f"denominator_prefix_total={prefix_finite+prefix_no_gap};"
        f"prefix_z2_cap={prefix_cap};prefix_bodies={prefix_bodies};"
        f"prefix_(E,z1)_rows={prefix_rows}",
        f"finite_z3_realizable_pairs={retained_finite};"
        f"no_positive_gap_denominator_pairs={retained_no_gap};"
        f"final_necessary_pairs={retained_finite+retained_no_gap};"
        f"retained_z2_cap={retained_cap};finite_z2_cap={finite_cap};"
        f"no_positive_gap_z2_cap={no_gap_cap};"
        f"retained_bodies={retained_bodies};"
        f"retained_(E,z1)_rows={retained_rows}",
        f"finite_prefix_without_physical_z3={prefix_finite-retained_finite}",
        f"prefix_digest={prefix_digest}",
        f"physical_realizability_digest={retained_digest}",
        "scope=exploratory necessary screen against an incoming,"
        "not-yet-canonical Lorenz+activity residual;not a realized-cover census",
        "all_exact_controls=PASS",
    ]
    output = "\n".join(lines) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
