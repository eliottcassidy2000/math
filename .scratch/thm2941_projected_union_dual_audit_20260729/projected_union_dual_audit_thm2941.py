#!/usr/bin/env python3
"""Independent direct-UNION audit of the THM-2941 projected closure.

The canonical five-aligned verifier computes a lower bound for projected-safe
mass by intersecting the two-drift danger unions over body-safe cells and then
taking a complement.  This audit rebuilds every one of its 236,985 terminal
pairs, but independently constructs

    P_prefix = union_j ([0,1] minus (E_z1(j) union E_z2(j))).

For every prefix it also constructs the common-danger intersection and checks
the exact finite De Morgan mass identity.  It then requires the direct-union
mass and stopping cell to agree with the canonical complement engine.

All arithmetic is rational.  Endpoint openness is immaterial to the measured
sets checked here; the theorem obtains its strict contradiction separately
from compact containment in a proper open aligned cover.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import math
import multiprocessing as mp
import os
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def repo_root(start: Path) -> Path:
    for parent in (start, *start.parents):
        if (parent / "04-computation").is_dir() and (parent / "01-canon").is_dir():
            return parent
    raise RuntimeError("cannot locate repository root")


def normalized_sha256(path: Path) -> str:
    payload = path.read_bytes()
    require(
        b"\r" not in payload.replace(b"\r\n", b""),
        f"{path.name}: lone carriage return",
    )
    return hashlib.sha256(payload.replace(b"\r\n", b"\n")).hexdigest()


HERE = Path(__file__).resolve().parent
ROOT = repo_root(HERE)
CANONICAL = (
    ROOT
    / "04-computation"
    / "lrc14_j7_five_aligned_two_drift_projected_closure_thm2941.py"
)
EXPECTED_CANONICAL_SHA256 = (
    "76f891edfcc029a08202481304a809e03e8bd81f247afaeabab685825c4d3662"
)
DEFAULT_OUTPUT = HERE / "projected_union_dual_audit_thm2941.out"

require(
    normalized_sha256(CANONICAL) == EXPECTED_CANONICAL_SHA256,
    "canonical five-aligned verifier changed",
)
spec = importlib.util.spec_from_file_location("thm2941_five_aligned", CANONICAL)
require(spec is not None and spec.loader is not None, "cannot load canonical verifier")
C = importlib.util.module_from_spec(spec)
spec.loader.exec_module(C)


Interval = tuple[F, F]
Intervals = tuple[Interval, ...]


def merge(rows: list[Interval] | tuple[Interval, ...]) -> Intervals:
    """Independent normalized finite-union implementation."""

    ordered = sorted((left, right) for left, right in rows if left < right)
    if not ordered:
        return ()
    out: list[Interval] = []
    left, right = ordered[0]
    for next_left, next_right in ordered[1:]:
        if next_left <= right:
            right = max(right, next_right)
        else:
            out.append((left, right))
            left, right = next_left, next_right
    out.append((left, right))
    return tuple(out)


def intersection(first: Intervals, second: Intervals) -> Intervals:
    """Independent two-pointer intersection of normalized finite unions."""

    out: list[Interval] = []
    i = 0
    j = 0
    while i < len(first) and j < len(second):
        left = max(first[i][0], second[j][0])
        right = min(first[i][1], second[j][1])
        if left < right:
            out.append((left, right))
        if first[i][1] <= second[j][1]:
            i += 1
        else:
            j += 1
    return tuple(out)


def complement(rows: Intervals) -> Intervals:
    """Complement a normalized finite union inside the unit interval."""

    out: list[Interval] = []
    cursor = F(0)
    for left, right in rows:
        require(F(0) <= left < right <= F(1), "interval left the unit line")
        if cursor < left:
            out.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < 1:
        out.append((cursor, F(1)))
    return tuple(out)


def mass(rows: Intervals) -> F:
    return sum((right - left for left, right in rows), F(0))


def floor_fraction(value: F) -> int:
    return value.numerator // value.denominator


def ceil_fraction(value: F) -> int:
    return -((-value.numerator) // value.denominator)


def cell_danger(cell: int, speed: int, canonical_l: int) -> Intervals:
    """Reconstruct E_speed(cell) from its tooth-center inequalities.

    The integer range is derived from

        m + 1/14 > speed*cell/L,
        m - 1/14 < speed*(cell+1)/L,

    rather than copied from the canonical phase routine.
    """

    x0 = F(speed * cell, canonical_l)
    x1 = F(speed * (cell + 1), canonical_l)
    radius = F(1, 14)
    first_tooth = floor_fraction(x0 - radius) + 1
    last_tooth = ceil_fraction(x1 + radius) - 1
    rows: list[Interval] = []
    for tooth in range(first_tooth, last_tooth + 1):
        left = max(
            F(0),
            F(canonical_l, speed) * (F(tooth) - radius) - cell,
        )
        right = min(
            F(1),
            F(canonical_l, speed) * (F(tooth) + radius) - cell,
        )
        if left < right:
            rows.append((left, right))
    return merge(rows)


def direct_projected_prefix(
    cells: tuple[int, ...],
    canonical_l: int,
    first: int,
    second: int,
) -> tuple[F, int]:
    """Return the direct safe-union certificate and its first stopping cell."""

    projected_safe: Intervals = ()
    common_danger: Intervals = ((F(0), F(1)),)
    for cells_used, cell in enumerate(cells, start=1):
        local_danger = merge(
            (
                *cell_danger(cell, first, canonical_l),
                *cell_danger(cell, second, canonical_l),
            )
        )
        local_safe = complement(local_danger)
        projected_safe = merge((*projected_safe, *local_safe))
        common_danger = intersection(common_danger, local_danger)

        direct_mass = mass(projected_safe)
        dual_mass = F(1) - mass(common_danger)
        require(
            direct_mass == dual_mass,
            (
                f"finite De Morgan failure: L={canonical_l},"
                f" z1={first}, z2={second}, prefix={cells_used}"
            ),
        )
        if direct_mass >= C.ALIGNED_UNION_CAP:
            return direct_mass, cells_used
    return mass(projected_safe), len(cells)


def update_minimum(
    current: tuple[F, tuple[object, ...]] | None,
    candidate: tuple[F, tuple[object, ...]],
) -> tuple[F, tuple[object, ...]]:
    if current is None or candidate < current:
        return candidate
    return current


def audit_body(body: tuple[int, ...]) -> dict[str, object]:
    primary = C.profile(body)
    carrier = C.A.carrier_for(body)
    canonical_l = 14 * math.lcm(*body)
    cells = C.body_cells(carrier, canonical_l)

    high_pairs = sorted(primary["high_pairs"])
    sub_pairs = sorted(primary["sub_pairs"])
    digest = hashlib.sha256()
    digest.update(b"LRC14/THM2941/direct-projected-union/body/v1\n")
    digest.update(repr(body).encode())
    minimum: tuple[F, tuple[object, ...]] | None = None
    maximum_prefix = 0

    for bank_name, pairs in (("high", high_pairs), ("subcritical", sub_pairs)):
        for pair in pairs:
            require(pair[0] == body, "pair escaped its body profile")
            first = pair[1]
            second = pair[2]
            require(pair[3] == canonical_l, "pair canonical period changed")
            direct_mass, direct_prefix = direct_projected_prefix(
                cells,
                canonical_l,
                first,
                second,
            )
            canonical_mass, canonical_prefix = C.projected_safe_lower_bound(
                cells,
                canonical_l,
                first,
                second,
            )
            require(
                (direct_mass, direct_prefix) == (canonical_mass, canonical_prefix),
                (
                    f"direct/complement mismatch: {body},"
                    f" z1={first}, z2={second}"
                ),
            )
            require(
                direct_mass >= C.ALIGNED_UNION_CAP,
                f"direct-union survivor: {body}, z1={first}, z2={second}",
            )
            certificate = (
                bank_name,
                body,
                first,
                second,
                canonical_l,
                direct_mass,
                direct_prefix,
            )
            digest.update(repr(certificate).encode())
            digest.update(b"\n")
            maximum_prefix = max(maximum_prefix, direct_prefix)
            minimum = update_minimum(
                minimum,
                (direct_mass - C.ALIGNED_UNION_CAP, certificate),
            )

    return {
        "body": body,
        "high_count": len(high_pairs),
        "sub_count": len(sub_pairs),
        "maximum_prefix": maximum_prefix,
        "minimum": minimum,
        "body_digest": digest.hexdigest(),
    }


def render(profiles: list[dict[str, object]]) -> str:
    profiles.sort(key=lambda row: row["body"])
    require(len(profiles) == 3_003, "six-body universe changed")
    high_count = sum(row["high_count"] for row in profiles)
    sub_count = sum(row["sub_count"] for row in profiles)
    maximum_prefix = max(row["maximum_prefix"] for row in profiles)
    minimum = min(
        (row["minimum"] for row in profiles if row["minimum"] is not None),
        default=None,
    )
    global_digest = hashlib.sha256(
        b"LRC14/THM2941/direct-projected-union/global/v1\n"
        + repr(
            tuple(
                (
                    row["body"],
                    row["high_count"],
                    row["sub_count"],
                    row["body_digest"],
                )
                for row in profiles
            )
        ).encode()
    ).hexdigest()

    require(high_count == 42_912, "high-excess pair count changed")
    require(sub_count == 194_073, "subcritical pair count changed")
    require(high_count + sub_count == 236_985, "combined pair count changed")
    require(maximum_prefix == 871, "maximum direct-union prefix changed")
    require(
        minimum
        == (
            F(1, 378_105),
            (
                "high",
                (1, 2, 3, 5, 9, 10),
                24,
                277,
                1_260,
                F(180, 277),
                1,
            ),
        ),
        "minimum direct-union margin changed",
    )
    require(
        global_digest
        == "e03500f6687b835ff1bf925b6731236caa8d3e2409d44a8cd62b353ede5d1dcf",
        "direct-union digest changed",
    )

    lines = [
        "THM-2941 direct projected-safe UNION dual audit",
        f"canonical_sha256={normalized_sha256(CANONICAL)}",
        "universe=(six_body_roots=3003,"
        "canonical_terminal_pairs=236985,"
        "all_prefix_arithmetic=exact_rational)",
        f"high_excess_pairs={high_count}",
        f"subcritical_pairs={sub_count}",
        f"combined_pairs={high_count + sub_count}",
        f"maximum_cell_prefix={maximum_prefix}",
        f"minimum_certificate={minimum}",
        f"direct_union_digest={global_digest}",
        "direct_union_equals_complement_engine=PASS",
        "finite_de_morgan_identity_every_prefix=PASS",
        "all_direct_union_certificates_close=PASS",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(8, mp.cpu_count() or 1))
    parser.add_argument("--hash-seed", type=int, default=0)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    require(args.hash_seed >= 0, "hash seed must be nonnegative")
    os.environ["PYTHONHASHSEED"] = str(args.hash_seed)

    roots = tuple(combinations(range(1, 15), 6))
    if args.workers == 1:
        profiles = [audit_body(body) for body in roots]
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            profiles = list(pool.imap(audit_body, roots, chunksize=1))
    output = render(profiles)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
