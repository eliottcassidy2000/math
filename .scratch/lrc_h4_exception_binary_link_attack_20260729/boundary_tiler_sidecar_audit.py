#!/usr/bin/env python3
"""Audit the ten surviving H4 flags as two weak-boundary LRC tilers.

The fresh-child-q1/H2 recursion leaves ten unordered H4 flags, all on

    E=(2,4,6,8,10,12,14),  apex=22.

This script independently checks that their parent-valid literal five-covers
are precisely the H4-pair views of two full thirteen-speed rows:

    2*{1,...,13}
    2*{1,...,11,13,24}.

The closed-danger interval model sees both rows as measure covers.  They are
not LRC counterexamples: after common-gcd normalization neither row contains a
multiple of 14, and the denominator-14 clock is a weak lonely witness.  In the
unscaled coordinates the witness time is 1/28.
"""

from __future__ import annotations

import hashlib
import importlib.util
from fractions import Fraction as F
from itertools import combinations
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
H4_PATH = (
    ROOT
    / "04-computation/lrc14_j6_paircap_exception_h4_membership_census_codex_20260729.py"
)
H4_SHA256 = "63a80908a6380a877345f0cc4aba7a5e0ef2bb3d59b1b10d58367444ed406b75"

BODY = (2, 4, 6, 8, 10, 12, 14)
APEX = 22
RAW_TAILS = {
    "AP": (16, 18, 20, 24, 26),
    "GW": (16, 18, 20, 26, 48),
}
NORMALIZED_ROWS = {
    "AP": tuple(range(1, 14)),
    "GW": (*range(1, 12), 13, 24),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_h4():
    require(file_sha256(H4_PATH) == H4_SHA256, "H4 membership engine changed")
    spec = importlib.util.spec_from_file_location("h4_boundary_tiler_base", H4_PATH)
    require(spec is not None and spec.loader is not None, "cannot load H4 engine")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


H4 = load_h4()
E = H4.E


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def circle_distance(value: F) -> F:
    floor_value = value.numerator // value.denominator
    residue = value - floor_value
    return min(residue, 1 - residue)


def main() -> None:
    matches = [
        H4.exact_row(fields)
        for fields in H4.parse_exceptions()
        if tuple(map(int, fields["E"].split(","))) == BODY
    ]
    require(len(matches) == 1, "special exception branch changed")
    row = matches[0]
    require(
        row["rank"] == 1 and row["apex"] == APEX and row["prefix"] == (APEX,),
        "special marked suffix changed",
    )
    carrier, components, mass = E.R.CORE.good_norm((*BODY, APEX))
    require(
        components == row["components"] and mass == row["mass"],
        "special parent carrier changed",
    )
    h4_core = frozenset(row["core"])
    semantic = hashlib.sha256(
        b"LRC14/j6/paircap-exception/H4-boundary-tiler-sidecar/v1\n"
    )
    pair_to_kinds: dict[tuple[int, int], list[str]] = {}
    incidence_lines: list[str] = []
    all_link_margins: list[F] = []

    for kind, tail in RAW_TAILS.items():
        full = tuple(sorted((*BODY, APEX, *tail)))
        require(len(full) == 13, f"{kind} full row is not distinct")
        common_gcd = 0
        for speed in full:
            common_gcd = gcd(common_gcd, speed)
        normalized = tuple(speed // common_gcd for speed in full)
        require(
            common_gcd == 2 and normalized == NORMALIZED_ROWS[kind],
            f"{kind} normalization changed",
        )
        require(
            all(speed % 14 for speed in normalized),
            f"{kind} normalized row contains a multiple of 14",
        )
        weak_time = F(1, 14 * common_gcd)
        distances = tuple(circle_distance(speed * weak_time) for speed in full)
        require(
            min(distances) == F(1, 14),
            f"{kind} denominator-clock boundary changed",
        )
        require(
            not E.R.subtract_local_multi(carrier, tail),
            f"{kind} tail no longer measure-covers the parent carrier",
        )
        require(
            all(label >= E.FIRST_EXTERNAL for label in tail)
            and not (set(tail) & set(row["prefix"])),
            f"{kind} violates the inherited ordered prefix",
        )

        retained = tuple(sorted(set(tail) & h4_core))
        require(
            retained == (
                (16, 18, 20, 26) if kind == "AP" else (16, 18, 20, 26, 48)
            ),
            f"{kind} H4 membership changed",
        )
        for hpair in combinations(retained, 2):
            pair_to_kinds.setdefault(hpair, []).append(kind)
            residual = E.R.subtract_local_multi(carrier, hpair)
            residual_mass = E.interval_mass(residual)
            link_threshold = residual_mass - row["q1"]
            triple = tuple(label for label in tail if label not in hpair)
            require(len(triple) == 3, "complementary triple changed")
            require(
                not E.R.subtract_local_multi(residual, triple),
                f"{kind}/{hpair} triple no longer covers the literal child",
            )
            link_margins = tuple(
                residual_mass
                - E.interval_mass(E.R.subtract_local_multi(residual, edge))
                - link_threshold
                for edge in combinations(triple, 2)
            )
            require(
                all(margin >= 0 for margin in link_margins),
                f"{kind}/{hpair} lost a parent binary-link edge",
            )
            all_link_margins.extend(link_margins)
            line = (
                f"kind={kind};L={hpair};T={triple};"
                f"hR={ftext(residual_mass)};"
                f"theta={ftext(link_threshold)};"
                f"link_margins={tuple(ftext(value) for value in link_margins)};"
                f"parent_cover=1;prefix_valid=1;gcd={common_gcd};"
                f"normalized={normalized};weak_t={ftext(weak_time)};"
                f"weak_min={ftext(min(distances))};exit=THM-369"
            )
            incidence_lines.append(line)
            semantic.update((line + "\n").encode())

    expected_pairs = tuple(combinations((16, 18, 20, 26, 48), 2))
    require(
        tuple(sorted(pair_to_kinds)) == tuple(sorted(expected_pairs)),
        "surviving unordered H4 flag set changed",
    )
    require(
        len(incidence_lines) == 16
        and sum(len(kinds) for kinds in pair_to_kinds.values()) == 16,
        "boundary tiler incidence count changed",
    )
    require(
        sum(len(kinds) == 2 for kinds in pair_to_kinds.values()) == 6
        and sum(len(kinds) == 1 for kinds in pair_to_kinds.values()) == 4,
        "AP/GW pair overlap changed",
    )

    print("LRC14 j6 H4 boundary-tiler sidecar audit")
    print("surviving_flags=" + repr(tuple(sorted(pair_to_kinds))))
    print(
        "pair_kind_partition="
        + repr(
            tuple(
                (pair, tuple(sorted(kinds)))
                for pair, kinds in sorted(pair_to_kinds.items())
            )
        )
    )
    print(
        f"incidences={len(incidence_lines)};"
        f"minimum_link_margin={ftext(min(all_link_margins))};"
        f"zero_link_margins={sum(value == 0 for value in all_link_margins)}"
    )
    for line in incidence_lines:
        print("BOUNDARY;" + line)
    print(f"semantic_digest={semantic.hexdigest()}")
    print(
        "scope=ten H4 flags on one exception branch;"
        "two closed-danger tilers;both weakly lonely after gcd normalization;"
        "the ordered prefix does not kill them"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
