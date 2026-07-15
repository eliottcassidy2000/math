#!/usr/bin/env python3
"""Exact replay for THM-825's positional-moment layer codec.

For a k-state word w on ordered positions 0,...,s-1, the degree-d carrier is
the per-state vector

    (sum_{i:w_i=c} i^j : c=0,...,k-1, j=0,...,d).

The script exhausts binary and four-state words through the first quadratic
failure, verifies the staircase B2 layer-size formulas, and records the
literal length-seven witness.  The proof in THM-825 is elementary and does
not depend on exhaustion; this is its finite exact replay.

Tournament Analysis treats M0, M0+M1, M0+M1+M2, and the literal word as
information carriers.  The pairwise observable is how many unordered word
pairs the carrier separates.  The gauges are retention and retention per
integer field.  These are diagnostics about carrier economy, not tournaments
whose vertices are runners.

Preservation boundary: the carrier reconstructs a literal mirror-layer word
through the proved size bound.  It does not preserve tournament node labels,
relations between different layers, LRC gaps, owner clocks, or continuation.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from itertools import product
from pathlib import Path

from tournament_tiling_metagraph_address_codex_S4 import carrier_tournament


OUTPUT = Path("05-knowledge/results/quadratic_positional_layer_codec_codex_S13_postjoin.out")
JSON_OUTPUT = Path("05-knowledge/results/quadratic_positional_layer_codec_codex_S13_postjoin.json")


def moment_key(word: tuple[int, ...], states: int, degree: int) -> tuple[int, ...]:
    return tuple(
        sum(position**power for position, state in enumerate(word) if state == colour)
        for colour in range(states)
        for power in range(degree + 1)
    )


def partition_census(states: int, size: int, degree: int) -> dict[str, int]:
    fibres = Counter(
        moment_key(word, states, degree)
        for word in product(range(states), repeat=size)
    )
    total = states**size
    return {
        "words": total,
        "cells": len(fibres),
        "collision_cells": sum(multiplicity > 1 for multiplicity in fibres.values()),
        "collision_excess": total - len(fibres),
        "collision_pairs": sum(
            multiplicity * (multiplicity - 1) // 2 for multiplicity in fibres.values()
        ),
        "max_fibre": max(fibres.values()),
    }


def staircase_layer_sizes(n: int) -> dict[str, object]:
    # A low-half position at clock tau is a tile (a,b) with a+b-1=tau,
    # a>=b+2.  Reflection supplies its high-half partner.  At tau=n the
    # positions are fixed; (n,1) is the canonical apex.
    nonfixed = [
        sum(
            a + b - 1 == tau
            for b in range(1, n - 1)
            for a in range(n, b + 1, -1)
        )
        for tau in range(3, n)
    ]
    # Only the low half tau<n appears above, so this direct count already has
    # one representative per reflection pair.
    fixed = sum(
        a + b - 1 == n
        for b in range(1, n - 1)
        for a in range(n, b + 1, -1)
    )
    expected_nonfixed = [(tau - 1) // 2 for tau in range(3, n)]
    expected_fixed = (n - 1) // 2
    assert nonfixed == expected_nonfixed
    assert fixed == expected_fixed
    return {
        "n": n,
        "nonfixed": nonfixed,
        "max_nonfixed": max(nonfixed, default=0),
        "max_nonfixed_formula": (n - 2) // 2,
        "fixed": fixed,
        "fixed_free_after_apex_zero": fixed - 1,
        "fixed_free_formula": (n - 1) // 2 - 1,
    }


def audit() -> dict:
    census = {
        str(states): {
            str(size): {
                str(degree): partition_census(states, size, degree)
                for degree in range(3)
            }
            for size in range(1, 8)
        }
        for states in (2, 4)
    }

    for states in (2, 4):
        for size in range(1, 7):
            assert census[str(states)][str(size)]["2"]["collision_excess"] == 0
        assert census[str(states)]["7"]["2"]["collision_excess"] > 0

    left = (0, 4, 5)
    right = (1, 2, 6)
    witness_moments = (
        len(left),
        sum(left),
        sum(position * position for position in left),
    )
    assert witness_moments == (
        len(right),
        sum(right),
        sum(position * position for position in right),
    ) == (3, 9, 41)
    binary_left = tuple(int(position in left) for position in range(7))
    binary_right = tuple(int(position in right) for position in range(7))
    assert moment_key(binary_left, 2, 2) == moment_key(binary_right, 2, 2)
    four_left = tuple(1 if position in left else 0 for position in range(7))
    four_right = tuple(1 if position in right else 0 for position in range(7))
    assert moment_key(four_left, 4, 2) == moment_key(four_right, 4, 2)

    geometry = [staircase_layer_sizes(n) for n in range(5, 18)]
    for row in geometry:
        assert row["max_nonfixed"] == row["max_nonfixed_formula"]
        assert row["fixed_free_after_apex_zero"] == row["fixed_free_formula"]
    assert all(
        max(row["max_nonfixed"], row["fixed_free_after_apex_zero"]) <= 6
        for row in geometry if row["n"] <= 15
    )
    assert next(
        row["n"] for row in geometry
        if max(row["max_nonfixed"], row["fixed_free_after_apex_zero"]) >= 7
    ) == 16

    # Four-state size-seven carrier comparison: this is the first failing
    # layer and therefore the informative economy boundary.
    s7 = census["4"]["7"]
    carrier_stats = {
        "M0": {"separated_pairs": 0, "cells": 4},
        "M0_M1": {"separated_pairs": s7["0"]["collision_pairs"] - s7["1"]["collision_pairs"], "cells": 8},
        "M0_M1_M2": {"separated_pairs": s7["0"]["collision_pairs"] - s7["2"]["collision_pairs"], "cells": 12},
        "literal_word": {"separated_pairs": s7["0"]["collision_pairs"], "cells": 7},
    }
    retention = carrier_tournament(carrier_stats, "retention")
    economy = carrier_tournament(carrier_stats, "economy")
    flips = sum(
        retention["adjacency"][i][j] != economy["adjacency"][i][j]
        for i in range(len(carrier_stats))
        for j in range(i + 1, len(carrier_stats))
    )

    return {
        "schema_version": 1,
        "theorem": "THM-825",
        "carrier_definition": "per-state power sums M_j(c), including M_0 counts",
        "census": census,
        "quadratic_exact_max_layer_size": 6,
        "staircase_exact_through_n": 15,
        "first_staircase_failure_size": 16,
        "length7_witness": {
            "left": list(left),
            "right": list(right),
            "common_M0_M1_M2": list(witness_moments),
            "binary_full_key_equal": True,
            "four_state_full_key_equal": True,
        },
        "geometry": geometry,
        "tournament_analysis": {
            "vertices": list(carrier_stats),
            "pairwise_observable": "number of unordered four-state length-seven word pairs separated relative to M0",
            "switches": ["separation retention", "separation retention per integer field"],
            "tie_hamiltonian_path": list(carrier_stats),
            "carrier_stats": carrier_stats,
            "retention": retention,
            "economy": economy,
            "edge_flips": flips,
        },
    }


def render(result: dict) -> str:
    lines = [
        "THM-825 QUADRATIC POSITIONAL-MOMENT LAYER CODEC",
        "=" * 72,
        "carrier: per-state M0=count, M1=sum(position), M2=sum(position^2)",
        "",
        "EXACT WORD CENSUS (cells/excess/max-fibre for degrees 0,1,2)",
    ]
    for states in (2, 4):
        lines.append(f"states={states}")
        for size in range(1, 8):
            entries = []
            for degree in range(3):
                row = result["census"][str(states)][str(size)][str(degree)]
                entries.append(
                    f"d{degree}:{row['cells']}/{row['collision_excess']}/{row['max_fibre']}"
                )
            lines.append(f"  size={size}: " + "  ".join(entries))
    witness = result["length7_witness"]
    lines.extend(
        [
            "",
            "SHARP BOUNDARY",
            f"  exact for every binary/four-state word of length <= {result['quadratic_exact_max_layer_size']}",
            f"  first witness left/right={witness['left']}/{witness['right']}",
            f"  common (M0,M1,M2)={witness['common_M0_M1_M2']}",
            "",
            "STAIRCASE B2 GEOMETRY",
        ]
    )
    for row in result["geometry"]:
        lines.append(
            f"  n={row['n']:2d}: max_nonfixed={row['max_nonfixed']} "
            f"fixed/free={row['fixed']}/{row['fixed_free_after_apex_zero']}"
        )
    lines.extend(
        [
            f"  S2+M1+M2 literal-exact through n={result['staircase_exact_through_n']}",
            f"  first layer-level failure at n={result['first_staircase_failure_size']}",
            "",
            "TOURNAMENT ANALYSIS (information carriers as planning vertices)",
            f"  vertices={tuple(result['tournament_analysis']['vertices'])}",
            f"  observable={result['tournament_analysis']['pairwise_observable']}",
            f"  switches={tuple(result['tournament_analysis']['switches'])}",
            f"  edge flips={result['tournament_analysis']['edge_flips']}",
            "",
            "Verdict: a quadratic positional sidecar is an unconditional literal",
            "mirror-layer repair throughout every staircase size relevant to LRC(14).",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--json", type=Path, default=JSON_OUTPUT)
    args = parser.parse_args()
    result = audit()
    text = render(result)
    print(text, end="")
    args.output.write_text(text)
    args.json.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
