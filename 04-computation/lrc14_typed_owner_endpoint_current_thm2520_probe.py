#!/usr/bin/env python3
"""Exact sparse endpoint-current probe for the canonical typed owner row.

This is a targeted positive control for THM-2520.  It imports the interval
constructor from ``lrc14_169_twist_variance_opus_20260726.py`` but does not
run that script's 169-twist bank.

The row is the canonical THM-2309/2334/2349 typed instance

    (H,q1,...,q5,c1,c2,c3)
      = (1,14,27,40,53,66,13,13^3,2*13^5),

with owner c1, word clock 13^2, and the three THM-2305 word strata.  The
source script explicitly does not assert that this numerical row is a
scalar-cover row.  Accordingly this companion verifies an exact lawful
typed-row signal, not a row exclusion or LRC(14).

For each word it constructs only the actual interval components of

    H_sigma = 1_E1 (1_Qsigma o T^2)

and the deep-summed integer response

    F_sigma = H_sigma * sum_r d(c3*x-r/13).

It then aggregates their endpoint jumps modulo the prime-to-13 part D0 of
the common grid D.  No array of D cells is allocated.
"""

from __future__ import annotations

from bisect import bisect_right
from collections import defaultdict
from fractions import Fraction
import importlib.util
from pathlib import Path
from typing import Iterator


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


ENGINE_PATH = Path(__file__).with_name(
    "lrc14_169_twist_variance_opus_20260726.py"
)
SPEC = importlib.util.spec_from_file_location("typed_owner_interval_engine", ENGINE_PATH)
require(SPEC is not None and SPEC.loader is not None, "cannot load interval engine")
engine = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(engine)


D = engine.NN
T_DEN = engine.T_DEN
RDIL = engine.RDIL
C3 = engine.W[engine.TB]


def split_13(number: int) -> tuple[int, int]:
    exponent = 0
    while number % 13 == 0:
        exponent += 1
        number //= 13
    return exponent, number


K, D0 = split_13(D)
require(D == 13**K * D0 and D0 % 13 != 0, "bad 13-primary split")
require(RDIL == 169 and D == RDIL * T_DEN, "unexpected word grid")
require(D % (182 * C3) == 0, "deep endpoints miss the common grid")


# A deep endpoint has c3*x=(14*r +/- 13)/182.  Multiplication by 169
# sends every such residue to +/-13 mod 182, exactly a phase-zero c3 word
# boundary.  Since every Q_sigma contains the c3 factor, no response
# component can cross a deep-multiplicity endpoint in its interior.
for root in range(13):
    for sign in (-1, 1):
        mapped = RDIL * (14 * root + sign * 13) % 182
        require(mapped in (13, 169), "deep/word endpoint alignment failed")


E_INTERVALS = engine.build_set(engine.PAT_E, engine.ZELL)
WORD_PATTERNS = (
    ("Qa", engine.PAT_QA),
    ("Qb", engine.PAT_QB),
    ("Qab", engine.PAT_QAB),
)


def response_components(pattern: dict[int, str]) -> Iterator[tuple[int, int]]:
    """Yield H_sigma components as integer half-open intervals on C_D."""

    require(pattern[engine.TB] in ("in", "out"), "word omits c3 boundary")
    word = engine.build_set(pattern, engine.ZELL)
    starts = [left for left, _ in word]
    previous_right = 0

    for slow_left, slow_right in E_INTERVALS:
        lifted_left = RDIL * slow_left
        phase_left = lifted_left % T_DEN
        phase_right = phase_left + RDIL * (slow_right - slow_left)
        require(phase_right - phase_left < T_DEN, "slow interval wraps twice")

        index = bisect_right(starts, phase_left) - 1
        offset = 0
        if index < 0:
            index = len(word) - 1
            offset = -T_DEN

        while True:
            word_left, word_right = word[index]
            word_left += offset
            word_right += offset
            if word_left >= phase_right:
                break
            if word_right > phase_left:
                local_left = max(phase_left, word_left)
                local_right = min(phase_right, word_right)
                if local_right > local_left:
                    left = lifted_left - phase_left + local_left
                    right = lifted_left - phase_left + local_right
                    require(0 <= left < right <= D, "component outside C_D")
                    require(left >= previous_right, "components overlap or reorder")
                    previous_right = right
                    yield left, right
            index += 1
            if index == len(word):
                index = 0
                offset += T_DEN


def deep_value_on_cell(cell: int) -> int:
    """Return sum_r d(c3*x-r/13) on the open D-grid cell ``cell``."""

    # Put x=(cell+1/2)/D.  The count is the number of integers in
    # (13*c3*x-13/14, 13*c3*x+13/14), evaluated by two exact floors.
    base = 14 * 13 * C3 * (2 * cell + 1)
    radius = 26 * D
    denominator = 28 * D
    value = (base + radius) // denominator - (base - radius) // denominator
    require(value in (1, 2), "unexpected deep multiplicity")
    return value


def prune(current: dict[int, int]) -> dict[int, int]:
    return {residue: value for residue, value in current.items() if value != 0}


def l1(current: dict[int, int]) -> int:
    return sum(abs(value) for value in current.values())


def lattice_verdict(mean: Fraction) -> str:
    return (
        "CERTIFIES_C_NONZERO"
        if (13**K * mean).denominator != 1
        else "INCONCLUSIVE"
    )


def denominator_signature(mean: Fraction) -> tuple[int, int]:
    return split_13(mean.denominator)


def analyse(pattern: dict[int, str]) -> dict[str, object]:
    pre_jumps: dict[int, int] = defaultdict(int)
    deep_jumps: dict[int, int] = defaultdict(int)
    pre_mass = 0
    deep_mass = 0
    components = 0

    # Consecutive deep endpoints are separated by at least
    # 2*D/(182*c3).  This plus equal endpoint-cell values independently
    # checks the endpoint-alignment argument above on every component.
    minimum_deep_gap = 2 * D // (182 * C3)

    for left, right in response_components(pattern):
        require(right - left < minimum_deep_gap, "component spans two deep jumps")
        deep_left = deep_value_on_cell(left)
        deep_right = deep_value_on_cell(right - 1)
        require(deep_left == deep_right, "component crosses a deep jump")

        length = right - left
        pre_mass += length
        deep_mass += deep_left * length
        pre_jumps[left % D0] += 1
        pre_jumps[right % D0] -= 1
        deep_jumps[left % D0] += deep_left
        deep_jumps[right % D0] -= deep_right
        components += 1

    pre = prune(pre_jumps)
    deep = prune(deep_jumps)
    require(sum(pre.values()) == sum(deep.values()) == 0, "jump current not cyclic")
    return {
        "components": components,
        "pre_mean": Fraction(pre_mass, D),
        "deep_mean": Fraction(deep_mass, D),
        "pre": pre,
        "deep": deep,
    }


RESULTS = {name: analyse(pattern) for name, pattern in WORD_PATTERNS}


EXPECTED_SUMMARIES = {
    "Qa": (
        188056,
        Fraction(21376087, 17907461390),
        Fraction(254882231, 116398499035),
        14,
        376112,
        14,
        694328,
    ),
    "Qb": (
        1933568,
        Fraction(35505957232, 16132831966251),
        Fraction(68054577200, 16132831966251),
        22,
        3867136,
        22,
        7412064,
    ),
    "Qab": (
        173550,
        Fraction(6675, 33787663),
        Fraction(23759, 62748517),
        4,
        347100,
        4,
        665252,
    ),
}


for name, expected in EXPECTED_SUMMARIES.items():
    result = RESULTS[name]
    actual = (
        result["components"],
        result["pre_mean"],
        result["deep_mean"],
        len(result["pre"]),
        l1(result["pre"]),
        len(result["deep"]),
        l1(result["deep"]),
    )
    require(actual == expected, f"{name} summary changed: {actual}")


def residue_table(
    entries: tuple[tuple[Fraction, int, int], ...]
) -> tuple[dict[int, int], dict[int, int]]:
    pre: dict[int, int] = {}
    deep: dict[int, int] = {}
    for residue, pre_value, deep_value in entries:
        scaled = residue * D0
        require(scaled.denominator == 1, "expected residue misses D0")
        index = scaled.numerator
        pre[index] = pre_value
        deep[index] = deep_value
    return pre, deep


QA_EXPECTED = residue_table(
    (
        (Fraction(1, 28), -86808, -159142),
        (Fraction(41, 742), -16, -16),
        (Fraction(1, 14), 7220, 14440),
        (Fraction(41, 560), -16, -32),
        (Fraction(127, 560), 16, 32),
        (Fraction(237, 742), -17, -34),
        (Fraction(13, 28), 93963, 173468),
        (Fraction(15, 28), -93963, -173468),
        (Fraction(505, 742), 17, 34),
        (Fraction(433, 560), -16, -32),
        (Fraction(519, 560), 16, 32),
        (Fraction(13, 14), -7220, -14440),
        (Fraction(701, 742), 16, 16),
        (Fraction(27, 28), 86808, 159142),
    )
)


QAB_EXPECTED = residue_table(
    (
        (Fraction(1, 28), 86808, 159142),
        (Fraction(13, 28), -86742, -173484),
        (Fraction(15, 28), 86742, 173484),
        (Fraction(27, 28), -86808, -159142),
    )
)


require(
    RESULTS["Qa"]["pre"] == QA_EXPECTED[0]
    and RESULTS["Qa"]["deep"] == QA_EXPECTED[1],
    "Qa sparse table changed",
)
require(
    RESULTS["Qab"]["pre"] == QAB_EXPECTED[0]
    and RESULTS["Qab"]["deep"] == QAB_EXPECTED[1],
    "Qab sparse table changed",
)


def aggregate_cell_jumps(values: tuple[int, ...], modulus: int) -> dict[int, int]:
    current: dict[int, int] = defaultdict(int)
    for index, value in enumerate(values):
        current[index % modulus] += value - values[index - 1]
    return prune(current)


# Sharp controls for the jump/multiplicity screen.
HOSTILE_VALUES = (1,) + (0,) * 12
HOSTILE_C = aggregate_cell_jumps(HOSTILE_VALUES, 1)
HOSTILE_MEAN = Fraction(sum(HOSTILE_VALUES), len(HOSTILE_VALUES))
require(HOSTILE_C == {} and HOSTILE_MEAN == Fraction(1, 13), "hostile failed")

POSITIVE_VALUES = (1, 0)
POSITIVE_C = aggregate_cell_jumps(POSITIVE_VALUES, 2)
POSITIVE_MEAN = Fraction(sum(POSITIVE_VALUES), len(POSITIVE_VALUES))
require(POSITIVE_C == {0: 1, 1: -1}, "positive control failed")
require(POSITIVE_MEAN == Fraction(1, 2), "positive mean failed")


def print_sparse(name: str) -> None:
    result = RESULTS[name]
    print(f"{name} sparse endpoint current: residue=r/D0 pre_jump deep_jump")
    for residue in sorted(set(result["pre"]) | set(result["deep"])):
        print(
            f"  {Fraction(residue, D0)} "
            f"{result['pre'].get(residue, 0):+d} "
            f"{result['deep'].get(residue, 0):+d}"
        )


def main() -> None:
    print("THM-2520 canonical typed-owner endpoint-current probe")
    print(f"engine={ENGINE_PATH.name}")
    print(f"row=(H,q1..q5,c1,c2,c3)={engine.W}")
    print(f"owner=c1={engine.W[engine.OWNER]}; word_clock={RDIL}; c3={C3}")
    print(f"T_DEN={T_DEN}")
    print(f"D={D}=13^{K}*{D0}")
    print("deep_endpoint_alignment: 26/26 residues map to the c3 word boundary")
    print("")

    for name, _ in WORD_PATTERNS:
        result = RESULTS[name]
        pre_v13, pre_non13 = denominator_signature(result["pre_mean"])
        deep_v13, deep_non13 = denominator_signature(result["deep_mean"])
        print(f"{name}: components={result['components']}")
        print(
            "  pre-deep:"
            f" mean={result['pre_mean']}"
            f" denominator=13^{pre_v13}*{pre_non13}"
            f" lattice_screen={lattice_verdict(result['pre_mean'])}"
            f" C_support={len(result['pre'])} C_L1={l1(result['pre'])}"
        )
        print(
            "  deep-summed:"
            f" mean={result['deep_mean']}"
            f" denominator=13^{deep_v13}*{deep_non13}"
            f" lattice_screen={lattice_verdict(result['deep_mean'])}"
            f" C_support={len(result['deep'])} C_L1={l1(result['deep'])}"
        )
    print("")

    print_sparse("Qa")
    print("")
    print_sparse("Qab")
    print("")

    print(
        "pure-13 hostile:"
        f" mean={HOSTILE_MEAN} C_support={len(HOSTILE_C)}"
        " P_13F=1/13 lattice_screen=INCONCLUSIVE"
    )
    print(
        "prime-to-13 positive:"
        f" mean={POSITIVE_MEAN} C={sorted(POSITIVE_C.items())}"
        " lattice_screen=CERTIFIES_C_NONZERO"
    )
    print(
        "scope=TYPED_CANONICAL_ROW_POSITIVE_CONTROL; "
        "SCALAR_COVER_IDENTIFICATION_OPEN; ROW_EXCLUSION=0; LRC14_OPEN"
    )
    print(
        "VERIFIED: all three canonical word strata have nonzero pre-deep and "
        "deep-summed prime-to-13 endpoint current; the denominator screen "
        "settles five of six cases and the Qab deep four-support current is "
        "the sharp endpoint-only survivor."
    )


if __name__ == "__main__":
    main()
