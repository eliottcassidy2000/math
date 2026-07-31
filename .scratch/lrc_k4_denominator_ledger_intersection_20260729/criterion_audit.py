#!/usr/bin/env python3
"""Exact on-demand intersection criterion for the physical k=4 bank.

THM-2928 stores denominator triples in sorted denominator order.  A physical
THM-2941 row is instead ordered by speed z1 < z2 < z3.  This module keeps
those coordinates distinct:

    delta_i = L/gcd(z_i,L)       (physical owner order),
    e = sorted(delta_1,delta_2,delta_3)  (ledger multiset order).

The final all-top relaxation is symmetric in the three denominators, so
sorting is lossless there.  It is not lawful to transfer the sorted positions
back to z1,z2,z3 for scalar-excess or residual recursion.
"""

from __future__ import annotations

import hashlib
import importlib.util
import math
from fractions import Fraction as F
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
THREE_SOURCE = (
    ROOT / "04-computation/lrc14_three_drift_body_projection_fiber_thm2928.py"
)
EXPECTED_THREE_SHA256 = (
    "42dc165781148c702dfcd3c6535f4d02aee516af60b5ddf602a19cb1d87695e4"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


require(
    hashlib.sha256(THREE_SOURCE.read_bytes()).hexdigest()
    == EXPECTED_THREE_SHA256,
    "THM-2928 three-drift source changed",
)
spec = importlib.util.spec_from_file_location("three_drift_referee", THREE_SOURCE)
require(spec is not None and spec.loader is not None, "cannot load referee")
T = importlib.util.module_from_spec(spec)
spec.loader.exec_module(T)
S = T.support_module


def physical_denominator(speed: int, canonical_l: int) -> int:
    return canonical_l // math.gcd(speed, canonical_l)


def ledger_accepts(
    body: tuple[int, ...],
    physical_denominators: tuple[int, int, int],
) -> bool:
    """Final THM-2928 necessary ledger predicate, including diagonal closure."""

    canonical_l, ranges = S.safe_cell_ranges(body)
    require(len(physical_denominators) == 3, "need three physical denominators")
    require(
        all(
            denominator > 1 and canonical_l % denominator == 0
            for denominator in physical_denominators
        ),
        "physical denominator left the genuine drift divisor universe",
    )
    resolving_d = math.lcm(*physical_denominators)
    require(canonical_l % resolving_d == 0, "resolving denominator left L")

    arcs = T.projected_support_arcs(resolving_d, ranges)
    support_count = sum(right - left for left, right in arcs)
    if F(support_count, resolving_d) > T.SUPPORT_CUTOFF:
        return False

    # The full diagonal sector is independently empty, not merely filtered.
    if len(set(physical_denominators)) == 1:
        return False

    top_capacity = 0
    for denominator in physical_denominators:
        width = (denominator + 6) // 7
        histogram = T.residue_load_histogram(arcs, denominator)
        top_capacity += T.top_class_load(histogram, width)
    return support_count <= top_capacity


def denominator_has_ordered_label(
    canonical_l: int,
    denominator: int,
    lower_speed: int,
    upper_speed: int | None,
) -> bool:
    """Whether an exact-denominator label occurs in the physical z3 interval."""

    scale = canonical_l // denominator
    first_c = lower_speed // scale + 1
    last_c = None if upper_speed is None else upper_speed // scale
    candidate = first_c
    while last_c is None or candidate <= last_c:
        if math.gcd(candidate, denominator) == 1:
            return True
        candidate += 1
    return False


def allowed_third_denominators(
    body: tuple[int, ...],
    first_speed: int,
    second_speed: int,
    third_lower: int,
    third_upper: int | None,
) -> tuple[int, ...]:
    """Project a physical ordered prefix into the symmetric denominator ledger."""

    canonical_l, _ranges = S.safe_cell_ranges(body)
    require(first_speed < second_speed < third_lower, "physical order changed")
    first_d = physical_denominator(first_speed, canonical_l)
    second_d = physical_denominator(second_speed, canonical_l)
    require(first_d > 1 and second_d > 1, "an alleged drift is aligned")
    result = []
    for third_d in S.divisors(canonical_l):
        if third_d <= 1:
            continue
        if not ledger_accepts(body, (first_d, second_d, third_d)):
            continue
        if not denominator_has_ordered_label(
            canonical_l,
            third_d,
            third_lower - 1,
            third_upper,
        ):
            continue
        result.append(third_d)
    return tuple(result)


def main() -> None:
    # Physical z-order and denominator order are independent.
    body = (1, 2, 3, 4, 5, 6)
    canonical_l, _ranges = S.safe_cell_ranges(body)
    require(canonical_l == 840, "ordering control ruler changed")
    require(17 < 20, "speed-order control changed")
    require(
        physical_denominator(17, canonical_l)
        == 840
        > physical_denominator(20, canonical_l)
        == 42,
        "denominator-order hostile changed",
    )

    # Repeated denominators retain multiplicity, but the all-equal case is
    # killed by the proved diagonal theorem.
    require(
        not ledger_accepts(body, (840, 840, 840)),
        "diagonal denominator triple revived",
    )

    # Every exact denominator class has arbitrarily large physical labels:
    # c can be increased within a reduced residue class modulo d.
    for denominator in S.divisors(canonical_l):
        if denominator > 1:
            require(
                denominator_has_ordered_label(
                    canonical_l,
                    denominator,
                    10_000,
                    None,
                ),
                f"unbounded denominator class {denominator} disappeared",
            )

    print("LRC14 k4 physical/denominator-ledger criterion audit")
    print(f"three_drift_source_sha256={EXPECTED_THREE_SHA256}")
    print(
        "physical_key=(E,z1,z2,z3);"
        "delta_i=L/gcd(z_i,L);"
        "ledger_key=(E,lcm(delta_i),sorted(delta_i))"
    )
    print("physical_z_order_implies_denominator_order=FALSE;(17,20)->(840,42)")
    print("denominator_multiplicity=RETAINED;diagonal_triples=KILLED")
    print("finite_z3_branch=requires_exact_coprime_label_in_interval")
    print("unbounded_z3_branch=every_d3>1_has_ordered_labels")
    print("criterion_controls=PASS")


if __name__ == "__main__":
    main()
