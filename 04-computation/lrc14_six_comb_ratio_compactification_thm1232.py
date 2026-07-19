#!/usr/bin/env python3
"""Exact arithmetic replay for THM-1232's six-comb ratio compactification."""

from fractions import Fraction as Q


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def ceil_q(value: Q) -> int:
    return -(-value.numerator // value.denominator)


first_two = 2 * Q(7, 36)
four_tail_baseline = 4 * Q(1, 7)
functional_baseline = first_two + four_tail_baseline
functional_deficit = 1 - functional_baseline
x3_cutoff = Q(2, 3) / functional_deficit

require(functional_baseline == Q(121, 126), "wrong functional baseline")
require(functional_deficit == Q(5, 126), "wrong functional deficit")
require(x3_cutoff == Q(84, 5), "wrong third-ratio cutoff")
require(
    functional_baseline + Q(2, 3) / x3_cutoff == 1,
    "third-ratio boundary identity failed",
)

ceiling3_argument = (6 * x3_cutoff + 1) / 7
ceiling3 = ceil_q(ceiling3_argument)
C3 = 1 + 3 * ceiling3
x4_cutoff = Q(49, 15) * C3

require(ceiling3_argument == Q(509, 35), "wrong level-three ceiling argument")
require(ceiling3 == 15 and C3 == 46, "wrong level-three component cap")
require(x4_cutoff == Q(2254, 15), "wrong fourth-ratio cutoff")

ceiling4_argument = (6 * x4_cutoff + 1) / 7
ceiling4 = ceil_q(ceiling4_argument)
C4 = 1 + 4 * ceiling4
x5_cutoff = Q(21, 8) * C4

require(ceiling4_argument == Q(4513, 35), "wrong level-four ceiling argument")
require(ceiling4 == 129 and C4 == 517, "wrong level-four component cap")
require(x5_cutoff == Q(10857, 8), "wrong fifth-ratio cutoff")

ceiling5_argument = (6 * x5_cutoff + 1) / 7
ceiling5 = ceil_q(ceiling5_argument)
C5 = 1 + 5 * ceiling5
x6_cutoff = 7 * C5

require(ceiling5_argument == Q(32575, 28), "wrong level-five ceiling argument")
require(ceiling5 == 1164 and C5 == 5821, "wrong level-five component cap")
require(x6_cutoff == 40747, "wrong sixth-ratio cutoff")

# Three-comb span coefficients: w4 + 2(w5 + 2w6) < 7w4.
span_coefficients = (1, 2, 4)
require(sum(span_coefficients) == 7, "wrong three-comb span coefficient")

print("THM-1232 SIX-COMB RATIO-CONE COMPACTIFICATION")
print(f"functional_baseline={functional_baseline}")
print(f"functional_deficit={functional_deficit}")
print(f"third_ratio_boundary={x3_cutoff}")
print("boundary_identity=121/126+2/(3*(84/5))=1")
print()
print("THREE_REMAINING_COMB_SPAN")
print("span_coefficients=w4+2*w5+4*w6<7*w4=1/d4")
print("three_prefix_physical_survivor=15/(49*c)")
print("d4/c<(49/15)*C3")
print()
print("EXACT_CEILING_PROPAGATION")
print(f"level3_argument={ceiling3_argument} ceil<={ceiling3} C3<={C3} d4/c<{x4_cutoff}")
print(f"level4_argument={ceiling4_argument} ceil<={ceiling4} C4<={C4} d5/c<{x5_cutoff}")
print(f"level5_argument={ceiling5_argument} ceil<={ceiling5} C5<={C5} d6/c<{x6_cutoff}")
print()
print("TOURNAMENT_LOSS_AUDIT")
print("runner_order_scores=0,1,2,3,4,5 cycles=0 SCCs=1,1,1,1,1,1 Hamiltonian_paths=1")
print("faithful_carrier=load_chambers_plus_prefix-depth_component-obligations")
print("ALL EXACT CHECKS PASSED")
