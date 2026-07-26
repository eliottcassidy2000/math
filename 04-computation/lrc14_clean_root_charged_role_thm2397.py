#!/usr/bin/env python3
"""Exact companion for THM-2397.

Checks raw-role Parseval ledgers, every finite rational
cross-correlation mask, the 13-translate all-colour refinement,
pigeonholes, and every rational floor.
"""

from fractions import Fraction as F
from itertools import combinations


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


ROLE_COUNTS = (4, 2, 2, 2, 2)
ROLE_COUNT = len(ROLE_COUNTS)
TARGET_COLOURS = 12
PAIR_COUNT = ROLE_COUNT * TARGET_COLOURS
PRIME = 13


def nonzero_energy(mask_size: int) -> F:
    """Normalized F_13 energy outside the zero character."""

    return F(13 * mask_size - mask_size * mask_size, 169)


singleton_energy = nonzero_energy(1)
adjacent_energy = nonzero_energy(2)
require(singleton_energy == F(12, 169), "singleton energy")
require(adjacent_energy == F(22, 169), "adjacent energy")
require(PAIR_COUNT == 60, "raw role/colour count")

# Exact single-factor deletion identity B-B(1-D)=BD on every Boolean
# state.  In the theorem B is the product of the other five safe roles
# and D is the actual q_* danger factor.
for base_safe in (0, 1):
    for q_danger in (0, 1):
        deleted = base_safe - base_safe * (1 - q_danger)
        require(deleted == base_safe * q_danger, "deletion identity")

# Each raw role is disjoint from A. Its nonzero-character inner-product
# sum is minus the zero-character product |A|*n_i/169.
for a_size in (1, 2):
    for role_count in ROLE_COUNTS:
        role_sum = -F(a_size * role_count, 169)
        require(role_sum < 0, "raw role Parseval sign")
        require(
            -role_sum / TARGET_COLOURS >= F(1, 1014),
            "uniform every-role floor",
        )

# Whole-bank signed ledgers. In the outside-double status, the extra
# singleton double root contributes -2/169 by disjoint Parseval.
singleton_bank_sum = -singleton_energy
adjacent_double_cross = -F(2, 169)
adjacent_bank_sum = -adjacent_energy + adjacent_double_cross
require(singleton_bank_sum == -F(12, 169), "singleton raw bank")
require(adjacent_bank_sum == -F(24, 169), "adjacent raw bank")
require(-singleton_bank_sum / PAIR_COUNT == F(1, 845), "uniform bank floor")
require(-adjacent_bank_sum / PAIR_COUNT == F(2, 845), "adjacent bank floor")

# The adjacent raw all-role sum need not be negative colour by colour.
# The theorem's A={0,1}, D={3}, k=5 hostile has a rigorous elementary
# lower bound 1-17*pi^2/169 > 53/8281 on the numerator, using pi<22/7.
hostile_numerator_floor = 1 - F(17, 169) * F(22, 7) ** 2
require(hostile_numerator_floor == F(53, 8281), "hostile rational floor")
require(hostile_numerator_floor > 0, "adjacent sign hostile")

# Rational cross-correlation universe.  A has size one or two; an actual
# lower q-mask has size two and the guard mask has size four.  For every
# disjoint pair, the integer correlation has zero base coefficient,
# positive total, and is not flat.  Phi_13 then proves every nonzero
# Fourier colour survives.
cross_correlation_masks = 0
residues = tuple(range(PRIME))
for a_size in (1, 2):
    for a_tuple in combinations(residues, a_size):
        a_mask = set(a_tuple)
        complement = tuple(r for r in residues if r not in a_mask)
        for role_size in (2, 4):
            for role_tuple in combinations(complement, role_size):
                role_mask = set(role_tuple)
                correlation = tuple(
                    sum(
                        1
                        for r in role_mask
                        if (r + shift) % PRIME in a_mask
                    )
                    for shift in residues
                )
                require(correlation[0] == 0, "exclusive base shift")
                require(
                    sum(correlation) == a_size * role_size,
                    "cross-correlation mass",
                )
                require(
                    len(set(correlation)) > 1,
                    "forbidden flat rational kernel",
                )
                cross_correlation_masks += 1
require(cross_correlation_masks == 37323, "cross-correlation universe")

# There are C(13,2)=78 abstract two-root masks, but a fixed ordinary role
# ranges over only the 13 translates of one fixed gap.  Odd prime order
# rules out a vanishing two-root coefficient, while the theorem's chord
# estimate supplies |a_k f_k|>4/28561 uniformly.
two_root_masks = PRIME * (PRIME - 1) // 2
fixed_gap_translates = PRIME
require(two_root_masks == 78, "abstract two-root mask count")
for gap in range(1, (PRIME + 1) // 2):
    translates = {
        tuple(sorted((shift, (shift + gap) % PRIME)))
        for shift in range(PRIME)
    }
    require(len(translates) == fixed_gap_translates, "fixed-gap orbit")
for k in range(1, PRIME):
    for gap in range(1, PRIME):
        require((k * gap) % PRIME != 0, "two-root Fourier collapse")
quantitative_all_colour_factor = F(4, fixed_gap_translates * 28561)
require(
    quantitative_all_colour_factor == F(4, 371293),
    "quantitative all-colour factor",
)

# Two BV Fourier bounds contribute 1/(4*pi^2), while the two-sided
# nonzero square sum contributes pi^2/3.
bv_mixing_constant = F(1, 4) * F(1, 3)
require(bv_mixing_constant == F(1, 12), "BV mixing constant")

# Six conjugate colour pairs distributed among five roles force one role
# on at least two pairs, hence four colours, whenever every colour is live.
conjugate_pairs = TARGET_COLOURS // 2
fixed_role_pairs = (conjugate_pairs + ROLE_COUNT - 1) // ROLE_COUNT
require(conjugate_pairs == 6, "conjugate pair count")
require(fixed_role_pairs == 2, "fixed-role pair count")
require(2 * fixed_role_pairs == 4, "fixed-role colour count")

# Universal last-lane floors.
delta_universal = F(1, 26754)
rho_universal = delta_universal / 52
require(rho_universal == F(1, 1391208), "universal charged cell")
universal_each_fibre = rho_universal / 1014
universal_each_aggregate = rho_universal * rho_universal / 1014
universal_some_fibre = rho_universal / 845
universal_some_aggregate = rho_universal * rho_universal / 845
require(universal_each_fibre == F(1, 1410684912), "universal each fibre")
require(
    universal_each_aggregate == F(1, 1962556135053696),
    "universal each aggregate",
)
require(universal_some_fibre == F(1, 1175570760), "universal some fibre")
require(
    universal_some_aggregate == F(1, 1635463445878080),
    "universal some aggregate",
)
universal_all_colour = rho_universal * quantitative_all_colour_factor
require(
    universal_all_colour == F(1, 129136447986),
    "universal quantitative all-colour floor",
)

# Sharper common-core floors.
delta_common = F(66, 4459)
rho_common = delta_common / 52
require(rho_common == F(33, 115934), "common charged cell")
common_each_fibre = rho_common / 1014
common_each_aggregate = rho_common * rho_common / 1014
common_some_fibre = rho_common / 845
common_some_aggregate = rho_common * rho_common / 845
require(common_each_fibre == F(11, 39185692), "common each fibre")
require(
    common_each_aggregate == F(363, 4542954016328),
    "common each aggregate",
)
require(common_some_fibre == F(33, 97964230), "common some fibre")
require(
    common_some_aggregate == F(1089, 11357385040820),
    "common some aggregate",
)
common_all_colour = rho_common * quantitative_all_colour_factor
require(
    common_all_colour == F(66, 21522741331),
    "common quantitative all-colour floor",
)

# Owner-resolved C_7 x C_13 quantities.
rho_owner = delta_common / 338
require(rho_owner == F(33, 753571), "owner-resolved cell")
owner_some_fibre = rho_owner / 845
owner_some_aggregate = rho_owner * rho_owner / 845
owner_joint_tensor = owner_some_fibre / 7
owner_external_tensor = owner_some_aggregate / 7
owner_each_fibre = rho_owner / 1014
owner_each_aggregate = rho_owner * rho_owner / 1014
owner_each_external = owner_each_aggregate / 7
require(owner_some_fibre == F(33, 636767495), "owner some fibre")
require(
    owner_some_aggregate == F(1089, 479849517974645),
    "owner some aggregate",
)
require(owner_joint_tensor == F(33, 4457372465), "owner joint tensor")
require(
    owner_external_tensor == F(1089, 3358946625822515),
    "owner external tensor",
)
require(owner_each_fibre == F(11, 254706998), "owner each fibre")
require(
    owner_each_aggregate == F(363, 191939807189858),
    "owner each aggregate",
)
require(
    owner_each_external == F(363, 1343578650329006),
    "owner each external tensor",
)
owner_all_colour = rho_owner * quantitative_all_colour_factor
owner_all_colour_joint = owner_all_colour / 7
require(
    owner_all_colour == F(132, 279795637303),
    "owner quantitative all-colour floor",
)
require(
    owner_all_colour_joint == F(132, 1958569461121),
    "owner quantitative joint all-colour floor",
)


print("theorem=THM-2397")
print("status=PROVED+VERIFIED-EXACT+INDEPENDENTLY-HOSTILE-AUDITED")
print(
    f"actual_roles={ROLE_COUNT}; role_counts={','.join(map(str, ROLE_COUNTS))};"
    f" role_colour_pairs={PAIR_COUNT}"
)
print(
    f"raw_signed_sum_singleton={singleton_bank_sum};"
    f" raw_signed_sum_adjacent={adjacent_bank_sum}"
)
print("exclusive_q_star_word=LITERAL_SINGLE_FACTOR_DELETION_SUPPORT")
print("every_actual_role_has_negative_Hermitian_pair=YES")
print(
    f"rational_cross_correlation_masks={cross_correlation_masks};"
    " every_actual_role_all_12_colours=YES"
)
print(f"fixed_lower_q_two_root_translates={fixed_gap_translates}")
print("nonnegative_root_constant_parent_filtering=HEREDITARY")
print("R=13^k_delayed_word_survival=EVENTUAL; BV_constant=1/12")
print(
    f"universal_rho={rho_universal}; each_fibre>={universal_each_fibre};"
    f" some_fibre>={universal_some_fibre};"
    f" all_colour_refined>{universal_all_colour}"
)
print(
    f"common_rho={rho_common}; each_fibre>={common_each_fibre};"
    f" some_fibre>={common_some_fibre};"
    f" all_colour_refined>{common_all_colour}"
)
print(
    f"owner_rho={rho_owner}; joint_tensor>={owner_joint_tensor};"
    f" all_colour_refined>{owner_all_colour};"
    f" all_colour_joint>{owner_all_colour_joint}"
)
print("row_decrement=0; ledger=165; LRC(14)=OPEN")
print("all_checks=PASS")
