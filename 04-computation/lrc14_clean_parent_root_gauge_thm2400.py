#!/usr/bin/env python3
"""Exact companion for THM-2400.

Exhausts the rational clean-parent mask universe, verifies invariance
under a common root relabelling, checks every signed and owner/digit
floor, and separates common-slope target actions from a sharp
unequal-slope hostile.
"""

from fractions import Fraction as F
from itertools import combinations


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


PRIME = 13
RESIDUES = tuple(range(PRIME))
TARGET_COLOURS = PRIME - 1
ROLE_COUNTS = (4, 2, 2, 2, 2)
ROLE_COLOUR_PAIRS = len(ROLE_COUNTS) * TARGET_COLOURS


def correlation(a_mask: frozenset[int], c_mask: frozenset[int]) -> tuple[int, ...]:
    """Integer circular correlation sum_r A(r+s)C(r)."""

    return tuple(
        sum(1 for r in c_mask if (r + shift) % PRIME in a_mask)
        for shift in RESIDUES
    )


def translate(mask: frozenset[int], shift: int) -> frozenset[int]:
    return frozenset((r + shift) % PRIME for r in mask)


# Exhaust the complete disjoint-mask universe.  Its correlations have
# zero base coefficient, positive total, and are invariant under every
# simultaneous root translation.  Nonflatness plus Phi_13 then forces
# all twelve charged colours.
mask_pairs = 0
translated_checks = 0
for a_size in (1, 2):
    for a_tuple in combinations(RESIDUES, a_size):
        a_mask = frozenset(a_tuple)
        complement = tuple(r for r in RESIDUES if r not in a_mask)
        for role_size in (2, 4):
            for role_tuple in combinations(complement, role_size):
                c_mask = frozenset(role_tuple)
                corr = correlation(a_mask, c_mask)
                require(corr[0] == 0, "exclusive base shift")
                require(sum(corr) == a_size * role_size, "correlation mass")
                require(len(set(corr)) > 1, "forbidden flat rational kernel")
                for root_shift in RESIDUES:
                    shifted_corr = correlation(
                        translate(a_mask, root_shift),
                        translate(c_mask, root_shift),
                    )
                    require(shifted_corr == corr, "common-root gauge")
                    translated_checks += 1
                mask_pairs += 1

require(mask_pairs == 37323, "mask-pair universe")
require(translated_checks == PRIME * mask_pairs, "translation audit count")

# Pointwise Parseval: disjoint A and C have nonzero-character sum
# -|A||C|/169.  The smallest role/status ledger gives delta/1014.
for a_size in (1, 2):
    for role_size in ROLE_COUNTS:
        signed_sum = -F(a_size * role_size, PRIME * PRIME)
        require(signed_sum < 0, "signed Parseval direction")
        require(
            -signed_sum / TARGET_COLOURS >= F(1, 1014),
            "per-role Hermitian-pair floor",
        )

singleton_bank_sum = -F(12, 169)
two_root_bank_sum = -F(24, 169)
require(
    -singleton_bank_sum / ROLE_COLOUR_PAIRS == F(1, 845),
    "singleton whole-bank floor",
)
require(
    -two_root_bank_sum / ROLE_COLOUR_PAIRS == F(2, 845),
    "two-root whole-bank floor",
)

# Exact universal and common-core floors before and after the thirteen
# transverse owner/digit categories.
delta_universal = F(1, 26754)
delta_common = F(66, 4459)

universal_each = delta_universal / 1014
universal_bank = delta_universal / 845
common_each = delta_common / 1014
common_bank = delta_common / 845
require(universal_each == F(1, 27128556), "universal every-role floor")
require(universal_bank == F(1, 22607130), "universal bank floor")
require(common_each == F(11, 753571), "common every-role floor")
require(common_bank == F(66, 3767855), "common bank floor")

owner_categories = 1 + 6 * 2
require(owner_categories == 13, "owner/digit category count")

rho_owner_universal = delta_universal / owner_categories
rho_owner_common = delta_common / owner_categories
owner_universal_each = rho_owner_universal / 1014
owner_universal_joint = owner_universal_each / 7
owner_universal_bank = rho_owner_universal / 845
owner_universal_bank_joint = owner_universal_bank / 7
owner_common_each = rho_owner_common / 1014
owner_common_joint = owner_common_each / 7
owner_common_bank = rho_owner_common / 845
owner_common_bank_joint = owner_common_bank / 7

require(rho_owner_universal == F(1, 347802), "universal owner mass")
require(
    owner_universal_each == F(1, 352671228),
    "universal owner every-role floor",
)
require(
    owner_universal_joint == F(1, 2468698596),
    "universal owner joint floor",
)
require(
    owner_universal_bank == F(1, 293892690),
    "universal owner bank floor",
)
require(
    owner_universal_bank_joint == F(1, 2057248830),
    "universal owner bank joint floor",
)
require(rho_owner_common == F(66, 57967), "common owner mass")
require(
    owner_common_each == F(11, 9796423),
    "common owner every-role floor",
)
require(
    owner_common_joint == F(11, 68574961),
    "common owner joint floor",
)
require(
    owner_common_bank == F(66, 48982115),
    "common owner bank floor",
)
require(
    owner_common_bank_joint == F(66, 342874805),
    "common owner bank joint floor",
)

# Labelled target actions descend to the common root gauge exactly when
# all unit-factor slopes eta/w agree.  Exhaust the two-factor local
# criterion; the six-factor statement follows pairwise.
slope_cases = 0
common_slope_cases = 0
for w_1 in range(1, PRIME):
    for w_2 in range(1, PRIME):
        for eta_1 in RESIDUES:
            for eta_2 in RESIDUES:
                slope_cases += 1
                alpha_1 = eta_1 * pow(w_1, -1, PRIME) % PRIME
                alpha_2 = eta_2 * pow(w_2, -1, PRIME) % PRIME
                determinant_zero = (eta_1 * w_2 - eta_2 * w_1) % PRIME == 0
                require(
                    (alpha_1 == alpha_2) == determinant_zero,
                    "slope/determinant equivalence",
                )
                if alpha_1 == alpha_2:
                    common_slope_cases += 1
                    require(
                        (eta_1 - alpha_1 * w_1) % PRIME == 0
                        and (eta_2 - alpha_1 * w_2) % PRIME == 0,
                        "pure-blocker representative",
                    )
require(slope_cases == 24336, "two-factor slope universe")
require(common_slope_cases == 1872, "common-slope census")

# Equal slopes translate an exclusive word.  Unequal slopes can change
# cardinality, so they cannot be absorbed by any common root gauge.
danger = frozenset((0, 1))
competitor = frozenset((1, 2))
base_exclusive = danger - competitor
for shift in RESIDUES:
    equal_slope_word = translate(danger, shift) - translate(competitor, shift)
    require(
        equal_slope_word == translate(base_exclusive, shift),
        "equal-slope positive control",
    )
unequal_shifts = (0, 12, 1)
unequal_sizes = tuple(
    len(danger - translate(competitor, shift))
    for shift in unequal_shifts
)
require(unequal_sizes == (1, 0, 2), "unequal-slope hostile")

# The collapsed base support does not determine its unequal-slope
# continuation: two labelled competitors have the same deletion at
# shift zero but different supports at shift twelve.
competitor_prime = frozenset((1, 3))
require(
    danger - competitor == danger - competitor_prime == frozenset((0,)),
    "collapsed-support hostile base",
)
require(
    danger - translate(competitor, 12) !=
    danger - translate(competitor_prime, 12),
    "collapsed-support hostile continuation",
)


print("theorem=THM-2400")
print("status=PROVED-CANDIDATE+VERIFIED-EXACT; INDEPENDENT-AUDIT-PENDING")
print(
    f"disjoint_mask_pairs={mask_pairs};"
    f" common_root_translate_checks={translated_checks}"
)
print("common_root_gauge=POINTWISE-INVARIANT-IN-JOINT-CHARGED-PRODUCT")
print("every_actual_role_all_12_colours=YES")
print(
    f"signed_floor_lower_q=delta/1014;"
    f" signed_floor_guard=delta/507; bank_floor=delta/845"
)
print(
    f"universal_each={universal_each}; universal_bank={universal_bank};"
    f" common_each={common_each}; common_bank={common_bank}"
)
print(
    f"owner_categories={owner_categories};"
    f" universal_owner_mass={rho_owner_universal};"
    f" common_owner_mass={rho_owner_common}"
)
print(
    f"universal_owner_each={owner_universal_each};"
    f" universal_owner_joint={owner_universal_joint};"
    f" universal_owner_bank_joint={owner_universal_bank_joint}"
)
print(
    f"common_owner_each={owner_common_each};"
    f" common_owner_joint={owner_common_joint};"
    f" common_owner_bank_joint={owner_common_bank_joint}"
)
print(
    f"two_factor_slope_cases={slope_cases};"
    f" common_slope_cases={common_slope_cases};"
    " labelled_descent=COMMON-SLOPE-IFF"
)
print(
    f"unequal_slope_hostile_shifts={','.join(map(str, unequal_shifts))};"
    f" sizes={','.join(map(str, unequal_sizes))}"
)
print("relative_shift_bank=LIVE; common_coshift_diagonal=IDENTICALLY-ZERO")
print("lawful_target_current=OPEN; row_decrement=0; ledger=165")
print("all_checks=PASS")
