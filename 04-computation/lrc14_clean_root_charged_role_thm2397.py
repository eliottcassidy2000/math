#!/usr/bin/env python3
"""Exact companion for THM-2397.

Checks the local least-role truth table, the root-energy and pigeonhole
ledgers, and every displayed rational mass/charged-product floor.
"""

from fractions import Fraction as F


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


ROLE_COUNT = 5
TARGET_COLOURS = 12
Q_STAR = ROLE_COUNT


# A root is represented by its nonempty active subset of the five other
# roles and q_*. The exclusive singleton {q_*} is exactly A. Every other
# root belongs to W and must have a least active role outside q_*.
activity_subsets = []
exclusive_q_star = 0
assigned_counts = [0] * ROLE_COUNT
for bits in range(1, 1 << (ROLE_COUNT + 1)):
    active = tuple(i for i in range(ROLE_COUNT + 1) if bits & (1 << i))
    activity_subsets.append(active)
    if active == (Q_STAR,):
        exclusive_q_star += 1
        continue
    other = tuple(i for i in active if i != Q_STAR)
    require(other, "a W-root had no lower role")
    assigned_counts[min(other)] += 1

require(len(activity_subsets) == 63, "local truth-table size")
require(exclusive_q_star == 1, "exclusive q_* state count")
require(sum(assigned_counts) == 62, "least-role partition count")


def nonzero_energy(mask_size: int) -> F:
    """Normalized sum over the twelve nonzero F_13 characters."""

    return F(13 * mask_size - mask_size * mask_size, 169)


singleton_energy = nonzero_energy(1)
adjacent_energy = nonzero_energy(2)
require(singleton_energy == F(12, 169), "singleton energy")
require(adjacent_energy == F(22, 169), "adjacent energy")

# An adjacent two-root transform cannot vanish at a thirteenth root:
# zeta^(k*eta)=-1 would give a nontrivial element of order two.
for k in range(1, 13):
    for eta in range(1, 13):
        require((2 * k * eta) % 13 != 0, "adjacent transform zero")

pair_count = ROLE_COUNT * TARGET_COLOURS
require(pair_count == 60, "role/colour pair count")
require(singleton_energy / pair_count == F(1, 845), "uniform pair floor")
require(adjacent_energy / pair_count == F(11, 5070), "adjacent pair floor")
colour_pairs = TARGET_COLOURS // 2
fixed_role_pairs = (colour_pairs + ROLE_COUNT - 1) // ROLE_COUNT
require(colour_pairs == 6, "conjugate colour pairs")
require(fixed_role_pairs == 2, "fixed-role pair pigeonhole")
fixed_role_colours = 2 * fixed_role_pairs
require(fixed_role_colours == 4, "fixed-role colour support")

# Universal last-septimal-lane floors.
delta_universal = F(1, 26754)
rho_universal = delta_universal / 52
require(rho_universal == F(1, 1391208), "universal charged cell")
universal_fibre = rho_universal / 845
universal_aggregate = rho_universal * rho_universal / 845
require(universal_fibre == F(1, 1175570760), "universal fibrewise floor")
require(
    universal_aggregate == F(1, 1635463445878080),
    "universal aggregate floor",
)

# Sharper common-core floors from THM-2396.
delta_common = F(66, 4459)
rho_common = delta_common / 52
require(rho_common == F(33, 115934), "common charged cell")
common_fibre = rho_common / 845
common_aggregate = rho_common * rho_common / 845
require(common_fibre == F(33, 97964230), "common fibrewise floor")
require(
    common_aggregate == F(1089, 11357385040820),
    "common aggregate floor",
)

# Owner-resolved C_7 x C_13 tensor.
rho_owner = delta_common / 338
require(rho_owner == F(33, 753571), "owner-resolved cell")
owner_fibre = rho_owner / 845
owner_aggregate = rho_owner * rho_owner / 845
owner_tensor = owner_aggregate / 7
require(owner_fibre == F(33, 636767495), "owner fibrewise floor")
require(
    owner_aggregate == F(1089, 479849517974645),
    "owner aggregate floor",
)
require(
    owner_tensor == F(1089, 3358946625822515),
    "owner septimal tensor floor",
)


print("theorem=THM-2397")
print("status=PROVED-CANDIDATE+VERIFIED-EXACT; independent-audit=PENDING")
print(
    f"local_activity_states={len(activity_subsets)};"
    f" exclusive_qstar={exclusive_q_star}; assigned_W={sum(assigned_counts)}"
)
print(
    f"singleton_energy={singleton_energy}; adjacent_energy={adjacent_energy};"
    f" role_colour_pairs={pair_count}"
)
print("every_nonzero_colour_has_lower_role=YES; one_role_colours>=4")
print(
    f"universal_rho={rho_universal}; fibrewise>={universal_fibre};"
    f" aggregate>={universal_aggregate}"
)
print(
    f"common_rho={rho_common}; fibrewise>={common_fibre};"
    f" aggregate>={common_aggregate}"
)
print(
    f"owner_rho={rho_owner}; fibrewise>={owner_fibre};"
    f" aggregate>={owner_aggregate}; septimal_tensor>={owner_tensor}"
)
print("row_decrement=0; ledger=165; LRC(14)=OPEN")
print("all_checks=PASS")
