#!/usr/bin/env python3
"""Bounded exact audit of the THM-2809/2819/2818 label interface.

This deliberately audits only consequences that follow from the displayed
canonical lift, support tables, and already exact selected-cell sidecars.  It
does not identify the slope-seven chart label with the THM-2771 endpoint-
dipole label as a typed physical object.
"""

from collections import Counter
from fractions import Fraction


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


P = 13
R = P**6
SOURCE_CARRY = 12
TARGET_CARRY = 6


def lift(label):
    """Canonical representative lift tau_delta=7 delta/R, 0<=delta<13."""
    return Fraction(7 * label, R)


def carry(base, label):
    return (base + 7 * label) % P


source_face = tuple(range(2, P))
transported_face = tuple((label - 1) % P for label in source_face)

# THM-2771/2818 nonzero-cell support, reconstructed from its three printed
# primitive target polynomials rather than copied as a 28-cell tuple.
clock_support = {
    1: tuple(range(3, 13)),
    2: (*range(2, 10), 12),
    3: (*range(2, 8), 10, 11, 12),
}
cofiber_cells = tuple(
    (clock, label)
    for clock, labels in clock_support.items()
    for label in labels
)
cofiber_projection = tuple(sorted({label for _clock, label in cofiber_cells}))
projection_multiplicity = Counter(label for _clock, label in cofiber_cells)


# For non-wrapping labels the proposed endpoint conjugacy is a literal
# equality of translations and of carry strata:
#
#   T_tau1 T_-taudelta P_(12+7delta)
#       = T_-tau(delta-1) P_(6+7(delta-1)).
#
# At delta=0, canonical reduction delta-1=12 introduces the odometer kernel.
nonwrap_rows = []
for label in range(1, P):
    target_label = label - 1
    left_translation = lift(1) - lift(label)
    right_translation = -lift(target_label)
    require(left_translation == right_translation,
            f"nonwrap translation failed at label {label}")
    require(carry(SOURCE_CARRY, label)
            == carry(TARGET_CARRY, target_label),
            f"nonwrap carry failed at label {label}")
    nonwrap_rows.append((
        label,
        target_label,
        carry(SOURCE_CARRY, label),
        left_translation,
    ))

wrap_left = lift(1) - lift(0)
wrap_right = -lift(P - 1)
wrap_defect = wrap_left - wrap_right
require(
    wrap_defect == Fraction(91, R) == Fraction(7, P**5),
    "wrap defect stopped being the odometer kernel",
)
require(
    carry(SOURCE_CARRY, 0) == carry(TARGET_CARRY, P - 1) == 12,
    "wrap carries no longer agree",
)

# THM-2818's exceptional chain lattice has physical step 1/(2*13^5).
# The conjugacy defect is fourteen such steps, so it preserves the 1010...
# parity shadow.  This is compatibility, not a physical identification.
half_step = Fraction(1, 2 * P**5)
require(wrap_defect == 14 * half_step,
        "wrap defect/exceptional-chain scale relation changed")
require(14 % 2 == 0, "wrap defect stopped preserving chain parity")

require(len(source_face) == 11, "THM-2809 source face changed")
require(
    transported_face == tuple(range(1, 12)),
    "nonwrap positive-face transfer changed",
)
require(len(cofiber_cells) == 28, "THM-2771 support size changed")
require(
    cofiber_projection == source_face,
    "cofiber projection stopped matching the source eleven-spine",
)
require(
    tuple(sorted(
        label for label, count in projection_multiplicity.items()
        if count == 2
    )) == (2, 8, 9, 10, 11),
    "two-clock spine fibres changed",
)
require(
    tuple(sorted(
        label for label, count in projection_multiplicity.items()
        if count == 3
    )) == (3, 4, 5, 6, 7, 12),
    "three-clock spine fibres changed",
)
require(
    set(cofiber_projection) ^ set(transported_face) == {1, 12},
    "identity-coordinate endpoint mismatch changed",
)
require(
    tuple(sorted((label - 1) % P for label in cofiber_projection))
    == transported_face,
    "affine shifted cofiber shadow stopped matching transported face",
)

# The entire nonlinear copy-count anomaly is at the canonical wrap label 12.
raw_at_wrap = (241, 528, 506)
live_at_wrap = (121, 265, 254)
dead_at_wrap = (120, 263, 252)
exceptional_block_lengths = {
    1: (145, 96),
    2: (143, 289, 96),
    3: (143, 289, 74),
}
require(
    tuple(live + dead for live, dead in zip(live_at_wrap, dead_at_wrap))
    == raw_at_wrap,
    "exceptional raw/live/dead partition changed",
)
require(sum(live_at_wrap) == 640 and sum(dead_at_wrap) == 635,
        "exceptional aggregate changed")
require(
    all(label != 1 for _clock, label in cofiber_cells)
    and all((clock, 12) in cofiber_cells for clock in (1, 2, 3)),
    "endpoint hostile labels changed",
)

# The independently reconstructed gap between consecutive blocks is 50h.
# Since kappa=14h, no piece can jump from one discovered block to another.
# Within a block, translation by kappa maps slot i to slot i+14.  It
# therefore leaves L-14 raw slots in common and removes exactly seven live
# and seven zero slots at the outgoing boundary of every block.
holonomy_overlap = {}
for clock, lengths in exceptional_block_lengths.items():
    raw_overlap = sum(length - 14 for length in lengths)
    live_overlap = sum((length - 14 + 1) // 2 for length in lengths)
    dead_overlap = raw_overlap - live_overlap
    block_count = len(lengths)
    require(
        live_at_wrap[clock - 1] - live_overlap == 7 * block_count
        and dead_at_wrap[clock - 1] - dead_overlap == 7 * block_count,
        f"holonomy boundary flux changed at clock {clock}",
    )
    holonomy_overlap[clock] = (
        raw_overlap,
        live_overlap,
        dead_overlap,
        7 * block_count,
        7 * block_count,
    )
require(
    holonomy_overlap
    == {
        1: (213, 107, 106, 14, 14),
        2: (486, 244, 242, 21, 21),
        3: (464, 233, 231, 21, 21),
    },
    "exceptional holonomy overlap table changed",
)

# Selected-cell hostile from the independently exact THM-2818 bridge.
# Both coefficient-effective copies fail native E3, one also native c2;
# their source carrier masks are empty.  Thus a "two copies = two edge
# choices" lift cannot preserve the marked-source packet at (clock,t)=(1,4).
factor_order = ("E3", "clock", "q1", "q2", "c2", "c3")
selected_native_masks = (
    (False, True, True, True, True, True),
    (False, True, True, True, False, True),
)
selected_source_carrier_counts = (0, 0)
require(
    all(not mask[factor_order.index("E3")]
        for mask in selected_native_masks),
    "selected copies regained native E3",
)
require(
    not selected_native_masks[1][factor_order.index("c2")],
    "selected plus copy regained native c2",
)
require(
    selected_source_carrier_counts == (0, 0),
    "selected copies regained source carrier support",
)

print("ENDPOINT CONJUGACY / ELEVEN-SPINE / POSITIVE-COPY AUDIT")
print("status=FINITE-EXACT scratch; no typed chart-to-cofiber map")
print(f"nonwrap_conjugacy_rows={tuple(nonwrap_rows)}")
print(
    f"wrap=(left={wrap_left},right={wrap_right},"
    f"defect={wrap_defect}=91/R=7/13^5)"
)
print(
    f"exceptional_chain_half_step={half_step};"
    "wrap_defect_in_half_steps=14;parity_shadow=preserved"
)
print(f"thm2809_source_face={source_face}")
print(f"nonwrap_transported_face={transported_face}")
print(f"thm2818_clock_support={clock_support}")
print(f"thm2818_projection={cofiber_projection}")
print(
    "projection_fibres="
    f"{tuple(sorted(projection_multiplicity.items()))};"
    "six_triple+five_double=28"
)
print(
    "identity_coordinate_symmetric_difference=(1,12);"
    "shifted_shadow_match=t->t-1"
)
print(
    f"wrap_column_raw_live_dead="
    f"{tuple(zip(raw_at_wrap, live_at_wrap, dead_at_wrap))};"
    "aggregate=1275=640_live+635_zero"
)
print(
    f"wrap_holonomy_block_overlap={holonomy_overlap};"
    "interblock_gap=50_half_steps;outgoing_flux_per_block=7_live+7_zero"
)
print(
    "selected_cell_hostile=(clock,t)=(1,4):"
    "both copies fail native E3;plus fails native c2;"
    "both source carrier masks empty"
)
print(
    "verdict=exact conjugacy on the twelve nonwrap labels and exact positive "
    "eleven-face transfer; cyclic conjugacy fails by odometer holonomy; "
    "cofiber support is the same eleven-spine only as a label shadow"
)
print("ALL EXACT CHECKS PASSED")
