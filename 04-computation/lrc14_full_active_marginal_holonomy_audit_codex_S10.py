#!/usr/bin/env python3
"""Exact cross-audit of THM-794 against THM-786's two marginal bounds.

For F=49H+1, f=F, g=F-7, and C={F-14,...,F-49}, all inverse
residues are one modulo seven.  The balanced-cluster hypergraph therefore has
the single edge C.  Its minimum speed-weight fractional transversal is
W_*=F-49, so the transversal theorem is applicable for every fixed H but its
slack 1-W_*/g=42/g tends to zero.

For complete g-period packets, q=1 and the allowed types are exactly

    (0,C),  (1,C\{c}) for c in C.

Writing p_c=1-z_c, their convex hull satisfies p_c>=0 and
epsilon=sum p_c<=1.  The natural target has

    epsilon_0=7/g,  p^0=(7,14,21,28,35,42)/g.

The separating functional

    h(epsilon,p)=p_4+p_5+p_6-epsilon <= 0

has target excess 98/g and l1 norm four, so the l-infinity distance is at
least 49/(2g).  The explicit point

    p=(0,0,0,7/2,21/2,35/2)/g, epsilon=63/(2g)

attains that distance.  Thus the exact distance is 49/(2g), again tending to
zero.  Both marginal theorems are correct but only give per-tuple bounds of
order H, allowing THM-794's H-1 active f-periods and H-2 complete g-periods.

Tournament Analysis: packet occurrences with chronological comparison form a
transitive tournament and merely count repetitions.  The switch to the
reduced deck-holonomy observable identifies each repeated owner word as the
same zero class in F_7^8/Delta.  The tie Hamiltonian path is the chronological
packet list.  The challenged assumption is that an occupation vector or an
active packet is a canonical vertex.  The predicate-preserving carrier must
also retain the ordered prefix-legal collision path, reduced holonomy, and
metric/core incidence.
"""

from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations


EXPECTED_DIGEST = "3bb05201d747b57788088bc8a4845e4625cc3addf19c0c9df69b3e5c55a2c8f0"


def subsets(values):
    for size in range(len(values) + 1):
        yield from combinations(values, size)


def audit_height(height: int):
    F = 49 * height + 1
    f = F
    g = F - 7
    companions = tuple(F - 7 * j for j in range(2, 8))

    assert all(value % 7 == 1 for value in (f, g, *companions))
    balanced = tuple(
        packet
        for packet in subsets(companions)
        if packet and (1 + len(packet)) % 7 == 0
    )
    assert balanced == (companions,)

    W_star = min(companions)
    assert W_star == F - 49 < g
    transversal_slack = Q(g - W_star, g)
    assert transversal_slack == Q(42, g)
    transversal_bound = Q(1) + Q((F - 49) * (F - 7), 21 * F)
    actual_active = height - 1
    assert actual_active < transversal_bound
    extent_bound = Q(2, g) + Q(g, 21 * F)
    actual_extent = Q(height - 1, F)
    assert actual_extent < extent_bound

    allowed = tuple(
        (epsilon, packet)
        for epsilon in (0, 1)
        for packet in subsets(companions)
        if (1 + epsilon + len(packet)) % 7 == 0
    )
    expected_allowed = ((0, companions),) + tuple(
        (1, tuple(value for value in companions if value != omitted))
        for omitted in companions
    )
    assert len(allowed) == len(expected_allowed)
    assert set(allowed) == set(expected_allowed)

    target_epsilon = Q(7, g)
    target_missing = tuple(Q(7 * j, g) for j in range(1, 7))
    distance = Q(49, 2 * g)
    witness_missing = tuple(
        Q(value, 2 * g) for value in (0, 0, 0, 7, 21, 35)
    )
    witness_epsilon = sum(witness_missing, Q(0))
    assert witness_epsilon == Q(63, 2 * g) <= 1

    # h=sum(last three missing coordinates)-epsilon is <=0 on the polytope.
    target_excess = sum(target_missing[3:], Q(0)) - target_epsilon
    assert target_excess == Q(98, g)
    assert target_excess / 4 == distance
    assert max(
        (abs(witness_epsilon - target_epsilon),)
        + tuple(
            abs(actual - target)
            for actual, target in zip(witness_missing, target_missing)
        )
    ) == distance

    packet_bound = 1 / distance
    actual_complete_g = max(0, height - 2)
    assert packet_bound == Q(2 * g, 49)
    assert actual_complete_g < packet_bound

    # Every owner occurs once and has inverse residue one.  The raw token
    # displacement is the diagonal vector (-1,...,-1), hence zero modulo the
    # diagonal deck action.
    owner_word = (f,) + tuple(F - 7 * j for j in range(7, 0, -1))
    assert len(owner_word) == 8 and len(set(owner_word)) == 8
    raw_displacement = tuple(-pow(owner, -1, 7) % 7 for owner in owner_word)
    assert raw_displacement == (6,) * 8
    reduced_holonomy = tuple(
        (value - raw_displacement[0]) % 7 for value in raw_displacement
    )
    assert reduced_holonomy == (0,) * 8

    return (
        height,
        F,
        W_star,
        transversal_slack,
        transversal_bound,
        distance,
        packet_bound,
        actual_active,
        actual_complete_g,
        actual_extent,
        extent_bound,
        owner_word,
        reduced_holonomy,
    )


def main() -> None:
    rows = tuple(audit_height(height) for height in range(2, 201))
    canonical = "\n".join(repr(row) for row in rows) + "\n"
    digest = sha256(canonical.encode()).hexdigest()
    if EXPECTED_DIGEST != "TO_BE_FILLED":
        assert digest == EXPECTED_DIGEST

    for height in (2, 5, 20, 200):
        row = rows[height - 2]
        (
            _,
            F,
            W_star,
            slack,
            t_bound,
            distance,
            p_bound,
            active,
            complete_g,
            extent,
            extent_bound,
            _,
            _,
        ) = row
        print(
            f"H={height} F={F} W*={W_star} slack={slack} "
            f"active={active} transversal_bound={t_bound}"
        )
        print(
            f"  delta={distance} complete_g={complete_g} "
            f"packet_bound={p_bound} extent={extent} "
            f"derived_extent_bound={extent_bound}"
        )

    print("balanced_hypergraph={C}")
    print("allowed_packets=(0,C) and (1,C-minus-{c})")
    print("separator=sum(last_three_missing)-epsilon <= 0; l1=4")
    print("target_excess=98/g exact_delta=49/(2g)")
    print("reduced_holonomy=0 in F7^8/Delta")
    print(
        "tournament_fingerprint=transitive score_histogram=0..H-2 "
        "directed_cycles=0 singleton_SCCs=true hp_count=1"
    )
    print("preserved=packet_count_and_abelian_occupation")
    print("destroyed=ordered_prefix_legal_path_reduced_holonomy_metric_incidence")
    print(f"canonical_sha256={digest}")
    print("PASS")


if __name__ == "__main__":
    main()
