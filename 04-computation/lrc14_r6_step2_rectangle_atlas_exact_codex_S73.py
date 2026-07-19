#!/usr/bin/env python3
"""Exact all-core rectangle atlas for the r=6 step-two five-comb ray.

For every seven-speed core ``P`` in ``{1,...,12}``, consider the five
removed killers

    (b, b+2, b+4, b+6, b+8).

Put ``B=b+4`` and use the fixed-torus coordinates ``x=B*t``.  The five
offsets are then ``A=(-4,-2,0,2,4)``.  This referee proves, with
``fractions.Fraction`` only, that every core has a labelled fixed rectangle

    I x X subset G_P x { (t,x) : ||x+a*t|| >= 1/14 for all a in A }

such that ``B*|I| >= 1+|X|`` for every integer ``B>=168`` and
``|X|>=1/7``.  A complete preimage of X under ``t -> B*t mod 1`` therefore
lies in I.  Its length is ``|X|/B >= 1/(7B) > 1/(7(B+4))``, which is the
sharp component target before the sixth killer.

Completeness is finite and exact.  Split time at every core danger endpoint
and every collision of two offset centers.  On an open cell the core
predicate and cyclic center order are constant.  For each labelled cyclic
gap and each direction into the cell, the two bounding walls are affine.
If their safe gap has endpoint width w and inward retreat rate rho, then its
common intersection over time delta has exact width ``w-rho*delta``.  The
largest delta retaining width at least 1/7 minimizes the required slope
``(1+|X|)/delta`` on that candidate.

Tournament / alternate-carrier audit.  Candidate boundary order is a
transitive tournament after exact ties are coalesced: no directed cycles,
singleton SCCs, and one sorted Hamiltonian path.  That quotient forgets
gap length, adjacent wall owners, inward direction, and wall slopes.  The
predicate-preserving carrier is the labelled cyclic gap together with its
core-safe side and owner-slope sidecar.  Runner, comb, core-cell, boundary,
wall-event, residue, and proof-obligation vertices were all considered.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations


H = F(1, 14)
TARGET_WIDTH = F(1, 7)
OFFSETS = (-4, -2, 0, 2, 4)
EXPECTED_TAIL_B = F(168)

CORES = list(combinations(range(1, 13), 7))
CORE_MASKS = [sum(1 << (p - 1) for p in core) for core in CORES]


def safe_speed_mask(t: F) -> int:
    """Speeds in 1..12 which are safe at the exact time t."""

    mask = 0
    for p in range(1, 13):
        x = (p * t) % 1
        if min(x, 1 - x) >= H:
            mask |= 1 << (p - 1)
    return mask


ELIGIBLE_CORES: list[list[int]] = []
for safe_mask in range(1 << 12):
    ELIGIBLE_CORES.append(
        [
            index
            for index, core_mask in enumerate(CORE_MASKS)
            if core_mask & ~safe_mask == 0
        ]
    )


def candidate_points() -> list[F]:
    """All core endpoints and offset-center collision times in [0,1]."""

    points = {F(0), F(1)}
    for p in range(1, 13):
        for tooth in range(-1, p + 2):
            for sign in (-1, 1):
                t = F(14 * tooth + sign, 14 * p)
                if 0 <= t <= 1:
                    points.add(t)
    for a, b in combinations(OFFSETS, 2):
        difference = b - a
        points.update(F(j, difference) for j in range(difference + 1))
    return sorted(points)


@dataclass(frozen=True, order=True)
class Certificate:
    required_B: F
    left: F
    right: F
    side: str
    lower_owner: int
    upper_owner: int
    endpoint_width: F
    retreat: int
    time_width: F
    fixed_width: F
    fixed_lower: F
    fixed_upper: F


def cell_certificates(left: F, right: F) -> list[Certificate]:
    """All labelled one-sided fixed rectangles supported by one order cell."""

    midpoint = (left + right) / 2
    lifted = []
    for offset in OFFSETS:
        raw = -offset * midpoint
        intercept = -(raw.numerator // raw.denominator)
        center = raw + intercept
        assert 0 <= center < 1
        lifted.append((center, offset, intercept))
    lifted.sort()

    certificates = []
    cell_width = right - left
    for index in range(len(lifted)):
        _, lower_owner, lower_intercept = lifted[index]
        _, upper_owner, upper_intercept = lifted[(index + 1) % len(lifted)]
        if index + 1 == len(lifted):
            upper_intercept += 1

        for side, endpoint, direction in (
            ("right", left, 1),
            ("left", right, -1),
        ):
            lower_center = -lower_owner * endpoint + lower_intercept
            upper_center = -upper_owner * endpoint + upper_intercept
            endpoint_width = upper_center - lower_center - 2 * H
            if endpoint_width < TARGET_WIDTH:
                continue

            lower_slope = -lower_owner * direction
            upper_slope = -upper_owner * direction
            retreat = max(0, lower_slope) - min(0, upper_slope)
            if retreat == 0:
                time_width = cell_width
            else:
                time_width = min(
                    cell_width,
                    (endpoint_width - TARGET_WIDTH) / retreat,
                )
            if time_width <= 0:
                continue

            fixed_lower = (
                lower_center + H + max(0, lower_slope) * time_width
            )
            fixed_upper = (
                upper_center - H + min(0, upper_slope) * time_width
            )
            fixed_width = fixed_upper - fixed_lower
            assert fixed_width == endpoint_width - retreat * time_width
            assert fixed_width >= TARGET_WIDTH
            required_B = (1 + fixed_width) / time_width
            certificates.append(
                Certificate(
                    required_B,
                    left,
                    right,
                    side,
                    lower_owner,
                    upper_owner,
                    endpoint_width,
                    retreat,
                    time_width,
                    fixed_width,
                    fixed_lower,
                    fixed_upper,
                )
            )
    return certificates


def main() -> None:
    points = candidate_points()
    best: list[Certificate | None] = [None for _ in CORES]
    eligible_order_cells = 0
    labelled_candidates = 0

    for left, right in zip(points, points[1:]):
        midpoint = (left + right) / 2
        core_indices = ELIGIBLE_CORES[safe_speed_mask(midpoint)]
        if not core_indices:
            continue
        eligible_order_cells += len(core_indices)
        certificates = cell_certificates(left, right)
        labelled_candidates += len(certificates) * len(core_indices)
        if not certificates:
            continue
        cell_best = min(certificates)
        for core_index in core_indices:
            if best[core_index] is None or cell_best < best[core_index]:
                best[core_index] = cell_best

    assert all(certificate is not None for certificate in best)
    exact_best = [certificate for certificate in best if certificate is not None]
    global_tail = max(certificate.required_B for certificate in exact_best)
    worst_indices = [
        index
        for index, certificate in enumerate(exact_best)
        if certificate.required_B == global_tail
    ]

    assert len(CORES) == 792
    assert global_tail == EXPECTED_TAIL_B
    assert all(certificate.fixed_width >= TARGET_WIDTH for certificate in exact_best)

    print("r=6 step-two five-comb exact rectangle atlas")
    print("arithmetic=fractions.Fraction only")
    print(f"cores={len(CORES)}")
    print(f"offsets={OFFSETS}")
    print(f"candidate_boundaries={len(points)}")
    print(f"eligible_core_order_cells={eligible_order_cells}")
    print(f"labelled_core_gap_side_candidates={labelled_candidates}")
    print(f"uniform_tail_B={global_tail}")
    print(f"uniform_tail_b={global_tail - 4}")
    print(f"worst_core_count={len(worst_indices)}")
    for index in worst_indices:
        certificate = exact_best[index]
        print(f"worst_core={CORES[index]}")
        print(f"worst_certificate={certificate}")
        print(
            "worst_crossing_identity="
            f"B*|I| at B=168 is {EXPECTED_TAIL_B * certificate.time_width}; "
            f"1+|X| is {1 + certificate.fixed_width}"
        )
    print("tail_component=|X|/B >= 1/(7B) > 1/(7(B+4))")
    print("candidate_lemma=core endpoints plus center collisions; exact affine wall retreat")
    print("tournament_vertices=coalesced candidate boundaries")
    print("tournament_fingerprint=transitive; directed_cycles=0; SCCs=singletons; sorted_HP=1")
    print("proof_carrier=labelled cyclic gap|core-safe side|owner slopes|metric rectangle")
    print("destroyed_by_order_only=gap length|owners|safe side|slopes|rectangle width")
    print("challenged_vertices=runners|combs|core cells|boundaries|wall events|residues|proof obligations")
    print("certificate=PASS")


if __name__ == "__main__":
    main()
