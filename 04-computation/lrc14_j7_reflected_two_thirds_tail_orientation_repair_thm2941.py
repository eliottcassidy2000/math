#!/usr/bin/env python3
r"""Repair the assignment cover in the reflected cap-5/2 tail.

The old low/high split for ``H=(1,2,3,4,6,12)`` reversed both the channel
and the ordered label pair.  Its two lanes therefore encoded the same
physical inequality ``q_0 <= q_1`` and did not cover the reverse ordering.
The analogous two lanes for ``H2`` had acyclic simultaneous failures.

This referee preserves that correction lineage and supplies a replacement:

* five Farey lanes, all on the same ordered physical pair (0,1), cover H;
* one whole-interval lane on (0,1) covers H2;
* singleton debt is coupled to the pair scale through
  ``max(q_e) <= (5/2) min(q_e)``;
* every finite reduced channel is checked directly until the transported
  primitive bracket becomes monotone and uniform;
* the predecessor's conclusion is not inherited: repeated-level closure,
  the CSP census, both orientations of every constrained located policy,
  and all ten unaffected whole-interval tails are recomputed here.

Together with the pinned exact CSP and valid policies in the predecessor,
this repairs the reflected ``D>=6, 3m>=2D`` sufficient-family theorem.  It is
not a physical-survivor census and not a proof of LRC(14).
"""

from __future__ import annotations

import argparse
import hashlib
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import ceil, gcd, lcm
from pathlib import Path


HERE = Path(__file__).resolve()
ROOT = next(parent for parent in HERE.parents if (parent / "00-navigation").is_dir())
OLD = ROOT / "04-computation/lrc14_j7_reflected_two_thirds_cone_closure_thm2941.py"
OLD_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_two_thirds_cone_closure_thm2941.out"
OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_two_thirds_tail_orientation_repair_thm2941.out"

EXPECTED_OLD_SHA256 = "e6e64c909c6bfcc776eaf6bf2ad210f75675a6a32a46ead70d8d29a37f607eb3"
EXPECTED_OLD_OUTPUT_SHA256 = "e77929c87f9d9f8fb7ce3a347c48522e87f63fa7c7085eb9e2cd8fe0bb4e4a90"
EXPECTED_OLD_SEMANTIC = "b34db3e25b9e4c81c1549f1ec7c7ab78e935ec778610daf617934b66b3a47304"
EXPECTED_TRAP_DIGEST = "b8e1f964d610af641a229d06951c4805167af4310429740501809034e9b2a716"
EXPECTED_VERDICT_DIGEST = "30650baad0b872b3df4a7278320f6ae4dbc7e3214c6aef16fb90e3a569b800c6"
EXPECTED_LOCATED_DIGEST = "e0a6c267f4f77017ea56229f3f2c80ba1b1a9c00b065ad9ab0540397ce238a75"
EXPECTED_HEAD_DIGEST = "0253483de002db1e2fa32a3feca1811f28b355b43d6f2426d19ba6fdfe58ffdc"
EXPECTED_SEMANTIC = "525c8d2d89cf11c9fcc6260dbf81e57d0e47c999049c89b9148ae731dbac4730"

CAP = F(5, 2)
MIN_LEVEL = 4
ETA_GUARD = F(29, 165)
H = (1, 2, 3, 4, 6, 12)
H2 = (1, 3, 4, 6, 8, 12)

# All five lanes use the same ordered pair (0,1).  Therefore their closed
# intervals cover the ratio q_1/q_0 for both physical orderings without an
# extra combinatorial assumption.
H_LANES = (
    (F(2, 5), F(3, 7), 48, F(16, 5)),
    (F(3, 7), F(1, 2), 38, F(22, 7)),
    (F(1, 2), F(2, 3), 26, F(3)),
    (F(2, 3), F(11, 9), 16, F(8, 3)),
    (F(11, 9), F(5, 2), 7, F(14, 9)),
)
H2_LANE = (F(2, 5), F(5, 2), 14, F(26, 5))


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def qtext(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def import_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("cannot import", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


require(sha256(OLD) == EXPECTED_OLD_SHA256, ("old source changed", sha256(OLD)))
require(sha256(OLD_OUTPUT) == EXPECTED_OLD_OUTPUT_SHA256,
        ("old output changed", sha256(OLD_OUTPUT)))
require(f"semantic_sha256={EXPECTED_OLD_SEMANTIC}" in OLD_OUTPUT.read_text(),
        "old semantic token missing")
M = import_module("two_thirds_orientation_repair_old", OLD)
T = M.T
T.MIN_LEVEL = MIN_LEVEL
T.RATIO_CAP = CAP


def interval_channels(left: F, right: F, start: int, oversize: int = 3):
    return tuple(
        (P, Q)
        for P in range(1, oversize * start + 1)
        for Q in range(1, oversize * start + 1)
        if P != Q and gcd(P, Q) == 1 and min(P, Q) < start
        and 1 / CAP <= F(Q, P) <= CAP and left <= F(Q, P) <= right
    )


def singleton_debt(body, ruler: int, level: int) -> F:
    return T.C2.singleton_debt(body, ruler, level)


def singleton_debt_drop(body, ruler: int, first: int, second: int) -> F:
    """Exact monotonicity gate for a (possibly multi-step) level increase."""
    require(MIN_LEVEL <= first <= second, (body, first, second))
    drop = singleton_debt(body, ruler, first) - singleton_debt(body, ruler, second)
    require((drop == 0) if first == second else (drop > 0),
            (body, first, second, drop))
    # Termwise this is e/[7(qL-e)]; every denominator increases with q.
    require(all(first * ruler > e for e in body), (body, ruler, first))
    return drop


def coupled_tail_envelope(body, pair, scale: int, numerator: F) -> tuple[F, int, F]:
    """Tail bound using the global minimum forced by the pair scale.

    If one pair has minimum level at least ``scale``, then its maximum level
    is at least ``scale``.  The cone gives ``q_max <= CAP*m``, hence
    ``m >= ceil(scale/CAP)``.  Keeping only this weakest implication is safe.
    """
    ruler = 14 * lcm(*body)
    a = body[pair[0]]
    debt_level = max(MIN_LEVEL, ceil(F(scale, CAP)))
    debt = singleton_debt(body, ruler, debt_level)
    value = (F(1, 49) - F(12, 49 * scale**2)
             - numerator / (ruler - F(a, scale)) - debt)
    return value, debt_level, debt


def analytic_tail_monotonicity(body, pair, scale: int, numerator: F,
                               coupled: bool) -> tuple:
    """A one-row certificate for identities valid at every later scale.

    Replacing ``scale`` by any integer ``s>=scale`` in the displayed exact
    identities leaves all denominators positive.  Thus this is an algebraic
    infinite-tail proof, not a finite sample of later scales.
    """
    ruler = 14 * lcm(*body)
    a = body[pair[0]]
    next_scale = scale + 1
    first_level = (max(MIN_LEVEL, ceil(F(scale, CAP)))
                   if coupled else MIN_LEVEL)
    second_level = (max(MIN_LEVEL, ceil(F(next_scale, CAP)))
                    if coupled else MIN_LEVEL)
    phase_gain = F(12, 49) * (F(1, scale**2) - F(1, next_scale**2))
    transport_gain = numerator * (
        F(scale, ruler * scale - a)
        - F(next_scale, ruler * next_scale - a)
    )
    require(phase_gain > 0, (body, scale, phase_gain))
    require(transport_gain == F(numerator * a,
            (ruler * scale - a) * (ruler * next_scale - a)) > 0,
            (body, scale, transport_gain))
    debt_gain = singleton_debt_drop(body, ruler, first_level, second_level)
    current = (coupled_tail_envelope(body, pair, scale, numerator)[0]
               if coupled else T.tail_envelope(body, pair, scale, numerator))
    later = (coupled_tail_envelope(body, pair, next_scale, numerator)[0]
             if coupled else T.tail_envelope(body, pair, next_scale, numerator))
    require(later - current == phase_gain + transport_gain + debt_gain > 0,
            (body, pair, scale, current, later, phase_gain, transport_gain, debt_gain))
    return (body, pair, scale, numerator, first_level, second_level,
            phase_gain, transport_gain, debt_gain, later - current)


def scale_monotonicity(body, pair, channel, scale: int, coupled: bool) -> tuple:
    """Certify the homotopy losses and ceiling debt for every later g.

    The identities checked below are rational identities in an arbitrary
    positive integer ``g``.  Their denominators remain positive after replacing
    this witness scale by any larger integer.  The ceiling level is monotone
    because it is the ceiling of a positive linear function of ``g``.
    """
    ruler = 14 * lcm(*body)
    a, b = body[pair[0]], body[pair[1]]
    P, Q = channel
    require(P > 0 and Q > 0 and P * scale * ruler > a,
            (body, pair, channel, scale))
    mismatch = abs(Q * a - P * b)
    eta_now = F(scale * mismatch, P * scale * ruler - a)
    eta_next = F((scale + 1) * mismatch, P * (scale + 1) * ruler - a)
    phase_drop = F(mismatch * a,
                   (P * scale * ruler - a)
                   * (P * (scale + 1) * ruler - a))
    require(eta_now - eta_next == phase_drop,
            (body, pair, channel, scale, eta_now, eta_next, phase_drop))
    require((phase_drop == 0) if mismatch == 0 else (phase_drop > 0),
            (body, pair, channel, mismatch, phase_drop))

    # c^{-1}-1=a/(PgL-a): the transport displacement falls strictly.  Once
    # the primitive bracket is positive we discard this favourable factor.
    transport_now = F(a, P * scale * ruler - a)
    transport_next = F(a, P * (scale + 1) * ruler - a)
    transport_drop = F(a * P * ruler,
                       (P * scale * ruler - a)
                       * (P * (scale + 1) * ruler - a))
    require(transport_now - transport_next == transport_drop > 0,
            (body, pair, channel, scale, transport_drop))

    if coupled:
        first_level = max(MIN_LEVEL, ceil(F(scale * max(P, Q), CAP)))
        second_level = max(MIN_LEVEL, ceil(F((scale + 1) * max(P, Q), CAP)))
    else:
        first_level = second_level = MIN_LEVEL
    require(second_level >= first_level, (channel, scale, first_level, second_level))
    debt_drop = singleton_debt_drop(body, ruler, first_level, second_level)
    require(min(P, Q) * scale >= MIN_LEVEL and eta_now <= ETA_GUARD < 1,
            (body, pair, channel, scale, eta_now))
    return (body, pair, channel, scale, coupled, mismatch, phase_drop,
            transport_drop, first_level, second_level, debt_drop)


def physical_predicate(side: str, pair: tuple[int, int]) -> tuple[int, int]:
    """Return the weak order encoded by a lane as ``q_left <= q_right``."""
    i, j = pair
    require(side in {"low", "high"}, side)
    return (j, i) if side == "low" else (i, j)


def old_lane_applies(side: str, pair: tuple[int, int], levels: tuple[int, ...]) -> bool:
    """Evaluate the old lane on physical packet levels, without relabelling."""
    i, j = pair
    ratio = F(levels[j], levels[i])
    require(1 / CAP <= ratio <= CAP, (side, pair, levels, ratio))
    if side == "all":
        return True
    require(side in {"low", "high"}, side)
    return ratio < 1 if side == "low" else ratio > 1


def closed_lane_indices(levels: tuple[int, int]) -> tuple[int, ...]:
    """H-lane indices containing q_1/q_0 on the fixed physical pair (0,1)."""
    ratio = F(levels[1], levels[0])
    require(1 / CAP <= ratio <= CAP, (levels, ratio))
    return tuple(
        index for index, (left, right, _, _) in enumerate(H_LANES, start=1)
        if left <= ratio <= right
    )


def recheck_unaffected_predecessor() -> dict[str, object]:
    """Rebuild only the predecessor pieces that survive the correction.

    In particular this routine never imports the predecessor's conclusion.
    It reruns the finite CSP and the proof obligations downstream from its
    valid branches, while deleting all four low/high rows before tail use.
    """
    bodies = tuple(T.residual_bodies())
    require(len(bodies) == 561 and len(set(bodies)) == 561, len(bodies))

    # P=Q is intentionally absent from every primitive/head channel bank.
    # Equal packet levels are already closed by the independently pinned
    # universal repeated-level theorem; neither of its two exceptions is in
    # this 561-body residual universe.
    repeated_exceptions = tuple(row[0] for row in T.C2.UNIVERSAL.EXPECTED_EXCEPTIONS)
    repeated_intersection = tuple(sorted(set(bodies) & set(repeated_exceptions)))
    require(len(repeated_exceptions) == 2 and not repeated_intersection,
            (repeated_exceptions, repeated_intersection))

    threshold_rows = []
    for body in bodies:
        ruler = 14 * lcm(*body)
        for pair in combinations(range(6), 2):
            threshold = M.coarse_threshold(body, ruler, pair)
            if threshold < T.PHASE_LIMIT:
                threshold_rows.append((threshold, body, pair))
    require(len(threshold_rows) == 6358, len(threshold_rows))
    maximum_bound = max(
        (T.product_bound(threshold), body, pair, threshold)
        for threshold, body, pair in threshold_rows
    )
    require(maximum_bound ==
            (F(483806351847017447820, 8426333817800417),
             (2, 4, 5, 6, 9, 12), (2, 5),
             F(5758395664776236224, 282220371910760177895)), maximum_bound)
    primitive_rows = M.THREE.ONE.C4.primitive_universe(maximum_bound[0], CAP)
    require(len(primitive_rows) == 15995, len(primitive_rows))
    require(all(P != Q for _, _, _, P, Q in primitive_rows),
            "P=Q leaked into primitive bank")

    closed = []
    traps = []
    constraints_by_body = {}
    verdict_digest = hashlib.sha256()
    for body in bodies:
        constraints, _ = M.body_constraints(body, 14 * lcm(*body), primitive_rows)
        witness, reason = M.solve_constraints(body, constraints)
        reverse, reverse_reason = M.solve_constraints(body, constraints, reverse=True)
        require((witness is None) == (reverse is None),
                ("CSP search-order mismatch", body, reason, reverse_reason))
        verdict_digest.update(repr((body, witness is not None,
                                    reason, reverse_reason)).encode())
        if witness is None:
            closed.append(body)
        else:
            traps.append(body)
            constraints_by_body[body] = constraints
    trap_digest = hashlib.sha256(repr(tuple(traps)).encode()).hexdigest()
    require((len(closed), len(traps), trap_digest, verdict_digest.hexdigest()) ==
            (389, 172, EXPECTED_TRAP_DIGEST, EXPECTED_VERDICT_DIGEST),
            (len(closed), len(traps), trap_digest, verdict_digest.hexdigest()))

    old_tail_bodies = {row[0] for row in M.TAIL_POLICIES}
    split_rows = tuple(row for row in M.TAIL_POLICIES if row[1] != "all")
    whole_rows = tuple(row for row in M.TAIL_POLICIES if row[1] == "all")
    require(len(old_tail_bodies) == 12 and len(whole_rows) == 10
            and len({row[0] for row in whole_rows}) == 10, (whole_rows, split_rows))
    require(split_rows == (
        (H, "low", (1, 0), 10, F(2)),
        (H, "high", (0, 1), 10, F(2)),
        (H2, "low", (2, 1), 6, F(14, 5)),
        (H2, "high", (0, 1), 8, F(4)),
    ), split_rows)
    require({row[0] for row in split_rows} == {H, H2}, split_rows)
    require(all(body in constraints_by_body and not constraints_by_body[body]
                for body in old_tail_bodies), old_tail_bodies)

    located = []
    policies = []
    located_monotonicity = []
    for body in traps:
        if body in old_tail_bodies:
            continue
        constraints = constraints_by_body[body]
        require(constraints, ("unexpected empty located graph", body))
        pair, allowed = min(constraints.items(), key=lambda row: (len(row[1]), row[0]))
        require(all(ratio > 1 for ratio in allowed), (body, pair, allowed))
        channels = tuple(
            channel for ratio in sorted(allowed)
            for channel in ((ratio.denominator, ratio.numerator),
                            (ratio.numerator, ratio.denominator))
        )
        require(len(channels) == 2 * len(allowed)
                and all(P != Q for P, Q in channels), (body, pair, channels))
        rows = tuple(M.located_best(body, pair, channel) for channel in channels)
        require(all(row[0] > 0 for row in rows), ("located policy", body, pair, rows))
        located.extend(rows)
        located_monotonicity.extend(
            scale_monotonicity(body, pair, channel, row[4], coupled=False)
            for channel, row in zip(channels, rows)
        )
        policies.append((body, pair, len(rows), min(rows)))
    located_digest = hashlib.sha256(repr(tuple(located)).encode()).hexdigest()
    require((len(policies), len(located), located_digest) ==
            (160, 1566, EXPECTED_LOCATED_DIGEST),
            (len(policies), len(located), located_digest))

    whole_tail_rows = []
    whole_tail_heads = []
    whole_tail_monotonicity = []
    for body, side, pair, start, numerator in whole_rows:
        require(side == "all" and not constraints_by_body[body],
                (body, side, constraints_by_body[body]))
        a, b = body[pair[0]], body[pair[1]]
        require(numerator == 2 * max(abs(F(2, 5) * a - b),
                                     abs(F(5, 2) * a - b)),
                (body, pair, numerator))
        margin = T.tail_envelope(body, pair, start, numerator)
        previous = T.tail_envelope(body, pair, start - 1, numerator)
        step = T.tail_step_gain(body, pair, start, numerator)
        tail_eta = numerator / 2 / (14 * lcm(*body) - F(a, start))
        require(margin > 0 >= previous and step > 0 and tail_eta < 1,
                ("whole tail", body, pair, start, previous, margin, step, tail_eta))
        whole_tail_monotonicity.append(
            analytic_tail_monotonicity(body, pair, start, numerator, coupled=False)
        )
        channels = M.oriented_channels_below("all", start)
        require(channels == M.oriented_channels_below("all", start, oversize=5)
                and all(P != Q for P, Q in channels),
                ("whole head exhaustion", body, len(channels)))
        rows = tuple(M.located_best(body, pair, channel) for channel in channels)
        require(all(row[0] > 0 for row in rows), ("whole head", body, min(rows)))
        located_monotonicity.extend(
            scale_monotonicity(body, pair, channel, row[4], coupled=False)
            for channel, row in zip(channels, rows)
        )
        whole_tail_heads.extend(rows)
        whole_tail_rows.append((body, pair, start, numerator, margin, previous,
                                step, tail_eta, len(channels), min(rows)))
    require(len(whole_tail_rows) == 10 and len(whole_tail_heads) == 300,
            (len(whole_tail_rows), len(whole_tail_heads)))

    return {
        "bodies": bodies,
        "repeated_exceptions": repeated_exceptions,
        "threshold_count": len(threshold_rows),
        "maximum_bound": maximum_bound,
        "primitive_count": len(primitive_rows),
        "closed": tuple(closed),
        "traps": tuple(traps),
        "trap_digest": trap_digest,
        "verdict_digest": verdict_digest.hexdigest(),
        "policies": tuple(policies),
        "located": tuple(located),
        "located_digest": located_digest,
        "located_monotonicity_digest": hashlib.sha256(
            repr(tuple(located_monotonicity)).encode()).hexdigest(),
        "whole_rows": whole_rows,
        "split_rows": split_rows,
        "whole_tail_rows": tuple(whole_tail_rows),
        "whole_tail_heads": tuple(whole_tail_heads),
        "whole_tail_head_digest": hashlib.sha256(
            repr(tuple(whole_tail_heads)).encode()).hexdigest(),
        "whole_tail_monotonicity": tuple(whole_tail_monotonicity),
    }


def multi_scale_head(body, pair, channel):
    """Direct finite prefix followed by a uniform transported tail in g."""
    ruler, safe_ranges = T.R.safe_cell_ranges(body)
    a, b = body[pair[0]], body[pair[1]]
    P, Q = channel
    reduced_min = min(P, Q)
    g0 = ceil(F(MIN_LEVEL, reduced_min))

    skeleton_rows = []
    for left, right in safe_ranges:
        for cell in range(left, right):
            skeleton = T.intersection_mass(
                T.primitive_arcs(P, F(a * cell, ruler) % 1),
                T.primitive_arcs(Q, F(b * cell, ruler) % 1),
            )
            skeleton_rows.append((skeleton, cell))
    skeleton, skeleton_cell = max(skeleton_rows)

    tail = None
    for g in range(g0, 1001):
        # The pair maximum is an actual packet level, so q_max<=CAP*m gives
        # the sharper coupled bound m>=ceil(g*max(P,Q)/CAP).
        debt_level = max(MIN_LEVEL, ceil(F(g * max(P, Q), CAP)))
        debt = singleton_debt(body, ruler, debt_level)
        eta = F(g * (Q * a - P * b), P * g * ruler - a)
        bracket = skeleton - 2 * abs(eta)
        if bracket <= debt:
            continue
        c = 1 - F(a, P * g * ruler)
        transported = bracket / c
        actual = T.intersection_mass(
            T.R.reflected_level_arcs(ruler, a, P * g, skeleton_cell),
            T.R.reflected_level_arcs(ruler, b, Q * g, skeleton_cell),
        )
        require(abs(eta) <= ETA_GUARD < 1 and 0 < c < 1,
                ("tail homotopy gate", body, pair, channel, g, eta, c))
        require(actual >= transported > bracket > debt,
                ("tail transport", body, pair, channel, g, actual,
                 transported, bracket, debt))
        # The direct check at this first tail scale is the positive control for
        # the transport lemma.  The returned algebraic certificate proves that
        # its hypotheses and the weaker bracket>debt inequality persist for
        # every later g; no finite sampling is used for that implication.
        monotonicity = scale_monotonicity(body, pair, channel, g, coupled=True)
        tail = (bracket - debt, g, debt_level, skeleton_cell,
                skeleton, eta, debt, actual, transported, monotonicity)
        break
    require(tail is not None, ("no common-scale tail", body, pair, channel))

    direct_rows = []
    for g in range(g0, tail[1]):
        debt_level = max(MIN_LEVEL, ceil(F(g * max(P, Q), CAP)))
        debt = singleton_debt(body, ruler, debt_level)
        best = max(
            (T.intersection_mass(
                T.R.reflected_level_arcs(ruler, a, P * g, cell),
                T.R.reflected_level_arcs(ruler, b, Q * g, cell),
             ) - debt, cell)
            for left, right in safe_ranges for cell in range(left, right)
        )
        require(best[0] > 0,
                ("finite direct scale", body, pair, channel, g, debt_level, best))
        direct_rows.append((g, debt_level, best))
    return (tail[0], body, pair, channel, g0, tail, tuple(direct_rows))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    surviving = recheck_unaffected_predecessor()

    # Diagnose the old quantifier failure explicitly.  Both H lanes encode
    # q_0<=q_1; H2's simultaneous failure order q_2<q_1<q_0 is consistent.
    old_by_body = {}
    for body, side, pair, start, numerator in M.TAIL_POLICIES:
        old_by_body.setdefault(body, []).append((side, pair, start, numerator))
    require(len(old_by_body) == 12, old_by_body)
    h_old = tuple(old_by_body[H])
    h2_old = tuple(old_by_body[H2])
    require(tuple(physical_predicate(side, pair) for side, pair, _, _ in h_old)
            == ((0, 1), (0, 1)), h_old)
    require(tuple(physical_predicate(side, pair) for side, pair, _, _ in h2_old)
            == ((1, 2), (0, 1)), h2_old)
    require(all(any(row[0] == "all" for row in rows)
                for body, rows in old_by_body.items() if body not in {H, H2}), old_by_body)

    # These are admissible physical assignments, not renamed channels.  The
    # first witnesses the missed reverse order for H.  The second lies on the
    # cap boundary and simultaneously falsifies H2's two acyclic predicates.
    h_hostile = (2, 1)
    h2_hostile = (5, 4, 2)
    h_hostile_truth = tuple(old_lane_applies(side, pair, h_hostile)
                            for side, pair, _, _ in h_old)
    h2_hostile_truth = tuple(old_lane_applies(side, pair, h2_hostile)
                             for side, pair, _, _ in h2_old)
    require(h_hostile_truth == (False, False), (h_hostile, h_hostile_truth))
    require(h2_hostile_truth == (False, False), (h2_hostile, h2_hostile_truth))

    # Exact interval audit.  Consecutive closed intervals meet without gaps;
    # a shared endpoint belongs to both adjacent lanes, which is harmless.
    h_intervals = tuple((left, right) for left, right, _, _ in H_LANES)
    require(h_intervals[0][0] == 1 / CAP and h_intervals[-1][1] == CAP,
            h_intervals)
    require(all(left < right for left, right in h_intervals), h_intervals)
    require(all(h_intervals[index][1] == h_intervals[index + 1][0]
                for index in range(len(h_intervals) - 1)), h_intervals)
    audited_endpoints = ((h_intervals[0][0],)
                         + tuple(right for _, right in h_intervals))
    endpoint_memberships = tuple(
        (endpoint, tuple(index for index, (left, right) in enumerate(h_intervals, start=1)
                         if left <= endpoint <= right))
        for endpoint in audited_endpoints
    )
    require(tuple(len(indices) for _, indices in endpoint_memberships)
            == (1, 2, 2, 2, 2, 1), endpoint_memberships)
    assignment_controls = tuple(
        (levels, F(levels[1], levels[0]), closed_lane_indices(levels))
        for levels in ((2, 1), (1, 2))
    )
    require(all(indices for _, _, indices in assignment_controls), assignment_controls)

    lanes = []
    lane_monotonicity = []
    for left, right, start, numerator in H_LANES:
        a, b = H[0], H[1]
        require(2 * max(abs(left * a - b), abs(right * a - b)) == numerator,
                ("H endpoint numerator", left, right, numerator))
        margin, debt_level, debt = coupled_tail_envelope(H, (0, 1), start, numerator)
        previous = coupled_tail_envelope(H, (0, 1), start - 1, numerator)[0]
        next_margin = coupled_tail_envelope(H, (0, 1), start + 1, numerator)[0]
        require(margin > 0 >= previous and next_margin > margin,
                ("H tail threshold", left, right, start, previous, margin, next_margin))
        lane_monotonicity.append(
            analytic_tail_monotonicity(H, (0, 1), start, numerator, coupled=True)
        )
        channels = interval_channels(left, right, start)
        require(channels == interval_channels(left, right, start, oversize=5),
                ("H head exhaustion", left, right, start))
        lanes.append((H, left, right, (0, 1), start, numerator,
                      margin, previous, debt_level, debt, channels))

    left, right, start, numerator = H2_LANE
    a, b = H2[0], H2[1]
    require(2 * max(abs(left * a - b), abs(right * a - b)) == numerator,
            ("H2 endpoint numerator", numerator))
    margin = T.tail_envelope(H2, (0, 1), start, numerator)
    previous = T.tail_envelope(H2, (0, 1), start - 1, numerator)
    require(margin > 0 >= previous,
            ("H2 tail threshold", previous, margin))
    lane_monotonicity.append(
        analytic_tail_monotonicity(H2, (0, 1), start, numerator, coupled=False)
    )
    channels = interval_channels(left, right, start)
    require(channels == interval_channels(left, right, start, oversize=5),
            "H2 head exhaustion")
    lanes.append((H2, left, right, (0, 1), start, numerator,
                  margin, previous, MIN_LEVEL,
                  singleton_debt(H2, 14 * lcm(*H2), MIN_LEVEL), channels))

    require(tuple(row[1] for row in lanes[:5]) ==
            (F(2, 5), F(3, 7), F(1, 2), F(2, 3), F(11, 9)), lanes)
    require(tuple(row[2] for row in lanes[:5]) ==
            (F(3, 7), F(1, 2), F(2, 3), F(11, 9), F(5, 2)), lanes)
    require(all(row[3] == (0, 1) for row in lanes), lanes)
    require(all(row[8] == max(MIN_LEVEL, ceil(F(row[4], CAP)))
                and row[9] == singleton_debt(row[0], 14 * lcm(*row[0]), row[8])
                for row in lanes[:5]), lanes[:5])
    require(lanes[5][8] == MIN_LEVEL
            and lanes[5][9] == singleton_debt(H2, 14 * lcm(*H2), MIN_LEVEL),
            lanes[5])
    require(sum(len(row[10]) for row in lanes) == 605,
            tuple(len(row[10]) for row in lanes))

    head_rows = []
    for body, _, _, pair, _, _, _, _, _, _, channels in lanes:
        head_rows.extend(multi_scale_head(body, pair, channel) for channel in channels)
    require(len(head_rows) == 605, len(head_rows))
    require(max(len(row[6]) for row in head_rows) == 4,
            max(len(row[6]) for row in head_rows))
    for _, body, pair, channel, _, tail, direct_rows in head_rows:
        P, Q = channel
        tail_level = max(MIN_LEVEL, ceil(F(tail[1] * max(P, Q), CAP)))
        require(tail[2] == tail_level
                and tail[6] == singleton_debt(body, 14 * lcm(*body), tail_level),
                (body, pair, channel, tail))
        require(all(debt_level == max(MIN_LEVEL, ceil(F(g * max(P, Q), CAP)))
                    for g, debt_level, _ in direct_rows),
                (body, pair, channel, direct_rows))
    head_digest = hashlib.sha256(repr(tuple(head_rows)).encode()).hexdigest()
    require(head_digest == EXPECTED_HEAD_DIGEST, head_digest)

    lane_rows = tuple(
        (body, left, right, pair, start, numerator, margin, previous,
         debt_level, debt, len(channels))
        for body, left, right, pair, start, numerator, margin, previous,
        debt_level, debt, channels in lanes
    )
    lane_digest = hashlib.sha256(repr(lane_rows).encode()).hexdigest()
    surviving_image = (
        surviving["bodies"], surviving["repeated_exceptions"],
        surviving["threshold_count"], surviving["maximum_bound"],
        surviving["primitive_count"], surviving["closed"], surviving["traps"],
        surviving["trap_digest"], surviving["verdict_digest"],
        surviving["policies"], surviving["located"], surviving["located_digest"],
        surviving["located_monotonicity_digest"], surviving["whole_rows"],
        surviving["split_rows"], surviving["whole_tail_rows"],
        surviving["whole_tail_heads"], surviving["whole_tail_head_digest"],
        surviving["whole_tail_monotonicity"],
    )
    semantic_image = (
        EXPECTED_OLD_SHA256, EXPECTED_OLD_OUTPUT_SHA256, EXPECTED_OLD_SEMANTIC,
        surviving_image,
        tuple((body, tuple(rows)) for body, rows in old_by_body.items()),
        h_old, h2_old, h_hostile, h_hostile_truth, h2_hostile, h2_hostile_truth,
        h_intervals, endpoint_memberships, assignment_controls,
        lane_rows, tuple(lane_monotonicity), lane_digest, tuple(head_rows), head_digest,
    )
    semantic = hashlib.sha256(repr(semantic_image).encode()).hexdigest()
    require(semantic == EXPECTED_SEMANTIC,
            ("semantic image changed", semantic, EXPECTED_SEMANTIC))

    weakest = min(head_rows)
    longest = max(head_rows, key=lambda row: len(row[6]))
    lines = [
        "LRC14 reflected cap-5/2 tail orientation repair",
        "predecessor_conclusion=RETRACTED_AND_NOT_INHERITED;only pinned code/data are reused and every surviving downstream obligation is recomputed",
        f"equal_level_gate=residual_bodies:{len(surviving['bodies'])};P=Q excluded only after universal repeated-level recheck;exception_bodies:{surviving['repeated_exceptions']};intersection:EMPTY",
        f"unaffected_csp=closed:{len(surviving['closed'])};traps:{len(surviving['traps'])};reverse_order_agrees;trap_digest:{surviving['trap_digest']};verdict_digest:{surviving['verdict_digest']}",
        f"unaffected_located=policies:{len(surviving['policies'])};both_orientation_controls:{len(surviving['located'])};digest:{surviving['located_digest']};monotonicity_digest:{surviving['located_monotonicity_digest']}",
        f"unaffected_whole_tails=bodies:{len(surviving['whole_rows'])};head_controls:{len(surviving['whole_tail_heads'])};head_digest:{surviving['whole_tail_head_digest']};discarded_split_bodies={(H,H2)}",
        "diagnosis=old H low(1,0) and high(0,1) both encode q_0<=q_1;old H2 simultaneous failure q_2<q_1<q_0 is consistent",
        f"hostile_old_policy_controls=H:q={h_hostile},truth={h_hostile_truth};H2:q={h2_hostile},truth={h2_hostile_truth}",
        "affected_tail_bodies=2;unaffected_whole-interval_tail_bodies=10;old constrained policies and CSP remain pinned exact data",
        "repair=H uses five closed Farey lanes on one fixed ordered pair (0,1);H2 uses one whole interval on fixed pair (0,1)",
        f"interval_endpoint_audit={endpoint_memberships};assignment_controls={assignment_controls};closed_endpoint_overlaps=HARMLESS",
        "coupled_debt=q_max<=5m/2 and pair_max=g*max(P,Q) imply m>=ceil(2g*max(P,Q)/5)",
        "singleton_debt_levels=H analytic tails use max(4,ceil(2*start/5));H2 whole tail retains the safe base m>=4;every direct and transported head scale uses max(4,ceil(2g*max(P,Q)/5))",
        "infinite_tail_gate=for every later scale,primitive skeleton is fixed;phase loss is strictly decreasing unless mismatch=0;transport displacement a/(PgL-a) strictly decreases and its favourable inverse factor is discarded;ceiling level is nondecreasing and singleton debt is nonincreasing",
        *(f"lane_{index+1}=body:{row[0]};interval:[{qtext(row[1])},{qtext(row[2])}];pair:{row[3]};start:{row[4]};numerator:{qtext(row[5])};margin:{qtext(row[6])};previous:{qtext(row[7])};debt_level:{row[8]};heads:{row[10]}" for index, row in enumerate(lane_rows)),
        f"replacement_lanes={len(lane_rows)};head_channels={len(head_rows)};weakest_tail_margin={qtext(weakest[0])}@body={weakest[1]},channel={weakest[3]},g_tail={weakest[5][1]};maximum_direct_prefix={len(longest[6])}@body={longest[1]},channel={longest[3]},g0={longest[4]},g_tail={longest[5][1]}",
        f"lane_digest={lane_digest};head_digest={head_digest}",
        "coverage=the five H intervals cover [2/5,5/2] on the same physical pair;the H2 lane covers [2/5,5/2];no reversed-pair disjunction remains",
        "conclusion=pinned cap-5/2 CSP+located data, ten valid old whole-interval tails, and the repaired H/H2 tails close D>=6,3m>=2D",
        "corollary=reflected sufficient-certificate failure is again confined to D>=6,1<=m<2D/3",
        "scope=reflected THM-2941 sufficient family only;physical LRC14 remains open",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"old_source_sha256={EXPECTED_OLD_SHA256}",
        f"old_output_sha256={EXPECTED_OLD_OUTPUT_SHA256}",
        f"old_semantic_sha256={EXPECTED_OLD_SEMANTIC}",
        f"source_sha256={sha256(HERE)}",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ]
    text = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
