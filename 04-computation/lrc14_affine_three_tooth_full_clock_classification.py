#!/usr/bin/env python3
"""Exact affine-handoff classification on the THM-2684 rail envelope.

Put

    D(x)={13x},                 T_(s,beta)(x)={s*13x+beta},

with ``s in {+1,-1}``, and let ``B`` be the exact full THM-2584 rail
envelope

    [0,1/28) union [13/28,15/28) union [27/28,1).

An inherited rail atom at state ``u`` has intrinsic stored edge

    c_7(Du) -> c_7(D^2u).

Consequently the exact three-event envelope/clock support is

    K intersect T^(-1)K intersect T^(-2)K,

where ``K={u in B:c_7(Du)!=c_7(D^2u)}``.  This script proves, by exact
rational polyhedral projection in the ``(x,beta)`` plane, that this support
has positive length for both signs exactly when

    beta notin {0,1/2} mod 1.

The proof is constructive.  Each choice of three components of ``K`` and
the two integer wrap branches gives six affine inequalities in ``x,beta``.
Fourier--Motzkin elimination of ``x`` gives one exact open beta corridor.
A deterministic largest-gap procedure selects enough such corridors to
cover ``(0,1/2) union (1/2,1)``; every apparent corridor seam is then healed
by a fresh positive witness, while beta 0 and 1/2 are checked directly empty.
Thus the parameter conclusion is a finite exact cover, not a grid sample.

The script then keeps the typing distinctions which the envelope alone
forgets.

* If ``y=T_(s,beta)x``, then the following shallow phase is
  ``s D^2x+13 beta``.  This descends to a permutation of the seven clock
  labels iff ``13 beta=k/7 mod 1``.  The inverse label transport is
  ``a -> s(a-k)``.
* At beta=1/13 that transport is identity for ``s=+1`` and reflection for
  ``s=-1``.  It also gives the inherited digit covariance
  ``j(Tx)=h(x)`` or ``12-h(x)``.  Exact THM-2584 reconstruction supplies an
  explicit positive rail-labelled triple for each sign, with all three
  intrinsic edges nonconstant.
* At beta=1/14 one has ``T_(+,beta)^2=D^2`` and positive intrinsic support,
  but ``13 beta=13/14`` is a half-clock-cell shift, not a label permutation.
  A phase-side sidecar is mandatory before calling this an inherited typed
  handoff.
* Every nontrivial THM-2657 physical lift ``beta=k/13^6`` has
  ``k=7 delta mod 13`` for a nonzero ``delta``, hence is a 13-adic unit.  It
  therefore avoids both exceptional beta values and has positive intrinsic
  support for both signs.  On the other hand global clock covariance would
  require ``7k/13^5`` integral, hence ``13^5`` divides ``k``, impossible for
  a 13-adic unit.  The script exhausts all residues modulo ``13^6`` as a
  finite control for this universal support-versus-typing tradeoff.
* The fixed odometer twist beta=14/13^6 has a positive endpoint/mixed
  intrinsic corridor even though the central orbit-clock scout sees only one
  flip.  Conversely the alternating minimal lifts -14,+14 have a cosmetic
  central state-clock cycle whose intrinsic stored edges are diagonal.
* The alternating lawful lifts +/-(13^5+1) instead have the central
  intrinsic edges 4->3 and 3->4 and exact owner-to-next-shallow gluing, as in
  the independently frozen odometer scout.
* Two selected cycle pairs from the independent THM-2640 coefficient atlas
  are rebuilt from their weighted physical rails and satisfy the exact unit
  relations ``Y_minus(ell)=10 Y_plus(-ell)`` and
  ``Y_minus(ell)=4 Y_plus(-ell) mod 13``.  These are frozen only as
  coefficient-level signals and next affine compatibility tests: they do not
  identify or transport the physical packets supporting the rows.

Positive envelope support distributes over the finite exact rail union, so
it guarantees some raw rail-labelled product.  It does not insert present
factors, delayed words, private units, guard semantics, or a row exclusion.
"""

from bisect import bisect_right
from fractions import Fraction
from math import ceil, floor
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_alternate_arrival_physical_rail_handoff as rail
import lrc14_predecessor_carry_private_root_atlas_thm2640 as carry_atlas


F = Fraction
P = 13
Q = 7
R = P**6
S = P**5
TOOTH_NAMES = (0, 6, 12)
TEETH = (
    (F(0), F(1, 28)),
    (F(13, 28), F(15, 28)),
    (F(27, 28), F(1)),
)
CLOCK_BOUNDARIES = tuple(F(2 * ell + 1, 2 * Q) for ell in range(Q))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def floor_fraction(value):
    return value.numerator // value.denominator


def ceil_fraction(value):
    return -floor_fraction(-value)


def clock(value):
    value %= 1
    shifted = Q * value + F(1, 2)
    return floor_fraction(shifted) % Q


def shallow(value):
    return clock(P * value)


def owner(value):
    return clock(P * P * value)


def digit_j(value):
    return floor_fraction(R * (value % 1)) % P


def digit_h(value):
    return floor_fraction(P * ((R * value) % 1))


def carry(value):
    return floor_fraction(R * (value % 1)) % P


def future_digit(value):
    return floor_fraction(P * ((R * value) % 1))


def tooth(value):
    value %= 1
    for index, (left, right) in enumerate(TEETH):
        if left <= value < right:
            return TOOTH_NAMES[index]
    return None


def merge_sorted(intervals, touching=True):
    """Merge an already sorted interval stream.

    ``touching=False`` treats open intervals which merely share an endpoint
    as different components.  That is essential in the beta seam audit.
    """
    out = []
    for left, right in intervals:
        if left >= right:
            continue
        joins = out and (left <= out[-1][1] if touching
                         else left < out[-1][1])
        if joins:
            out[-1] = (out[-1][0], max(out[-1][1], right))
        else:
            out.append((left, right))
    return out


def merge(intervals, touching=True):
    return merge_sorted(sorted(intervals), touching=touching)


def intersect(left, right):
    out = []
    i = j = 0
    while i < len(left) and j < len(right):
        lo = max(left[i][0], right[j][0])
        hi = min(left[i][1], right[j][1])
        if lo < hi:
            out.append((lo, hi))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return out


def measure(intervals):
    return sum((right - left for left, right in intervals), F(0))


def reflect_support(intervals):
    # Canonical half-open reflection; null endpoint assignments are irrelevant
    # to every positive-length statement in this script.
    return merge((1 - right, 1 - left) for left, right in intervals)


def pullback_positive(intervals, slope, shift):
    """Positive-length pullback under x -> {slope*x+shift}."""
    out = []
    first = floor_fraction(shift) - 1
    last = ceil_fraction(shift) + slope + 1
    # Branch-major order and sorted target intervals produce sorted output.
    for branch in range(first, last + 1):
        for left, right in intervals:
            lo = max(F(0), (branch + left - shift) / slope)
            hi = min(F(1), (branch + right - shift) / slope)
            if lo < hi:
                out.append((lo, hi))
    return merge_sorted(out)


def pullback(intervals, slope, shift):
    if slope > 0:
        return pullback_positive(intervals, slope, shift)
    # {-a*x+shift} in I iff {a*x-shift} in rho(I).
    return pullback_positive(reflect_support(intervals), -slope, -shift)


def build_intrinsic_legal_set():
    """Partition B exactly where its intrinsic stored edge is nonconstant."""
    out = []
    for tooth_index, (left, right) in enumerate(TEETH):
        cuts = {left, right}
        for slope in (P, P * P):
            for boundary in CLOCK_BOUNDARIES:
                first = floor_fraction(slope * left - boundary) - 1
                last = ceil_fraction(slope * right - boundary) + 1
                for branch in range(first, last + 1):
                    point = (branch + boundary) / slope
                    if left < point < right:
                        cuts.add(point)
        cuts = sorted(cuts)
        for lo, hi in zip(cuts, cuts[1:]):
            midpoint = (lo + hi) / 2
            if shallow(midpoint) != owner(midpoint):
                out.append((lo, hi, tooth_index,
                            (shallow(midpoint), owner(midpoint))))
    return tuple(out)


LEGAL_DATA = build_intrinsic_legal_set()
K = tuple((left, right) for left, right, _, _ in LEGAL_DATA)
K_STARTS = tuple(left for left, _ in K)


def locate_k(value):
    index = bisect_right(K_STARTS, value) - 1
    if index >= 0 and K[index][0] < value < K[index][1]:
        return K[index]
    return None


def affine_states(x, sign, beta):
    return (
        x % 1,
        (sign * P * x + beta) % 1,
        (P * P * x + (sign * P + 1) * beta) % 1,
    )


CHAIN_CACHE = {}


def chain_components(beta, sign, carrier=K):
    """States x,T x,T^2 x all lie in the supplied carrier."""
    if carrier is K:
        key = (sign, beta)
        if key in CHAIN_CACHE:
            return CHAIN_CACHE[key]
    first = pullback(carrier, sign * P, beta)
    second = pullback(carrier, P * P, (sign * P + 1) * beta)
    result = tuple(intersect(intersect(carrier, first), second))
    if carrier is K:
        CHAIN_CACHE[sign, beta] = result
    return result


def find_corridor_witness(beta, sign):
    """Return three K components and wrap branches at a positive beta."""
    for left, right in chain_components(beta, sign):
        for fraction in (F(1, 3), F(2, 5), F(1, 2), F(3, 5)):
            x = left + fraction * (right - left)
            x0, x1, x2 = affine_states(x, sign, beta)
            j0, j1, j2 = locate_k(x0), locate_k(x1), locate_k(x2)
            if not all((j0, j1, j2)):
                continue
            first_branch = sign * P * x + beta - x1
            second_branch = P * P * x + (sign * P + 1) * beta - x2
            require(first_branch.denominator == 1
                    and second_branch.denominator == 1,
                    "affine wrap branch stopped being integral")
            return (j0, j1, j2, int(first_branch),
                    int(second_branch), x)
    return None


def corridor_projection(witness, sign):
    """Fourier--Motzkin projection of one exact (x,beta) corridor.

    Every lower x-bound has the form a*beta+b and every upper bound has the
    same form.  Positive x-width is equivalent to all lower bounds being
    strictly below all upper bounds; eliminating x therefore leaves one open
    rational beta interval.
    """
    j0, j1, j2, first_branch, second_branch, _ = witness
    coefficient = sign * P + 1
    lowers = [(F(0), j0[0])]
    uppers = [(F(0), j0[1])]
    if sign == 1:
        lowers.append((-F(1, P), (first_branch + j1[0]) / P))
        uppers.append((-F(1, P), (first_branch + j1[1]) / P))
    else:
        lowers.append((F(1, P), -(first_branch + j1[1]) / P))
        uppers.append((F(1, P), -(first_branch + j1[0]) / P))
    lowers.append((-F(coefficient, P * P),
                   (second_branch + j2[0]) / (P * P)))
    uppers.append((-F(coefficient, P * P),
                   (second_branch + j2[1]) / (P * P)))

    beta_lo = F(0)
    beta_hi = F(1)
    for lower_slope, lower_offset in lowers:
        for upper_slope, upper_offset in uppers:
            slope = lower_slope - upper_slope
            offset = upper_offset - lower_offset
            if slope == 0:
                require(offset > 0, "corridor has incompatible parallel walls")
            elif slope > 0:
                beta_hi = min(beta_hi, offset / slope)
            else:
                beta_lo = max(beta_lo, offset / slope)
    beta_lo = max(F(0), beta_lo)
    beta_hi = min(F(1), beta_hi)
    require(beta_lo < beta_hi,
            "positive affine witness projected to an empty beta corridor")
    return beta_lo, beta_hi


def subtract_interval(gaps, interval):
    left, right = interval
    out = []
    for lo, hi in gaps:
        if right <= lo or hi <= left:
            out.append((lo, hi))
            continue
        if lo < left:
            out.append((lo, left))
        if right < hi:
            out.append((right, hi))
    return [(lo, hi) for lo, hi in out if lo < hi]


def strict_union(intervals):
    return tuple(merge(intervals, touching=False))


def cover_beta_parameter(sign):
    """Construct and verify a finite exact corridor cover."""
    target = ((F(0), F(1, 2)), (F(1, 2), F(1)))
    gaps = list(target)
    corridors = []
    while gaps:
        index = max(range(len(gaps)),
                    key=lambda i: gaps[i][1] - gaps[i][0])
        lo, hi = gaps[index]
        beta = (lo + hi) / 2
        witness = find_corridor_witness(beta, sign)
        require(witness is not None,
                f"parameter cover found an open gap at sign={sign}, beta={beta}")
        corridor = corridor_projection(witness, sign)
        require(corridor[0] < beta < corridor[1],
                "chosen beta left its projected corridor")
        corridors.append(corridor)
        gaps = subtract_interval(gaps, corridor)
    greedy_count = len(corridors)

    # The open greedy intervals may merely touch at finitely many seams.
    # Positive support at such a seam persists to an open beta neighbourhood;
    # explicitly derive that additional corridor until only 0 and 1/2 remain.
    seam_repairs = 0
    while True:
        union = strict_union(corridors)
        seams = tuple(
            left[1]
            for left, right in zip(union, union[1:])
            if left[1] == right[0]
            and left[1] not in (F(0), F(1, 2), F(1))
        )
        if not seams:
            break
        for beta in seams:
            witness = find_corridor_witness(beta, sign)
            require(witness is not None,
                    "an interior beta seam had no three-event support")
            corridor = corridor_projection(witness, sign)
            require(corridor[0] < beta < corridor[1],
                    "seam witness did not heal a neighbourhood of the seam")
            corridors.append(corridor)
            seam_repairs += 1

    final_union = strict_union(corridors)
    require(final_union == target,
            f"affine beta cover changed at sign={sign}: {final_union}")
    require(not chain_components(F(0), sign)
            and not chain_components(F(1, 2), sign),
            "one of the two exceptional betas acquired legal support")
    return greedy_count, seam_repairs, len(corridors), final_union


def support_from_bank(bank, key):
    denominator = rail.old.T
    return tuple(merge(
        (F(left, denominator), F(right, denominator))
        for left, right, weight in bank[key]
        if weight
    ))


def labelled_joint(bank, sign, beta, keys):
    base, following, third = tuple(support_from_bank(bank, key)
                                   for key in keys)
    return tuple(intersect(
        intersect(
            intersect(chain_components(beta, sign), base),
            pullback(following, sign * P, beta),
        ),
        pullback(third, P * P, (sign * P + 1) * beta),
    ))


def component_signature(component, sign, beta):
    x = (component[0] + component[1]) / 2
    states = affine_states(x, sign, beta)
    return {
        "component": component,
        "length": component[1] - component[0],
        "arrivals": tuple(tooth(value) for value in states),
        "state_clocks": tuple(clock(value) for value in states),
        "stored_edges": tuple((shallow(value), owner(value))
                              for value in states),
    }


def raw_special(beta, sign):
    signatures = []
    total = F(0)
    for first, base_tooth in enumerate(TEETH):
        for second, next_tooth in enumerate(TEETH):
            for third, last_tooth in enumerate(TEETH):
                raw = intersect(
                    intersect(
                        (base_tooth,),
                        pullback((next_tooth,), sign * P, beta),
                    ),
                    pullback(
                        (last_tooth,), P * P,
                        (sign * P + 1) * beta,
                    ),
                )
                if not raw:
                    continue
                # At both exceptional parameters the first stored edge is
                # diagonal on every raw triple.  This is stronger than merely
                # knowing that at least one of the three events fails.
                require(not intersect(raw, K),
                        "an exceptional raw triple acquired a nonconstant first edge")
                mass = measure(raw)
                signatures.append((
                    (TOOTH_NAMES[first], TOOTH_NAMES[second],
                     TOOTH_NAMES[third]),
                    mass,
                ))
                total += mass
    return tuple(signatures), total


def rebuild_thm2640_unit_vector(assets, row_type):
    """Rebuild one selected normalized THM-2640 coefficient row exactly.

    ``row_type`` is ``(rail_meta,c,h,kappa,sector,edge,root)``.  In
    particular, the returned vector is not inferred from the claimed
    reflection: its seven entries are recomputed from weighted support,
    present-section intersections, the selected deep half-edge, delayed
    half-digit prefixes, and predecessor-carry descent.
    """
    module, prefixes, rails, present, starts, caches = assets
    meta, c, h, kappa, sector, edge, root = row_type
    matches = tuple(pieces for s, ell, theta, pieces in rails
                    if (s, ell, theta) == meta)
    require(len(matches) == 1,
            f"selected THM-2640 rail metadata is not unique: {meta}")
    pieces = matches[0]
    predicted_root = (2 * c + (2 * h + kappa) // P
                      + (edge == 0)) % P
    require(root == predicted_root != 0,
            "selected THM-2640 carry/half-edge root typing changed")
    values = []
    for ell5 in range(Q):
        overlap = carry_atlas.old.intersect_weighted_union(
            pieces, present[ell5, (-h) % P], starts[ell5, (-h) % P]
        )
        if edge == 0:
            lower, upper = 14 * root - 13, 14 * root
        else:
            lower, upper = 14 * root, 14 * root + 13
        half = carry_atlas.old.intersect_weighted_comb(
            overlap, module.C3, 182, lower, upper
        )
        pair = carry_atlas.delayed_carry_pair(
            half, prefixes[sector][ell5][h],
            caches.setdefault((sector, ell5, h), {}),
        )
        value = pair[c][kappa]
        require(value % carry_atlas.core.GLOBAL_CONTENT == 0,
                "selected THM-2640 entry lost global-content divisibility")
        values.append(
            (value // carry_atlas.core.GLOBAL_CONTENT)
            * pow(root, -1, P) % P
        )
    values = tuple(values)
    require(carry_atlas.is_unit(
        tuple(
            value * carry_atlas.core.GLOBAL_CONTENT * root
            for value in values
        ), root, carry_atlas.core.GLOBAL_CONTENT
    ), "selected normalized THM-2640 coefficient row stopped being a unit")
    return values


def main():
    require(len(K) == 148 and measure(K) == F(146, 1183),
            "intrinsic-legal three-tooth partition changed")
    require(all(shallow((left + right) / 2)
                != owner((left + right) / 2)
                for left, right in K),
            "K retained a diagonal stored edge")

    covers = {}
    for sign in (1, -1):
        covers[sign] = cover_beta_parameter(sign)

    # Exact finite/general audit for the physical lift scale.  For
    # beta=k/R, the label-permutation criterion becomes
    #
    #       Q*P*beta = 7k/13^5 in Z  iff  13^5 divides k.
    #
    # Exhaust the full residue universe as a control.  Every THM-2657
    # nonzero lift is in the complementary unit locus because
    # k=7*delta mod 13 with delta nonzero.  Since R is odd, no integral k/R
    # is the other exceptional parameter 1/2.
    covariant_lift_residues = tuple(
        k for k in range(R) if (Q * k) % S == 0
    )
    require(covariant_lift_residues == tuple(j * S for j in range(P)),
            "global seven-clock physical-lift residue atlas changed")
    thm2657_residue_map = tuple((delta, (Q * delta) % P)
                                for delta in range(1, P))
    require(tuple(sorted(residue for _, residue in thm2657_residue_map))
            == tuple(range(1, P)),
            "THM-2657 nonzero carry/root lift stopped being a 13-unit")
    covariant_set = set(covariant_lift_residues)
    unit_lift_count = 0
    tradeoff_violations = []
    for k in range(1, R):
        if k % P == 0:
            continue
        unit_lift_count += 1
        if k in covariant_set or 2 * k == R:
            tradeoff_violations.append(k)
    require(unit_lift_count == (P - 1) * S
            and not tradeoff_violations,
            "physical unit-lift support/covariance tradeoff changed")

    plus_zero, plus_zero_mass = raw_special(F(0), 1)
    plus_half, plus_half_mass = raw_special(F(1, 2), 1)
    minus_zero, minus_zero_mass = raw_special(F(0), -1)
    minus_half, minus_half_mass = raw_special(F(1, 2), -1)
    require(all(mass == F(1, 1183) for mass in (
        plus_zero_mass, plus_half_mass, minus_zero_mass, minus_half_mass
    )), "exceptional raw triple mass changed")
    # beta=1/13: strongest fixed typed representatives and actual rail keys.
    beta_typed = F(1, 13)
    require((P * beta_typed) % 1 == 0
            and (R * beta_typed).denominator == 1
            and int(R * beta_typed) % P == 0,
            "beta=1/13 lost clock/digit covariance")
    bank = rail.build_rail_bank()
    require(len(bank) == 324, "THM-2584 positive rail bank changed")
    plus_keys = ((12, 1, 3, 12), (0, 1, 2, 0), (6, 1, 3, 12))
    minus_keys = ((12, 1, 3, 12), (6, 1, 5, 0), (6, 1, 3, 12))
    plus_joint = labelled_joint(bank, 1, beta_typed, plus_keys)
    minus_joint = labelled_joint(bank, -1, beta_typed, minus_keys)
    require(plus_joint and minus_joint,
            "beta=1/13 lost a rail-labelled three-event product")
    plus_signature = component_signature(plus_joint[0], 1, beta_typed)
    minus_signature = component_signature(minus_joint[0], -1, beta_typed)
    require(all(left != right for left, right in plus_signature["stored_edges"]
                + minus_signature["stored_edges"]),
            "typed rail witness acquired a diagonal intrinsic edge")
    for signature, sign in ((plus_signature, 1), (minus_signature, -1)):
        midpoint = sum(signature["component"], F(0)) / 2
        _, following, third = affine_states(midpoint, sign, beta_typed)
        if sign == 1:
            require(owner(midpoint) == shallow(following)
                    and owner(following) == shallow(third),
                    "positive typed clock interfaces stopped gluing")
            require(digit_j(following) == digit_h(midpoint)
                    and digit_j(third) == digit_h(following),
                    "positive beta=1/13 digit covariance changed")
        else:
            require(owner(midpoint) == (-shallow(following)) % Q
                    and owner(following) == (-shallow(third)) % Q,
                    "negative reflected clock interfaces stopped gluing")
            require(digit_j(following)
                    == (P - 1 - digit_h(midpoint)) % P
                    and digit_j(third)
                    == (P - 1 - digit_h(following)) % P,
                    "negative beta=1/13 digit covariance changed")

    # beta=1/14: intrinsic support is real, but label-only covariance is not.
    beta_half_clock = F(1, 14)
    half_clock_components = chain_components(beta_half_clock, 1)
    require((P + 1) * beta_half_clock == 1
            and len(half_clock_components) == 78
            and measure(half_clock_components) == F(891, 399854),
            "beta=1/14 intrinsic corridor changed")
    require((Q * P * beta_half_clock).denominator == 2,
            "beta=1/14 unexpectedly became a seven-clock permutation")

    # A fixed minimal odometer twist: the central state-flip monotonicity does
    # not see this endpoint/mixed intrinsic corridor.
    beta_small = F(14, R)
    small_components = chain_components(beta_small, 1)
    require(len(small_components) == 2
            and measure(small_components) == F(28, 815730721),
            "fixed minimal odometer corridor changed")
    small_signature = component_signature(small_components[0], 1, beta_small)
    require(small_signature["arrivals"] == (0, 0, 6)
            and small_signature["stored_edges"]
            == ((0, 3), (3, 0), (0, 3)),
            "fixed odometer endpoint/mixed signature changed")

    # Alternating minimal lifts: orbit clocks alternate, but the stored edges
    # and existing owner/shallow interface both fail.
    small_plus = F(1, 2) + F(1, R)
    small_minus = F(1, 2) - F(1, R)
    require((P * small_plus - F(14, R)) % 1 == small_minus
            and (P * small_minus + F(14, R)) % 1 == small_plus,
            "minimal alternating odometer state cycle changed")
    require((clock(small_plus), clock(small_minus)) == (4, 3)
            and (shallow(small_plus), owner(small_plus)) == (4, 4)
            and (shallow(small_minus), owner(small_minus)) == (3, 3)
            and owner(small_plus) != shallow(small_minus)
            and owner(small_minus) != shallow(small_plus),
            "minimal alternating state/internal/interface hostile changed")

    # Relation to the independently frozen lawful alternating central scout.
    lift = S + 1
    amplitude = F(lift, (P + 1) * R)
    lawful_plus = F(1, 2) + amplitude
    lawful_minus = F(1, 2) - amplitude
    require((P * lawful_plus - F(lift, R)) % 1 == lawful_minus
            and (P * lawful_minus + F(lift, R)) % 1 == lawful_plus,
            "lawful alternating odometer cycle changed")
    lawful_edges = (
        (shallow(lawful_plus), owner(lawful_plus)),
        (shallow(lawful_minus), owner(lawful_minus)),
    )
    require(lawful_edges == ((4, 3), (3, 4))
            and owner(lawful_plus) == shallow(lawful_minus)
            and owner(lawful_minus) == shallow(lawful_plus)
            and (carry(lawful_plus), carry(lawful_minus)) == (7, 5)
            and (future_digit(lawful_plus), future_digit(lawful_minus))
            == (6, 6)
            and ((-2 * lift) % P, (2 * lift) % P) == (11, 2),
            "lawful alternating odometer typing changed")
    require(lift % P != 0 and lift not in covariant_set
            and (R - lift) not in covariant_set,
            "lawful alternating local cycle became globally clock-covariant")

    # Coefficient-level connection to four selected THM-2640 unit rows.
    # Reflection rho acts on the seven coefficient coordinates by ell -> -ell.
    # Rebuild every vector from exact rail data before checking the relations;
    # the retained metadata make explicit that this does not yet map supports.
    carrier = carry_atlas.core.build_carrier_data()
    coefficient_assets = (
        carrier[0], carry_atlas.build_pair_prefixes(carrier[0]), carrier[4],
        carrier[5], carrier[6], {},
    )
    positive_plus_type = (
        (1, 3, 0),  # rail metadata (s,ell,theta)
        7, 6, 1, 0, 0, 3,  # carry,h,kappa,sector,edge,root
    )
    positive_minus_type = ((1, 4, 12), 5, 6, 1, 0, 0, 12)
    negative_plus_type = ((1, 3, 0), 7, 6, 0, 0, 1, 1)
    negative_minus_type = ((1, 4, 12), 5, 6, 0, 0, 1, 10)
    positive_plus = rebuild_thm2640_unit_vector(
        coefficient_assets, positive_plus_type
    )
    positive_minus = rebuild_thm2640_unit_vector(
        coefficient_assets, positive_minus_type
    )
    negative_plus = rebuild_thm2640_unit_vector(
        coefficient_assets, negative_plus_type
    )
    negative_minus = rebuild_thm2640_unit_vector(
        coefficient_assets, negative_minus_type
    )
    require(positive_plus == (0, 0, 0, 0, 8, 8, 7)
            and positive_minus == (0, 5, 2, 2, 0, 0, 0)
            and negative_plus == (0, 0, 0, 0, 11, 11, 8)
            and negative_minus == (0, 6, 5, 5, 0, 0, 0),
            "selected THM-2640 normalized coefficient vectors changed")
    positive_reflection = tuple(positive_plus[(-ell) % Q]
                                for ell in range(Q))
    negative_reflection = tuple(negative_plus[(-ell) % Q]
                                for ell in range(Q))
    require(tuple(10 * value % P for value in positive_reflection)
            == positive_minus
            and tuple(4 * value % P for value in negative_reflection)
            == negative_minus,
            "selected THM-2640 reflected unit relations changed")

    print("LRC14 affine three-tooth full intrinsic-clock classification")
    print("scope=positive envelope/intrinsic-clock support; exact rail witnesses at beta=1/13; downstream full packets untested")
    print(f"K_components={len(K)} K_measure={measure(K)}")
    for sign in (1, -1):
        greedy, seams, total, union = covers[sign]
        print(
            f"sign={sign:+d}:greedy_corridors={greedy}:"
            f"seam_repairs={seams}:total_corridors={total}:"
            f"positive_beta_union={union}:exceptions=(0,1/2)"
        )
    print(
        f"physical_lift_residue_audit=modulus={R}:"
        f"THM2657_delta_to_kmod13={thm2657_residue_map}:"
        f"unit_lifts={unit_lift_count}:tradeoff_violations="
        f"{tuple(tradeoff_violations)}"
    )
    print(
        f"global_7clock_covariant_kmodR={covariant_lift_residues}:"
        "criterion=13^5_divides_k"
    )
    print(
        "universal_physical_lift_tradeoff="
        "every_nontrivial_THM2657_unit_lift_has_positive_intrinsic_"
        "three_event_support_for_both_signs_and_no_global_7clock_permutation"
    )
    print(f"exception_beta0_raw_mass={plus_zero_mass}:plus_words={plus_zero}:minus_words={minus_zero}")
    print(f"exception_beta_half_raw_mass={plus_half_mass}:plus_words={plus_half}:minus_words={minus_half}")
    print(
        f"beta_1_over_13_plus_keys={plus_keys}:component="
        f"{plus_signature['component']}:length={plus_signature['length']}:"
        f"arrivals={plus_signature['arrivals']}:"
        f"stored_edges={plus_signature['stored_edges']}:"
        "clock_transport=identity:digit_transport=j_next=h"
    )
    print(
        f"beta_1_over_13_minus_keys={minus_keys}:component="
        f"{minus_signature['component']}:length={minus_signature['length']}:"
        f"arrivals={minus_signature['arrivals']}:"
        f"stored_edges={minus_signature['stored_edges']}:"
        "clock_transport=reflection:digit_transport=j_next=12-h"
    )
    print(
        f"beta_1_over_14=intrinsic_positive:components="
        f"{len(half_clock_components)}:mass={measure(half_clock_components)}:"
        "T_squared=D_squared:13beta=13/14:not_a_7_clock_permutation:"
        "phase_sidecar_required=True"
    )
    print(
        f"fixed_beta_14_over_R=intrinsic_positive:components="
        f"{len(small_components)}:mass={measure(small_components)}:"
        f"first_component={small_signature['component']}:"
        f"arrivals={small_signature['arrivals']}:"
        f"stored_edges={small_signature['stored_edges']}"
    )
    print(
        "alternating_minimal_lifts=(-14,14):state_clocks=(4,3):"
        "stored_edges=((4,4),(3,3)):existing_interface=False:"
        "verdict=cosmetic_state_cycle"
    )
    print(
        f"alternating_lawful_lifts=({-lift},{lift}):amplitude={amplitude}:"
        f"stored_edges={lawful_edges}:existing_interface=True:"
        "carries=(7,5):future_digits=(6,6):quotient_root_steps=(11,2):"
        "global_7clock_permutation=False:scope=phase_local_cycle"
    )
    print(
        "coefficient_reflection_positive_cycle="
        f"plus_type={positive_plus_type}:Yplus={positive_plus}:"
        f"minus_type={positive_minus_type}:Yminus={positive_minus}:"
        "law=Yminus(ell)=10*Yplus(-ell)_mod13"
    )
    print(
        "coefficient_reflection_negative_cycle="
        f"plus_type={negative_plus_type}:Yplus={negative_plus}:"
        f"minus_type={negative_minus_type}:Yminus={negative_minus}:"
        "law=Yminus(ell)=4*Yplus(-ell)_mod13"
    )
    print(
        "coefficient_reflection_scope="
        "independently_rebuilt_selected_THM2640_coefficient_units_only:"
        "physical_packet_transport=UNTESTED:"
        "next_test=weighted_support_intersection_under_affine_reflection"
    )
    print("typing_law=following_shallow_phase=s*owner_phase+13beta; label permutation iff 7*13*beta is integral")
    print("verdict=PASS: fixed affine rail-clock support is positive exactly off beta=0,1/2; only typed sidecars can promote it to a full handoff")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
