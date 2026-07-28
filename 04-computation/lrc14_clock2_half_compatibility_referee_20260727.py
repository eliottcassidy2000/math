#!/usr/bin/env python3
"""Exact referee for the THM-2625/THM-2635 common-atom question.

The chronological THM-2635 unit uses

    (h, epsilon, r) = (3, 1, 9),

where the literal left half is ``{c3*x} in [113/182,126/182)``.  The
THM-2625 clock-two word Q_a contains the target-a danger condition

    {c2*(169*x)} in [-1/14,1/14) mod 1.

On the canonical row ``c3=2*169*c2``.  If ``c3*x=n+z``, ``0<=z<1``, then

    n even: {c2*169*x}=z/2,       so z in [0,1/7),
    n odd:  {c2*169*x}=1/2+z/2,   so z in [6/7,1).

The second lower endpoint is included because ``13/14`` is the included
left endpoint of ``[-1/14,1/14) mod 1``.  Thus the projection of the
clock-two target-a condition to z is exactly

    [0,26/182) union [156/182,182/182).

It is disjoint from ``[113/182,126/182)`` with a strict gap.  Hence the
h=3 common atom is empty before imposing E^ell, the middle-rail weight, or
the clock-six word.  This is an atomwise obstruction, not Fourier
cancellation.

The same calculation identifies the unique compatible canonical left-half
unit: h=10, r=2, whose half ``[15/182,28/182)`` meets the low target-a
window.  The companion constructs one exact positive two-clock atom for
that reroute and computes its frequency-13 endpoint numerator in a
Lucas-certified finite-field specialization of Z[zeta_N], where

    T = 297836897838480,
    N = lcm(169*T,13^6*T) = 13^6*T.

The ambient order N is safe for every endpoint, relation-DFT, and septimal
phase.  The endpoint numerator alone uses zeta_N^13 and therefore has
minimal phase order N/13.

Finally, the script audits the representative gauge.  A fixed physical
rail is positive for only eight of the thirteen representatives ell=sW and
is empty for five, so it does not descend to L_13^perp/<W>.  Translating the
carrier equivariantly,

    H_s(x) = H_0(x+s/13),

repairs the gauge: both clocks are invariant under this translation, and
X=13 makes the endpoint phase invariant.  Orbit saturation is also checked,
but the thirteen orbit labels are retained when reporting the source
seven-clock unit.  Summing those labels modulo 13 would erase the vector and
is deliberately not used as a nonlinear unit test.

This is a finite exact typed audit.  It is not a semantic two-root gluing,
an endpoint-current transport theorem, a row exclusion, or LRC(14).
"""

from bisect import bisect_left, bisect_right
from collections import Counter
from math import gcd

import lrc14_canonical_endpoint_current_thm2625 as current
import lrc14_cross_time_target_future_diagonal_thm2616 as cross
import lrc14_old_wall_successor_sector_thm2630 as wall


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


P = 13
T = current.T_DEN
R2 = P**2
R6 = P**6
COMMON_FACTOR = R6 // R2
N = R6 * T
X = current.X

require(T == cross.T == 297836897838480, "base-grid drift")
require(R2 == current.RDIL == 169 and R6 % R2 == 0, "clock divisibility")
require(N == 1437601819018855810320, "common cyclotomic order drift")
require(N // P == 110584755309142754640, "minimal phase order drift")
require(current.W[current.TB] == 2 * R2 * current.W[current.TA],
        "c3=2*169*c2 identity drift")


# ---------------------------------------------------------------------------
# Structural half compatibility with the clock-two target-a danger window
# ---------------------------------------------------------------------------
CLOCK2_TARGET_PROJECTION = ((0, 26), (156, 182))


def positive_interval_overlap(left, right):
    return max(left[0], right[0]) < min(left[1], right[1])


def half_interval(epsilon, probe):
    require(epsilon in (0, 1) and 1 <= probe < P, "invalid half label")
    if epsilon == 0:
        return 14 * probe, 14 * probe + 13
    return 14 * probe - 13, 14 * probe


def compatible_probes(epsilon):
    return tuple(
        probe for probe in range(1, P)
        if any(positive_interval_overlap(half_interval(epsilon, probe), window)
               for window in CLOCK2_TARGET_PROJECTION)
    )


COMPATIBLE = {0: compatible_probes(0), 1: compatible_probes(1)}
require(COMPATIBLE == {0: (1, 11, 12), 1: (1, 2, 12)},
        "global half compatibility changed")
require(half_interval(1, 9) == (113, 126), "h=3 half changed")
require(not any(positive_interval_overlap((113, 126), window)
                for window in CLOCK2_TARGET_PROJECTION),
        "h=3 unexpectedly meets the clock-two target window")
require(half_interval(1, 2) == (15, 28), "h=10 half changed")
require(positive_interval_overlap((15, 28), (0, 26)),
        "h=10 lost its clock-two-compatible subinterval")


# ---------------------------------------------------------------------------
# Lucas prime and exact common-order root
# ---------------------------------------------------------------------------
FIELD_PRIME = 5 * N + 1
PRIMITIVE_WITNESS = 58
ROOT = pow(PRIMITIVE_WITNESS, 5, FIELD_PRIME)
P_MINUS_ONE_FACTORS = {2: 4, 3: 3, 5: 2, 7: 2, 11: 1, 13: 12, 53: 1}


def factor_product(factors):
    product = 1
    for prime, exponent in factors.items():
        product *= prime**exponent
    return product


require(FIELD_PRIME == 7188009095094279051601, "field prime drift")
require(factor_product(P_MINUS_ONE_FACTORS) == FIELD_PRIME - 1,
        "incomplete p-1 factorization")
require(pow(PRIMITIVE_WITNESS, FIELD_PRIME - 1, FIELD_PRIME) == 1,
        "Lucas Fermat condition")
for prime_factor in P_MINUS_ONE_FACTORS:
    require(
        gcd(pow(PRIMITIVE_WITNESS,
                (FIELD_PRIME - 1) // prime_factor,
                FIELD_PRIME) - 1, FIELD_PRIME) == 1,
        f"Lucas gcd condition failed at {prime_factor}",
    )
require(ROOT == 656356768, "common-order root drift")
require(pow(ROOT, N, FIELD_PRIME) == 1, "root lacks N-power one")
require(all(pow(ROOT, N // prime_factor, FIELD_PRIME) != 1
            for prime_factor in P_MINUS_ONE_FACTORS),
        "root lacks exact order N")


# ---------------------------------------------------------------------------
# Common-grid interval and twisted endpoint-prefix machinery
# ---------------------------------------------------------------------------
def translate_unweighted(intervals, shift):
    output = []
    for left, right in intervals:
        new_left = (left + shift) % T
        new_right = new_left + right - left
        if new_right <= T:
            output.append((new_left, new_right))
        else:
            output.append((new_left, T))
            output.append((0, new_right - T))
    return sorted(output)


def translate_weighted(pieces, shift):
    output = []
    for left, right, weight in pieces:
        new_left = (left + shift) % T
        new_right = new_left + right - left
        if new_right <= T:
            output.append((new_left, new_right, weight))
        else:
            output.append((new_left, T, weight))
            output.append((0, new_right - T, weight))
    return sorted(output)


def restrict_base(pieces, intervals):
    return cross.old.intersect_weighted_union(
        pieces, intervals, [left for left, _ in intervals]
    )


QA = current.build_set(current.PAT_QA, current.ZERO_ELL)
QA_STARTS = [left for left, _ in QA]


def restrict_clock_two(pieces):
    """Restrict weighted T-grid pieces by Q_a(169*x), output on N-grid."""
    output = []
    count_q = len(QA)
    for left, right, weight in pieces:
        lifted_left = R2 * left
        lifted_right = R2 * right
        period = lifted_left // T
        residue = lifted_left - period * T
        index = bisect_right(QA_STARTS, residue) - 1
        if index < 0:
            index = count_q - 1
            period -= 1
        while True:
            q_left = period * T + QA[index][0]
            q_right = period * T + QA[index][1]
            if q_left >= lifted_right:
                break
            overlap_left = max(lifted_left, q_left)
            overlap_right = min(lifted_right, q_right)
            if overlap_left < overlap_right:
                output.append((
                    COMMON_FACTOR * overlap_left,
                    COMMON_FACTOR * overlap_right,
                    weight,
                ))
            index += 1
            if index == count_q:
                index = 0
                period += 1
    return output


def merge_intervals(intervals):
    merged = []
    for left, right in intervals:
        require(0 <= left < right <= T, "invalid delayed-word interval")
        if merged and left <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], right))
        else:
            merged.append((left, right))
    return merged


class DelayedEndpointEvaluator:
    """Mass and frequency-X endpoint prefixes for 1_Q(13^6*x)."""

    def __init__(self, prefix):
        starts, lengths, _ = prefix
        intervals = merge_intervals([
            (left, left + length) for left, length in zip(starts, lengths)
        ])
        require(intervals and intervals[0][0] > 0 and intervals[-1][1] < T,
                "selected future digit must avoid the circle boundary")
        self.intervals = intervals
        self.starts = [left for left, _ in intervals]
        self.lengths = [right - left for left, right in intervals]
        self.length_prefix = [0]
        for length in self.lengths:
            self.length_prefix.append(self.length_prefix[-1] + length)

        self.z = pow(ROOT, (-X) % N, FIELD_PRIME)
        self.rho = pow(self.z, T, FIELD_PRIME)
        require(self.rho != 1 and pow(self.rho, R6, FIELD_PRIME) == 1,
                "twisted period ratio changed")

        events = {}
        for left, right in intervals:
            events[left] = events.get(left, 0) + 1
            events[right] = events.get(right, 0) - 1
        boundaries = sorted((position, jump)
                            for position, jump in events.items() if jump)
        self.boundary_positions = [position for position, _ in boundaries]
        self.boundary_prefix = [0]
        for position, jump in boundaries:
            term = jump * pow(self.z, position, FIELD_PRIME)
            self.boundary_prefix.append(
                (self.boundary_prefix[-1] + term) % FIELD_PRIME
            )
        self.full_boundary = self.boundary_prefix[-1]
        self.inverse_one_minus_rho = pow(
            (1 - self.rho) % FIELD_PRIME, -1, FIELD_PRIME
        )

    def left_inside(self, coordinate):
        if coordinate == 0:
            return False
        index = bisect_right(self.starts, coordinate) - 1
        return (
            index >= 0
            and self.starts[index] < coordinate <= self.intervals[index][1]
        )

    def length_phi(self, coordinate):
        index = bisect_right(self.starts, coordinate) - 1
        if index < 0:
            return 0
        return (
            self.length_prefix[index]
            + min(coordinate - self.starts[index], self.lengths[index])
        )

    def mass_prefix(self, coordinate):
        period, residue = divmod(coordinate, T)
        return period * self.length_prefix[-1] + self.length_phi(residue)

    def endpoint_prefix(self, coordinate):
        """Endpoint sum of [0,coordinate) meet T^-6 Q on the N-grid.

        For coordinate=k*T+u, repeated internal Q-boundaries contribute

          B_Q * (1-rho^k)/(1-rho) + rho^k B_Q(<u).

        The final ``-z^coordinate`` term is present exactly when the left
        limit at coordinate lies in Q.  This half-open prefix convention
        handles a clock boundary at the truncation point without ambiguity.
        """
        period, residue = divmod(coordinate, T)
        geometric = (
            (1 - pow(self.rho, period, FIELD_PRIME))
            * self.inverse_one_minus_rho
        ) % FIELD_PRIME
        value = self.full_boundary * geometric % FIELD_PRIME
        boundary_index = bisect_left(self.boundary_positions, residue)
        value = (
            value
            + pow(self.rho, period, FIELD_PRIME)
            * self.boundary_prefix[boundary_index]
        ) % FIELD_PRIME
        if self.left_inside(residue):
            value = (value - pow(self.z, coordinate, FIELD_PRIME)) % FIELD_PRIME
        return value

    def weighted_mass_and_endpoint(self, pieces):
        mass = 0
        endpoint = 0
        for left, right, weight in pieces:
            mass += weight * (
                self.mass_prefix(right) - self.mass_prefix(left)
            )
            endpoint = (
                endpoint
                + weight * (
                    self.endpoint_prefix(right) - self.endpoint_prefix(left)
                )
            ) % FIELD_PRIME
        return mass, endpoint


# ---------------------------------------------------------------------------
# Selected THM-2635 rail, h=3 no-go, and h=10 reroute
# ---------------------------------------------------------------------------
(MODULE, PREFIXES, _WHOLE_PREFIXES, _DIGIT_MASSES, RAILS,
 PRESENT, PRESENT_STARTS) = cross.build_carrier_data()


def selected_rail(s, ell4):
    theta = wall.selected_theta(s, ell4)
    matches = [
        row for row in RAILS
        if row[0] == s and row[1] == ell4 and (row[2] - 12) % P == theta
    ]
    require(len(matches) == 1, "selected rail is not unique")
    return matches[0]


RAIL = selected_rail(1, 3)


def selected_left_half(h, ell5):
    probe = (-h - 1) % P
    require(probe != 0, "punctured graph probe")
    shift = (-h) % P
    pieces = cross.old.intersect_weighted_union(
        RAIL[3], PRESENT[ell5, shift], PRESENT_STARTS[ell5, shift]
    )
    left, right = half_interval(1, probe)
    return cross.old.intersect_weighted_comb(
        pieces, MODULE.C3, 182, left, right
    )


E_ZERO = current.build_set(current.PAT_E, current.ZERO_ELL)
H3 = selected_left_half(3, 0)
H3_E = restrict_base(H3, E_ZERO)
require((len(H3), len(H3_E), len(restrict_clock_two(H3_E)))
        == (157, 145, 0),
        "selected positive h=3 hostile control changed")

H10_BY_ELL5 = tuple(selected_left_half(10, ell5) for ell5 in range(7))
RAW_H10_VECTOR = tuple(
    cross.delayed_all_digits(pieces, PREFIXES[ell5], {})[10]
    for ell5, pieces in enumerate(H10_BY_ELL5)
)
require(all(value % cross.GLOBAL_CONTENT == 0 for value in RAW_H10_VECTOR),
        "h=10 vector left the proved primitive lattice")
NORMALIZED_H10_VECTOR = tuple(
    (value // cross.GLOBAL_CONTENT) * pow(2, -1, P) % P
    for value in RAW_H10_VECTOR
)
require(NORMALIZED_H10_VECTOR == (0, 0, 0, 1, 0, 0, 0),
        "selected h=10 normalized seven-clock vector changed")
REDUCED_H10_VECTOR = tuple(
    (NORMALIZED_H10_VECTOR[index] - NORMALIZED_H10_VECTOR[-1]) % P
    for index in range(6)
)
H10_UNIT_DETERMINANT = cross.old.sat.multiplication_determinant_7(
    REDUCED_H10_VECTOR
)
require(H10_UNIT_DETERMINANT == 1, "selected h=10 source unit changed")

# Type the reroute in THM-2640's predecessor-carry notation.  Here ``b=1``
# is THM-2635's coarse ``floor(2u)`` quotient.  It is not THM-2640's fine
# half digit ``kappa=floor(26y)-2h``: both fine-kappa values give b=1 at
# h=10, and this common-atom computation does not select between them.
# Epsilon is the half label, and the inverse of two is seven in F_13.
H10_H, H10_B, H10_EPSILON, H10_ROOT = (10, 1, 1, 2)
H10_CARRY = 7 * (H10_ROOT - H10_B - H10_EPSILON) % P
H10_OUTGOING = (-H10_ROOT) % P
H10_ORIGIN_OFFSET = (H10_OUTGOING - H10_CARRY) % P
require(tuple((2 * H10_H + fine_kappa) // P for fine_kappa in (0, 1))
        == (H10_B, H10_B),
        "h=10 coarse quotient changed across the fine half digits")
require((H10_CARRY, H10_B, H10_EPSILON, H10_ROOT,
         H10_OUTGOING, H10_ORIGIN_OFFSET) == (0, 1, 1, 2, 11, 11),
        "h=10 private-root/origin labels changed")

DELAYED_H10 = DelayedEndpointEvaluator(PREFIXES[0][10])


def endpoint_atom(carrier, ell):
    present = restrict_base(
        carrier, current.build_set(current.PAT_E, ell)
    )
    clock_two = restrict_clock_two(present)
    mass, endpoint = DELAYED_H10.weighted_mass_and_endpoint(clock_two)
    return len(present), len(clock_two), mass, endpoint


H10 = H10_BY_ELL5[0]
WMOD = tuple(value % P for value in current.W)
GAUGE_ELLS = tuple(
    tuple(s * value % P for value in WMOD) for s in range(P)
)
FIXED_GAUGE_RESULTS = tuple(endpoint_atom(H10, ell) for ell in GAUGE_ELLS)
FIXED_GAUGE_MASSES = tuple(row[2] for row in FIXED_GAUGE_RESULTS)
FIXED_GAUGE_ENDPOINTS = tuple(row[3] for row in FIXED_GAUGE_RESULTS)

WITNESS_MASS = 320917389308122335557400
WITNESS_ENDPOINT = 1441723002435168223705
require(FIXED_GAUGE_RESULTS[0]
        == (97, 187, WITNESS_MASS, WITNESS_ENDPOINT),
        "h=10 two-clock endpoint witness changed")
require(FIXED_GAUGE_MASSES
        == (WITNESS_MASS,) * 4 + (0,) * 5 + (WITNESS_MASS,) * 4,
        "fixed-carrier thirteen-gauge mass pattern changed")
require(FIXED_GAUGE_ENDPOINTS
        == (WITNESS_ENDPOINT,) * 4 + (0,) * 5 + (WITNESS_ENDPOINT,) * 4,
        "fixed-carrier thirteen-gauge endpoint pattern changed")

FIXED_GAUGE_SUPPORT = tuple(
    index for index, mass in enumerate(FIXED_GAUGE_MASSES) if mass
)
MASK_MASS = len(FIXED_GAUGE_SUPPORT)
MASK_ENERGY = MASK_MASS  # the gauge mask is an indicator
MASK_DEFECT = MASK_MASS**2 - MASK_ENERGY
MASK_RETURN = len(set(FIXED_GAUGE_SUPPORT).intersection(
    {(-index) % P for index in FIXED_GAUGE_SUPPORT}
))
require((FIXED_GAUGE_SUPPORT, MASK_MASS, MASK_ENERGY,
         MASK_DEFECT, MASK_RETURN)
        == ((0, 1, 2, 3, 9, 10, 11, 12), 8, 8, 56, 7),
        "fixed-gauge purity/return hostile changed")


# ---------------------------------------------------------------------------
# Lawful equivariant translation and orbit saturation
# ---------------------------------------------------------------------------
GAUGE_SHIFT = T // P
require(all(
    current.build_set(current.PAT_E, GAUGE_ELLS[s])
    == translate_unweighted(E_ZERO, -s * GAUGE_SHIFT)
    for s in range(P)
), "E^(ell+sW) translation law changed")

TRANSLATED_H10 = tuple(
    translate_weighted(H10, -s * GAUGE_SHIFT) for s in range(P)
)
EQUIVARIANT_RESULTS = tuple(
    endpoint_atom(TRANSLATED_H10[s], GAUGE_ELLS[s]) for s in range(P)
)
require(all(row[2:] == (WITNESS_MASS, WITNESS_ENDPOINT)
            for row in EQUIVARIANT_RESULTS),
        "equivariantly translated carrier failed to preserve the witness")
require(current.W[current.TB] % P == 0 and R6 % P == 0,
        "gauge translation changed the deep/future labels")

# Translation by 1/13 is invisible to the clock-six digit and preserves each
# of the seven raw source coefficients.  The orbit label is kept explicitly.
for s in range(P):
    for ell5, pieces in enumerate(H10_BY_ELL5):
        translated = translate_weighted(pieces, -s * GAUGE_SHIFT)
        translated_value = cross.delayed_all_digits(
            translated, PREFIXES[ell5], {}
        )[10]
        require(translated_value == RAW_H10_VECTOR[ell5],
                "translated orbit copy changed the seven-clock vector")

ORBIT_CARRIERS = TRANSLATED_H10
ORBIT_SATURATED = sorted(
    piece for carrier in ORBIT_CARRIERS for piece in carrier
)
ORBIT_RESULTS = tuple(endpoint_atom(ORBIT_SATURATED, ell)
                      for ell in GAUGE_ELLS)
require(len({row[2:] for row in ORBIT_RESULTS}) == 1,
        "orbit-saturated carrier did not descend")
ORBIT_MASS, ORBIT_ENDPOINT = ORBIT_RESULTS[0][2:]
require(ORBIT_MASS == 2567339114464978684459200,
        "orbit-saturated mass changed")
require(ORBIT_ENDPOINT == 4345774924387066738039,
        "orbit-saturated endpoint changed")

UNLABELLED_ORBIT_VECTOR = tuple(
    P * value % P for value in NORMALIZED_H10_VECTOR
)
require(UNLABELLED_ORBIT_VECTOR == (0,) * 7,
        "thirteen-copy unlabelled hostile control changed")


def main():
    print("THM-2625/2635 clock-two half-compatibility referee")
    print(f"T={T}; N6={N}; minimal_endpoint_phase_order={N // P}")
    print(f"clock2_target_projection_182={CLOCK2_TARGET_PROJECTION}")
    print(f"compatible_probes_eps0={COMPATIBLE[0]}; "
          f"compatible_probes_eps1={COMPATIBLE[1]}")
    print("h3_eps1_r9_common_atom=EMPTY: "
          "[113,126) lies strictly between 26 and 156")
    print(f"h10_source_normalized_mod13={NORMALIZED_H10_VECTOR}; "
          f"unit_determinant={H10_UNIT_DETERMINANT}")
    print("h10_private_labels(c,b,epsilon,r,outgoing,offset)="
          f"({H10_CARRY},{H10_B},{H10_EPSILON},{H10_ROOT},"
          f"{H10_OUTGOING},{H10_ORIGIN_OFFSET})")
    print(f"field_prime={FIELD_PRIME}; exact_N6_root={ROOT}; "
          "Lucas/order_gate=PASS")
    print(f"h10_witness_mass={WITNESS_MASS}; "
          f"endpoint_mod_p={WITNESS_ENDPOINT}")
    print(f"fixed_gauge_mass_hist={tuple(sorted(Counter(FIXED_GAUGE_MASSES).items()))}; "
          f"positive_gauges={(0, 1, 2, 3, 9, 10, 11, 12)}")
    print(f"fixed_mask(M,E,delta,R)=({MASK_MASS},{MASK_ENERGY},"
          f"{MASK_DEFECT},{MASK_RETURN}); purity_return_gate=NO")
    print(f"equivariant_gauges=13/13; mass={WITNESS_MASS}; "
          f"endpoint_mod_p={WITNESS_ENDPOINT}")
    print(f"orbit_saturated_mass={ORBIT_MASS}; "
          f"endpoint_mod_p={ORBIT_ENDPOINT}; descended_gauges=13/13")
    print(f"orbit_labelled_unit_copies=13; "
          f"unlabelled_mod13_vector={UNLABELLED_ORBIT_VECTOR}")
    print("verdict=PASS: h=3 is atomwise incompatible with clock two; "
          "h=10 has a nonzero common endpoint atom, but a fixed carrier "
          "does not descend; equivariant translation/orbit labelling repairs "
          "the quotient gauge")
    print("SCOPE: exact canonical-control referee; not semantic two-root "
          "gluing, current transport, a row exclusion, or LRC(14)")


if __name__ == "__main__":
    main()
