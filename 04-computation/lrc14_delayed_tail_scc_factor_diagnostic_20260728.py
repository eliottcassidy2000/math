#!/usr/bin/env python3
"""Exact delayed-tail SCC and factor-minimality diagnostic for LRC(14).

The affine odometer coordinate of THM-2657 is

    R*x = N+y,       R=13^6,
    y'={13y},        N'=13N+k+floor(13y) (mod R).

Thus recurrence of an inherited delayed word must be decided in the tail
coordinate before choosing the integer lift ``k``.  This companion compares
the two natural delayed target teeth

    target a: D_(13^3),       target b: D_(2*13^5),

and both guard sectors, while retaining the five ordinary safe factors, the
other deep safe factor, and the seven owner-clock refinements.

It proves the following exact diagnostic facts.

* The raw recurrent SCC of a repeated target tooth consists only of its
  fixed tooth centres: ``{0}`` for target a and ``{0,1/2}`` for target b.
  The safe speed 14 deletes every one of these centres.  Changing the guard
  sector therefore cannot restore an SCC.
* The target-a word has nilpotence index four in either guard sector.  The
  target-b guard-safe word has the sharp index six: a strict five-event word
  is exhibited, while two target-b hits followed by the state-five speed-14
  safe bit force zero.  The guard-danger target-b word is also nilpotent by
  depth six (a strict four-event control is retained; sharpness is not
  claimed for that sector).
* Among single factor-type deletions, only deleting the repeated target tooth
  opens recurrence.  Deleting it gives fixed tails 5/12 (guard-safe) and
  1/12 (guard-danger).  Each fixed tail admits an exact two-state odometer
  fibre lift with nonconstant glued clock edges and the dynamically typed
  THM-2640 present factor.  This is a positive carrier control, not a rail,
  private-unit, endpoint, row-exclusion, or LRC(14) certificate.
* Refining by intrinsic guard/source edge debt rather than deleting a target
  succeeds at the exact carrier level.  The resulting ``B_debt,A_debt``
  automaton has the exact period-four tail orbit
  ``7,91,163,79`` divided by 170.  Its guard sectors are respectively
  danger,safe,danger,safe.  It has a lawful odometer lift with glued clocks
  ``0->1->0->1->0`` and strictly positive aligned THM-2584 rail, delayed,
  present, pointwise private-root, and primitive-unit factors at every event.
  Here ``a,b`` name delayed-tooth roles, not THM-2305 terminal words, and all
  four rails retain fixed ``PAT_QB``.  The orbit is a recurrent SCC of this
  newly defined guard-cospan edge automaton, but not an inherited physical or
  semantic SCC.  Literal target-word alternation and endpoint transport are
  not supplied.

All interval arithmetic is exact on the canonical T grid.  The finite-horizon
controls are midpoints of grid atoms.  The mixed rational orbit is checked
strictly against every retained interval, so no half-open endpoint convention
is used as a positive witness.
"""

from bisect import bisect_right
from fractions import Fraction

import lrc14_cross_time_target_future_diagonal_thm2616 as cross
import lrc14_alternate_arrival_physical_rail_handoff as rail


P = 13
R = P**6
Q7 = 7


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


(MODULE, PREFIXES, WHOLE_PREFIXES, DIGIT_MASSES, RAILS,
 PRESENT, PRESENT_STARTS) = cross.build_carrier_data()
T = cross.T

require(MODULE.C1 == P and MODULE.C2 == P**3
        and MODULE.C3 == 2 * P**5, "canonical delayed scales changed")
require(T % R == 0, "canonical grid lost the odometer scale")


def merge_support(intervals):
    """Union sorted/unsorted half-open interval lists."""
    out = []
    for left, right in sorted(intervals):
        if out and left <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], right))
        else:
            out.append((left, right))
    return out


def intersect_sorted(left, right):
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


def build_word(target, guard_mode, retain_target=True):
    """Build one target/guard delayed word before its owner-clock cut."""
    require(target in ("a", "b"), "unknown target")
    require(guard_mode in ("safe", "danger"), "unknown guard sector")
    target_speed = MODULE.C2 if target == "a" else MODULE.C3
    other_speed = MODULE.C3 if target == "a" else MODULE.C2
    intervals = [(0, T)]
    if retain_target:
        intervals = MODULE.intersect_comb(
            intervals, target_speed, 182, -13, 13
        )
    if guard_mode == "safe":
        intervals = MODULE.subtract_comb(
            intervals, MODULE.W[MODULE.GUARD], 91, -13, 13
        )
    else:
        intervals = MODULE.intersect_comb(
            intervals, MODULE.W[MODULE.GUARD], 91, -13, 13
        )
    for index in MODULE.UNIT_IDX:
        intervals = MODULE.subtract_comb(
            intervals, MODULE.W[index], 182, -13, 13
        )
    intervals = MODULE.subtract_comb(
        intervals, other_speed, 182, -13, 13
    )
    return intervals


def owner_words(word):
    return tuple(
        MODULE.subtract_comb(
            word, MODULE.C1, 182, 26 * ell - 13, 26 * ell + 13
        )
        for ell in range(Q7)
    )


def strict_hit(intervals, value):
    """Return a strict normalized margin, or None."""
    coordinate = value * T
    starts = [left for left, _right in intervals]
    index = bisect_right(starts, coordinate) - 1
    if index < 0:
        return None
    left, right = intervals[index]
    if not left < coordinate < right:
        return None
    return min(coordinate - left, right - coordinate) / T


def strict_labels(words, value):
    return tuple(
        (ell, margin)
        for ell, intervals in enumerate(words)
        if (margin := strict_hit(intervals, value)) is not None
    )


def tail_step(value):
    return (P * value) % 1


EXPECTED_WORD_CENSUS = {
    ("a", "safe"): (32_578, 10_786_297_480_556),
    ("a", "danger"): (14_906, 4_936_568_463_918),
    ("b", "safe"): (188_634, 10_807_572_370_920),
    ("b", "danger"): (86_950, 4_981_755_322_960),
}

WORDS = {}
OWNER_WORDS = {}
for key, expected in EXPECTED_WORD_CENSUS.items():
    word = build_word(*key)
    census = (len(word), sum(right - left for left, right in word))
    require(census == expected, f"{key} word census changed")
    refinements = owner_words(word)
    require(
        merge_support(piece for row in refinements for piece in row) == word,
        f"{key} owner-clock refinements stopped covering the base word",
    )
    WORDS[key] = word
    OWNER_WORDS[key] = refinements

require(WORDS["a", "safe"] == MODULE.build_word_Ta(),
        "target-a guard-safe word no longer matches the canonical builder")

# The two guard sectors partition the word with the guard factor deleted.
for target in ("a", "b"):
    without_guard = [(0, T)]
    target_speed = MODULE.C2 if target == "a" else MODULE.C3
    other_speed = MODULE.C3 if target == "a" else MODULE.C2
    without_guard = MODULE.intersect_comb(
        without_guard, target_speed, 182, -13, 13
    )
    for index in MODULE.UNIT_IDX:
        without_guard = MODULE.subtract_comb(
            without_guard, MODULE.W[index], 182, -13, 13
        )
    without_guard = MODULE.subtract_comb(
        without_guard, other_speed, 182, -13, 13
    )
    require(
        merge_support(WORDS[target, "safe"] + WORDS[target, "danger"])
        == without_guard,
        f"{target} guard sectors stopped partitioning the raw word",
    )
    require(
        not intersect_sorted(WORDS[target, "safe"],
                             WORDS[target, "danger"]),
        f"{target} guard sectors acquired positive overlap",
    )


# -------------------------------------------------------------------------
# Strict finite-horizon controls
# -------------------------------------------------------------------------

# Each stored integer z represents the strict midpoint (2z+1)/(2T).
HORIZON_WITNESSES = {
    ("a", "safe"): (53_548_336_995_378, 3),
    ("a", "danger"): (1_897_955_067_622, 3),
    ("b", "safe"): (240_416_984_448_501, 5),
    ("b", "danger"): (295_924_946_234_062, 4),
}

HORIZON_ROWS = {}
for key, (grid_atom, horizon) in HORIZON_WITNESSES.items():
    value = Fraction(2 * grid_atom + 1, 2 * T)
    path = []
    margins = []
    for _step in range(horizon):
        labels = strict_labels(OWNER_WORDS[key], value)
        require(len(labels) == 6, f"{key} strict owner-label census changed")
        path.append(tuple(ell for ell, _margin in labels))
        margins.extend(margin for _ell, margin in labels)
        value = tail_step(value)
    require(not strict_labels(OWNER_WORDS[key], value),
            f"{key} stored transient unexpectedly extended")
    HORIZON_ROWS[key] = (
        Fraction(2 * grid_atom + 1, 2 * T), horizon,
        tuple(path), min(margins), value,
    )

# Analytic nilpotence invoices.  For target a, three consecutive target hits
# shrink z={13^3 y} to radius 1/(14*13^2), so speed 14 is dangerous at state
# three.  For target b, two consecutive target hits shrink the deviation of
# z={13^5 y} from {0,1/2} to radius 1/(28*13), so speed 14 is dangerous at
# state five.  These are the (p,r,a,q)=(13,3,1,14) and (13,5,2,14) cases of
# the sparse nested-danger mechanism.
require(Fraction(14, 14 * P**2) < Fraction(1, 14),
        "target-a contraction no longer beats speed-14 safety")
require(Fraction(14, 28 * P) < Fraction(1, 14),
        "target-b contraction no longer beats speed-14 safety")

NILPOTENCE_UPPER = {
    ("a", "safe"): 4,
    ("a", "danger"): 4,
    ("b", "safe"): 6,
    ("b", "danger"): 6,
}

# The target-a bounds and target-b guard-safe bound are sharp by the strict
# witnesses.  Only a lower/upper bracket is asserted for target-b danger.
require(HORIZON_ROWS["a", "safe"][1] + 1 == 4
        and HORIZON_ROWS["a", "danger"][1] + 1 == 4,
        "target-a sharp nilpotence controls changed")
require(HORIZON_ROWS["b", "safe"][1] + 1 == 6,
        "target-b guard-safe sharp nilpotence control changed")


# -------------------------------------------------------------------------
# The recurrent centre automata and the unique useful factor deletion
# -------------------------------------------------------------------------

CENTRES = {
    "a": (Fraction(0),),
    "b": (Fraction(0), Fraction(1, 2)),
}


def circular_distance(value):
    residue = value % 1
    return min(residue, 1 - residue)


def target_danger(target, value):
    multiplier = 1 if target == "a" else 2
    return circular_distance(multiplier * value) < Fraction(1, 14)


def speed_safe(speed, value):
    return circular_distance(speed * value) >= Fraction(1, 14)


for target, centres in CENTRES.items():
    require(all(tail_step(centre) == centre for centre in centres),
            f"{target} target centres stopped being fixed")
    require(all(target_danger(target, centre) for centre in centres),
            f"{target} target centre left its danger tooth")
    require(all(not speed_safe(14, centre) for centre in centres),
            f"{target} recurrent centre escaped the speed-14 killer")

# Half-open endpoint controls: every included left endpoint maps to the
# excluded right endpoint, so no endpoint-only recurrent SCC was hidden by
# the strict-centre argument.
ENDPOINT_IMAGES = {
    "a": ((Fraction(13, 14), Fraction(1, 14)),),
    "b": (
        (Fraction(27, 28), Fraction(15, 28)),
        (Fraction(13, 28), Fraction(1, 28)),
    ),
}
for target, rows in ENDPOINT_IMAGES.items():
    for left_endpoint, excluded_right in rows:
        require(tail_step(left_endpoint) == excluded_right,
                f"{target} endpoint transition changed")
        require(not target_danger(target, excluded_right),
                f"{target} excluded endpoint entered the strict target tooth")

# If the repeated target tooth is deleted, exact fixed tails survive every
# remaining factor.  Both targets share the same control in a given guard
# sector because the other deep speed is safe at that control.
FIXED_TAILS = {
    "safe": (Fraction(5, 12), Fraction(7, 12)),
    "danger": (Fraction(1, 12), Fraction(11, 12)),
}

DELETED_WORDS = {}
DELETED_OWNER_WORDS = {}
for target in ("a", "b"):
    for guard_mode in ("safe", "danger"):
        key = (target, guard_mode)
        word = build_word(target, guard_mode, retain_target=False)
        words = owner_words(word)
        DELETED_WORDS[key] = word
        DELETED_OWNER_WORDS[key] = words
        for value in FIXED_TAILS[guard_mode]:
            require(tail_step(value) == value,
                    "target-deleted tail stopped being fixed")
            labels = strict_labels(words, value)
            require(len(labels) == 6,
                    f"{key} target-deleted fixed tail lost owner support")


def clock(value):
    return int((Q7 * (value % 1) + Fraction(1, 2)) // 1) % Q7


def shallow(value):
    return clock(P * value)


def owner(value):
    return clock(P * P * value)


def handoff(value, lift):
    return (P * value + Fraction(lift, R)) % 1


def present_margin(value, ell, future_digit):
    return strict_hit(PRESENT[ell, (-future_digit) % P], value)


# The first cycle uses guard-safe y=5/12 and clocks 0->1->0.  The second
# uses guard-danger y=1/12 and clocks 0->2->0.  All N values are chosen so
# that the two unique fibre-closing lifts are thirteenth-units.
FIBRE_CYCLES = {
    "safe": (Fraction(5, 12), (744_626, 769_107),
             (742_582, 399_848), ((0, 1), (1, 0))),
    "danger": (Fraction(1, 12), (748_707, 826_230),
               (746_656, 4_488_143), ((0, 2), (2, 0))),
}

FIBRE_ROWS = {}
for guard_mode, (fixed_tail, coarse_states, lifts,
                 expected_edges) in FIBRE_CYCLES.items():
    states = tuple((Fraction(N) + fixed_tail) / R for N in coarse_states)
    future_digit = int(P * fixed_tail)
    require(tuple(handoff(states[i], lifts[i])
                  for i in range(2)) == (states[1], states[0]),
            f"{guard_mode} target-deleted fibre cycle failed to close")
    require(all((R * state) % 1 == fixed_tail for state in states),
            f"{guard_mode} fibre cycle changed its fixed base tail")
    require(all(lift % P for lift in lifts),
            f"{guard_mode} fibre cycle acquired a quotient-zero lift")
    edges = tuple((shallow(state), owner(state)) for state in states)
    require(edges == expected_edges and all(a != b for a, b in edges),
            f"{guard_mode} intrinsic clock cycle changed")
    require(edges[0][1] == edges[1][0]
            and edges[1][1] == edges[0][0],
            f"{guard_mode} owner-to-next-shallow gluing failed")
    present_margins = tuple(
        present_margin(state, edges[index][0], future_digit)
        for index, state in enumerate(states)
    )
    require(all(margin is not None and margin > 0
                for margin in present_margins),
            f"{guard_mode} fibre cycle left the dynamic present packet")

    delayed_margins = {}
    for target in ("a", "b"):
        words = DELETED_OWNER_WORDS[target, guard_mode]
        margins = tuple(
            strict_hit(words[edges[index][0]], fixed_tail)
            for index in range(2)
        )
        require(all(margin is not None and margin > 0 for margin in margins),
                f"{target, guard_mode} aligned deleted-word clock failed")
        delayed_margins[target] = margins

    require((14 * fixed_tail) % 1 != 0,
            f"{guard_mode} positive cycle became speed-14 resonant")
    FIBRE_ROWS[guard_mode] = (
        fixed_tail, coarse_states, states, lifts, edges,
        present_margins, delayed_margins,
    )


# -------------------------------------------------------------------------
# A positive recurrent SCC in the intrinsic guard-cospan debt automaton
# -------------------------------------------------------------------------

# The two homogeneous obstruction patterns are aaa (three consecutive
# target-a hits before state-three speed-14 safety) and bb (two consecutive
# target-b hits before state-five speed-14 safety).  The alternating delayed-
# tooth schedule avoids both.  The schedule is only a search coordinate: the
# final A_debt/B_debt predicates below recover its labels intrinsically.  They
# remain cospan labels, not terminal-word or rail-target labels.  The exact
# tail and fibre data below were found by first
# searching the periodic base lattice 1/(13^4-1), and only then solving the
# odometer fibre, as required by the skew-product typing.
MIXED_TARGETS = ("b", "a", "b", "a")
MIXED_GUARDS = ("danger", "safe", "danger", "safe")
MIXED_TAILS = (
    Fraction(7, 170), Fraction(91, 170),
    Fraction(163, 170), Fraction(79, 170),
)
MIXED_COARSE = (2_572_844, 2_254_279, 2_572_531, 2_254_279)
MIXED_LIFTS = (2_594_970, 2_227_752, 2_599_027, 2_228_065)
MIXED_EDGES = ((0, 1), (1, 0), (0, 1), (1, 0))
MIXED_FUTURE_DIGITS = (0, 6, 12, 6)
MIXED_RAIL_KEYS = (
    (6, 1, 1, 0), (6, 1, 0, 12),
    (6, 1, 1, 0), (6, 1, 0, 12),
)
MIXED_RAIL_WEIGHTS = (
    27_581_135_604, 27_582_102_210,
    27_581_135_604, 27_582_102_210,
)
MIXED_HALF_EDGES = (1, 1, 0, 0)  # 0=left, 1=right
MIXED_PRIVATE_ROOTS = (2, 3, 2, 3)
MIXED_D_ROOTS = (1, 1, 12, 12)
MIXED_BROAD_UNIT_VECTORS = (
    (9, 0, 0, 0, 0, 1, 0),
    (0, 0, 0, 2, 0, 0, 0),
    (9, 0, 0, 0, 10, 10, 0),
    (0, 0, 0, 11, 0, 0, 0),
)
MIXED_BROAD_UNIT_DETERMINANTS = (12, 12, 10, 12)
MIXED_DEBT_UNIT_VECTORS = (
    (9, 0, 0, 0, 0, 12, 0),
    (0, 0, 0, 2, 0, 0, 0),
    (9, 0, 0, 0, 2, 4, 0),
    (0, 0, 0, 11, 0, 0, 0),
)
MIXED_DEBT_UNIT_DETERMINANTS = (12, 12, 8, 12)
MIXED_ANNULUS_B_VECTORS = (
    (0, 0, 0, 0, 0, 2, 0),
    (0, 0, 0, 0, 8, 6, 0),
)
require(tuple((narrow + annulus) % P
              for narrow, annulus in zip(
                  MIXED_DEBT_UNIT_VECTORS[0],
                  MIXED_ANNULUS_B_VECTORS[0]))
        == MIXED_BROAD_UNIT_VECTORS[0]
        and tuple((narrow + annulus) % P
                  for narrow, annulus in zip(
                      MIXED_DEBT_UNIT_VECTORS[2],
                      MIXED_ANNULUS_B_VECTORS[1]))
        == MIXED_BROAD_UNIT_VECTORS[2],
        "mixed narrow/annulus unit vectors stopped summing to broad danger")

# Complete period-four tail-lattice audit.  Every fixed point of B^4 is
# j/(13^4-1).  The physical delayed BABA teeth, five ordinary safeties, and
# opposite deep safeties leave only the two b-phase rotations of the displayed
# orbit.  (The teeth alone also contain the centre j=0, killed by speed 14.)
# At both survivors, the b vertices fail unshifted guard safety and the a
# vertices fail unshifted c1 safety.  Hence no other period-four tail can
# repair the missing terminal-word factors while retaining the delayed core.
PERIOD4_DENOMINATOR = P**4 - 1


def period4_path(index):
    return tuple(
        Fraction(index * P**event % PERIOD4_DENOMINATOR,
                 PERIOD4_DENOMINATOR)
        for event in range(4)
    )


def mixed_delayed_core(target, value):
    target_speed = MODULE.C2 if target == "a" else MODULE.C3
    other_speed = MODULE.C3 if target == "a" else MODULE.C2
    return (
        circular_distance(target_speed * value) < Fraction(1, 14)
        and all(speed_safe(MODULE.W[index], value)
                for index in MODULE.UNIT_IDX)
        and speed_safe(other_speed, value)
    )


BABA_PERIOD4_INDICES = tuple(
    index for index in range(PERIOD4_DENOMINATOR)
    if all(mixed_delayed_core(target, value)
           for target, value in zip(MIXED_TARGETS, period4_path(index)))
)
require(BABA_PERIOD4_INDICES == (1_176, 27_384),
        "complete BABA period-four tail lattice changed")
require(tuple(Fraction(index, PERIOD4_DENOMINATOR)
              for index in BABA_PERIOD4_INDICES)
        == (Fraction(7, 170), Fraction(163, 170)),
        "BABA lattice stopped matching the displayed orbit")
BABA_WITH_B_GUARD_SAFE = tuple(
    index for index in BABA_PERIOD4_INDICES
    if all(circular_distance(value) >= Fraction(1, 7)
           for target, value in zip(MIXED_TARGETS, period4_path(index))
           if target == "b")
)
BABA_WITH_A_C1_SAFE = tuple(
    index for index in BABA_PERIOD4_INDICES
    if all(speed_safe(P, value)
           for target, value in zip(MIXED_TARGETS, period4_path(index))
           if target == "a")
)
require(BABA_WITH_B_GUARD_SAFE == (),
        "a BABA period-four tail acquired b-vertex guard safety")
require(BABA_WITH_A_C1_SAFE == (),
        "a BABA period-four tail acquired a-vertex unshifted-c1 safety")


def broad_guard_safe(value):
    """The canonical guard safety used in the delayed cospan."""
    return circular_distance(value) >= Fraction(1, 7)


def narrow_guard_danger(value):
    return circular_distance(value) < Fraction(1, 14)


def broad_guard_annulus(value):
    distance = circular_distance(value)
    return Fraction(1, 14) <= distance < Fraction(1, 7)


def unshifted_c1_safe(value):
    return speed_safe(P, value)


def debt_sector_labels(value, ell):
    """Intrinsic labels of the new guard-cospan edge automaton.

    A_debt retains the guard-safe target-a owner word but records danger of
    the unshifted c1 factor.  B_debt retains the guard-danger target-b owner
    word, requires the current guard to lie in its narrow 1/14 tooth (hence
    also in the broad 1/7 tooth), and records unshifted c1 safety.  Membership
    in either owner word includes safety of the selected *shifted* c1
    complement.
    """
    labels = []
    if (strict_hit(OWNER_WORDS["a", "safe"][ell], value) is not None
            and broad_guard_safe(value)
            and not unshifted_c1_safe(value)):
        labels.append("A_debt")
    if (strict_hit(OWNER_WORDS["b", "danger"][ell], value) is not None
            and narrow_guard_danger(value)
            and unshifted_c1_safe(value)):
        labels.append("B_debt")
    return tuple(labels)

require(tuple(tail_step(value) for value in MIXED_TAILS)
        == MIXED_TAILS[1:] + MIXED_TAILS[:1],
        "mixed target tail orbit failed to close")
require(len(set(MIXED_TAILS)) == 4,
        "mixed tail orbit lost exact period four")
require(tuple(int(P * value) for value in MIXED_TAILS)
        == MIXED_FUTURE_DIGITS,
        "mixed tail future digits changed")

MIXED_STATES = tuple(
    (Fraction(N) + value) / R
    for N, value in zip(MIXED_COARSE, MIXED_TAILS)
)
require(tuple(handoff(MIXED_STATES[index], MIXED_LIFTS[index])
              for index in range(4))
        == MIXED_STATES[1:] + MIXED_STATES[:1],
        "mixed target odometer lift failed to close")
require(all(lift % P for lift in MIXED_LIFTS),
        "mixed target cycle acquired a quotient-zero lift")
require(tuple((shallow(state), owner(state)) for state in MIXED_STATES)
        == MIXED_EDGES,
        "mixed target intrinsic clock edges changed")
require(all(MIXED_EDGES[index][1]
            == MIXED_EDGES[(index + 1) % 4][0]
            for index in range(4)),
        "mixed target owner-to-next-shallow gluing failed")

MIXED_DEBT_SECTORS = tuple(
    debt_sector_labels(value, edge[0])
    for value, edge in zip(MIXED_TAILS, MIXED_EDGES)
)
require(MIXED_DEBT_SECTORS
        == (("B_debt",), ("A_debt",), ("B_debt",), ("A_debt",)),
        "mixed orbit lost its intrinsic debt-sector labels")
MIXED_EDGE_DEBT = tuple(
    int(narrow_guard_danger(P * value))
    for value in MIXED_TAILS
)
require(MIXED_EDGE_DEBT == (0, 1, 0, 1),
        "mixed orbit edge-debt word changed")
require(tuple(circular_distance(P * value) for value in MIXED_TAILS)
        == tuple(circular_distance(value)
                 for value in MIXED_TAILS[1:] + MIXED_TAILS[:1]),
        "source-to-next-guard debt covariance changed")
require(all((not broad_guard_safe(value))
            == (narrow_guard_danger(value) or broad_guard_annulus(value))
            and not (narrow_guard_danger(value)
                     and broad_guard_annulus(value))
            for value in MIXED_TAILS),
        "broad guard danger stopped splitting into narrow tooth and annulus")


def iterate_tail(value, steps):
    for _step in range(steps):
        value = tail_step(value)
    return value


# Global identities behind the proposed debt queue are c1=13 and c2=13^3:
# their narrow danger bits at y are the narrow guard bit one and three tail
# steps later.  The exact period-four orbit makes both debts land in the same
# B_debt class.
require(all(
    narrow_guard_danger(P * value)
    == narrow_guard_danger(iterate_tail(value, 1))
    and narrow_guard_danger(MODULE.C2 * value)
    == narrow_guard_danger(iterate_tail(value, 3))
    for value in MIXED_TAILS
), "mixed c1/c2 debt-queue identities changed")
MIXED_DEBT_QUEUE_1_3 = tuple(
    (int(narrow_guard_danger(P * value)),
     int(narrow_guard_danger(MODULE.C2 * value)))
    for value in MIXED_TAILS
)
require(MIXED_DEBT_QUEUE_1_3 == ((0, 0), (1, 1), (0, 0), (1, 1)),
        "mixed +1/+3 debt queue stopped identifying the B class")

MIXED_DELAYED_MARGINS = []
MIXED_PRESENT_MARGINS = []
for index, (target, guard_mode, value, state, edge, future_digit) in enumerate(
        zip(MIXED_TARGETS, MIXED_GUARDS, MIXED_TAILS, MIXED_STATES,
            MIXED_EDGES, MIXED_FUTURE_DIGITS)):
    ell = edge[0]
    delayed = strict_hit(OWNER_WORDS[target, guard_mode][ell], value)
    present = present_margin(state, ell, future_digit)
    require(delayed is not None and delayed > 0,
            f"mixed target event {index} left its aligned delayed word")
    require(present is not None and present > 0,
            f"mixed target event {index} left its dynamic present packet")
    require((14 * value) % 1 != 0,
            f"mixed target event {index} became speed-14 resonant")
    MIXED_DELAYED_MARGINS.append(delayed)
    MIXED_PRESENT_MARGINS.append(present)
MIXED_DELAYED_MARGINS = tuple(MIXED_DELAYED_MARGINS)
MIXED_PRESENT_MARGINS = tuple(MIXED_PRESENT_MARGINS)


def strict_weighted_hit(pieces, value):
    coordinate = value * T
    starts = [left for left, _right, _weight in pieces]
    index = bisect_right(starts, coordinate) - 1
    if index < 0:
        return None
    left, right, weight = pieces[index]
    if not left < coordinate < right:
        return None
    return weight, min(coordinate - left, right - coordinate) / T


require(rail.old.T == T, "mixed rail bank changed its exact grid")
MIXED_RAIL_BANK = rail.build_rail_bank()
MIXED_RAIL_ROWS = tuple(
    strict_weighted_hit(MIXED_RAIL_BANK[key], state)
    for key, state in zip(MIXED_RAIL_KEYS, MIXED_STATES)
)
require(all(row is not None for row in MIXED_RAIL_ROWS),
        "mixed target cycle left its selected THM-2584 rail")
require(tuple(row[0] for row in MIXED_RAIL_ROWS) == MIXED_RAIL_WEIGHTS,
        "mixed target selected rail weights changed")
require(all(key[2] == edge[1]
            for key, edge in zip(MIXED_RAIL_KEYS, MIXED_EDGES)),
        "mixed target rail owner stopped matching the intrinsic edge")


def root_memberships(value):
    """Literal THM-2640 left/right C3 half-root memberships."""
    phase = (MODULE.C3 * value) % 1
    out = []
    for edge in (0, 1):
        for root in range(1, P):
            low = 14 * root - 13 if edge == 0 else 14 * root
            high = low + 13
            if Fraction(low, 182) <= phase < Fraction(high, 182):
                out.append((edge, root))
    return tuple(out)


MIXED_DELTAS = tuple(2 * lift % P for lift in MIXED_LIFTS)
for index, state in enumerate(MIXED_STATES):
    require((MIXED_HALF_EDGES[index], MIXED_PRIVATE_ROOTS[index])
            in root_memberships(state),
            f"mixed target event {index} lost its private half/root")
    d_state = (P * state) % 1
    require(any(root == MIXED_D_ROOTS[index]
                for _edge, root in root_memberships(d_state)),
            f"mixed target event {index} intermediate D-root changed")
    require((MIXED_D_ROOTS[index] + MIXED_DELTAS[index]) % P
            == MIXED_PRIVATE_ROOTS[(index + 1) % 4],
            f"mixed target event {index} private-root handoff failed")

# The complete coefficient normalization is recomputed in
# lrc14_target_b_successor_content_mixed_unit_20260728.py.  Here we retain a
# cheap independent use of the frozen /26 vectors: normalize by the selected
# private root, reduce modulo Phi_7, and require the exact nonzero determinant.
mixed_unit_determinants = []
for vector, root in zip(MIXED_DEBT_UNIT_VECTORS, MIXED_PRIVATE_ROOTS):
    scaled = tuple(value * pow(root, -1, P) % P for value in vector)
    reduced = tuple((scaled[i] - scaled[-1]) % P for i in range(Q7 - 1))
    mixed_unit_determinants.append(
        cross.old.sat.multiplication_determinant_7(reduced)
    )
require(tuple(mixed_unit_determinants) == MIXED_DEBT_UNIT_DETERMINANTS,
        "mixed target primitive-unit determinants changed")
broad_unit_determinants = []
for vector, root in zip(MIXED_BROAD_UNIT_VECTORS, MIXED_PRIVATE_ROOTS):
    scaled = tuple(value * pow(root, -1, P) % P for value in vector)
    reduced = tuple((scaled[i] - scaled[-1]) % P for i in range(Q7 - 1))
    broad_unit_determinants.append(
        cross.old.sat.multiplication_determinant_7(reduced)
    )
require(tuple(broad_unit_determinants) == MIXED_BROAD_UNIT_DETERMINANTS,
        "mixed broad guard unit determinants changed")


def quotient_multiply(left, right):
    """Multiply in F_13[z]/(1+z+...+z^6), degree at most five."""
    coefficients = [0] * 11
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            coefficients[i + j] = (coefficients[i + j] + a * b) % P
    for degree in range(10, 5, -1):
        coefficient = coefficients[degree] % P
        for offset in range(6):
            target = degree - 6 + offset
            coefficients[target] = (coefficients[target] - coefficient) % P
        coefficients[degree] = 0
    return tuple(coefficients[:6])


MIXED_UNIT_POLYNOMIALS = tuple(
    tuple(
        (value * pow(root, -1, P)
         - vector[-1] * pow(root, -1, P)) % P
        for value in vector[:-1]
    )
    for vector, root in zip(MIXED_DEBT_UNIT_VECTORS, MIXED_PRIVATE_ROOTS)
)
require(MIXED_UNIT_POLYNOMIALS == (
    (11, 0, 0, 0, 0, 6),
    (0, 0, 0, 5, 0, 0),
    (11, 0, 0, 0, 1, 2),
    (0, 0, 0, 8, 0, 0),
), "mixed normalized unit polynomials changed")
MIXED_NAIVE_HOLONOMY = (1, 0, 0, 0, 0, 0)
for polynomial in MIXED_UNIT_POLYNOMIALS:
    MIXED_NAIVE_HOLONOMY = quotient_multiply(
        MIXED_NAIVE_HOLONOMY, polynomial
    )
require(MIXED_NAIVE_HOLONOMY == (9, 2, 8, 7, 6, 9),
        "mixed naive multiplicative holonomy changed")
MIXED_NAIVE_HOLONOMY_DETERMINANT = (
    cross.old.sat.multiplication_determinant_7(MIXED_NAIVE_HOLONOMY)
)
require(MIXED_NAIVE_HOLONOMY_DETERMINANT == 5,
        "mixed naive holonomy stopped being a unit")
power = (1, 0, 0, 0, 0, 0)
MIXED_NAIVE_HOLONOMY_ORDER = None
for exponent in range(1, P**6):
    power = quotient_multiply(power, MIXED_NAIVE_HOLONOMY)
    if power == (1, 0, 0, 0, 0, 0):
        MIXED_NAIVE_HOLONOMY_ORDER = exponent
        break
require(MIXED_NAIVE_HOLONOMY_ORDER == 168,
        "mixed naive holonomy order changed")


def main():
    print("LRC14 delayed-tail recurrent-SCC factor diagnostic")
    print("status=VERIFIED-EXACT DIAGNOSTIC; canonical typed row")
    print(f"p={P}; R={R}; T={T}; tail_map=y->{P}y mod1")
    print("universe=targets(a:13^3,b:2*13^5) x "
          "guards(safe,danger) x seven owner-clock complements")
    for key in (("a", "safe"), ("a", "danger"),
                ("b", "safe"), ("b", "danger")):
        intervals, mass = EXPECTED_WORD_CENSUS[key]
        row = HORIZON_ROWS[key]
        upper = NILPOTENCE_UPPER[key]
        sharp = row[1] + 1 == upper
        print(
            f"word={key}:intervals={intervals}:grid_mass={mass}:"
            f"strict_transient={row[1]}:nilpotence_upper={upper}:"
            f"sharp={sharp}:witness_y0={row[0]}:"
            f"min_margin={row[3]}:clock_label_words={row[2]}"
        )
    print("analytic_a_bits=D_13^3@[0,1,2]+safe_14@3 -> empty; "
          "guard-independent")
    print("analytic_b_bits=D_(2*13^5)@[0,1]+safe_14@5 -> empty; "
          "guard-independent")
    print("raw_target_recurrent_SCCs="
          "a:{0->0};b:{0->0,1/2->1/2}; safe_14 deletes every vertex")
    print("endpoint_only_SCCs=0; included left endpoints map to excluded right "
          "endpoints")
    print("single_factor_deletion=only deleting the repeated target tooth "
          "opens recurrence; every non-target deletion leaves a centre killer")
    for guard_mode in ("safe", "danger"):
        (fixed_tail, coarse, states, lifts, edges,
         present_margins, delayed_margins) = FIBRE_ROWS[guard_mode]
        print(
            f"target_deleted_cycle={guard_mode}:tail={fixed_tail}:"
            f"coarse={coarse}:x={states}:lifts={lifts}:"
            f"lift_mod13={tuple(lift % P for lift in lifts)}:edges={edges}:"
            f"present_margins={present_margins}:"
            f"delayed_margins={delayed_margins}:speed14_resonant=False"
        )
    print(
        f"guard_cospan_edge_recurrent_SCC=delayed_teeth:{MIXED_TARGETS}:"
        f"guards:{MIXED_GUARDS}:tails:{MIXED_TAILS}:"
        f"future_digits:{MIXED_FUTURE_DIGITS}:coarse:{MIXED_COARSE}:"
        f"x:{MIXED_STATES}:lifts:{MIXED_LIFTS}:"
        f"lift_mod13:{tuple(lift % P for lift in MIXED_LIFTS)}:"
        f"edges:{MIXED_EDGES}:rail_keys:{MIXED_RAIL_KEYS}:"
        f"rail_rows:{MIXED_RAIL_ROWS}:"
        f"half_edges:{MIXED_HALF_EDGES}:private_roots:{MIXED_PRIVATE_ROOTS}:"
        f"D_roots:{MIXED_D_ROOTS}:deltas:{MIXED_DELTAS}:"
        f"broad_unit_vectors_div26_mod13:{MIXED_BROAD_UNIT_VECTORS}:"
        f"broad_unit_determinants:{MIXED_BROAD_UNIT_DETERMINANTS}:"
        f"narrow_debt_unit_vectors_div26_mod13:{MIXED_DEBT_UNIT_VECTORS}:"
        f"narrow_debt_unit_determinants:{MIXED_DEBT_UNIT_DETERMINANTS}:"
        f"annulus_b_vectors_div26_mod13:{MIXED_ANNULUS_B_VECTORS}:"
        f"debt_sectors:{MIXED_DEBT_SECTORS}:edge_debt:{MIXED_EDGE_DEBT}:"
        f"debt_queue_plus1_plus3:{MIXED_DEBT_QUEUE_1_3}:"
        f"conditional_naive_unit_holonomy:{MIXED_NAIVE_HOLONOMY}:"
        f"holonomy_determinant:{MIXED_NAIVE_HOLONOMY_DETERMINANT}:"
        f"holonomy_order:{MIXED_NAIVE_HOLONOMY_ORDER}:"
        f"delayed_margins:{MIXED_DELAYED_MARGINS}:"
        f"present_margins:{MIXED_PRESENT_MARGINS}:speed14_resonant=False"
    )
    print(
        f"complete_period4_BABA_lattice=denominator:{PERIOD4_DENOMINATOR}:"
        f"indices:{BABA_PERIOD4_INDICES}:"
        f"tails:{tuple(Fraction(index, PERIOD4_DENOMINATOR) for index in BABA_PERIOD4_INDICES)}:"
        f"with_b_vertex_guard_safe:{BABA_WITH_B_GUARD_SAFE}:"
        f"with_a_vertex_unshifted_c1_safe:{BABA_WITH_A_C1_SAFE}"
    )
    print("verdict=homogeneous repeated targets are nilpotent, but intrinsic "
          "A_debt/B_debt sectors open an exact period-four recurrent SCC of "
          "the new guard-cospan edge automaton")
    print("next_carrier=use the SCC as a signed mapping-cone/Cech guard-debt "
          "control or change chronology; literal positive endpoint attachment "
          "is not supplied and naive unit gluing has order-168 holonomy")
    print("SCOPE: recurrent SCC only in the newly defined debt automaton; a,b "
          "are delayed-tooth labels, not terminal words; fixed-PAT_QB rails; "
          "no unit transport, semantic endpoint transport, row exclusion, or "
          "LRC14 conclusion")


if __name__ == "__main__":
    main()
