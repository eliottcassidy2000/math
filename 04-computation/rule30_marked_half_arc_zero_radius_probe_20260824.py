#!/usr/bin/env python3
"""Exact controls for the Rule 30 marked half-arc and zero-radius probe.

This is an unnumbered, deterministic, dependency-free companion.  It checks
the physical terminal-line identities used in the accompanying research note
by two independent Rule 30 implementations through DIRECT_HORIZON.  It then
enumerates the nearest backward physical zero for the single-seed evolution
on the explicit finite universe 2 <= k <= HORIZON.

The finite enumeration is not a density theorem and proves none of the Rule
30 prizes.  In particular, the observed geometric-looking histogram is only
a finite signal, while the first radius-nine witness is an exact hostile to a
universal radius-eight bound.
"""

from collections import Counter, deque
from hashlib import sha256
import json
from pathlib import Path


DIRECT_HORIZON = 512
HORIZON = (1 << 18) - 1
SEARCH_CAP = 64
DYADIC_EXPONENTS = tuple(range(8, 19))
HAAR_EXHAUSTIVE_RADIUS = 8


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_bytes(path):
    return Path(path).read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")


def lf_sha256(path):
    return sha256(lf_bytes(path)).hexdigest()


def next_packed(row):
    """Rule 30 with bit j at physical site x=t-j in the time-t row."""
    return row ^ ((row << 1) | (row << 2))


def packed_rows(horizon):
    rows = [1]
    for _ in range(horizon):
        rows.append(next_packed(rows[-1]))
    return tuple(rows)


def sparse_rows(horizon):
    """Independent literal update on sets of physical integer sites."""
    black = {0}
    rows = []
    for time in range(horizon + 1):
        rows.append(frozenset(black))
        black = {
            site
            for site in range(-time - 1, time + 2)
            if ((site - 1 in black) ^ ((site in black) or (site + 1 in black)))
        }
    return tuple(rows)


def right_characteristic_all_ones(word, radius):
    """Literal Rule 30 on a finite word, read on the right-moving ray."""
    offset = radius
    row = word << offset
    mask = (1 << (4 * radius + 1)) - 1
    for step in range(radius):
        if ((row >> (offset + step)) & 1) == 0:
            return False
        row = ((row << 1) ^ (row | (row >> 1))) & mask
    return True


def packed_black_sites(row, time):
    return frozenset(
        time - bit
        for bit in range(row.bit_length())
        if (row >> bit) & 1
    )


def terminal(rows, depth, phase):
    """T_k(h)=a_(k+h)(h); packed exponent is identically k."""
    time = depth + phase
    require(0 <= time < len(rows), ("terminal time outside stored rows", depth, phase))
    return (rows[time] >> depth) & 1


def current_by_difference(rows, depth, phase):
    return terminal(rows, depth, phase + 1) ^ terminal(rows, depth, phase)


def current_by_parents(rows, depth, phase):
    """Q_k(h) as the adjacent-parent OR in the literal Rule 30 update."""
    time = depth + phase
    row = rows[time]
    return ((row >> (depth - 1)) & 1) | ((row >> (depth - 2)) & 1)


repository_root = Path(__file__).resolve().parents[1]
script_hash = lf_sha256(Path(__file__))

# Independent direct controls on complete finite Rule 30 rows.
direct_packed = packed_rows(DIRECT_HORIZON)
direct_sparse = sparse_rows(DIRECT_HORIZON)
full_row_checks = 0
for time, row in enumerate(direct_packed):
    require(
        packed_black_sites(row, time) == direct_sparse[time],
        ("packed and sparse Rule 30 rows disagree", time),
    )
    require(row.bit_length() <= 2 * time + 1, ("light-cone support failed", time))
    require((row >> (2 * time)) & 1, ("extreme-left ray failed", time))
    full_row_checks += 1

prefix20 = "".join(str((direct_packed[time] >> time) & 1) for time in range(20))
require(prefix20 == "11011100110001011001", ("center prefix changed", prefix20))

# On a Haar-iid base row, an r-cell all-one right characteristic has exactly
# one compatible (2r-1)-bit dependency word: 1 followed by 2r-2 zeros.  This
# is exhaustively checked here; the all-r induction is given in the note.
haar_characteristic_rows = []
haar_words_tested = 0
for radius in range(1, HAAR_EXHAUSTIVE_RADIUS + 1):
    width = 2 * radius - 1
    survivors = []
    for word in range(1 << width):
        haar_words_tested += 1
        if right_characteristic_all_ones(word, radius):
            survivors.append(word)
    require(survivors == [1], ("Haar characteristic preimage changed", radius, survivors))
    haar_characteristic_rows.append((radius, width, len(survivors), survivors[0]))

# Direct current, telescope, universal half-arc, and sharp boundary controls.
current_checks = 0
full_arc_checks = 0
half_arc_checks = 0
sharp_even_boundary_checks = 0
marked_cylinder_checks = 0
direct_zero_radii = {}
for depth in range(2, DIRECT_HORIZON + 1):
    for phase in range(-depth, 0):
        q_difference = current_by_difference(direct_packed, depth, phase)
        q_parents = current_by_parents(direct_packed, depth, phase)
        require(q_difference == q_parents, ("current implementations disagree", depth, phase))
        current_checks += 1

    center = terminal(direct_packed, depth, 0)
    full_arc = 0
    for phase in range(-depth, 0):
        full_arc ^= current_by_parents(direct_packed, depth, phase)
    require(full_arc == center, ("full terminal arc failed", depth))
    full_arc_checks += 1

    half_length = depth // 2 + 1
    require(
        terminal(direct_packed, depth, -half_length) == 0,
        ("universal midpoint-exterior zero failed", depth, half_length),
    )
    half_arc = 0
    for phase in range(-half_length, 0):
        half_arc ^= current_by_parents(direct_packed, depth, phase)
    require(half_arc == center, ("shortened half arc failed", depth, half_length))
    half_arc_checks += 1

    if depth % 2 == 0:
        boundary_length = depth // 2
        require(
            terminal(direct_packed, depth, -boundary_length) == 1,
            ("extreme-left midpoint boundary changed", depth),
        )
        boundary_arc = 0
        for phase in range(-boundary_length, 0):
            boundary_arc ^= current_by_parents(direct_packed, depth, phase)
        require(boundary_arc == (center ^ 1), ("sharp boundary correction failed", depth))
        sharp_even_boundary_checks += 1

    zero_radius = next(
        radius
        for radius in range(1, half_length + 1)
        if terminal(direct_packed, depth, -radius) == 0
    )
    nearest_arc = 0
    for phase in range(-zero_radius, 0):
        nearest_arc ^= current_by_parents(direct_packed, depth, phase)
    require(nearest_arc == center, ("nearest-zero arc failed", depth, zero_radius))
    direct_zero_radii[depth] = zero_radius

    terminal_run_is_one = True
    current_zero_run = terminal(direct_packed, depth, -1) == 1
    for radius in range(1, half_length + 1):
        terminal_run_is_one = terminal_run_is_one and (
            terminal(direct_packed, depth, -radius) == 1
        )
        if radius >= 2:
            current_zero_run = current_zero_run and (
                current_by_parents(direct_packed, depth, -radius) == 0
            )
        base_row = direct_packed[depth - radius]
        zero_width = 2 * radius - 2
        zero_mask = ((1 << zero_width) - 1) << (depth - 2 * radius + 2)
        earlier_cylinder = ((base_row >> depth) & 1) == 1 and (
            base_row & zero_mask
        ) == 0
        require(
            terminal_run_is_one == current_zero_run == earlier_cylinder,
            ("marked characteristic cylinder equivalence failed", depth, radius),
        )
        marked_cylinder_checks += 1

# Large finite exact enumeration.  Only SEARCH_CAP recent rows are retained;
# failure to find a zero within that cap is an explicit failure, not truncation.
recent_rows = deque([1], maxlen=SEARCH_CAP + 2)
current_row = 1
radius_histogram = Counter()
radius_by_center = Counter()
record_maxima = []
large_radius_witnesses = []
dyadic_histograms = []
maximum_radius = 0
for depth in range(1, HORIZON + 1):
    current_row = next_packed(current_row)
    recent_rows.append(current_row)
    if depth < 2:
        continue

    center = (current_row >> depth) & 1
    universal_length = depth // 2 + 1
    search_limit = min(SEARCH_CAP, universal_length)
    zero_radius = None
    for radius in range(1, search_limit + 1):
        earlier_row = recent_rows[-radius - 1]
        if ((earlier_row >> depth) & 1) == 0:
            zero_radius = radius
            break
    require(
        zero_radius is not None,
        ("nearest zero exceeded explicit search cap", depth, SEARCH_CAP),
    )
    require(zero_radius <= universal_length, ("universal zero bound failed", depth, zero_radius))
    if depth <= DIRECT_HORIZON:
        require(
            zero_radius == direct_zero_radii[depth],
            ("streaming and direct zero radii disagree", depth, zero_radius),
        )

    radius_histogram[zero_radius] += 1
    radius_by_center[(zero_radius, center)] += 1
    if zero_radius > maximum_radius:
        maximum_radius = zero_radius
        record_maxima.append((depth, zero_radius, center))
    if zero_radius >= 8:
        large_radius_witnesses.append((depth, zero_radius, center))
    if depth + 1 in {1 << exponent for exponent in DYADIC_EXPONENTS}:
        dyadic_histograms.append(
            (
                depth,
                tuple(sorted(radius_histogram.items())),
                maximum_radius,
            )
        )

require(sum(radius_histogram.values()) == HORIZON - 1, "finite universe count changed")
require(maximum_radius < SEARCH_CAP, "search cap was reached")
require(any(radius == 9 for _, radius, _ in large_radius_witnesses), "radius-nine hostile vanished")
first_radius_nine = next(item for item in large_radius_witnesses if item[1] == 9)
require(first_radius_nine == (79883, 9, 0), ("first radius-nine witness changed", first_radius_nine))

# Compare the deterministic temporal histogram with the exact Haar law:
# p(1)=1/2 and p(r)=3/2^(2r-1) for r>=2.  The scaled discrepancies are exact
# finite comparisons; they are not a test of asymptotic convergence.
universe_size = HORIZON - 1
geometric_scaled_discrepancies = []
for radius, count in sorted(radius_histogram.items()):
    if radius == 1:
        numerator, denominator = 1, 2
    else:
        numerator, denominator = 3, 1 << (2 * radius - 1)
    scaled_discrepancy = denominator * count - numerator * universe_size
    geometric_scaled_discrepancies.append(
        (radius, numerator, denominator, count, scaled_discrepancy)
    )
require(
    any(row[-1] != 0 for row in geometric_scaled_discrepancies),
    "finite temporal histogram unexpectedly equals Haar proportions",
)

histogram = tuple(sorted(radius_histogram.items()))
center_split = tuple(
    (radius, radius_by_center[(radius, 0)], radius_by_center[(radius, 1)])
    for radius, _ in histogram
)
radius_nine_witnesses = tuple(item for item in large_radius_witnesses if item[1] == 9)

semantic_object = {
    "direct_horizon": DIRECT_HORIZON,
    "horizon": HORIZON,
    "prefix20": prefix20,
    "haar_exhaustive_radius": HAAR_EXHAUSTIVE_RADIUS,
    "haar_characteristic_rows": tuple(haar_characteristic_rows),
    "haar_words_tested": haar_words_tested,
    "full_row_checks": full_row_checks,
    "current_checks": current_checks,
    "full_arc_checks": full_arc_checks,
    "half_arc_checks": half_arc_checks,
    "sharp_even_boundary_checks": sharp_even_boundary_checks,
    "marked_cylinder_checks": marked_cylinder_checks,
    "search_cap": SEARCH_CAP,
    "histogram": histogram,
    "center_split": center_split,
    "record_maxima": tuple(record_maxima),
    "large_radius_witnesses": tuple(large_radius_witnesses),
    "dyadic_histograms": tuple(dyadic_histograms),
    "geometric_scaled_discrepancies": tuple(geometric_scaled_discrepancies),
}
semantic_hash = sha256(
    json.dumps(semantic_object, sort_keys=True, separators=(",", ":")).encode("ascii")
).hexdigest()

print("UNNUMBERED Rule30 marked half-arc and nearest-zero-radius probe")
print(f"script_sha256_lf={script_hash}")
print(f"rule30_center_prefix20={prefix20}")
print(f"independent_complete_rule30_rows_match_through={DIRECT_HORIZON}")
print(f"complete_row_checks={full_row_checks}")
print(f"haar_characteristic_unique_preimage_checks={tuple(haar_characteristic_rows)}")
print(f"haar_characteristic_words_exhausted={haar_words_tested}")
print("haar_tail_law=P(Z>r)=2^{-(2r-1)}")
print("haar_mass_law=P(Z=1)=1/2;P(Z=r)=3*2^{-(2r-1)}_for_r>=2")
print("haar_center_given_radius=P(c=1|Z=1)=3/4;P(c=1|Z=r)=1/4_for_r>=2")
print("haar_joint_law=P(c=1,Z=1)=3/8;P(c=1,Z=r)=3*2^{-(2r+1)}_for_r>=2")
print(f"direct_transition_current_checks={current_checks}")
print(f"direct_full_terminal_arc_checks={full_arc_checks}")
print(f"direct_shortened_half_arc_checks={half_arc_checks}")
print(f"direct_sharp_even_midpoint_boundary_checks={sharp_even_boundary_checks}")
print(f"direct_marked_current_zero_run_and_earlier_cylinder_checks={marked_cylinder_checks}")
print("universal_identity_target=c_k=XOR_{h=-floor(k/2)-1}^{-1} Q_k(h)")
print("universal_boundary_target=for_even_k_T_k(-k/2)=1")
print(f"finite_universe=single_seed_rule30;2<=k<={HORIZON}")
print(f"finite_universe_size={universe_size}")
print(f"nearest_zero_search_cap={SEARCH_CAP}")
print("nearest_zero_unresolved_count=0")
print(f"nearest_zero_histogram={histogram}")
print(f"nearest_zero_center_split_radius_zero_one={center_split}")
print(f"nearest_zero_record_maxima={tuple(record_maxima)}")
print(f"nearest_zero_radius_at_least_8_witnesses={tuple(large_radius_witnesses)}")
print(f"nearest_zero_radius_9_witnesses={radius_nine_witnesses}")
print(f"first_radius_9_witness={first_radius_nine}")
print("universal_radius_at_most_8=REFUTED_FINITE_EXACT")
print(f"dyadic_prefix_histograms={tuple(dyadic_histograms)}")
print(f"single_seed_vs_haar_scaled_discrepancies={tuple(geometric_scaled_discrepancies)}")
print("single_seed_finite_temporal_histogram_equals_haar_proportions=NO")
print(f"semantic_sha256={semantic_hash}")
print("status=UNNUMBERED;UNIVERSAL_IDENTITIES_DIRECTLY_CHECKED_FINITE;RADIUS_DATA_FINITE_EXACT")
print("rule30_density_prize=OPEN")
print("ALL CHECKS PASSED")
