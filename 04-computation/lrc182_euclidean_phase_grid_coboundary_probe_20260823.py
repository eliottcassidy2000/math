#!/usr/bin/env python3
"""Exact Euclidean-address / coefficient-cocycle probe on THM-3710.

This is deliberately scoped to the thirty typed THM-3672 control charts used
by THM-3710.  It does not assert that a genuine cover has this skeleton.

For each of the two 13-cosets

    theta_(sigma,j) = (14*j+sigma)/182 (mod 1),

we retain the canonical continued-fraction word (the THM-778 Euclidean
address), the transported midpoint phase trace, and the integer cyclic
difference cocycle of the signed endpoint multiplicity.  One integer
basepoint per sign fibre reconstructs the coefficient profile exactly.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
import importlib.util
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
THM3710_PATH = ROOT / "04-computation/lrc_successor_endpoint_182_bad_prime_thm3710.py"
P = 13


def require(condition: bool, payload: object) -> None:
    if condition is not True:
        raise RuntimeError(payload)


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("cannot load", path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def canonical_cf(value: Fraction) -> tuple[int, ...]:
    """Canonical digits after the initial zero for 0<value<1."""

    require(0 < value < 1, value)
    numerator, denominator = value.numerator, value.denominator
    digits: list[int] = []
    while numerator:
        digit, remainder = divmod(denominator, numerator)
        digits.append(digit)
        denominator, numerator = numerator, remainder
    require(digits[-1] > 1 or len(digits) == 1, digits)
    return tuple(digits)


def fraction_from_cf(digits: tuple[int, ...]) -> Fraction:
    value = Fraction()
    for digit in reversed(digits):
        value = Fraction(1, digit + value)
    return value


def midpoint_trace(digits: tuple[int, ...]) -> tuple[int, ...]:
    """THM-778's s'=(a-s) mod 2, with s_0=1."""

    phase = 1
    trace = [phase]
    for digit in digits:
        phase = (digit - phase) % 2
        trace.append(phase)
    return tuple(trace)


def cyclic_difference(values: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(values[(index + 1) % len(values)] - values[index] for index in range(len(values)))


def reconstruct_from_base(basepoint: int, differences: tuple[int, ...]) -> tuple[int, ...]:
    values = [basepoint]
    for difference in differences[:-1]:
        values.append(values[-1] + difference)
    require(values[0] - values[-1] == differences[-1], (basepoint, differences, values))
    return tuple(values)


def exact_mask(values: tuple[int, ...], sign: int) -> tuple[int, ...]:
    """Delete the denominator-14 corner from one full 13-coset."""

    corner = 12 if sign == 1 else 1
    return tuple(0 if index == corner else value for index, value in enumerate(values))


def semantic_hash(lines: list[str]) -> str:
    return sha256(("\n".join(lines) + "\n").encode("utf-8")).hexdigest()


thm3710 = load_module(THM3710_PATH, "thm3710_bridge_parent")
source = thm3710.load_source()
ref = source.load_referee()

require(ref.R == 169, ref.R)
require(ref.W[ref.TARGET_B] == 2 * P**5, ref.W)
successor_speed = P * ref.W[ref.TARGET_B]
period = ref.NN // successor_speed
require((successor_speed, period) == (9653618, 5214048840), (successor_speed, period))
require(period % 182 == 0, period)
period_scale = period // 182

# Rebuild exactly the endpoint counters used by THM-3710.
zero = (0,) * 9
word = ref.build_boolean_set(ref.PATTERN_QA, zero)
interval_sets = {
    "zero": source.marked_intervals(
        ref,
        ref.build_boolean_set(ref.PATTERN_E, zero),
        word,
    )
}
for graft in source.UNITS:
    interval_sets[f"a{graft}"] = source.marked_intervals(
        ref,
        ref.build_boolean_set(
            ref.PATTERN_E,
            source.negative_target_dipole(ref.TARGET_A, graft),
        ),
        word,
    )
    interval_sets[f"b{graft}"] = source.marked_intervals(
        ref,
        ref.build_boolean_set(
            ref.PATTERN_E,
            source.negative_target_dipole(ref.TARGET_B, graft),
        ),
        word,
    )

labels = (
    ("zero",)
    + tuple(f"a{k}" for k in source.UNITS)
    + tuple(f"b{k}" for k in source.UNITS)
)
counters: dict[str, dict[int, int]] = {}
for label in labels:
    counter: defaultdict[int, int] = defaultdict(int)
    thm3710.add_endpoints(counter, interval_sets[label], 1, period)
    counters[label] = {residue: coefficient for residue, coefficient in counter.items() if coefficient}


def chart_counter(k: int, ell: int) -> dict[int, int]:
    counter: defaultdict[int, int] = defaultdict(int)
    for residue, coefficient in counters["zero"].items():
        counter[residue] += coefficient
    for residue, coefficient in counters[f"a{k}"].items():
        counter[residue] += coefficient
    for residue, coefficient in counters[f"b{ell}"].items():
        counter[residue] -= 2 * coefficient
    return {residue: coefficient for residue, coefficient in counter.items() if coefficient}


charts = tuple((k, ell) for k in source.UNITS for ell in source.UNITS if k != ell)
require(len(charts) == 30, charts)
chart_counters = {chart: chart_counter(*chart) for chart in charts}
fibres = {
    (k, ell, sign): thm3710.kgrid_fiber(chart_counters[(k, ell)], period, sign)
    for k, ell in charts
    for sign in (1, -1)
}

# The exact ordered Euclidean address on all 26 full-grid vertices.
phase_rows: list[tuple[int, int, Fraction, tuple[int, ...], tuple[int, ...]]] = []
for sign in (1, -1):
    for index in range(P):
        numerator = (14 * index + sign) % 182
        phase = Fraction(numerator, 182)
        digits = canonical_cf(phase)
        trace = midpoint_trace(digits)
        require(fraction_from_cf(digits) == phase, (phase, digits))
        phase_rows.append((sign, index, phase, digits, trace))

require(len({row[2] for row in phase_rows}) == 26, "phase collision")
require(len({row[3] for row in phase_rows}) == 26, "ordered word collision")
corner_rows = tuple((sign, index, phase) for sign, index, phase, *_ in phase_rows if phase.denominator == 14)
require(corner_rows == ((1, 12, Fraction(13, 14)), (-1, 1, Fraction(1, 14))), corner_rows)

# Exact reconstruction from one basepoint and the integer edge cocycle.
reconstruction_checks = 0
cycle_checks = 0
hasse_checks = 0
for key, values in fibres.items():
    differences = cyclic_difference(values)
    require(sum(differences) == 0, (key, differences))
    require(reconstruct_from_base(values[0], differences) == values, key)
    sign = key[2]
    reconstructed_exact = exact_mask(reconstruct_from_base(values[0], differences), sign)
    require(reconstructed_exact == exact_mask(values, sign), key)
    require(
        thm3710.radical_order(tuple(value % P for value in reconstructed_exact))
        == thm3710.radical_order(tuple(value % P for value in exact_mask(values, sign))),
        key,
    )
    reconstruction_checks += P
    cycle_checks += 1
    hasse_checks += P


def exact_fibre_residual(values: tuple[int, ...], sign: int) -> Fraction:
    masked = exact_mask(values, sign)
    return sum(
        coefficient
        * thm3710.primitive(((14 * index + sign) % 182) * period_scale, period)
        for index, coefficient in enumerate(masked)
    ) / successor_speed


def direct_exact_residual(counter: dict[int, int]) -> Fraction:
    return sum(
        coefficient * thm3710.primitive(residue, period)
        for residue, coefficient in counter.items()
        if Fraction(residue, period).denominator == 182
    ) / successor_speed


chart_exact_residuals: dict[tuple[int, int], Fraction] = {}
chart_full_residuals: dict[tuple[int, int], Fraction] = {}
exact_orders: list[tuple[int, int, int, int]] = []
for chart in charts:
    reconstructed = sum(
        exact_fibre_residual(fibres[chart + (sign,)], sign)
        for sign in (1, -1)
    )
    direct = direct_exact_residual(chart_counters[chart])
    require(reconstructed == direct, (chart, reconstructed, direct))
    require(direct > 0, (chart, direct))
    chart_exact_residuals[chart] = direct
    chart_full_residuals[chart] = thm3710.endpoint_residual(
        chart_counters[chart], period, successor_speed
    )
    for sign in (1, -1):
        order = thm3710.radical_order(
            tuple(value % P for value in exact_mask(fibres[chart + (sign,)], sign))
        )
        exact_orders.append((*chart, sign, order))

order_histogram = {
    order: sum(row[3] == order for row in exact_orders)
    for order in sorted({row[3] for row in exact_orders})
}
require(order_histogram == {0: 57, 2: 3}, order_histogram)

# The 30 charts collapse to exactly 13 endpoint-profile classes.  This is the
# quotient {0,1,2}->A, with 3,4,5 retained separately, on each chart slot.
def graft_class(graft: int) -> str:
    return "A" if graft in (0, 1, 2) else str(graft)


integer_profile_groups: defaultdict[
    tuple[tuple[int, ...], tuple[int, ...]], list[tuple[int, int]]
] = defaultdict(list)
for chart in charts:
    integer_profile_groups[
        tuple(fibres[chart + (sign,)] for sign in (1, -1))
    ].append(chart)
require(len(integer_profile_groups) == 13, len(integer_profile_groups))

profile_classes = tuple(
    sorted((tuple(group) for group in integer_profile_groups.values()), key=lambda group: group[0])
)
profile_by_graft_class: dict[tuple[str, str], tuple[tuple[int, ...], tuple[int, ...]]] = {}
for profile, grouped_charts in integer_profile_groups.items():
    quotient_keys = {(graft_class(k), graft_class(ell)) for k, ell in grouped_charts}
    require(len(quotient_keys) == 1, (grouped_charts, quotient_keys))
    quotient_key = next(iter(quotient_keys))
    require(quotient_key not in profile_by_graft_class, quotient_key)
    profile_by_graft_class[quotient_key] = profile
    first_counter = chart_counters[grouped_charts[0]]
    require(
        all(chart_counters[chart] == first_counter for chart in grouped_charts),
        ("full counter does not factor", grouped_charts),
    )
    require(
        all(chart_exact_residuals[chart] == chart_exact_residuals[grouped_charts[0]] for chart in grouped_charts),
        ("exact residual does not factor", grouped_charts),
    )
    require(
        all(chart_full_residuals[chart] == chart_full_residuals[grouped_charts[0]] for chart in grouped_charts),
        ("full residual does not factor", grouped_charts),
    )
require(len(profile_by_graft_class) == 13, profile_by_graft_class)

# The quotient forgets absolute interval lifts even though it preserves every
# endpoint residue and all residuals.  Pin one literal period-translation.
a0_sorted = sorted(interval_sets["a0"])
a1_sorted = sorted(interval_sets["a1"])
lift_hostile = next(
    (first, second)
    for first, second in zip(a0_sorted, a1_sorted)
    if first != second
    and first[0] % period == second[0] % period
    and first[1] % period == second[1] % period
)
lift_translation = lift_hostile[1][0] - lift_hostile[0][0]
require(
    lift_hostile[1][1] - lift_hostile[0][1] == lift_translation
    and lift_translation % period == 0,
    (lift_hostile, lift_translation),
)
lift_translation_periods = lift_translation // period

# Hostile 1: an edge cocycle without a basepoint cannot preserve the
# corner-deleted Hasse order.  Use THM-3710's exceptional (4,0,minus) fibre.
base_hostile_key = (4, 0, -1)
base_original = fibres[base_hostile_key]
base_shifted = tuple(value + 1 for value in base_original)
require(cyclic_difference(base_original) == cyclic_difference(base_shifted), "shift changes cocycle")
base_original_order = thm3710.radical_order(tuple(value % P for value in exact_mask(base_original, -1)))
base_shifted_order = thm3710.radical_order(tuple(value % P for value in exact_mask(base_shifted, -1)))
require((base_original_order, base_shifted_order) == (2, 0), (base_original_order, base_shifted_order))
base_residual_delta = exact_fibre_residual(base_shifted, -1) - exact_fibre_residual(base_original, -1)
require(base_residual_delta != 0, base_residual_delta)

# Hostile 2: no chart-independent additive digit cocycle with only one scalar
# chart basepoint can produce the actual fibres.  Such a model would make all
# within-fibre coefficient differences independent of the chart.
additive_hostile = None
for sign in (1, -1):
    for first_chart, second_chart in combinations(charts, 2):
        first = fibres[first_chart + (sign,)]
        second = fibres[second_chart + (sign,)]
        for left, right in combinations(range(P), 2):
            first_difference = first[right] - first[left]
            second_difference = second[right] - second[left]
            if first_difference != second_difference:
                additive_hostile = (
                    sign,
                    first_chart,
                    second_chart,
                    left,
                    right,
                    first_difference,
                    second_difference,
                )
                break
        if additive_hostile is not None:
            break
    if additive_hostile is not None:
        break
require(additive_hostile is not None, "unexpected chart-independent additive cocycle")

word_only_hostile = (
    Fraction(1, 182),
    (0, 1),
    fibres[(0, 1, 1)][0],
    (0, 3),
    fibres[(0, 3, 1)][0],
)
require(word_only_hostile[2] != word_only_hostile[4], word_only_hostile)

# Hostile 3: characteristic-13 profiles do not retain the rational endpoint
# residual.  Seek an actual pair of typed charts, not a synthetic lift.
mod_profile_groups: defaultdict[tuple[tuple[int, ...], tuple[int, ...]], list[tuple[int, int]]] = defaultdict(list)
for chart in charts:
    profile = tuple(
        tuple(value % P for value in fibres[chart + (sign,)])
        for sign in (1, -1)
    )
    mod_profile_groups[profile].append(chart)

mod13_hostile = None
for grouped_charts in mod_profile_groups.values():
    for first_chart, second_chart in combinations(grouped_charts, 2):
        if chart_exact_residuals[first_chart] != chart_exact_residuals[second_chart]:
            mod13_hostile = (
                first_chart,
                second_chart,
                chart_exact_residuals[first_chart],
                chart_exact_residuals[second_chart],
            )
            break
    if mod13_hostile is not None:
        break

# In this fixed bank the complete mod-13 profile happens to be injective on
# the 13 realized integer profile classes.  This is a finite lookup fact, not
# an all-ledger theorem.  A one-coordinate +13 lift is the minimal ambient
# hostile to exporting it.
require(mod13_hostile is None, mod13_hostile)
mod13_profile_bank_lookup_injective = len(
    {
        tuple(tuple(value % P for value in sign_fibre) for sign_fibre in profile)
        for profile in integer_profile_groups
    }
) == len(integer_profile_groups)
require(mod13_profile_bank_lookup_injective, "mod-13 bank lookup collision")

mod13_lift_original = fibres[(1, 2, 1)]
mod13_lift_shifted = tuple(
    value + (P if index == 0 else -P if index == 1 else 0)
    for index, value in enumerate(mod13_lift_original)
)
require(
    tuple(value % P for value in mod13_lift_original)
    == tuple(value % P for value in mod13_lift_shifted),
    "mod-13 lift changed residue profile",
)
require(sum(mod13_lift_original) == sum(mod13_lift_shifted), "lift changed augmentation")
mod13_lift_residual_delta = (
    exact_fibre_residual(mod13_lift_shifted, 1)
    - exact_fibre_residual(mod13_lift_original, 1)
)
require(mod13_lift_residual_delta != 0, mod13_lift_residual_delta)

# Ordered-word positive control / orderless-word hostile from the reflection.
phase_lookup = {phase: (sign, index, digits, trace) for sign, index, phase, digits, trace in phase_rows}
first_phase = Fraction(43, 182)
second_phase = Fraction(55, 182)
first_address = phase_lookup[first_phase]
second_address = phase_lookup[second_phase]
require(first_address[2] == (4, 4, 3, 3), first_address)
require(second_address[2] == (3, 3, 4, 4), second_address)
require(sorted(first_address[2]) == sorted(second_address[2]), "digit multisets differ")
require(midpoint_trace(first_address[2])[-1] == midpoint_trace(second_address[2])[-1] == 1, "phase mismatch")
canonical_chart = (1, 2)
first_coefficient = fibres[canonical_chart + (first_address[0],)][first_address[1]]
second_coefficient = fibres[canonical_chart + (second_address[0],)][second_address[1]]
require(first_coefficient != second_coefficient, (first_coefficient, second_coefficient))

# Literal Euclidean stripping leaves the 182 phase skeleton, so THM-3710 has
# no coefficient fibre at the target.  This is the typing obstruction to
# calling the grid-difference cocycle an Euclidean-descent cocycle.
def first_euclidean_step(value: Fraction) -> tuple[int, Fraction]:
    quotient, remainder = divmod(value.denominator, value.numerator)
    return quotient, Fraction(remainder, value.numerator)


phase_set = {row[2] for row in phase_rows}
euclidean_exit_controls = (
    (first_phase, first_euclidean_step(first_phase)),
    (second_phase, first_euclidean_step(second_phase)),
)
require(
    euclidean_exit_controls
    == (
        (Fraction(43, 182), (4, Fraction(10, 43))),
        (Fraction(55, 182), (3, Fraction(17, 55))),
    ),
    euclidean_exit_controls,
)
require(
    all(step[1] not in phase_set for _, step in euclidean_exit_controls),
    euclidean_exit_controls,
)
require(
    all(
        all(
            not any(Fraction(residue, period) == step[1] for residue in counter)
            for counter in chart_counters.values()
        )
        for _, step in euclidean_exit_controls
    ),
    ("Euclidean target unexpectedly present in a full chart counter", euclidean_exit_controls),
)

integer_profile_count = len(
    {
        tuple(fibres[chart + (sign,)] for sign in (1, -1))
        for chart in charts
    }
)
mod13_profile_count = len(mod_profile_groups)
require(mod13_profile_count == 13, mod13_profile_count)

integer_delta_groups: defaultdict[tuple[int, ...], list[tuple[tuple[int, int, int], int]]] = defaultdict(list)
for key, values in fibres.items():
    integer_delta_groups[cyclic_difference(values)].append((key, values[0]))
integer_sign_profile_count = len(set(fibres.values()))
require(integer_sign_profile_count == 26, integer_sign_profile_count)
require(len(integer_delta_groups) == integer_sign_profile_count, len(integer_delta_groups))
require(
    all(len({basepoint for _, basepoint in group}) == 1 for group in integer_delta_groups.values()),
    "realized delta class has multiple basepoints",
)

summary_lines = [
    "lrc182 ordered Euclidean-word / coefficient-cocycle probe",
    f"parent_thm3710_sha256={sha256(THM3710_PATH.read_bytes()).hexdigest()}",
    f"typed_chart_count={len(charts)}",
    f"sign_fibre_count={len(fibres)}",
    f"full_grid_phase_count={len(phase_rows)}",
    f"primitive_phase_count={sum(row[2].denominator == 182 for row in phase_rows)}",
    f"corner_rows={corner_rows}",
    f"ordered_word_count={len({row[3] for row in phase_rows})}",
    f"euclidean_digit_step_count={sum(len(row[3]) for row in phase_rows)}",
    f"euclidean_depth_range={min(len(row[3]) for row in phase_rows)}..{max(len(row[3]) for row in phase_rows)}",
    f"reconstruction_checks={reconstruction_checks}",
    f"cycle_closure_checks={cycle_checks}",
    f"hasse_reconstruction_checks={hasse_checks}",
    f"exact_residual_reconstruction_checks={len(charts)}",
    f"exact_order_histogram={order_histogram}",
    f"integer_profile_count={integer_profile_count}",
    f"mod13_profile_count={mod13_profile_count}",
    f"profile_classes={profile_classes}",
    f"profile_graft_quotient_keys={tuple(sorted(profile_by_graft_class))}",
    f"integer_sign_profile_count={integer_sign_profile_count}",
    f"integer_delta_class_count={len(integer_delta_groups)}",
    f"realized_delta_lookup_recovers_basepoint=True",
    f"basepoint_hostile_key={base_hostile_key}",
    f"basepoint_hostile_orders={base_original_order}->{base_shifted_order}",
    f"basepoint_hostile_exact_residual_delta={base_residual_delta}",
    f"additive_digit_cocycle_hostile={additive_hostile}",
    f"word_only_chart_hostile={word_only_hostile}",
    f"mod13_profile_bank_lookup_injective={mod13_profile_bank_lookup_injective}",
    f"ambient_mod13_lift_residual_delta={mod13_lift_residual_delta}",
    f"absolute_lift_hostile={lift_hostile}",
    f"absolute_lift_translation_periods={lift_translation_periods}",
    "orderless_word_hostile=(43/182,(4,4,3,3),55/182,(3,3,4,4))",
    f"orderless_word_hostile_canonical_coefficients={first_coefficient},{second_coefficient}",
    f"euclidean_stripping_exit_controls={euclidean_exit_controls}",
    "euclidean_stripping_targets_in_full_chart_counters=0,0",
    "scope=THM-3710 typed non-cover control only; no owner/root/arrival semantics",
]

for line in summary_lines:
    print(line)
print(f"semantic_sha256={semantic_hash(summary_lines)}")
print("ALL CHECKS PASSED;ARTIFACT=LRC182_EUCLIDEAN_PHASE_GRID_COBOUNDARY")
