#!/usr/bin/env python3
"""Finite-exact companions for THM-3713's deep-offset detector.

The first control is a bounded diagonal-zero tensor with constant successor
marginal and positive target drift.  The second retains the deep-offset
colour inside all thirty typed THM-3672 charts.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys


sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
P = 13
SOURCE = ROOT / "04-computation/lrc_successor_mass_all_pair_swap_control_thm3672.py"
REFEREE = ROOT / "04-computation/lrc14_169_twist_two_twist_referee_thm2334.py"

PINS = {
    "04-computation/lrc_successor_mass_all_pair_swap_control_thm3672.py":
        "a191f934b20494b98e878fc5504a328c47d4cb6a92ea332246f782fac80c01c8",
    "05-knowledge/results/lrc_successor_mass_all_pair_swap_control_thm3672.out":
        "0d72ff662c908d9d9934976edce67895135c17a1eb0ba409e0abc382563b52c6",
    "04-computation/lrc14_169_twist_two_twist_referee_thm2334.py":
        "0e4a9e181263647e13d2a6738b6996c45df901d9d2b37d4d589dfddfbdd91480",
}

# Filled after the first independent exact reconstruction, then checked on
# every run.  The digest covers all thirty aggregates and all 390 coloured
# defects, not merely their sign pattern.
EXPECTED_TYPED_DIGEST = "fee602c81a407062806886a4d8d2945e86b706e99b8a9d90f6275cbe83ab5a27"


def require(condition, payload):
    if condition is not True:
        raise RuntimeError(payload)


def digest(value):
    encoded = json.dumps(value, separators=(",", ":")).encode("ascii")
    return sha256(encoded).hexdigest()


def dft_mod_79(values):
    """Exact C_13 transform in F_79; 8 has multiplicative order 13."""
    require(pow(8, P, 79) == 1 and 8 % 79 != 1, "bad order-13 root")
    return tuple(
        sum(value * pow(8, frequency * u, 79) for u, value in enumerate(values)) % 79
        for frequency in range(P)
    )


def balanced_first_moment(values):
    lifts = tuple(u if u <= 6 else u - P for u in range(P))
    return sum(lifts[u] * value for u, value in enumerate(values))


def load_source():
    spec = importlib.util.spec_from_file_location("thm3672_source", SOURCE)
    require(spec is not None and spec.loader is not None, "cannot load THM-3672")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def hostile_f(s):
    s %= P
    if s == 0:
        return 1
    if s == 1:
        return -1
    return 0


def hostile_g(u):
    u %= P
    if u == 1:
        return 1
    if u == 2:
        return -1
    return 0


def hostile_H(r, s, t):
    u = (r - t) % P
    if u == 0:
        return 0
    return 2 + hostile_f(s) * hostile_g(u)


def hostile_h(u, s, t):
    return hostile_H((t + u) % P, s, t)


def hostile_successor(s, t):
    return sum(hostile_H(r, s, t) for r in range(P))


def hostile_orbit_defect(u, s, t):
    """The plus-oriented mask h+h(s+1)-2h(t+1)."""
    return (
        hostile_h(u, s, t)
        + hostile_h(u, s + 1, t)
        - 2 * hostile_h(u, s, t + 1)
    )


def hostile_successor_defect(s, t):
    return (
        hostile_successor(s, t)
        + hostile_successor(s + 1, t)
        - 2 * hostile_successor(s, t + 1)
    )


def audit_bounded_hostile():
    entries = [
        hostile_H(r, s, t)
        for r in range(P)
        for s in range(P)
        for t in range(P)
    ]
    require((min(entries), max(entries)) == (0, 3), (min(entries), max(entries)))
    require(
        all(hostile_H(t, s, t) == 0 for s in range(P) for t in range(P)),
        "hostile diagonal zero",
    )

    successor_values = {
        hostile_successor(s, t) for s in range(P) for t in range(P)
    }
    require(successor_values == {24}, successor_values)
    successor_mean = Fraction(
        sum(hostile_successor(s, t) for s in range(P) for t in range(P)),
        P * P,
    )
    successor_variance = Fraction(
        sum(
            (Fraction(hostile_successor(s, t)) - successor_mean) ** 2
            for s in range(P)
            for t in range(P)
        ),
        P * P,
    )
    require(successor_variance == 0, successor_variance)
    require(
        all(hostile_successor_defect(s, t) == 0 for s in range(P) for t in range(P)),
        "hostile successor defect",
    )

    drift_numerator = Fraction()
    for u in range(P):
        values = [hostile_h(u, s, t) for s in range(P) for t in range(P)]
        mean = Fraction(sum(values), P * P)
        for value in values:
            drift_numerator += (Fraction(value) - mean) ** 2
    drift = drift_numerator / P**3
    require(drift == Fraction(4, 169), drift)
    deep_target_energy = drift - successor_variance / P**2
    require(deep_target_energy == drift, (deep_target_energy, drift))

    defects = tuple(
        hostile_orbit_defect(u, s, t)
        for u in range(P)
        for s in range(P)
        for t in range(P)
    )
    nonzero = tuple(value for value in defects if value)
    require(len(nonzero) == 78, len(nonzero))
    require(set(nonzero) == {-2, -1, 1, 2}, set(nonzero))
    nonzero_profiles = tuple(
        tuple(hostile_orbit_defect(u, s, t) for u in range(P))
        for s in range(P)
        for t in range(P)
        if any(hostile_orbit_defect(u, s, t) for u in range(P))
    )
    require(len(nonzero_profiles) == 39, len(nonzero_profiles))
    require(all(profile[0] == 0 for profile in nonzero_profiles), "anchored zero")
    hostile_spectra = tuple(dft_mod_79(profile) for profile in nonzero_profiles)
    require(
        all(spectrum[0] == 0 and all(spectrum[a] for a in range(1, P)) for spectrum in hostile_spectra),
        "hostile deep-character saturation",
    )
    require(
        all(
            hostile_successor_defect(s, t)
            == sum(hostile_orbit_defect(u, s, t) for u in range(P))
            for s in range(P)
            for t in range(P)
        ),
        "deep-offset quotient identity",
    )

    edge_energy = Fraction(
        sum(
            (hostile_h(u, s + 1, t) - hostile_h(u, s, t)) ** 2
            + (hostile_h(u, s, t + 1) - hostile_h(u, s, t)) ** 2
            for u in range(P)
            for s in range(P)
            for t in range(P)
        ),
        P**3,
    )
    defect_energy = Fraction(sum(value * value for value in defects), P**3)
    require(edge_energy == Fraction(12, 169), edge_energy)
    require(defect_energy == edge_energy, (defect_energy, edge_energy))

    largest = max(abs(value) for value in nonzero)
    three_site_tariff = Fraction(largest**2, 6 * P**3)
    edge_tariff = Fraction(largest**2, 2 * P**3)
    require(three_site_tariff == Fraction(2, 6591), three_site_tariff)
    require(edge_tariff == Fraction(2, 2197), edge_tariff)
    require(drift >= edge_tariff > three_site_tariff, (drift, edge_tariff, three_site_tariff))

    return {
        "entry_range": (min(entries), max(entries)),
        "successor_value": next(iter(successor_values)),
        "successor_variance": successor_variance,
        "nonzero_defects": len(nonzero),
        "nonzero_profiles": len(nonzero_profiles),
        "defect_values": tuple(sorted(set(nonzero))),
        "drift": drift,
        "deep_target_energy": deep_target_energy,
        "edge_energy": edge_energy,
        "three_site_energy": defect_energy,
        "largest_defect": largest,
        "three_site_tariff": three_site_tariff,
        "edge_tariff": edge_tariff,
    }


def periodic_arc_prefix(x, period, start, length):
    """Measure of [0,x) in one periodic half-open arc."""
    quotient, residue = divmod(x, period)
    result = quotient * length
    end = start + length
    if end <= period:
        pieces = ((start, end),)
    else:
        pieces = ((start, period), (0, end - period))
    for left, right in pieces:
        result += max(0, min(residue, right) - left)
    return result


def audit_typed_control():
    for relative, expected in PINS.items():
        actual = sha256((ROOT / relative).read_bytes()).hexdigest()
        require(actual == expected, ("parent hash drift", relative, actual))

    source = load_source()
    ref = source.load_referee()
    require(ref.W[ref.TARGET_B] == 742586 and ref.NN == 50334435734703120, (ref.W, ref.NN))

    zero = (0,) * 9
    word = ref.build_boolean_set(ref.PATTERN_QA, zero)

    def marked(shift):
        return source.marked_intervals(
            ref,
            ref.build_boolean_set(ref.PATTERN_E, shift),
            word,
        )

    marked_zero = marked(zero)
    marked_a = {
        k: marked(source.negative_target_dipole(ref.TARGET_A, k))
        for k in source.UNITS
    }
    marked_b = {
        k: marked(source.negative_target_dipole(ref.TARGET_B, k))
        for k in source.UNITS
    }

    deep_speed = ref.W[ref.TARGET_B]
    period = ref.NN // deep_speed
    half_width = ref.NN // (14 * deep_speed)
    phase_unit = ref.NN // (P * deep_speed)
    require(period == P * phase_unit, (period, phase_unit))

    def colour_mass(intervals, r):
        center = (r * phase_unit) % period
        start = (center - half_width) % period
        length = 2 * half_width
        return sum(
            periodic_arc_prefix(right, period, start, length)
            - periodic_arc_prefix(left, period, start, length)
            for left, right in intervals
        )

    h_zero = tuple(colour_mass(marked_zero, r) for r in range(P))
    h_a = {
        k: tuple(colour_mass(marked_a[k], r) for r in range(P))
        for k in source.UNITS
    }
    h_b = {
        k: tuple(colour_mass(marked_b[k], r) for r in range(P))
        for k in source.UNITS
    }

    expected = dict(source.EXPECTED_VALUES)
    require(sum(h_zero) == expected["zero"][0], "zero successor sum")
    for k in source.UNITS:
        require(sum(h_a[k]) == expected[f"a{k}"][0], ("A successor sum", k))
        require(sum(h_b[k]) == expected[f"b{k}"][0], ("B successor sum", k))

    rows = []
    for k in source.UNITS:
        for ell in source.UNITS:
            if k == ell:
                continue
            # h_u(0,0)+h_u(1,0)-2h_u(0,1).  At t=1 the fixed
            # offset u=r-t samples the raw deep colour r=u+1.
            defects = tuple(
                h_zero[u] + h_a[k][u] - 2 * h_b[ell][(u + 1) % P]
                for u in range(P)
            )
            aggregate = sum(defects)
            expected_aggregate = (
                expected["zero"][0]
                + expected[f"a{k}"][0]
                - 2 * expected[f"b{ell}"][0]
            )
            require(aggregate == expected_aggregate, (k, ell, aggregate, expected_aggregate))
            rows.append((k, ell, aggregate, defects))

    require(len(rows) == 30, len(rows))
    row_digest = digest(rows)
    require(row_digest == EXPECTED_TYPED_DIGEST, ("typed digest", row_digest))

    chart_support = {
        tuple(u for u, value in enumerate(defects) if value)
        for _, _, _, defects in rows
    }
    require(chart_support == {(1, 2, 10, 11, 12)}, chart_support)
    for _, _, _, defects in rows:
        require(defects[0] == 0, defects)
        require(all(defects[u] > 0 for u in (1, 2)), defects)
        require(all(defects[u] < 0 for u in (10, 11, 12)), defects)
        require(all(defects[u] == 0 for u in (0, 3, 4, 5, 6, 7, 8, 9)), defects)

    typed_spectra = tuple(dft_mod_79(defects) for _, _, _, defects in rows)
    require(
        all(all(spectrum[a] for a in range(P)) for spectrum in typed_spectra),
        "typed deep-character saturation",
    )
    balanced_moments = tuple(
        balanced_first_moment(defects) for _, _, _, defects in rows
    )
    require(all(value > 0 for value in balanced_moments), balanced_moments)
    require(
        (min(balanced_moments), max(balanced_moments))
        == (430999196213378, 441421418113351),
        (min(balanced_moments), max(balanced_moments)),
    )

    a_edges = {
        k: tuple(h_a[k][u] - h_zero[u] for u in range(P))
        for k in source.UNITS
    }
    b_edges = {
        ell: tuple(h_b[ell][(u + 1) % P] - h_zero[u] for u in range(P))
        for ell in source.UNITS
    }
    for values in a_edges.values():
        require(all(values[u] < 0 for u in (1, 2, 11, 12)), values)
        require(all(values[u] == 0 for u in (0, 3, 4, 5, 6, 7, 8, 9, 10)), values)
    for values in b_edges.values():
        require(all(values[u] < 0 for u in (1, 2)), values)
        require(all(values[u] > 0 for u in (10, 11, 12)), values)
        require(all(values[u] == 0 for u in (0, 3, 4, 5, 6, 7, 8, 9)), values)

    max_abs = max(
        abs(value)
        for _, _, _, defects in rows
        for value in defects
    )
    strongest_cases = tuple(
        (k, ell, u, value)
        for k, ell, _, defects in rows
        for u, value in enumerate(defects)
        if abs(value) == max_abs
    )
    require(
        strongest_cases
        == (
            (5, 0, 11, -66828200140260),
            (5, 1, 11, -66828200140260),
            (5, 2, 11, -66828200140260),
        ),
        strongest_cases,
    )
    strongest_k, strongest_ell, strongest_u, strongest_value = strongest_cases[-1]
    normalized_strongest = Fraction(strongest_value, ref.NN)
    drift_tariff = normalized_strongest**2 / (6 * P**3)
    deep_target_tariff = normalized_strongest**2 / (6 * P**4)

    return {
        "parent_sha256": PINS[str(SOURCE.relative_to(ROOT)).replace("\\", "/")],
        "deep_speed": deep_speed,
        "period": period,
        "half_width": half_width,
        "chart_count": len(rows),
        "support": next(iter(chart_support)),
        "positive_colours": (1, 2),
        "negative_colours": (10, 11, 12),
        "all_13_deep_characters_live_mod79": True,
        "balanced_first_moment_range": (min(balanced_moments), max(balanced_moments)),
        "a_negative_colours": (1, 2, 11, 12),
        "b_negative_colours": (1, 2),
        "b_positive_colours": (10, 11, 12),
        "strongest_cases": strongest_cases,
        "strongest": (strongest_k, strongest_ell, strongest_u, strongest_value),
        "normalized_strongest": normalized_strongest,
        "drift_tariff": drift_tariff,
        "deep_target_tariff": deep_target_tariff,
        "digest": row_digest,
    }


def main():
    hostile = audit_bounded_hostile()
    typed = audit_typed_control()

    print(f"p={P}")
    print(f"hostile_entry_range={hostile['entry_range']}")
    print("hostile_diagonal_zero=True")
    print(f"hostile_successor_value={hostile['successor_value']}")
    print(f"hostile_successor_variance={hostile['successor_variance']}")
    print("hostile_all_169_successor_defects_zero=True")
    print(f"hostile_orbitwise_nonzero_defects={hostile['nonzero_defects']}")
    print(f"hostile_nonzero_orbit_profiles={hostile['nonzero_profiles']}")
    print(f"hostile_orbitwise_defect_values={hostile['defect_values']}")
    print(f"hostile_target_drift={hostile['drift']}")
    print(f"hostile_deep_target_energy={hostile['deep_target_energy']}")
    print(f"hostile_edge_energy={hostile['edge_energy']}")
    print(f"hostile_three_site_energy={hostile['three_site_energy']}")
    print(f"hostile_largest_defect={hostile['largest_defect']}")
    print(f"hostile_edge_tariff={hostile['edge_tariff']}")
    print(f"hostile_three_site_tariff={hostile['three_site_tariff']}")
    print("hostile_all_12_nontrivial_deep_characters_live_mod79=True")
    print("successor_defect_is_deep_offset_sum=True")
    print(f"parent_sha256={typed['parent_sha256']}")
    print(f"typed_deep_speed={typed['deep_speed']}")
    print(f"typed_period={typed['period']}")
    print(f"typed_half_width={typed['half_width']}")
    print("typed_successor_colour_sum_identity=True")
    print(f"typed_chart_count={typed['chart_count']}")
    print(f"typed_common_nonzero_colours={typed['support']}")
    print(f"typed_positive_colours={typed['positive_colours']}")
    print(f"typed_negative_colours={typed['negative_colours']}")
    print(f"typed_all_13_deep_characters_live_mod79={typed['all_13_deep_characters_live_mod79']}")
    print(f"typed_balanced_first_moment_range={typed['balanced_first_moment_range']}")
    print(f"typed_A_negative_colours={typed['a_negative_colours']}")
    print(f"typed_B_negative_colours={typed['b_negative_colours']}")
    print(f"typed_B_positive_colours={typed['b_positive_colours']}")
    print(f"typed_strongest_cases={typed['strongest_cases']}")
    print(f"typed_strongest_k_l_u_numerator={typed['strongest']}")
    print(f"typed_strongest_reduced={typed['normalized_strongest']}")
    print(f"typed_strongest_drift_tariff={typed['drift_tariff']}")
    print(f"typed_strongest_deep_target_tariff={typed['deep_target_tariff']}")
    print(f"typed_defect_digest={typed['digest']}")
    print("PASS")


if __name__ == "__main__":
    main()
