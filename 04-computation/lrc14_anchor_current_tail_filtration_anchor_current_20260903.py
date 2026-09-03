#!/usr/bin/env python3
"""Exact audits for the LRC(14) sharp-current tail/7-adic fibre lane.

All integrals below use ordinary ``dt`` on ``[0,1/2)``.  For an odd tail
speed ``w`` put

    sigma_w(t) = 1_{||wt||<1/14} - 1_{||w(t+1/2)||<1/14}.

The script checks four claims.

1.  The sharp current rebate has an exact binomial/layer-cake form.
2.  At the critical 7-adic height of an even anchor, every runner is a unit
    directed edge on a seven-point fibre.  A dynamic program independently
    checks the resulting sharp relaxed six-safe-vertex transport envelope.
3.  Exact rational wall sweeps test the anchored filtration inequality on
    structured and deterministic pseudo-random twelve-tail families.
4.  The counterexample-local 4/7 coefficient is attained by an exact
    five-tail row even after retaining residue-labelled edge types, and the
    same fibre embeds into a primitive twelve-tail 7/2/3 valuation profile.
5.  Common-dilation controls check the Euclidean residual-wall formula for
    the nonlinear rebate itself, including the exact D=7 independence case.

This is a lane-specific research audit, not an LRC(14) proof.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from math import ceil, comb, floor, gcd
from random import Random


HALF = Fraction(1, 2)
RADIUS = Fraction(1, 14)


def v7(number: int) -> int:
    exponent = 0
    while number % 7 == 0:
        number //= 7
        exponent += 1
    return exponent


def circle_distance(value: Fraction) -> Fraction:
    residue = value % 1
    return min(residue, 1 - residue)


def danger(speed: int, time: Fraction, shift: Fraction = Fraction()) -> int:
    return int(circle_distance(speed * (time + shift)) < RADIUS)


def sigma(speed: int, time: Fraction) -> int:
    return danger(speed, time) - danger(speed, time, HALF)


def rebate(depth: int) -> Fraction:
    """The sharp current-only envelope g(depth), for integer depth >= 0."""
    assert depth >= 0
    if depth <= 2:
        return Fraction()
    if depth == 3:
        return Fraction(1, 2)
    return Fraction(depth * (depth * depth - 6 * depth + 11), 12)


def binomial_rebate(depth: int) -> Fraction:
    """Equivalent form 1/2(1_{d>=3}+C(d-1,3))."""
    return Fraction(int(depth >= 3) + (comb(depth - 1, 3) if depth >= 4 else 0), 2)


def rebate_tail_weights(maximum: int = 12) -> tuple[Fraction, ...]:
    """w_j=g(j)-g(j-1), so g(d)=sum_j w_j 1_{d>=j}."""
    return tuple(rebate(j) - rebate(j - 1) for j in range(1, maximum + 1))


def safe_transport_envelope(level: int, budget: int) -> Fraction:
    """Sharp relaxed cost on six safe vertices.

    Six vertices initially carry the common integer current ``level``.  One
    unit edge can lower at most one safe vertex by one while its positive end
    is hidden at the deleted anchor vertex.  Convexity makes equalized
    spending optimal.  Spending beyond level two has no benefit.
    """
    assert 0 <= level <= 12 and budget >= 0
    demand = 6 * max(level - 2, 0)
    spend = min(budget, demand)
    quotient, remainder = divmod(spend, 6)
    answer = (6 - remainder) * rebate(level - quotient)
    if remainder:
        answer += remainder * rebate(level - quotient - 1)
    return answer


def brute_safe_transport_envelope(level: int, budget: int) -> Fraction:
    """Independent DP over six reduction allocations with total <= budget."""
    costs = {(0, 0): Fraction()}
    for vertex in range(6):
        updated: dict[tuple[int, int], Fraction] = {}
        for (_, used), value in costs.items():
            for reduction in range(budget - used + 1):
                key = (vertex + 1, used + reduction)
                residual = max(level - reduction, 0)
                candidate = value + rebate(residual)
                if key not in updated or candidate < updated[key]:
                    updated[key] = candidate
        costs = updated
    return min(value for (vertices, _), value in costs.items() if vertices == 6)


def universal_twelve_tail_envelope(level: int) -> Fraction:
    """Worst transport envelope using n_a <= 12-|Z|."""
    return safe_transport_envelope(level, 12 - level)


def clipped_arcs(
    speed: int,
    shift: Fraction,
    radius: Fraction = RADIUS,
    left: Fraction = Fraction(),
    right: Fraction = HALF,
) -> tuple[tuple[Fraction, Fraction], ...]:
    lo = floor(speed * left + shift - radius) - 1
    hi = ceil(speed * right + shift + radius) + 1
    answer: list[tuple[Fraction, Fraction]] = []
    for index in range(lo, hi + 1):
        arc_left = (index - shift - radius) / speed
        arc_right = (index - shift + radius) / speed
        a = max(left, arc_left)
        b = min(right, arc_right)
        if a < b:
            answer.append((a, b))
    return tuple(answer)


def add_current_events(
    events: dict[Fraction, list[int]], speeds: tuple[int, ...], coordinate: int
) -> None:
    for speed in speeds:
        for a, b in clipped_arcs(speed, Fraction()):
            events[a][coordinate] += 1
            events[b][coordinate] -= 1
        # Since speed is odd, speed/2 is 1/2 modulo one.  These arcs enter
        # sigma with a minus sign.
        for a, b in clipped_arcs(speed, Fraction(speed, 2)):
            events[a][coordinate] -= 1
            events[b][coordinate] += 1


def exact_profile(h: int, speeds: tuple[int, ...]) -> dict[str, object]:
    """Exact current-rebate profile and the critical-height fibre bounds."""
    assert h > 0
    assert 1 <= len(speeds) <= 12 and len(speeds) == len(set(speeds))
    assert all(speed > 0 and speed % 2 for speed in speeds)
    anchor_height = v7(h)
    exact_shell = tuple(w for w in speeds if v7(w) == anchor_height)
    high_shells = tuple(w for w in speeds if v7(w) > anchor_height)

    # Coordinates: total current, high-shell current Z, anchor danger.
    events: dict[Fraction, list[int]] = defaultdict(lambda: [0, 0, 0])
    add_current_events(events, speeds, 0)
    add_current_events(events, high_shells, 1)
    for a, b in clipped_arcs(2 * h, Fraction()):
        events[a][2] += 1
        events[b][2] -= 1
    events[Fraction()]
    events[HALF]

    current = high = anchor = 0
    full_rebate = core_rebate = strip_rebate = Fraction()
    high_rebate = high_ge4_rebate = Fraction()
    exact_fibre_bound = uniform_fibre_bound = Fraction()
    current_tails = [Fraction() for _ in range(13)]
    core_current_tails = [Fraction() for _ in range(13)]
    strip_current_tails = [Fraction() for _ in range(13)]
    high_tails = [Fraction() for _ in range(13)]

    walls = sorted(events)
    for wall, next_wall in zip(walls, walls[1:]):
        dc, dz, da = events[wall]
        current += dc
        high += dz
        anchor += da
        width = next_wall - wall
        assert anchor in (0, 1)
        d = abs(current)
        z = abs(high)
        value = rebate(d)
        full_rebate += width * value
        high_rebate += width * rebate(z)
        if z >= 4:
            high_ge4_rebate += width * rebate(z)
        exact_fibre_bound += width * safe_transport_envelope(z, len(exact_shell)) / 7
        uniform_fibre_bound += width * universal_twelve_tail_envelope(z) / 7
        for threshold in range(1, 13):
            if d >= threshold:
                current_tails[threshold] += width
            if z >= threshold:
                high_tails[threshold] += width
        if anchor:
            strip_rebate += width * value
            for threshold in range(1, 13):
                if d >= threshold:
                    strip_current_tails[threshold] += width
        else:
            core_rebate += width * value
            for threshold in range(1, 13):
                if d >= threshold:
                    core_current_tails[threshold] += width

    assert full_rebate == core_rebate + strip_rebate
    assert core_rebate >= exact_fibre_bound >= uniform_fibre_bound
    assert exact_fibre_bound >= Fraction(max(6 - len(exact_shell), 0), 7) * high_rebate
    assert uniform_fibre_bound >= Fraction(2, 7) * high_ge4_rebate

    # Independent layer-cake reconstruction of both the full and core costs.
    weights = rebate_tail_weights()
    assert full_rebate == sum(
        (weights[j - 1] * current_tails[j] for j in range(1, 13)),
        start=Fraction(),
    )
    assert core_rebate == sum(
        (weights[j - 1] * core_current_tails[j] for j in range(1, 13)),
        start=Fraction(),
    )

    return {
        "h": h,
        "speeds": speeds,
        "anchor_height": anchor_height,
        "exact_shell_count": len(exact_shell),
        "high_shell_count": len(high_shells),
        "full_rebate": full_rebate,
        "core_rebate": core_rebate,
        "strip_rebate": strip_rebate,
        "high_rebate": high_rebate,
        "high_level3_mass": high_tails[3] - high_tails[4],
        "high_ge4_rebate": high_ge4_rebate,
        "exact_fibre_bound": exact_fibre_bound,
        "uniform_fibre_bound": uniform_fibre_bound,
        "full_current_tails": tuple(current_tails),
        "core_current_tails": tuple(core_current_tails),
        "strip_current_tails": tuple(strip_current_tails),
    }


def residual_wall_rebate(
    base: tuple[int, ...], wall_speed: int, wall_shift: Fraction, wall_radius: Fraction
) -> tuple[Fraction, Fraction, tuple[Fraction, ...], tuple[Fraction, ...]]:
    """Return full rebate and its restriction to one residual wall."""
    events: dict[Fraction, list[int]] = defaultdict(lambda: [0, 0])
    add_current_events(events, base, 0)
    if wall_radius:
        for a, b in clipped_arcs(wall_speed, wall_shift, wall_radius):
            events[a][1] += 1
            events[b][1] -= 1
    events[Fraction()]
    events[HALF]
    current = wall_depth = 0
    full = restricted = Fraction()
    full_tails = [Fraction() for _ in range(13)]
    restricted_tails = [Fraction() for _ in range(13)]
    walls = sorted(events)
    for wall, next_wall in zip(walls, walls[1:]):
        dc, dw = events[wall]
        current += dc
        wall_depth += dw
        width = next_wall - wall
        value = rebate(abs(current))
        full += width * value
        restricted += width * value * wall_depth
        for threshold in range(1, 13):
            if abs(current) >= threshold:
                full_tails[threshold] += width
                restricted_tails[threshold] += width * wall_depth
        assert wall_depth in (0, 1)
    return full, restricted, tuple(full_tails), tuple(restricted_tails)


def audit_common_dilation(h: int, base: tuple[int, ...], multiplier: int) -> dict[str, object]:
    scaled = tuple(multiplier * speed for speed in base)
    profile = exact_profile(h, scaled)
    common = 1
    # multiplier is odd, so gcd(multiplier,2h)=gcd(multiplier,h).
    common = gcd(multiplier, h)
    degree = multiplier // common
    reduced_anchor = 2 * h // common
    quotient, remainder = divmod(degree, 7)
    full, residual, full_tails, residual_tails = residual_wall_rebate(
        base, reduced_anchor, Fraction(quotient, 2), Fraction(remainder, 14)
    )
    predicted_strip = (quotient * full + residual) / degree
    assert profile["full_rebate"] == full
    assert profile["strip_rebate"] == predicted_strip
    assert profile["full_current_tails"] == full_tails
    for threshold in range(1, 13):
        predicted_tail = (
            quotient * full_tails[threshold] + residual_tails[threshold]
        ) / degree
        assert profile["strip_current_tails"][threshold] == predicted_tail
    if remainder == 0:
        assert profile["strip_rebate"] == full / 7
        assert profile["core_rebate"] == 6 * full / 7
    return {
        "multiplier": multiplier,
        "degree": degree,
        "quotient": quotient,
        "remainder": remainder,
        "full": full,
        "residual": residual,
        "strip": profile["strip_rebate"],
        "core": profile["core_rebate"],
        "fibre_bound": profile["exact_fibre_bound"],
    }


def fibre_snapshot(h: int, speeds: tuple[int, ...], time: Fraction) -> dict[str, object]:
    anchor_height = v7(h)
    step = Fraction(1, 7 ** (anchor_height + 1))
    high_shells = tuple(w for w in speeds if v7(w) > anchor_height)
    exact_shell = tuple(w for w in speeds if v7(w) == anchor_height)
    points = tuple(time + j * step for j in range(7))
    anchor = tuple(danger(2 * h, point) for point in points)
    high = tuple(sum(sigma(w, point) for w in high_shells) for point in points)
    shell = tuple(sum(sigma(w, point) for w in exact_shell) for point in points)
    total = tuple(sum(sigma(w, point) for w in speeds) for point in points)
    assert sum(anchor) == 1
    assert len(set(high)) == 1
    assert sum(shell) == 0
    return {
        "anchor": anchor,
        "high": high,
        "exact_shell": shell,
        "total": total,
        "safe_cost": sum(rebate(abs(total[j])) for j in range(7) if not anchor[j]),
    }


def audit_random_profiles() -> tuple[int, Fraction]:
    """Deterministic exact hostile sweep; this is evidence, not exhaustiveness."""
    rng = Random(20260903)
    universe = tuple(range(1, 100, 2))
    anchors = (1, 7, 14, 49, 98, 343)
    checked = 0
    minimum_slack: Fraction | None = None
    for h in anchors:
        structured = [
            tuple(range(1, 24, 2)),
            tuple(7 * w for w in range(1, 24, 2)),
            tuple(49 * w for w in range(1, 24, 2)),
        ]
        families = structured + [tuple(sorted(rng.sample(universe, 12))) for _ in range(12)]
        for speeds in families:
            profile = exact_profile(h, speeds)
            slack = profile["core_rebate"] - profile["exact_fibre_bound"]
            assert isinstance(slack, Fraction)
            minimum_slack = slack if minimum_slack is None else min(minimum_slack, slack)
            checked += 1
    assert minimum_slack is not None
    return checked, minimum_slack


def main() -> None:
    for depth in range(13):
        assert rebate(depth) == binomial_rebate(depth)
    weights = rebate_tail_weights()
    assert weights[:4] == (Fraction(), Fraction(), Fraction(1, 2), Fraction(1, 2))
    for threshold in range(5, 13):
        assert weights[threshold - 1] == Fraction((threshold - 2) * (threshold - 3), 4)

    for level in range(13):
        for budget in range(13):
            assert safe_transport_envelope(level, budget) == brute_safe_transport_envelope(
                level, budget
            )
    for budget in range(13):
        values = [safe_transport_envelope(level, budget) for level in range(13)]
        slopes = [values[level] - values[level - 1] for level in range(1, 13)]
        assert all(left <= right for left, right in zip(slopes, slopes[1:]))

    universal = tuple(universal_twelve_tail_envelope(level) for level in range(13))
    expected = (
        Fraction(), Fraction(), Fraction(), Fraction(), Fraction(2), Fraction(11, 2),
        Fraction(15), Fraction(38), Fraction(78), Fraction(279, 2), Fraction(227),
        Fraction(345), Fraction(498),
    )
    assert universal == expected
    universal_tail_weights = tuple(
        universal[level] - universal[level - 1] for level in range(1, 13)
    )
    assert universal_tail_weights == (
        Fraction(), Fraction(), Fraction(), Fraction(2), Fraction(7, 2),
        Fraction(19, 2), Fraction(23), Fraction(40), Fraction(123, 2),
        Fraction(175, 2), Fraction(118), Fraction(153),
    )
    for level in range(4, 13):
        assert universal[level] >= 2 * rebate(level)

    # At a tiny positive time, all listed runners have their positive fibre
    # endpoint at the anchor-danger vertex.  The first family has exactly six
    # height-zero runners, whose residues send one negative endpoint to each
    # safe vertex and erase a constant high current of three completely.
    critical_six_family = (1, 3, 5, 7, 9, 11, 13, 21, 35)
    critical_six_snapshot = fibre_snapshot(
        1, critical_six_family, Fraction(1, 10000)
    )
    assert critical_six_snapshot["high"] == (3,) * 7
    assert critical_six_snapshot["exact_shell"] == (6, -1, -1, -1, -1, -1, -1)
    assert critical_six_snapshot["total"] == (9, 2, 2, 2, 2, 2, 2)
    assert critical_six_snapshot["safe_cost"] == 0

    # Padding by three more height-zero runners preserves zero safe cost and
    # gives the twelve-tail control used by the integrated sweep below.
    level3_family = (1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 35)
    level3_snapshot = fibre_snapshot(1, level3_family, Fraction(1, 10000))
    assert level3_snapshot["high"] == (3,) * 7
    assert level3_snapshot["safe_cost"] == 0

    level4_family = (1, 3, 5, 7, 9, 11, 13, 17, 19, 21, 35, 49)
    level4_snapshot = fibre_snapshot(1, level4_family, Fraction(1, 10000))
    assert level4_snapshot["high"] == (4,) * 7
    assert level4_snapshot["safe_cost"] == universal[4] == 2

    # The 4/7 coefficient forced by the THM-4370 five-tail upper cone is
    # already sharp for physical residue-labelled edges.  Both critical
    # runners put their positive endpoint at the deleted vertex, and their
    # negative endpoints at two distinct safe vertices.
    sharp_four_sevenths_family = (1, 3, 7, 21, 35)
    sharp_four_sevenths_time = Fraction(1, 10000)
    sharp_four_sevenths_snapshot = fibre_snapshot(
        1, sharp_four_sevenths_family, sharp_four_sevenths_time
    )
    sharp_points = tuple(sharp_four_sevenths_time + Fraction(j, 7) for j in range(7))
    edge_w1 = tuple(sigma(1, point) for point in sharp_points)
    edge_w3 = tuple(sigma(3, point) for point in sharp_points)
    assert edge_w1 == (1, 0, 0, -1, 0, 0, 0)
    assert edge_w3 == (1, -1, 0, 0, 0, 0, 0)
    assert 1 * (3 - 0) % 7 == 3 and 3 * (1 - 0) % 7 == 3
    assert sharp_four_sevenths_snapshot["anchor"] == (1, 0, 0, 0, 0, 0, 0)
    assert sharp_four_sevenths_snapshot["high"] == (3,) * 7
    assert sharp_four_sevenths_snapshot["exact_shell"] == (
        2, -1, 0, -1, 0, 0, 0
    )
    assert sharp_four_sevenths_snapshot["total"] == (5, 2, 3, 2, 3, 3, 3)
    assert sharp_four_sevenths_snapshot["safe_cost"] == 2
    sharp_four_sevenths_profile = exact_profile(1, sharp_four_sevenths_family)
    assert sharp_four_sevenths_profile["high_rebate"] == Fraction(1, 70)
    assert sharp_four_sevenths_profile["core_rebate"] == Fraction(2, 245)
    assert sharp_four_sevenths_profile["exact_fibre_bound"] == Fraction(2, 245)
    assert sharp_four_sevenths_profile["core_rebate"] == Fraction(
        4, 7
    ) * sharp_four_sevenths_profile["high_rebate"]

    # Embed the same upper projection at anchor height one into the exact
    # 7-lower / 2-critical / 3-higher profile forced by THM-4370.  This is a
    # primitive local control, not a claimed strict counterexample.
    primitive_lower = (5, 9, 11, 15, 17, 29, 31)
    primitive_critical = (7, 21)
    primitive_higher = (49, 147, 245)
    primitive_family = tuple(sorted(primitive_lower + primitive_critical + primitive_higher))
    primitive_time = Fraction(1, 1_000_000)
    primitive_snapshot = fibre_snapshot(7, primitive_family, primitive_time)
    primitive_points = tuple(primitive_time + Fraction(j, 49) for j in range(7))
    primitive_lower_current = tuple(
        sum(sigma(w, point) for w in primitive_lower) for point in primitive_points
    )
    primitive_upper_projection = tuple(
        high + shell
        for high, shell in zip(
            primitive_snapshot["high"], primitive_snapshot["exact_shell"]
        )
    )
    assert len(set(primitive_family)) == 12
    assert gcd(14, *primitive_family) == 1
    assert sum(v7(w) < 1 for w in primitive_family) == 7
    assert sum(v7(w) == 1 for w in primitive_family) == 2
    assert sum(v7(w) > 1 for w in primitive_family) == 3
    assert primitive_snapshot["anchor"] == (1, 0, 0, 0, 0, 0, 0)
    assert primitive_snapshot["high"] == (3,) * 7
    assert primitive_snapshot["exact_shell"] == (2, -1, 0, -1, 0, 0, 0)
    assert primitive_upper_projection == (5, 2, 3, 2, 3, 3, 3)
    assert primitive_lower_current == (7, 0, -1, 0, -1, -1, -1)
    assert primitive_snapshot["total"] == (12, 2, 2, 2, 2, 2, 2)
    assert primitive_snapshot["safe_cost"] == 0

    # Exact integrated controls for the same two physical twelve-tail rows.
    level3_profile = exact_profile(1, level3_family)
    level4_profile = exact_profile(1, level4_family)

    # Minimal tail-count failure of the tempting ``6/7 of the high-shell
    # rebate survives the core`` assertion.  Three high runners are necessary
    # before g can be nonzero; with only those three exact independence holds.
    # Adding the single smallest height-zero runner makes the inequality fail.
    high_only = exact_profile(1, (7, 21, 35))
    minimal_failure = exact_profile(1, (1, 7, 21, 35))
    assert high_only["core_rebate"] == Fraction(3, 245)
    assert minimal_failure["high_rebate"] == Fraction(1, 70)
    assert minimal_failure["core_rebate"] == Fraction(1, 98)
    assert minimal_failure["core_rebate"] - Fraction(6, 7) * minimal_failure[
        "high_rebate"
    ] == -Fraction(1, 490)
    assert minimal_failure["core_rebate"] == minimal_failure["exact_fibre_bound"]

    # A twelve-tail, smallest-possible-high-triple control shows that padding
    # does not repair the failed implication.  Its maximum speed 35 is forced
    # by the need for three distinct positive odd multiples of seven.
    twelve_failure_family = (1, 5, 7, 9, 11, 13, 17, 19, 21, 29, 31, 35)
    twelve_failure = exact_profile(1, twelve_failure_family)
    assert twelve_failure["core_rebate"] == Fraction(23423, 8369690)
    assert twelve_failure["core_rebate"] - Fraction(6, 7) * twelve_failure[
        "high_rebate"
    ] == -Fraction(79063, 8369690)

    base = (1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 45)
    dilation_rows = tuple(audit_common_dilation(420, base, m) for m in (1, 13, 49, 127))
    checked, minimum_slack = audit_random_profiles()

    print("sharp_current_binomial_and_tail_identity")
    print(" g(d)=1/2*(1_{d>=3}+binom(d-1,3))")
    print(" E[g(|C|)]=1/2*P(|C|>=3)+1/2*sum_{j=4}^12 binom(j-2,2)*P(|C|>=j)")
    print(" tail_weights=" + ",".join(str(value) for value in weights))
    print()
    print("critical_septimal_fibre_transport")
    print(" universal_Lambda_12=" + ",".join(str(value) for value in universal))
    print(
        " universal_tail_weights="
        + ",".join(str(value) for value in universal_tail_weights)
    )
    print(" exact_bound=core_g(C)>=1/7*integral_half L_{N_a}(|C_{>a}|)")
    print(" sparse_shell_bound>=(6-N_a)_+/7*integral_half g(|C_{>a}|)")
    print(" uniform_bound>=2/7*integral_half g(|C_{>a}|)*1_{|C_{>a}|>=4}")
    print(f" dynamic_program_rows={13 * 13}")
    print(" all_L_n_even_extensions_convex=True")
    print()
    print("physical_fibre_witnesses")
    print(f" critical_six_family={critical_six_family}")
    print(f" critical_six_snapshot={critical_six_snapshot}")
    print(" exactly_six_edges_erase_level3_on_all_six_safe_vertices")
    print(f" level3_family={level3_family}")
    print(f" level3_snapshot={level3_snapshot}")
    print(" level3_mean_current=3_but_all_six_safe_vertices_have_zero_rebate")
    print(f" level4_family={level4_family}")
    print(f" level4_snapshot={level4_snapshot}")
    print(" level4_safe_cost=2=sharp_uniform_fibre_envelope")
    print()
    print("counterexample_locus_four_sevenths_sharpness")
    print(f" upper_five_family={sharp_four_sevenths_family}")
    print(f" upper_snapshot={sharp_four_sevenths_snapshot}")
    print(f" edge_w1={edge_w1} endpoints=(0,3) residue_type=3")
    print(f" edge_w3={edge_w3} endpoints=(0,1) residue_type=3")
    print(
        f" global_core={sharp_four_sevenths_profile['core_rebate']}"
        f" high={sharp_four_sevenths_profile['high_rebate']}"
        f" exact_bound={sharp_four_sevenths_profile['exact_fibre_bound']}"
        " ratio=4/7"
    )
    print(f" primitive_7_2_3_family={primitive_family}")
    print(" primitive_valuation_counts=(7,2,3) gcd_with_anchor=1")
    print(f" primitive_upper_projection={primitive_upper_projection}")
    print(f" primitive_lower_current={primitive_lower_current}")
    print(f" primitive_total={primitive_snapshot['total']}")
    print(" primitive_safe_cost=0; local_control_not_strict_counterexample")
    print()
    print("minimal_exact_failure_of_naive_high_shell_independence")
    print(" high_only=(7,21,35): core_g=3/245=(6/7)*high_g")
    print(
        " add_speed_1: core_g="
        f"{minimal_failure['core_rebate']} high_g={minimal_failure['high_rebate']}"
        " naive_slack=-1/490 exact_transport_bound=1/98"
    )
    print(f" twelve_tail_max35_family={twelve_failure_family}")
    print(
        f" twelve_tail_core={twelve_failure['core_rebate']}"
        f" naive_slack={twelve_failure['core_rebate'] - Fraction(6, 7) * twelve_failure['high_rebate']}"
    )
    print()
    print("integrated_physical_controls")
    for name, profile in (("level3", level3_profile), ("level4", level4_profile)):
        print(
            f" {name}: core={profile['core_rebate']} high={profile['high_rebate']}"
            f" high_level3_mass={profile['high_level3_mass']}"
            f" exact_bound={profile['exact_fibre_bound']}"
            f" uniform_bound={profile['uniform_fibre_bound']}"
        )
    print()
    print("common_dilation_residual_wall_controls")
    print(" every_tail_event_obeys T_j(strip)=(q*T_j(full)+T_j(residual_wall))/D")
    for row in dilation_rows:
        print(
            f" m={row['multiplier']} D={row['degree']}=7*{row['quotient']}+{row['remainder']}"
            f" full={row['full']} residual={row['residual']}"
            f" strip={row['strip']} core={row['core']} fibre_bound={row['fibre_bound']}"
        )
    print()
    print("deterministic_exact_profile_sweep")
    print(f" profiles={checked} minimum_exact_fibre_slack={minimum_slack}")
    print()
    print("canonical_hostile_scope")
    print(" W_9=(1,9,...,9^11) occupies one 7-adic shell; for h=1 C_{>a}=0,")
    print(" and for 7|h the projection M_a is already zero.  Hence this shell/transducer")
    print(" lower bound is exactly silent on the THM-4346 geometric-chain hostile.")


if __name__ == "__main__":
    main()
