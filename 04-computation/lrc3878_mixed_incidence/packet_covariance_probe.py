#!/usr/bin/env python3
"""Scratch-only mixed packet-incidence probe for THM-3878.

For a body-safe set G and the primitive pair-danger set A_{p,q}, the full
scale-one safe mass is

    |G intersect (A_{p,q}(t .))^c|
      = |G|(1-|A|) - Cov(1_G, 1_A(t .)).

The covariance is evaluated in two independent ways: direct exact interval
intersection and a signed endpoint/Bernoulli-B2 formula.  We also test the
Cauchy--Schwarz certificate using the exact t-grid energy of G.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction as Q
from hashlib import sha256
import json
from math import isqrt
import sys


sys.stdout.reconfigure(newline="\n")

DELTA = Q(1, 14)
ONE = Q(1)
EXPECTED_SEMANTIC_SHA256 = "9eafb7bb0fc15266386b2e5cfa2f155ef402a63bea8a5fa045de3124de390a88"

SCALE1 = (
    (1, 3), (1, 4), (1, 9), (1, 10),
    (2, 3), (2, 9), (2, 15), (2, 21), (2, 23),
    (3, 7), (3, 8), (3, 14), (3, 17), (3, 19), (3, 20),
    (3, 22), (3, 26), (3, 31), (3, 38),
    (4, 7), (4, 13), (4, 19), (4, 21), (4, 25), (4, 37),
    (4, 43), (4, 49), (4, 51),
    (5, 6), (5, 12), (5, 17), (5, 18), (5, 24), (5, 29),
    (5, 36), (5, 39), (5, 41), (5, 42), (5, 48), (5, 53),
    (5, 54), (5, 63),
    (6, 11), (6, 17), (6, 19), (6, 23), (6, 41), (6, 47),
    (6, 53), (6, 65),
    (7, 10), (7, 13), (7, 15), (7, 22),
    (8, 9), (8, 21), (9, 11),
)


def require(ok: bool, label: str) -> None:
    if not ok:
        raise RuntimeError(label)


def qfmt(x: Q) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def merge(intervals: list[tuple[Q, Q]]) -> list[tuple[Q, Q]]:
    out: list[list[Q]] = []
    for a, b in sorted(intervals):
        a, b = Q(a), Q(b)
        if a >= b:
            continue
        if not out or a > out[-1][1]:
            out.append([a, b])
        elif b > out[-1][1]:
            out[-1][1] = b
    return [(a, b) for a, b in out]


def danger(speed: int) -> list[tuple[Q, Q]]:
    radius = Q(1, 14 * speed)
    pieces: list[tuple[Q, Q]] = []
    for k in range(speed):
        c = Q(k, speed)
        a, b = c - radius, c + radius
        if a < 0:
            pieces.extend(((0, b), (a + 1, 1)))
        elif b > 1:
            pieces.extend(((a, 1), (0, b - 1)))
        else:
            pieces.append((a, b))
    return merge(pieces)


def union_danger(speeds: tuple[int, ...]) -> list[tuple[Q, Q]]:
    return merge(sum((danger(w) for w in speeds), []))


def complement(pieces: list[tuple[Q, Q]]) -> list[tuple[Q, Q]]:
    out = []
    at = Q(0)
    for a, b in pieces:
        if at < a:
            out.append((at, a))
        at = max(at, b)
    if at < 1:
        out.append((at, Q(1)))
    return out


def measure(pieces: list[tuple[Q, Q]]) -> Q:
    return sum((b - a for a, b in pieces), Q(0))


def intersection_measure(
    left: list[tuple[Q, Q]], right: list[tuple[Q, Q]]
) -> Q:
    i = j = 0
    total = Q(0)
    while i < len(left) and j < len(right):
        a = max(left[i][0], right[j][0])
        b = min(left[i][1], right[j][1])
        if a < b:
            total += b - a
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return total


def pullback(pieces: list[tuple[Q, Q]], t: int) -> list[tuple[Q, Q]]:
    return merge([
        (Q(k, t) + a / t, Q(k, t) + b / t)
        for k in range(t) for a, b in pieces
    ])


def endpoints(pieces: list[tuple[Q, Q]]) -> tuple[tuple[Q, int], ...]:
    """Periodic distributional endpoints; cancel the identified 0=1 wall."""
    weights: defaultdict[Q, int] = defaultdict(int)
    for a, b in pieces:
        weights[a % 1] += 1
        weights[b % 1] -= 1
    return tuple(sorted((x, s) for x, s in weights.items() if s))


def b2(x: Q) -> Q:
    x %= 1
    return x * x - x + Q(1, 6)


def circle_distance(x: Q) -> Q:
    x %= 1
    return min(x, 1 - x)


def covariance_endpoint(
    g_pieces: list[tuple[Q, Q]],
    a_pieces: list[tuple[Q, Q]],
    t: int,
) -> Q:
    eg = endpoints(g_pieces)
    ea = endpoints(a_pieces)
    return Q(1, 2 * t) * sum(
        sg * sa * b2(f - t * e)
        for e, sg in eg for f, sa in ea
    )


def grid_energy(g_pieces: list[tuple[Q, Q]], t: int) -> Q:
    """disc_t(G)=sum_{m != 0}|g-hat(mt)|^2, exact endpoint B2 form."""
    eg = endpoints(g_pieces)
    return Q(1, 2 * t * t) * sum(
        se * sf * b2(t * (e - f))
        for e, se in eg for f, sf in eg
    )


def circle_component_lengths(pieces: list[tuple[Q, Q]]) -> tuple[Q, ...]:
    pieces = merge(pieces)
    lengths = [b - a for a, b in pieces]
    if len(pieces) > 1 and pieces[0][0] == 0 and pieces[-1][1] == 1:
        lengths = [pieces[0][1] + 1 - pieces[-1][0]] + [
            b - a for a, b in pieces[1:-1]
        ]
    return tuple(sorted(lengths))


def cs_positive(m: Q, alpha: Q, disc: Q) -> bool:
    # m(1-alpha) > sqrt(alpha(1-alpha)disc)
    return m * m * (1 - alpha) > alpha * disc


def packet_tax_positive(m: Q, alpha: Q, disc: Q, t: int) -> bool:
    """Contradict the THM-2048 abstract occupancy floor for G subset A(t.)."""
    theta = (Q(t) * m / alpha) % 1
    necessary = (1 - alpha) * m * m / alpha + alpha * theta * (1 - theta) / (t * t)
    return disc < necessary


def main() -> None:
    body = tuple(range(1, 12))
    g = complement(union_danger(body))
    m = measure(g)
    require(m == Q(10931, 194040), "AP11 safe mass")
    require(len(endpoints(g)) == 28, "AP11 positive-arc endpoint count")

    # The t=15 scalar/topological liar: identical packet mass and identical
    # body energy, but the signed cross phase may differ.
    liar = []
    for p, q in ((3, 14), (7, 10)):
        a = union_danger((p, q))
        alpha = measure(a)
        at = pullback(a, 15)
        overlap = intersection_measure(g, at)
        cov_direct = overlap - m * alpha
        cov_b2 = covariance_endpoint(g, a, 15)
        require(cov_direct == cov_b2, f"liar B2 covariance {(p, q)}")
        full_mass = m - overlap
        liar.append((p, q, alpha, cov_b2, full_mass))
    require(liar[0][2] == liar[1][2] == Q(13, 49), "liar equal mass")
    require(liar[0][3] != liar[1][3], "liar cross-phase separation")

    # At t=11 the same liar is stronger: equal packet mass, pair sum,
    # component counts, pullback count, inversion, and body energy, but the
    # signed cross covariance even has opposite sign.
    liar11 = []
    isolated = tuple(Q(r, 14) for r in (3, 5, 9, 11))
    for p, q in ((3, 14), (7, 10)):
        pair_danger = union_danger((p, q))
        alpha = measure(pair_danger)
        cov = covariance_endpoint(g, pair_danger, 11)
        overlap = intersection_measure(g, pullback(pair_danger, 11))
        require(cov == overlap - m * alpha, f"liar11 covariance {(p,q)}")
        require(len(circle_component_lengths(pullback(pair_danger, 11))) == 154,
                f"liar11 pullback count {(p,q)}")
        liar11.append((p, q, alpha, len(circle_component_lengths(pair_danger)), cov, m - overlap))
    require(liar11[0][4] < 0 < liar11[1][4], "liar11 opposite covariance signs")
    require(not any(
        circle_distance(11 * 3 * y) >= DELTA
        and circle_distance(11 * 14 * y) >= DELTA
        for y in isolated
    ), "liar11 first pair kills isolated walls")
    require(all(
        circle_distance(11 * 7 * y) >= DELTA
        and circle_distance(11 * 10 * y) >= DELTA
        for y in isolated
    ), "liar11 second pair preserves isolated walls")

    # A complete AP11 family certificate over all t>=U=11.  The universal
    # endpoint estimate disc_t(G)<=r^2/(3t^2), with r=14, closes the tail once
    # it closes the largest pair-danger mass.  Below that threshold use exact
    # grid energy; only its resonant failures need the signed cross phase.
    alphas = {
        (p, q): measure(union_danger((p, q))) for p, q in SCALE1
    }
    worst_pair, worst_alpha = max(alphas.items(), key=lambda item: item[1])
    require((worst_pair, worst_alpha) == ((1, 10), Q(19, 70)), "worst pair mass")
    tail_square = Q(14 * 14) * worst_alpha / (
        3 * m * m * (1 - worst_alpha)
    )
    tail_start = isqrt(tail_square.numerator // tail_square.denominator)
    while Q(tail_start * tail_start) <= tail_square:
        tail_start += 1
    require(tail_start == 88, "AP packet-energy tail threshold")

    global_cells = 0
    global_cs = 0
    global_tax = 0
    global_phase = []
    for t in range(11, tail_start):
        disc = grid_energy(g, t)
        for (p, q), alpha in alphas.items():
            tails = (t * p, t * q)
            if tails[0] in body or tails[1] in body or tails[0] == tails[1]:
                continue
            global_cells += 1
            if cs_positive(m, alpha, disc):
                global_cs += 1
                continue
            if packet_tax_positive(m, alpha, disc, t):
                global_tax += 1
                continue
            pair_danger = union_danger((p, q))
            cov = covariance_endpoint(g, pair_danger, t)
            direct_overlap = intersection_measure(g, pullback(pair_danger, t))
            require(cov == direct_overlap - m * alpha,
                    f"global B2 covariance {(p,q,t)}")
            full_mass = m - direct_overlap
            require(full_mass > 0, f"global resonant phase closure {(p,q,t)}")
            global_phase.append((p, q, t, disc, cov, full_mass))
    require(global_cells == 4385, "global finite physical cell count")
    require(global_tax == 32, "global packet-tax closure count")
    require(len(global_phase) == 69, "global phase-residual count")
    require({t for _, _, t, *_ in global_phase} == {13, 26},
            "global energy resonance support")
    phase_min = min(global_phase, key=lambda row: row[-1])

    # Exhaust the AP11 body's actually relevant integer t cells in the
    # conditional cyclic-width window for all 57 scale-one survivors.
    rows = []
    actual_zero = []
    cs_closed = []
    cs_open = []
    tax_closed = []
    crude_closed = []
    covariance_signs = Counter()
    by_pair = {}
    for p, q in SCALE1:
        a = union_danger((p, q))
        alpha = measure(a)
        beta = max(circle_component_lengths(a))
        # Necessary THM-3878 window: U<=t<42 beta U, with U=11.
        t_max = (42 * beta * 11).numerator // (42 * beta * 11).denominator
        if Q(t_max) == 42 * beta * 11:
            t_max -= 1
        pair_stats = Counter()
        for t in range(11, t_max + 1):
            # Physical 13-speed rows require distinct speeds.
            tails = (t * p, t * q)
            if tails[0] in body or tails[1] in body or tails[0] == tails[1]:
                continue
            at = pullback(a, t)
            overlap = intersection_measure(g, at)
            cov = overlap - m * alpha
            require(cov == covariance_endpoint(g, a, t), f"B2 covariance {(p,q,t)}")
            full_mass = m - overlap
            if full_mass == 0:
                actual_zero.append((p, q, t))
            disc = grid_energy(g, t)
            require(disc >= 0, f"grid energy {(p,q,t)}")
            csp = cs_positive(m, alpha, disc)
            if csp:
                cs_closed.append((p, q, t))
            else:
                cs_open.append((p, q, t, disc, cov, full_mass))
                if packet_tax_positive(m, alpha, disc, t):
                    tax_closed.append((p, q, t))
            # THM-732 universal endpoint-energy bound.
            crude_disc = Q(14 * 14, 3 * t * t)
            if cs_positive(m, alpha, crude_disc):
                crude_closed.append((p, q, t))
            covariance_signs["negative" if cov < 0 else "zero" if cov == 0 else "positive"] += 1
            pair_stats["cells"] += 1
            pair_stats["cs"] += int(csp)
            pair_stats["positive"] += int(full_mass > 0)
            rows.append({
                "p": p, "q": q, "t": t,
                "alpha": qfmt(alpha), "disc": qfmt(disc),
                "cov": qfmt(cov), "safe": qfmt(full_mass), "cs": csp,
            })
        by_pair[(p, q)] = pair_stats

    require(not actual_zero, "THM-733 AP-body positivity in relevant cells")
    semantic = {
        "scope": "AP11 body; all 57 THM3878 scale-one survivors for t>=U=11",
        "rows": rows,
        "liar": [tuple(qfmt(x) if isinstance(x, Q) else x for x in row) for row in liar],
        "liar11": [tuple(qfmt(x) if isinstance(x, Q) else x for x in row) for row in liar11],
        "global_counts": [global_cells, global_cs, global_tax, len(global_phase), tail_start],
        "global_phase": [
            (p, q, t, qfmt(disc), qfmt(cov), qfmt(full_mass))
            for p, q, t, disc, cov, full_mass in global_phase
        ],
        "formula": "Cov=(1/(2t))*sum sigma_e tau_f B2({f-te})",
    }
    digest = sha256(json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")).hexdigest()
    if digest != EXPECTED_SEMANTIC_SHA256:
        raise RuntimeError("frozen semantic transcript")

    print("THM3878_MIXED_PACKET_COVARIANCE_SCRATCH_PROBE_20260823")
    print("scope=AP11_body;all_57_scale1_rows;t>=U=11;LRC14=OPEN")
    print(f"ap11_mass={qfmt(m)};positive_arcs=14;endpoint_count={len(endpoints(g))}")
    print(f"tested_physical_cells={len(rows)};actual_positive={len(rows)-len(actual_zero)};actual_zero={len(actual_zero)}")
    print(f"exact_energy_CS_closures={len(cs_closed)};crude_endpoint_CS_closures={len(crude_closed)}")
    print(f"conditional_packet_tax_extra_closures={len(tax_closed)};signed_phase_needed={len(cs_open)-len(tax_closed)}")
    print("exact_energy_CS_survivors=" + repr(tuple(
        (p, q, t, qfmt(disc), qfmt(cov), qfmt(full_mass))
        for p, q, t, disc, cov, full_mass in cs_open
    )))
    print("covariance_sign_histogram=" + repr(dict(covariance_signs)))
    print("pair_CS_counts=" + repr(tuple((p, q, s['cs'], s['cells']) for (p, q), s in by_pair.items())))
    for p, q, alpha, cov, full_mass in liar:
        print(f"liar_pair=({p},{q});alpha={qfmt(alpha)};disc15={qfmt(grid_energy(g,15))};cov={qfmt(cov)};safe_mass={qfmt(full_mass)}")
    print("liar_result=same_alpha_same_body_energy_different_signed_cross_phase")
    for p, q, alpha, count, cov, full_mass in liar11:
        print(f"strong_liar_t11_pair=({p},{q});sum={p+q};alpha={qfmt(alpha)};base_components={count};pullback_components={11*count};disc11={qfmt(grid_energy(g,11))};cov={qfmt(cov)};safe_mass={qfmt(full_mass)}")
    print("strong_liar_t11_result=same_sum_mass_components_inversion_body_energy;opposite_endpoint_incidence_and_covariance_sign")
    print(f"global_AP11_tail_start={tail_start};worst_alpha={qfmt(worst_alpha)}@{worst_pair};endpoint_bound_closes_all_t_ge_{tail_start}")
    print(f"global_AP11_finite_cells_t11_to_{tail_start-1}={global_cells};exact_energy_closures={global_cs};quantized_packet_tax_closures={global_tax};signed_phase_closures={len(global_phase)}")
    print("global_AP11_energy_failure_times=13,26;packet_tax_reduces_101_to_69;all_remaining_signed_cross_covariances_positive_safe")
    print("global_AP11_min_phase_safe=" + repr((
        phase_min[0], phase_min[1], phase_min[2], qfmt(phase_min[-1])
    )))
    print("global_AP11_conclusion=all_57_scale1_rows_with_physical_distinctness_and_t_ge_U_are_lonely;known_THM734_scope_but_new_packet_covariance_proof")
    print("semantic_sha256=" + digest)
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
