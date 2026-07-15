#!/usr/bin/env python3
"""Exact verifier for THM-831's small-half-frequency two-ball theorem."""

from __future__ import annotations

from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import gcd
from pathlib import Path


THETA = F(11, 13)
EPS = F(2, 13)


def circle_norm(x: F) -> F:
    x %= 1
    return min(x, 1 - x)


def folded_value(t: F, alpha: int, beta: int) -> F:
    return circle_norm(alpha * t) + circle_norm(beta * t)


def merge_intervals(intervals: list[tuple[F, F]]) -> list[tuple[F, F]]:
    out: list[list[F]] = []
    for lo, hi in sorted(intervals):
        if not out or lo > out[-1][1]:
            out.append([lo, hi])
        else:
            out[-1][1] = max(out[-1][1], hi)
    return [(lo, hi) for lo, hi in out]


def affine_sweep(alpha: int, beta: int) -> list[tuple[F, F]]:
    """Build {||alpha*t||+||beta*t|| >= 11/13} from all affine cells."""
    walls = {F(k, 2 * alpha) for k in range(2 * alpha + 1)}
    walls |= {F(k, 2 * beta) for k in range(2 * beta + 1)}
    walls |= {F(0), F(1)}
    base = sorted(walls)
    cuts = set(base)
    for lo, hi in zip(base, base[1:]):
        flo = folded_value(lo, alpha, beta)
        fhi = folded_value(hi, alpha, beta)
        if (flo - THETA) * (fhi - THETA) < 0:
            cross = lo + (THETA - flo) * (hi - lo) / (fhi - flo)
            cuts.add(cross)
    grid = sorted(cuts)
    good: list[tuple[F, F]] = []
    for lo, hi in zip(grid, grid[1:]):
        if folded_value((lo + hi) / 2, alpha, beta) >= THETA:
            good.append((lo, hi))
    return merge_intervals(good)


def odd_offsets(alpha: int) -> tuple[int, ...]:
    return tuple(
        m
        for m in range(-4 * alpha, 4 * alpha + 1)
        if m and m % 2 and 13 * abs(m) <= 4 * alpha
    )


def bezout_intervals(alpha: int, beta: int) -> list[tuple[F, F, int]]:
    """Intersect the two odd-speed teeth indexed by rp-sq=m."""
    r, s = alpha - beta, alpha + beta
    rows: list[tuple[F, F, int]] = []
    for m in odd_offsets(alpha):
        p = next(p for p in range(s) if (r * p - m) % s == 0)
        q = (r * p - m) // s
        assert (p - q) % 2 == 1
        xlo, xhi = F(p, s) - F(2, 13 * s), F(p, s) + F(2, 13 * s)
        ylo, yhi = F(q, r) - F(2, 13 * r), F(q, r) + F(2, 13 * r)
        lo, hi = max(xlo, ylo), min(xhi, yhi)
        assert F(0) <= lo < hi <= F(1)
        rows.append((lo, hi, m))
    return sorted(rows)


def formula_balls(alpha: int, beta: int) -> list[tuple[F, F, int]]:
    """Closed component formula, returned as (centre, radius, signed q)."""
    r, s = alpha - beta, alpha + beta
    u = pow(r, -1, s)
    n = (4 * alpha + 13) // 26
    rows: list[tuple[F, F, int]] = []
    for j in range(1, n + 1):
        q = 2 * j - 1
        p = (u * q) % s
        assert 1 <= p < s
        h = F(min(4 * r, 4 * alpha - 13 * q), 26 * r * s)
        c = (F(p, s) + F(min(0, 4 * beta - 13 * q), 26 * r * s)) % 1
        rows.extend(((c, h, q), ((-c) % 1, h, -q)))
    return sorted(rows)


def weighted_switch_slack(balls: list[tuple[F, F, int]]) -> F:
    slacks: list[F] = []
    for i, (ci, hi, _) in enumerate(balls):
        for j, (cj, hj, _) in enumerate(balls):
            for k, (ck, hk, _) in enumerate(balls):
                if i == j == k:
                    continue
                slacks.append(circle_norm(2 * ci - cj - ck) - (2 * hi + hj + hk))
    return min(slacks)


def in_intervals(t: F, intervals: list[tuple[F, F]]) -> bool:
    t %= 1
    return any(lo <= t <= hi for lo, hi in intervals)


def radius_to_centres(t: F, centres: tuple[F, F]) -> F:
    return min(circle_norm(t - c) for c in centres)


def direct_packet(
    e_set: tuple[F, ...], r_set: tuple[F, ...], intervals: list[tuple[F, F]]
) -> bool:
    return all(in_intervals(t + u, intervals) for t in e_set for u in r_set)


def tournament_edges(order: list[int]) -> set[tuple[int, int]]:
    rank = {v: i for i, v in enumerate(order)}
    return {
        (u, v) if rank[u] < rank[v] else (v, u)
        for u, v in combinations(order, 2)
    }


def main() -> None:
    admissible = [
        (alpha, beta)
        for alpha in range(4, 10)
        for beta in range(1, alpha)
        if gcd(alpha, beta) == 1 and (alpha - beta) % 2 == 1
    ]
    assert len(admissible) == 16

    print("THM831_SMALL_HALF_FREQUENCY_TWO_BALL_RADIUS_BUDGET")
    print(f"threshold={THETA} epsilon={EPS}")
    print(f"admissible_rows={len(admissible)} alpha_range=4..9")
    print(
        "columns=alpha,beta,r,s,centre,h,second_difference_gap,slack,"
        "full_s_tooth,intervals"
    )

    certificate_rows: list[str] = []
    summaries: list[dict[str, F | int | bool]] = []
    synthetic_cases = 0
    synthetic_mismatches = 0

    for index, (alpha, beta) in enumerate(admissible):
        r, s = alpha - beta, alpha + beta
        assert odd_offsets(alpha) == (-1, 1)
        swept = affine_sweep(alpha, beta)
        built = bezout_intervals(alpha, beta)
        literal = [(lo, hi) for lo, hi, _ in built]
        assert swept == literal
        assert len(literal) == 2
        (lo0, hi0), (lo1, hi1) = literal
        assert lo0 == 1 - hi1 and hi0 == 1 - lo1
        h0, h1 = (hi0 - lo0) / 2, (hi1 - lo1) / 2
        assert h0 == h1
        h = h0
        centres = ((lo0 + hi0) / 2, (lo1 + hi1) / 2)
        c = min(centres)
        assert centres == (c, 1 - c)

        formula_h = F(2, 13 * s) if beta >= 4 else F(4 * alpha - 13, 26 * r * s)
        assert h == formula_h
        full_s_tooth = beta >= 4

        gaps = []
        for i, ci in enumerate(centres):
            for j, cj in enumerate(centres):
                for k, ck in enumerate(centres):
                    if i == j == k:
                        continue
                    gaps.append(circle_norm(2 * ci - cj - ck))
        delta = min(gaps)
        assert delta == min(circle_norm(2 * c), circle_norm(4 * c))
        assert delta > 4 * h
        assert circle_norm(centres[0] - centres[1]) > 2 * h
        slack = delta - 4 * h

        # Exact finite packets exercise both directions of the radius theorem,
        # including return radii larger than h and centre-switch candidates.
        e_sets: list[tuple[F, ...]] = []
        for centre in centres:
            for z in range(-6, 7):
                e_sets.append(((centre + F(z, 4) * h) % 1,))
        for z in range(-4, 5):
            e_sets.append(
                (
                    (centres[0] + F(z, 4) * h) % 1,
                    (centres[1] - F(z, 4) * h) % 1,
                )
            )
        return_radii = {F(k, 4) * h for k in range(9)}
        return_radii |= {
            circle_norm(centres[0] - centres[1]),
            F(1, 2),
        }
        r_sets = [(F(0), u % 1, (-u) % 1) for u in sorted(return_radii)]
        for e_set in e_sets:
            for r_set in r_sets:
                direct = direct_packet(e_set, r_set, literal)
                budget = (
                    max(radius_to_centres(t, centres) for t in e_set)
                    + max(circle_norm(u) for u in r_set)
                    <= h
                )
                synthetic_cases += 1
                synthetic_mismatches += direct != budget
        assert synthetic_mismatches == 0

        intervals_word = ";".join(f"[{lo},{hi}]" for lo, hi in literal)
        row = (
            f"{alpha},{beta},{r},{s},{c},{h},{delta},{slack},"
            f"{int(full_s_tooth)},{intervals_word}"
        )
        certificate_rows.append(row)
        summaries.append(
            {
                "alpha": alpha,
                "beta": beta,
                "slack": slack,
                "relative": slack / h,
            }
        )
        print(f"row[{index:02d}]={row}")

    # Independent general component-formula and sharp-switch audit.  This
    # range contains 518 primitive opposite-parity pairs.
    general_rows: list[str] = []
    empty_rows = 0
    positive_switch_rows = 0
    nonpositive_switch_rows = 0
    equal_radius_law_rows = 0
    for alpha in range(2, 51):
        for beta in range(1, alpha):
            if gcd(alpha, beta) != 1 or (alpha - beta) % 2 != 1:
                continue
            swept = affine_sweep(alpha, beta)
            built = bezout_intervals(alpha, beta)
            balls = formula_balls(alpha, beta)
            formula_intervals = sorted((c - h, c + h) for c, h, _ in balls)
            assert all(F(0) <= lo < hi <= F(1) for lo, hi in formula_intervals)
            assert swept == [(lo, hi) for lo, hi, _ in built] == formula_intervals
            n = (4 * alpha + 13) // 26
            assert len(balls) == 2 * n
            if n == 0:
                empty_rows += 1
                general_rows.append(f"{alpha},{beta},0,0,empty")
                continue
            radii = {h for _, h, _ in balls}
            equal_radius = len(radii) == 1
            predicted_equal = n == 1 or 4 * beta >= 13 * (2 * n - 1)
            assert equal_radius == predicted_equal
            equal_radius_law_rows += 1
            slack = weighted_switch_slack(balls)
            if alpha <= 9:
                assert slack > 0
                positive_switch_rows += 1
            else:
                assert slack <= 0
                nonpositive_switch_rows += 1
                by_q = {q: (c, h) for c, h, q in balls}
                c1, h1 = by_q[1]
                c3, h3 = by_q[3]
                defect = circle_norm(c3 - 3 * c1)
                tax = 3 * h1 + h3
                assert defect <= tax
            general_rows.append(
                f"{alpha},{beta},{n},{len(balls)},{int(equal_radius)},{slack}"
            )
    assert len(general_rows) == 518
    assert empty_rows == 2
    assert positive_switch_rows == 16
    assert nonpositive_switch_rows == 500

    # Alpha=10 is the exact topological and no-switch boundary: every
    # admissible row has offsets -3,-1,1,3 and the q=(-1,-3,+1) switch.
    boundary_rows = []
    for beta in range(1, 10):
        if gcd(10, beta) != 1 or (10 - beta) % 2 != 1:
            continue
        offsets = odd_offsets(10)
        swept = affine_sweep(10, beta)
        built = bezout_intervals(10, beta)
        assert offsets == (-3, -1, 1, 3)
        assert swept == [(lo, hi) for lo, hi, _ in built]
        assert len(swept) == 4
        by_q = {q: (c, h) for c, h, q in formula_balls(10, beta)}
        c1, h1 = by_q[1]
        c3, h3 = by_q[3]
        defect = circle_norm(c3 - 3 * c1)
        tax = 3 * h1 + h3
        boundary_rows.append((beta, offsets, len(swept), defect, tax, defect - tax))
    assert len(boundary_rows) == 4

    raw_order = sorted(
        range(len(summaries)),
        key=lambda i: (summaries[i]["slack"], -summaries[i]["alpha"], -summaries[i]["beta"]),
        reverse=True,
    )
    relative_order = sorted(
        range(len(summaries)),
        key=lambda i: (summaries[i]["relative"], -summaries[i]["alpha"], -summaries[i]["beta"]),
        reverse=True,
    )
    raw_edges = tournament_edges(raw_order)
    relative_edges = tournament_edges(relative_order)
    edge_flips = len(raw_edges ^ relative_edges) // 2

    digest = sha256("\n".join(certificate_rows).encode()).hexdigest()
    formula_digest = sha256("\n".join(general_rows).encode()).hexdigest()
    source_digest = sha256(Path(__file__).read_bytes()).hexdigest()
    print("\nEXACT_CHECKS")
    print("affine_sweep_equals_bezout_intersections=16/16")
    print("two_antipodal_components=16/16")
    print("radius_formula=16/16")
    print("second_difference_gap_gt_4h=16/16")
    print(f"synthetic_compact_packets={synthetic_cases} mismatches={synthetic_mismatches}")
    print("primitive_formula_audit_alpha_le_50=518/518 including_2_empty")
    print("equal_radius_iff_n1_or_4beta_ge_13qmax=516/516_nonempty")
    print(
        "weighted_no_switch_frontier="
        f"positive_alpha4..9:{positive_switch_rows} "
        f"nonpositive_alpha10..50:{nonpositive_switch_rows}"
    )
    boundary_word = "[" + ", ".join(
        f"(beta={beta},offsets={offsets},components={count},"
        f"defect={defect},tax={tax},slack={slack})"
        for beta, offsets, count, defect, tax, slack in boundary_rows
    ) + "]"
    print(f"alpha10_boundary_rows={boundary_word}")
    print(f"certificate_sha256={digest}")
    print(f"formula_audit_sha256={formula_digest}")
    print(f"source_sha256={source_digest}")
    print("\nTOURNAMENT_TELEMETRY")
    print("vertices=16_reduced_ratio_types")
    print("pair_observable=no_switch_slack_difference")
    print("switch=raw_slack_vs_slack_per_target_radius")
    print("tie_path=lexicographic_alpha_beta")
    print(f"edge_flips={edge_flips}")
    print("both_tournaments=transitive score_histogram=0..15 cycles=0 SCCs=16_singletons HP=1")
    print("preserves=planning_order_only destroys=centres_endpoints_owner_transport_LRC_predicate")
    print("PASS: exact components, radii, no-switch gaps, packets, and boundary")


if __name__ == "__main__":
    main()
