#!/usr/bin/env python3
"""Exact companion for THM-2587: deep-root/co-shift incidence wall.

On the canonical THM-2584 ``sigma={b}``, depth-five packet, write

    2x = (z+r)/13 mod 1,     0 <= z < 1,

where ``r=floor(26x) mod 13`` is the absolute deepest-speed cell.  The
deep-coordinate translate with label ``tau`` is dangerous exactly when

    d_1(2x-tau/13) = 1.

The pointwise incidence is

    tau=r       if z<13/14,
    tau=r+1     if z>1/14.

Thus the middle wall has two incident labels.  The script computes the fine
``(collision displacement, absolute root, translate, owner cell)`` tensor by
two exact routes: root-chart wall integration and one-circle physical
integration.  It also computes the danger and safe translate marginals,
checks the theta-selector no-go, and exhausts their joint cyclotomic spectra.

Scope: this is the projection of the deep coordinate only.  A lawful target
co-shift is the paired dipole ``e_(c_3)-e_(k_b)``; the compensating graft
coordinate is absent here.  No THM-2365 target current or row exclusion is
claimed.
"""

from fractions import Fraction

import lrc14_base_only_bridge_opus_20260728 as base
import lrc14_b_r5_owner_clock_host_thm2581 as host
import lrc14_b_r5_theta_target_tensor_thm2584 as root_tensor


P = 13
QMOD = 7
GRID = base.T_DEN
RPACK = 13**2
DEPTH = 13**5
DEN = RPACK * DEPTH * DEPTH * GRID
ZONE_NAMES = ("low", "middle", "high")


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def danger_intervals(speed, shift):
    """Half-open model of ||speed*x-shift/13||<1/14 on the exact grid."""
    require(GRID % (182 * speed) == 0, "grid does not resolve danger comb")
    unit = GRID // (182 * speed)
    start_residue = (-13 + 14 * shift) % 182
    out = []
    for turn in range(speed):
        left = (start_residue + 182 * turn) * unit
        right = left + 26 * unit
        if right <= GRID:
            out.append((left, right))
        else:
            out.append((left, GRID))
            out.append((0, right - GRID))
    out.sort()
    return out


def wall_pieces():
    """(theta, zone, y-intervals), where z={2y}."""
    require(GRID % 28 == 0, "grid does not resolve the 1/14 wall")
    return (
        (0, 0, [(0, GRID // 28)]),
        (0, 1, [(GRID // 28, 13 * GRID // 28)]),
        (0, 2, [(13 * GRID // 28, GRID // 2)]),
        (1, 0, [(GRID // 2, 15 * GRID // 28)]),
        (1, 1, [(15 * GRID // 28, 27 * GRID // 28)]),
        (1, 2, [(27 * GRID // 28, GRID)]),
    )


def root_chart_route(E, Q):
    Ust, Uv, Vst, Vv = root_tensor.packet_profiles(E, Q)
    arrival = [base.extract_window(Ust, Uv, u, P, GRID) for u in range(P)]
    source = [base.extract_window(Vst, Vv, u, P, GRID) for u in range(P)]
    owner = base.clock_cells(P, QMOD, GRID, P)

    # wall[theta][zone][s][r][ell]
    wall = [
        [[[[0] * QMOD for _ in range(P)] for _ in range(P)]
         for _ in range(3)]
        for _ in range(2)
    ]
    # fine[s][r][tau][ell]
    fine = [[[[0] * QMOD for _ in range(P)] for _ in range(P)]
            for _ in range(P)]
    old = [[0] * QMOD for _ in range(P)]

    for u in range(P):
        source_starts, source_values = source[u]
        for s in range(P):
            v = (u + s) % P
            arrival_starts, arrival_values = arrival[v]
            ps, pv, pc, _ = base.product_cum(
                arrival_starts,
                arrival_values,
                source_starts,
                source_values,
                GRID,
            )
            for ell in range(QMOD):
                for theta, zone, interval in wall_pieces():
                    cell = root_tensor.intersect_lists(owner[ell], interval)
                    mass = root_tensor.integrate_profile(ps, pv, pc, cell)
                    if not mass:
                        continue
                    root = (2 * v + theta) % P
                    wall[theta][zone][s][root][ell] += mass
                    old[s][ell] += mass
                    offsets = (0,) if zone == 0 else ((0, 1) if zone == 1 else (1,))
                    for offset in offsets:
                        tau = (root + offset) % P
                        fine[s][root][tau][ell] += mass
    return (Ust, Uv, Vst, Vv), wall, fine, old


def physical_route(profiles):
    Ust, Uv, Vst, Vv = profiles
    owner_pull = base.clock_cells(P, QMOD, GRID, P * P)
    root_cells = root_tensor.deep_cells()
    translates = [danger_intervals(2, tau) for tau in range(P)]

    cells = [[[
        root_tensor.intersect_lists(
            root_tensor.intersect_lists(owner_pull[ell], root_cells[root]),
            translates[tau],
        )
        for tau in range(P)] for root in range(P)] for ell in range(QMOD)]

    fine = [[[[0] * QMOD for _ in range(P)] for _ in range(P)]
            for _ in range(P)]
    old = [[0] * QMOD for _ in range(P)]
    for s in range(P):
        rotated_starts, rotated_values = base.rotate_profile(
            Vst, Vv, s * (GRID // P), GRID
        )
        ps, pv, pc, _ = base.product_cum(
            Ust, Uv, rotated_starts, rotated_values, GRID
        )
        for ell in range(QMOD):
            old[s][ell] = P * root_tensor.integrate_profile(
                ps, pv, pc, owner_pull[ell]
            )
            for root in range(P):
                for tau in range(P):
                    fine[s][root][tau][ell] = P * root_tensor.integrate_profile(
                        ps, pv, pc, cells[ell][root][tau]
                    )
    return fine, old


def marginal(fine):
    return [
        [
            [sum(fine[s][root][tau][ell] for root in range(P))
             for ell in range(QMOD)]
            for tau in range(P)
        ]
        for s in range(P)
    ]


def bucket(table, k, h, ell=None):
    out = [0] * P
    owners = range(QMOD) if ell is None else (ell,)
    for s in range(P):
        for tau in range(P):
            out[(-k * s - h * tau) % P] += sum(
                table[s][tau][owner] for owner in owners
            )
    return out


def centered_bucket(table, k, h, ell):
    cell = bucket(table, k, h, ell)
    total = bucket(table, k, h)
    # Seven times the centred transform; multiplication by 7 is harmless.
    return [QMOD * x - y for x, y in zip(cell, total)]


def spectrum(table):
    global_buckets = []
    cell_buckets = []
    centered_buckets = []
    for k in range(P):
        for h in range(P):
            value = bucket(table, k, h)
            require(not base.is_zero_b(value), "global Fourier zero")
            global_buckets.append(value)
            for ell in range(QMOD):
                value = bucket(table, k, h, ell)
                require(not base.is_zero_b(value), "cellwise Fourier zero")
                cell_buckets.append(value)
                value = centered_bucket(table, k, h, ell)
                require(not base.is_zero_b(value), "centred Fourier zero")
                centered_buckets.append(value)
    return global_buckets, cell_buckets, centered_buckets


def eval_bucket(value, modulus, root):
    return sum(
        (coefficient % modulus) * pow(root, exponent, modulus)
        for exponent, coefficient in enumerate(value)
    ) % modulus


def modular_cover(tagged_buckets):
    witnesses = ((79, 18), (131, 107))
    remaining = set(range(len(tagged_buckets)))
    caught = []
    for modulus, root in witnesses:
        require(pow(root, P, modulus) == 1 and root != 1,
                "listed modular root is not primitive of order 13")
        require(DEN % modulus != 0, "modular witness divides denominator")
        hit = {
            index for index in remaining
            if eval_bucket(tagged_buckets[index], modulus, root)
        }
        caught.append((modulus, len(hit)))
        remaining -= hit
    require(not remaining, "fixed finite-field witnesses missed a transform")
    return caught


def main():
    print("== THM-2587: deep-root/co-shift incidence wall ==")
    print("row:", base.W)
    print("sigma={b}, K=2, r=5; deep coordinate C/d=2")

    E = base.build_set(base.PAT_E, base.ZELL)
    Q = base.build_set(host.PAT_QB, base.ZELL)
    profiles, wall, fine1, old1 = root_chart_route(E, Q)
    fine2, old2 = physical_route(profiles)
    require(fine1 == fine2, "root-chart/physical fine incidence mismatch")
    require(old1 == old2, "root-chart/physical old host mismatch")
    print("two exact routes agree on all 13^3x7 fine incidence cells: PASS")
    print("common rational denominator:", DEN)

    # Pointwise law and the actual four-edge graph.
    edges = sorted({
        (root, tau)
        for s in range(P) for root in range(P) for tau in range(P)
        for ell in range(QMOD)
        if fine1[s][root][tau][ell]
    })
    require(edges == [(0, 0), (0, 1), (12, 0), (12, 12)],
            "incidence graph is not the four-edge toothpick")
    require(all(
        sum(fine1[s][root][tau][ell] for s in range(P)) > 0
        for root, tau in edges for ell in range(QMOD)
    ), "an incidence edge disappeared in an owner cell")
    print("positive (absolute root, translated danger label) edges:", edges)
    print("every one of the four edges is positive in every owner cell: PASS")

    wall_counts = [[0] * 3 for _ in range(2)]
    wall_masses = [[0] * 3 for _ in range(2)]
    wall_roots = [[set() for _ in range(3)] for _ in range(2)]
    for theta in range(2):
        for zone in range(3):
            for s in range(P):
                for root in range(P):
                    for ell in range(QMOD):
                        value = wall[theta][zone][s][root][ell]
                        if value:
                            wall_counts[theta][zone] += 1
                            wall_masses[theta][zone] += value
                            wall_roots[theta][zone].add(root)
            require(wall_masses[theta][zone] > 0, "empty theta/wall state")
    require(wall_counts == [[48, 154, 48], [48, 154, 48]],
            "theta/wall support census changed")
    require([[sorted(values) for values in row] for row in wall_roots]
            == [[[0], [0, 12], [12]], [[0], [0, 12], [12]]],
            "theta/wall root support changed")
    print("positive cells by theta x (low,middle,high):", wall_counts)
    print("mass numerators by theta x wall state:", wall_masses)

    # Any diagonal-translation-equivariant selector has tau-root=c_theta.
    # Positive low mass forces c_theta=0; positive high mass forces c_theta=1.
    # The pure low and pure high packets are the two hostile controls.
    for theta in range(2):
        low_offsets = {0}
        high_offsets = {1}
        require(low_offsets != high_offsets and wall_masses[theta][0] > 0
                and wall_masses[theta][2] > 0,
                "theta-selector hostile lost one side")
    print("pure-low selector tau=r: PASS; pure-high selector tau=r+1: PASS")
    print("each theta fibre contains both controls: no equivariant tau(r,theta)")
    print("middle wall has two edges; one overlap-edge bit is losslessly minimal")

    danger = marginal(fine1)
    safe = [
        [[old1[s][ell] - danger[s][tau][ell] for ell in range(QMOD)]
         for tau in range(P)]
        for s in range(P)
    ]
    require(all(value >= 0 for plane in safe for row in plane for value in row),
            "safe complement became negative")
    require(all(danger[s][0][ell] == old1[s][ell]
                for s in range(P) for ell in range(QMOD)),
            "original translate does not recover THM-2581 host")
    require(all(safe[s][0][ell] == 0
                for s in range(P) for ell in range(QMOD)),
            "original translate gained a safe cell")

    danger_support = [
        sum(danger[s][tau][ell] > 0
            for s in range(P) for ell in range(QMOD))
        for tau in range(P)
    ]
    safe_support = [
        sum(safe[s][tau][ell] > 0
            for s in range(P) for ell in range(QMOD))
        for tau in range(P)
    ]
    require(danger_support == [84, 77] + [0] * 10 + [77],
            "danger translate support changed")
    require(safe_support == [0] + [84] * 12,
            "safe translate support changed")
    danger_mass = [
        sum(danger[s][tau][ell] for s in range(P) for ell in range(QMOD))
        for tau in range(P)
    ]
    safe_mass = [
        sum(safe[s][tau][ell] for s in range(P) for ell in range(QMOD))
        for tau in range(P)
    ]
    print("positive (s,ell) cells by danger translate:", danger_support)
    print("positive (s,ell) cells by safe translate:", safe_support)
    print("danger translate masses:", [str(Fraction(value, DEN)) for value in danger_mass])
    print("safe translate masses:", [str(Fraction(value, DEN)) for value in safe_mass])

    danger_spectrum = spectrum(danger)
    safe_spectrum = spectrum(safe)
    print("danger joint Fourier: global/cell/centred = 169/1183/1183 nonzero")
    print("safe joint Fourier: global/cell/centred = 169/1183/1183 nonzero")
    tagged = [
        value
        for family in (danger_spectrum, safe_spectrum)
        for group in family
        for value in group
    ]
    cover = modular_cover(tagged)
    print("finite-field independent cover (prime,new witnesses):", cover)
    print("total transform buckets certified:", len(tagged))

    print("\nsource: THM-2584 absolute root/displacement/owner tensor")
    print("target: projected c3 danger/safe translate incidence")
    print("preserved: one common positive b/r=5 ancestry and deep membership")
    print("lost: within-cell wall if roots alone; always lost: paired graft k_b")
    print("consequence: projected deep action is saturated, full target action remains open")
    print("\nall exact checks passed")


if __name__ == "__main__":
    main()
