#!/usr/bin/env python3
"""Exact probe of the endogenous r=5 deep-root tensor on the b-word packet.

The same-fibre THM-2581 host retains the collision displacement ``s`` and
owner-clock cell ``ell`` but sums the absolute deep root.  At depth d=13^5
and deepest speed C=2d, the THM-2471 stalk gives, for the arrival root v,

    t = 2v + floor(2y) (mod 13),    theta=t-2v=floor(2y) in {0,1}.

This companion retains t before integration:

    C_ell(s,t) = int_{H_ell} sum_u A(y,u+s)F(y,u)
                   1[t=2(u+s)+floor(2y)] dy.

It computes the positive 13 x 13 x 7 tensor by two exact routes.  The second
route uses x=(y+u+s)/13 and the identity

    C_ell(s,t)=13 int U(x)V(x-s/13)
       1[frac(169x) in win_ell] 1[floor(26x)=t (mod 13)] dx.

The script also retains the two theta halves separately and audits the full
two-dimensional cyclotomic transform in (s,t).  All decisions are exact.

Scope: canonical typed row, sigma={b}, K=2, r=5 only.  This is a deep-root
ancestry tensor, not yet an identification with the THM-2365 target shift or
a physical LRC current.
"""

from bisect import bisect_right
from fractions import Fraction

import lrc14_base_only_bridge_opus_20260728 as base
import lrc14_b_r5_owner_clock_host_thm2581 as host


P = 13
QMOD = 7
GRID = base.T_DEN
RPACK = 13**2
DEPTH = 13**5


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def intersect_lists(left, right):
    """Intersection of two sorted half-open interval unions."""
    out = []
    i = j = 0
    while i < len(left) and j < len(right):
        a = max(left[i][0], right[j][0])
        b = min(left[i][1], right[j][1])
        if a < b:
            out.append((a, b))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return out


def integrate_profile(ps, pv, pc, intervals):
    return sum(
        base.pquery(ps, pv, pc, b) - base.pquery(ps, pv, pc, a)
        for a, b in intervals
    )


def packet_profiles(E_iv, Q_iv):
    """Rebuild U=P_d(1_Q P_R 1_E), V=P_d1_E exactly."""
    n_st, n_v = base.weighted_fold(
        [(a, b, 1) for a, b in E_iv], RPACK, GRID
    )
    fp = []
    for qa, qb in Q_iv:
        i = bisect_right(n_st, qa) - 1
        while True:
            pa = n_st[i]
            pb = n_st[i + 1] if i + 1 < len(n_st) else GRID
            a = max(qa, pa)
            b = min(qb, pb)
            if a < b and n_v[i]:
                fp.append((a, b, n_v[i]))
            if pb >= qb:
                break
            i += 1
    Ust, Uv = base.weighted_fold(fp, DEPTH, GRID)
    Vst, Vv = base.weighted_fold(
        [(a, b, 1) for a, b in E_iv], DEPTH, GRID
    )
    return Ust, Uv, Vst, Vv


def deep_cells():
    """Cells floor(26x)=t mod 13; two intervals of length 1/26."""
    require(GRID % (2 * P) == 0, "grid does not resolve deep-root cells")
    unit = GRID // (2 * P)
    cells = []
    for t in range(P):
        cells.append([
            (t * unit, (t + 1) * unit),
            ((t + P) * unit, (t + P + 1) * unit),
        ])
    base.check_partition(cells, GRID)
    return cells


def exact_tensor(E_iv, Q_iv):
    Ust, Uv, Vst, Vv = packet_profiles(E_iv, Q_iv)
    owner = base.clock_cells(P, QMOD, GRID, P)
    owner_pull = base.clock_cells(P, QMOD, GRID, P * P)
    halves = [[(0, GRID // 2)], [(GRID // 2, GRID)]]
    owner_half = [
        [intersect_lists(owner[ell], halves[eps]) for eps in range(2)]
        for ell in range(QMOD)
    ]
    deep = deep_cells()
    arrival = [[(v * (GRID // P), (v + 1) * (GRID // P))]
               for v in range(P)]
    owner_deep = [
        [intersect_lists(owner_pull[ell], deep[t]) for t in range(P)]
        for ell in range(QMOD)
    ]
    owner_arrival_deep = [
        [
            [intersect_lists(owner_deep[ell][t], arrival[v])
             for t in range(P)]
            for v in range(P)
        ]
        for ell in range(QMOD)
    ]

    Aprof = [base.extract_window(Ust, Uv, u, P, GRID) for u in range(P)]
    Fprof = [base.extract_window(Vst, Vv, u, P, GRID) for u in range(P)]

    # Route 1: root charts in y, split by owner cell and theta half.
    C1 = [[[0] * QMOD for _ in range(P)] for _ in range(P)]
    K1 = [[[[0] * QMOD for _ in range(P)] for _ in range(P)]
          for _ in range(P)]
    CE = [[[0] * QMOD for _ in range(P)] for _ in range(2)]
    for u0 in range(P):
        sF, vF = Fprof[u0]
        for s in range(P):
            u1 = (u0 + s) % P
            sA, vA = Aprof[u1]
            ps, pv, pc, _ = base.product_cum(sA, vA, sF, vF, GRID)
            for ell in range(QMOD):
                for eps in range(2):
                    mass = integrate_profile(ps, pv, pc, owner_half[ell][eps])
                    CE[eps][s][ell] += mass
                    t = (2 * u1 + eps) % P
                    C1[s][t][ell] += mass
                    K1[s][u1][t][ell] += mass

    # Route 2: one physical x-coordinate, absolute deep-root digit retained.
    C2 = [[[0] * QMOD for _ in range(P)] for _ in range(P)]
    K2 = [[[[0] * QMOD for _ in range(P)] for _ in range(P)]
          for _ in range(P)]
    for s in range(P):
        rVst, rVv = base.rotate_profile(Vst, Vv, s * (GRID // P), GRID)
        ps, pv, pc, _ = base.product_cum(Ust, Uv, rVst, rVv, GRID)
        for ell in range(QMOD):
            for t in range(P):
                C2[s][t][ell] = P * integrate_profile(
                    ps, pv, pc, owner_deep[ell][t]
                )
                for v in range(P):
                    K2[s][v][t][ell] = P * integrate_profile(
                        ps, pv, pc, owner_arrival_deep[ell][v][t]
                    )
    require(C1 == C2, "two-route deep-root tensor mismatch")
    require(K1 == K2, "two-route arrival/deep-root tensor mismatch")
    return C1, CE, K1


def bucket_transform(C, k, h, ell=None):
    """Integer buckets for sum C(s,t,ell) zeta^(-ks-ht)."""
    bucket = [0] * P
    ells = range(QMOD) if ell is None else (ell,)
    for s in range(P):
        for t in range(P):
            exponent = (-k * s - h * t) % P
            bucket[exponent] += sum(C[s][t][e] for e in ells)
    return bucket


def vector_transform(D, k, ell=None):
    bucket = [0] * P
    ells = range(QMOD) if ell is None else (ell,)
    for s in range(P):
        bucket[(-k * s) % P] += sum(D[s][e] for e in ells)
    return bucket


def projective_zeros(C, ell=None):
    zeros = []
    for k in range(P):
        for h in range(P):
            if base.is_zero_b(bucket_transform(C, k, h, ell)):
                zeros.append((k, h))
    return zeros


def main():
    print("== THM-2584 probe: b/r=5 theta and absolute deep-root tensor ==")
    print("row:", base.W)
    print("sigma={b}, K=2, r=5, d=13^5, Cdeep=2d")
    print("theta=t-2v=floor(2y) in {0,1}; all restrictions are pre-marginal")

    E = base.build_set(base.PAT_E, base.ZELL)
    QB = base.build_set(host.PAT_QB, base.ZELL)
    C, CE, K = exact_tensor(E, QB)
    DENC = RPACK * DEPTH * DEPTH * GRID
    print("two exact routes agree on all 13x13x7 tensor entries: PASS")
    print("common denominator:", DENC)

    old = host.exact_host(E, QB, RPACK, DEPTH)[1]
    for s in range(P):
        for ell in range(QMOD):
            require(sum(C[s][t][ell] for t in range(P)) == old[s][ell],
                    "deep-root marginal does not recover THM-2581 host")
            require(CE[0][s][ell] + CE[1][s][ell] == old[s][ell],
                    "theta marginal does not recover THM-2581 host")
    require(all(C[0][t][ell] == 0 for t in range(P)
                for ell in range(QMOD)), "s=0 diagonal is not zero")
    print("t-marginal and theta-marginal recover THM-2581 C_ell(s): PASS")
    print("C_ell(0,t)=0 cellwise for all 91 (t,ell): PASS")

    nonzero = sum(C[s][t][ell] != 0 for s in range(P)
                  for t in range(P) for ell in range(QMOD))
    positive_per_t = [sum(C[s][t][ell] > 0 for s in range(P)
                          for ell in range(QMOD)) for t in range(P)]
    target_mass = [sum(C[s][t][ell] for s in range(P)
                       for ell in range(QMOD)) for t in range(P)]
    print(f"positive tensor cells: {nonzero}/1183")
    print("positive (s,ell) cells per deep root t:", positive_per_t)
    print("deep-root mass numerators:", target_mass)
    require(len(set(target_mass)) > 1, "absolute deep-root marginal is uniform")
    target_hat = []
    for h in range(P):
        b = [0] * P
        for t in range(P):
            b[(-h * t) % P] += target_mass[t]
        target_hat.append(b)
    require(all(not base.is_zero_b(target_hat[h]) for h in range(1, P)),
            "a nonzero absolute deep-root colour vanished")
    print("absolute deep-root marginal is nonuniform; all 12 colours survive: PASS")

    # The finer tensor keeps the arrival collision root v.  After the unit
    # reparametrization w=t/2=7t, theta=0 is the literal diagonal w=v and
    # theta=1 is the parallel rail w=v+7.
    edge_mass = [[sum(K[s][v][t][ell] for s in range(P)
                      for ell in range(QMOD))
                  for t in range(P)] for v in range(P)]
    edges = [(v, t, edge_mass[v][t]) for v in range(P) for t in range(P)
             if edge_mass[v][t]]
    print("positive arrival/deep-root edges (v,t,numerator):", edges)
    require(all(((7 * t - v) % P) in (0, 7) for v, t, _ in edges),
            "arrival/deep-root tensor escaped the two affine rails")
    rail_mass = [sum(m for v, t, m in edges if (7 * t - v) % P == shift)
                 for shift in (0, 7)]
    require(all(m > 0 for m in rail_mass), "one theta rail is empty")
    require(sum(rail_mass) == sum(target_mass), "rail masses lost service")
    print("in w=7t coordinates, support is w=v or w=v+7: PASS")
    print("theta=0/1 rail mass numerators:", rail_mass)
    for s in range(P):
        for t in range(P):
            for ell in range(QMOD):
                require(sum(K[s][v][t][ell] for v in range(P))
                        == C[s][t][ell], "arrival-root marginal lost C")
    print("two exact routes agree on the finer 13x13x13x7 tensor: PASS")

    expected_edges = {(0, 0), (6, 0), (6, 12), (12, 12)}
    hall_diagonal = []
    for ell in range(QMOD):
        edge_ell = {
            (v, t): sum(K[s][v][t][ell] for s in range(P))
            for v in range(P) for t in range(P)
            if sum(K[s][v][t][ell] for s in range(P))
        }
        require(set(edge_ell) == expected_edges,
                "an owner cell changed the four-edge support graph")
        # In w=7t coordinates the allowed edges are
        # 0->0, 6->0, 6->6, 12->6.  Deleting equal-label edges leaves
        # arrival vertex 0 with no neighbour, so weighted Hall forces its
        # entire positive margin onto the diagonal 0->0.
        hall_diagonal.append(edge_ell[(0, 0)])
    require(all(value > 0 for value in hall_diagonal),
            "an owner cell lost the forced diagonal")
    print("each owner cell has exactly the four-edge path")
    print("  v=0--w=0--v=6--w=6--v=12: PASS")
    print("after deleting v=w, S={0} has empty neighbourhood;")
    print("forced Hall diagonal numerators by owner cell:", hall_diagonal)

    global_zeros = projective_zeros(C)
    local_zeros = {ell: projective_zeros(C, ell) for ell in range(QMOD)}
    print("global 2D Fourier zeros (k,h):", global_zeros)
    print("per-owner-cell 2D Fourier zero counts:",
          [len(local_zeros[ell]) for ell in range(QMOD)])

    # Owner-centred two-dimensional transform, represented as 7B-B_global.
    centred_zero = []
    for k in range(P):
        for h in range(P):
            global_b = bucket_transform(C, k, h)
            for ell in range(QMOD):
                centred = base.bsub(
                    base.bscale(bucket_transform(C, k, h, ell), QMOD),
                    global_b,
                )
                if base.is_zero_b(centred):
                    centred_zero.append((k, h, ell))
    print(f"owner-centred 2D transform nonzero: {1183-len(centred_zero)}/1183")
    print("owner-centred zeros:", centred_zero)

    # Signed theta half current and its exact reflection law.
    D = [[CE[1][s][ell] - CE[0][s][ell] for ell in range(QMOD)]
         for s in range(P)]
    for s in range(P):
        for ell in range(QMOD):
            require(D[s][ell] == -D[(-s) % P][(-ell) % QMOD],
                    "signed theta reflection law failed")
    dzeros_global = [k for k in range(P)
                     if base.is_zero_b(vector_transform(D, k))]
    dzeros_local = [
        [k for k in range(P) if base.is_zero_b(vector_transform(D, k, ell))]
        for ell in range(QMOD)
    ]
    print("signed theta law D_ell(s)=-D_-ell(-s): PASS")
    print("signed theta global Fourier zeros k:", dzeros_global)
    print("signed theta per-cell Fourier zeros:", dzeros_local)
    print("signed theta numerator table (rows s, columns ell):")
    for s, row in enumerate(D):
        print(f"s={s:2d}: {row}")

    # Full tensor reflection: (s,t,ell)->(-s,-t-1,-ell).
    for s in range(P):
        for t in range(P):
            for ell in range(QMOD):
                require(C[s][t][ell]
                        == C[(-s) % P][(-t - 1) % P][(-ell) % QMOD],
                        "deep-root tensor reflection law failed")
    print("C_ell(s,t)=C_-ell(-s,-t-1): PASS")
    for s in range(P):
        for v in range(P):
            for t in range(P):
                for ell in range(QMOD):
                    require(K[s][v][t][ell]
                            == K[(-s) % P][(-v - 1) % P]
                                 [(-t - 1) % P][(-ell) % QMOD],
                            "fine tensor reflection law failed")
    print("K_ell(s,v,t)=K_-ell(-s,-v-1,-t-1): PASS")
    print("SCOPE: exact positive ancestry tensor on one typed b/K=2 packet;")
    print("no THM-2365 target identification, physical current, row exclusion,")
    print("or LRC(14) conclusion. All exact checks passed.")


if __name__ == "__main__":
    main()
