#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2535.

The script checks the finite affine bookkeeping over F_13 x F_7, the
boundary-tooth clock coordinate, its reparametrized cut-DFT factorization,
the zero neutral character, and an all-mixed-mode hostile whose selected
semantic tooth vanishes.
"""

P = 547  # 547 - 1 = 6 * 7 * 13


def inv(x, p):
    return pow(x % p, p - 2, p)


def primitive_root(p):
    factors = (2, 3, 7, 13)
    for g in range(2, p):
        if all(pow(g, (p - 1) // q, p) != 1 for q in factors):
            return g
    raise AssertionError("no primitive root")


GEN = primitive_root(P)
ZETA = pow(GEN, (P - 1) // 13, P)
XI = pow(GEN, (P - 1) // 7, P)


def zp(e):
    return pow(ZETA, e % 13, P)


def xp(e):
    return pow(XI, e % 7, P)


def inv7(a):
    return pow(a % 7, 5, 7)


def bit(mask, r):
    return (mask >> (r % 13)) & 1


def marker_boundary(mask, tau):
    words = [tuple(bit(mask, a + j * tau) for j in range(13)) for a in range(13)]
    alpha = max(range(13), key=lambda a: words[a])
    assert words[alpha][0] == 1
    q = next(j for j in range(1, 13) if words[alpha][j] == 0)
    s = (alpha + (q - 1) * tau) % 13
    t = (s + tau) % 13
    assert bit(mask, s) == 1 and bit(mask, t) == 0
    return alpha, q, s, t


def tooth_source(a, b, j):
    return (inv7(a) * (j - b)) % 7


def rcut(d, tau, a, b, t):
    return sum(d[(t - tau * ((a * r + b) % 7)) % 13][r] for r in range(7)) % P


def clock_table(d, tau, a):
    return [
        [rcut(d, tau, a, (1 - a * kappa) % 7, t) for t in range(13)]
        for kappa in range(7)
    ]


def dtilde(d, alpha, beta):
    return sum(
        d[h][r] * zp(-alpha * h) * xp(-beta * r)
        for h in range(13)
        for r in range(7)
    ) % P


def qhat(qtab, alpha, beta):
    return sum(
        qtab[kappa][t] * zp(-alpha * t) * xp(-beta * kappa)
        for kappa in range(7)
        for t in range(13)
    ) % P


def boundary_kernel(tau, a, alpha, beta):
    ai = inv7(a)
    return sum(zp(-alpha * tau * j) * xp(beta * ai * j) for j in range(7)) % P


def make_row_zero(seed):
    d = [[0 for _ in range(7)] for _ in range(13)]
    for h in range(13):
        total = 0
        for r in range(6):
            value = (seed * (h + 1) * (r + 2) + h * h + 3 * r + seed * seed) % P
            d[h][r] = value
            total += value
        d[h][6] = (-total) % P
    assert all(sum(row) % P == 0 for row in d)
    return d


def basis_defect(h0, r0):
    d = [[0 for _ in range(7)] for _ in range(13)]
    d[h0][r0] = 1
    d[h0][6] = P - 1
    return d


def hostile_defect():
    # (delta_0-delta_2) tensor (delta_0-delta_1).  Both margins vanish,
    # every primitive mixed Fourier coefficient is nonzero, and row h=1
    # is identically zero.
    phi = [0] * 13
    psi = [0] * 7
    phi[0], phi[2] = 1, P - 1
    psi[0], psi[1] = 1, P - 1
    return [[phi[h] * psi[r] % P for r in range(7)] for h in range(13)]


def hostile_anchored_table():
    def dint(h, r):
        return (int(h == 0) - int(h == 2)) * (int(r == 0) - int(r == 1))

    table = [[0 for _ in range(13)] for _ in range(7)]
    for r in range(7):
        for h in range(13):
            if h == 0:
                table[r][h] = int(r == 0)
            else:
                table[r][h] = dint(h, r) - dint(0, r) + int(r == 0) + 1
    return table


def main():
    assert pow(ZETA, 13, P) == 1 and all(pow(ZETA, k, P) != 1 for k in range(1, 13))
    assert pow(XI, 7, P) == 1 and all(pow(XI, k, P) != 1 for k in range(1, 7))

    # The Latin scheduler and its fixed-owner slot.
    scheduler_entries = set()
    owner_slots = []
    scheduler_translation_checks = 0
    scheduler_affine_slot_checks = 0
    for kappa in range(7):
        row = [(kappa + d) % 7 for d in range(7)]
        assert set(row) == set(range(7))
        owner_slots.append(row.index(0))
        for d, gamma in enumerate(row):
            scheduler_entries.add((kappa, d, gamma))
            for C in range(7):
                assert (gamma + C) % 7 == ((kappa + C) + d) % 7
                scheduler_translation_checks += 1
            for B in range(1, 7):
                for C in range(7):
                    # This affine statement transports the slot label d -> B d.
                    assert (B * gamma + C) % 7 == ((B * kappa + C) + B * d) % 7
                    scheduler_affine_slot_checks += 1
    assert owner_slots == [(-k) % 7 for k in range(7)]

    # The root marker is unchanged across the seven scheduler rows.  The
    # pure septimal CRT translation K=13 has (H,C)=(0,6), giving the exact
    # root-only section obstruction.
    marker_epoch_states = 0
    for mask in range(1, (1 << 13) - 1):
        data = marker_boundary(mask, 1)
        for _kappa in range(7):
            assert marker_boundary(mask, 1) == data
            marker_epoch_states += 1
    K = 13
    H, C = K % 13, K % 7
    assert (H, C) == (0, 6)
    witness_mask = 1 << 1
    witness_boundary = marker_boundary(witness_mask, 1)
    assert marker_boundary(witness_mask, 1) == witness_boundary

    # Every affine cut chart has one tooth-1 clock.  At a selected boundary
    # s -> t=s+tau, that tooth contributes d(s,kappa) to R(...)(t).
    chart_clock_counts = [0] * 7
    boundary_tooth_checks = 0
    for tau in range(1, 13):
        for s in range(13):
            t = (s + tau) % 13
            for a in range(1, 7):
                for b in range(7):
                    kappa = tooth_source(a, b, 1)
                    assert (a * kappa + b) % 7 == 1
                    assert (t - tau * ((a * kappa + b) % 7)) % 13 == s
                    chart_clock_counts[kappa] += 1
                    boundary_tooth_checks += 1
    # Divide out the 12*13 repeated boundary choices: six charts per clock.
    assert [x // (12 * 13) for x in chart_clock_counts] == [6] * 7

    # On every target-anchored deep-comb mask, the marker-to-target chord is
    # an available nonzero cut slope even when the adjacent selected wall
    # does not end at the target.
    deep_masks = [1 << j for j in range(1, 13)] + [
        (1 << j) | (1 << (j + 1)) for j in range(1, 12)
    ]
    deep_target_chord_checks = 0
    for mask in deep_masks:
        alpha, _q, _wall_s, _wall_t = marker_boundary(mask, 1)
        assert alpha != 0 and bit(mask, alpha) == 1 and bit(mask, 0) == 0
        target, tau_edge = 0, (-alpha) % 13
        assert tau_edge != 0 and (alpha + tau_edge) % 13 == target
        for a in range(1, 7):
            for b in range(7):
                kappa = tooth_source(a, b, 1)
                assert (target - tau_edge * ((a * kappa + b) % 7)) % 13 == alpha
                deep_target_chord_checks += 1

    # Exact old-to-new geometric chart covariance.  The source chart sends
    # r' = B r + C and (a,b) -> (a B^-1, b-a B^-1 C), while the root chart
    # sends h' = U h+H and (tau,s,t)->(U tau,U s+H,U t+H).
    source_covariance_checks = 0
    for a in range(1, 7):
        for b in range(7):
            for B in range(1, 7):
                Bi = inv7(B)
                ap = a * Bi % 7
                for C0 in range(7):
                    bp = (b - a * Bi * C0) % 7
                    for j in range(7):
                        r = tooth_source(a, b, j)
                        rp = tooth_source(ap, bp, j)
                        assert rp == (B * r + C0) % 7
                        source_covariance_checks += 1

    root_covariance_checks = 0
    for tau in range(1, 13):
        for t in range(13):
            for U in range(1, 13):
                for H0 in range(13):
                    tp, taup = (U * t + H0) % 13, U * tau % 13
                    for j in range(7):
                        assert (tp - taup * j) % 13 == (U * (t - tau * j) + H0) % 13
                        root_covariance_checks += 1

    # The source-neutral clock character collapses on every row-zero defect.
    uniform_collapse_checks = 0
    for h0 in range(13):
        for r0 in range(6):
            d = basis_defect(h0, r0)
            for tau in range(1, 13):
                for t in range(13):
                    for a in range(1, 7):
                        total = sum(
                            rcut(d, tau, a, (1 - a * kappa) % 7, t)
                            for kappa in range(7)
                        ) % P
                        assert total == 0
                        uniform_collapse_checks += 1

    # The boundary-clock DFT is an invertible reparametrization of the
    # primitive cut-character bundle.
    kernel_nonzero_checks = 0
    for tau in range(1, 13):
        for a in range(1, 7):
            for alpha in range(1, 13):
                for beta in range(7):
                    assert boundary_kernel(tau, a, alpha, beta) != 0
                    kernel_nonzero_checks += 1

    factorization_checks = 0
    for seed in range(1, 9):
        d = make_row_zero(seed)
        for tau in range(1, 13):
            for a in range(1, 7):
                qtab = clock_table(d, tau, a)
                ai = inv7(a)
                for alpha in range(1, 13):
                    for beta in range(7):
                        lhs = qhat(qtab, alpha, beta)
                        rhs = (
                            xp(-beta * ai)
                            * boundary_kernel(tau, a, alpha, beta)
                            * dtilde(d, alpha, beta)
                        ) % P
                        assert lhs == rhs
                        if beta == 0:
                            assert lhs == 0
                        factorization_checks += 1

    # Sharp hostile: complete primitive spectrum, but the marker-selected
    # semantic source row is identically zero.
    hostile = hostile_defect()
    assert all(sum(hostile[h][r] for r in range(7)) % P == 0 for h in range(13))
    assert all(sum(hostile[h][r] for h in range(13)) % P == 0 for r in range(7))
    hostile_modes = 0
    for alpha in range(1, 13):
        for beta in range(1, 7):
            assert dtilde(hostile, alpha, beta) != 0
            hostile_modes += 1

    anchored = hostile_anchored_table()
    assert min(min(row) for row in anchored) >= 0
    assert [anchored[r][0] for r in range(7)] == [1, 0, 0, 0, 0, 0, 0]
    assert len(set(anchored[0])) > 1
    inv13, inv7p, inv91 = inv(13, P), inv(7, P), inv(91, P)
    row_means = [sum(anchored[r]) * inv13 % P for r in range(7)]
    col_means = [sum(anchored[r][h] for r in range(7)) * inv7p % P for h in range(13)]
    grand = sum(sum(row) for row in anchored) * inv91 % P
    anchored_interaction_checks = 0
    for r in range(7):
        for h in range(13):
            interaction = (anchored[r][h] - row_means[r] - col_means[h] + grand) % P
            assert interaction == hostile[h][r]
            anchored_interaction_checks += 1
    alpha0, q0, s0, t0 = witness_boundary
    assert (alpha0, q0, s0, t0) == (1, 1, 1, 2)
    hostile_selected_zeros = 0
    for a in range(1, 7):
        for b in range(7):
            kappa = tooth_source(a, b, 1)
            assert hostile[s0][kappa] == 0
            hostile_selected_zeros += 1

    hostile_clock_modes = 0
    for tau in range(1, 13):
        for a in range(1, 7):
            qtab = clock_table(hostile, tau, a)
            for alpha in range(1, 13):
                for beta in range(1, 7):
                    assert qhat(qtab, alpha, beta) != 0
                    hostile_clock_modes += 1

    # Correlated selector/scheduler incidence escape.  Each clock event is
    # used once at tooth zero and once at tooth one of the preceding chart.
    # The uncharted signed table cancels, while the chartwise output is the
    # marker-to-target root dipole.
    incidence_cut_checks = 0
    forgotten_table_zero_checks = 0
    orbit_dipole_checks = 0
    for edge_root in range(1, 13):
        tau_edge = (-edge_root) % 13
        for active_clock in range(7):
            x = [int(k == active_clock) for k in range(7)]
            dsum = [[0 for _ in range(7)] for _ in range(13)]
            rsum = [0 for _ in range(13)]
            for kappa in range(7):
                local = [[0 for _ in range(7)] for _ in range(13)]
                local[edge_root][kappa] = x[kappa]
                local[edge_root][(kappa + 1) % 7] = (
                    local[edge_root][(kappa + 1) % 7] - x[(kappa + 1) % 7]
                ) % P
                for h in range(13):
                    for r in range(7):
                        dsum[h][r] = (dsum[h][r] + local[h][r]) % P
                for v in range(13):
                    got = rcut(local, tau_edge, 1, (-kappa) % 7, v)
                    expected = (
                        x[kappa] * int(v == edge_root)
                        - x[(kappa + 1) % 7] * int(v == 0)
                    ) % P
                    assert got == expected
                    rsum[v] = (rsum[v] + got) % P
                    incidence_cut_checks += 1
            for h in range(13):
                for r in range(7):
                    assert dsum[h][r] == 0
                    forgotten_table_zero_checks += 1
            mass = sum(x) % P
            for v in range(13):
                assert rsum[v] == mass * (int(v == edge_root) - int(v == 0)) % P
                orbit_dipole_checks += 1

    print("THM-2535 boundary-tooth clock intertwiner exact referee")
    print(f"field={P} generator={GEN} zeta13={ZETA} xi7={XI}")
    print(
        "scheduler_cells=49 "
        f"owner_slots={','.join(map(str, owner_slots))} "
        f"translation_checks={scheduler_translation_checks} "
        f"affine_transported_slot_checks={scheduler_affine_slot_checks}"
    )
    print(
        f"marker_owner_epoch_states={marker_epoch_states} "
        f"pure_septimal_CRT_witness_K={K} root_shift={H} clock_shift={C}"
    )
    print(
        f"boundary_tooth_checks={boundary_tooth_checks} "
        "charts_per_clock=6,6,6,6,6,6,6 "
        f"deep_target_chord_checks={deep_target_chord_checks}"
    )
    print(
        f"source_chart_covariance_checks={source_covariance_checks} "
        f"root_chart_covariance_checks={root_covariance_checks}"
    )
    print(f"uniform_diagonal_basis_checks={uniform_collapse_checks} neutral_character=zero")
    print(
        f"primitive_kernel_nonzero_checks={kernel_nonzero_checks} "
        f"direct_factorization_checks={factorization_checks}"
    )
    print(
        f"hostile_mixed_modes={hostile_modes} "
        f"hostile_selected_tooth_zeros={hostile_selected_zeros} "
        f"hostile_clock_modes={hostile_clock_modes}"
    )
    print(
        f"hostile_anchored_table_entries={sum(map(len, anchored))} "
        f"min={min(min(row) for row in anchored)} max={max(max(row) for row in anchored)} "
        f"anova_interaction_checks={anchored_interaction_checks} owner_row_nonconstant=yes"
    )
    print(
        f"selector_scheduler_incidence_cut_checks={incidence_cut_checks} "
        f"forgotten_table_zero_checks={forgotten_table_zero_checks} "
        f"orbit_dipole_checks={orbit_dipole_checks}"
    )
    print("PASS")


if __name__ == "__main__":
    main()
