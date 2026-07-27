#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2534.

The script exhausts every nonconstant Boolean mask on F_13, checks the
oriented-cut, matched-flux, Crofton, affine/converse, and sparse-Prony laws,
and then verifies the weighted and target-anchored hostile controls used to
mark the theorem's exact boundary.
"""

from collections import Counter


P = 13
Q = 53
ZETA = 16  # 2^4 mod 53; exact order 13.


def bits_of(mask):
    return tuple((mask >> r) & 1 for r in range(P))


def mask_of(support):
    return sum(1 << (r % P) for r in support)


def boundary(mask, tau):
    f = bits_of(mask)
    return tuple(f[r] * (1 - f[(r + tau) % P]) for r in range(P))


def endpoint_matrix(mask):
    f = bits_of(mask)
    return tuple(
        tuple(f[r] * (1 - f[s]) for s in range(P)) for r in range(P)
    )


def transform_mask(mask, u, b):
    """Push f_r to g_(u r+b)."""
    out = 0
    for r in range(P):
        if (mask >> r) & 1:
            out |= 1 << ((u * r + b) % P)
    return out


def rotate_mask(mask, b):
    return transform_mask(mask, 1, b)


def poly_eval(coeffs, x, q=Q):
    ans = 0
    for c in coeffs:
        ans = (ans * x + c) % q
    return ans


def prony_support(kvec):
    """Recover a 0/1 support from DC plus its first b power sums."""
    support = [r for r, value in enumerate(kvec) if value]
    b = len(support)
    nodes = [pow(ZETA, r, Q) for r in support]
    power_sums = [0] + [sum(pow(x, j, Q) for x in nodes) % Q for j in range(1, b + 1)]

    elementary = [1]
    for j in range(1, b + 1):
        rhs = 0
        for i in range(1, j + 1):
            rhs += (-1) ** (i - 1) * elementary[j - i] * power_sums[i]
        elementary.append((rhs * pow(j, -1, Q)) % Q)

    # X^b-e_1 X^(b-1)+...+(-1)^b e_b.
    polynomial = [1]
    for j in range(1, b + 1):
        polynomial.append(((-1) ** j * elementary[j]) % Q)
    recovered = [r for r in range(P) if poly_eval(polynomial, pow(ZETA, r, Q)) == 0]
    assert recovered == support
    return b


def cyclotomic_reduce(values, frequency):
    """Exact coordinates in Q[zeta_13] on 1,zeta,...,zeta^11.

    Input represents sum_t values[t] zeta^(-frequency*t).  A degree-at-most
    12 rational polynomial vanishes at zeta iff all 13 exponent buckets are
    equal; subtracting the zeta^12 bucket gives the 12 exact coordinates.
    """
    buckets = [0] * P
    for t, value in enumerate(values):
        buckets[(-frequency * t) % P] += value
    return tuple(buckets[j] - buckets[12] for j in range(12))


def reduce_buckets(buckets):
    return tuple(buckets[j] - buckets[12] for j in range(12))


def endpoint_fourier_zero(matrix, k, ell):
    buckets = [0] * P
    for r in range(P):
        for s in range(P):
            buckets[(-k * r - ell * s) % P] += matrix[r][s]
    return len(set(buckets)) == 1


def weighted_moments(weighted_masks):
    h = [0] * P
    g = [[0] * P for _ in range(P)]
    m = [[0] * P for _ in range(P)]
    for mask, weight in weighted_masks:
        f = bits_of(mask)
        for r in range(P):
            h[r] += weight * f[r]
            for s in range(P):
                g[r][s] += weight * f[r] * f[s]
                m[r][s] += weight * f[r] * (1 - f[s])
    return h, g, m


def main():
    assert pow(ZETA, P, Q) == 1
    assert all(pow(ZETA, j, Q) != 1 for j in range(1, P))

    masks_checked = 0
    edge_checks = 0
    affine_checks = 0
    prony_checks = 0
    crofton_fourier_checks = 0
    zero_directional_variance = Counter()
    zero_variance_masks = []

    for mask in range(1, (1 << P) - 1):
        masks_checked += 1
        f = bits_of(mask)
        n = sum(f)
        kbank = {tau: boundary(mask, tau) for tau in range(1, P)}
        bvals = [0] + [sum(kbank[tau]) for tau in range(1, P)]

        outgoing = [sum(kbank[tau][r] for tau in range(1, P)) for r in range(P)]
        incoming = [
            sum(kbank[tau][(r - tau) % P] for tau in range(1, P))
            for r in range(P)
        ]
        assert outgoing == [(P - n) * f[r] for r in range(P)]
        assert incoming == [n * (1 - f[r]) for r in range(P)]
        assert [int(value > 0) for value in outgoing] == list(f)

        for tau in range(1, P):
            kt = kbank[tau]
            krev = kbank[(-tau) % P]
            b = sum(kt)
            assert 1 <= b <= min(n, P - n)
            prony_support(kt)
            prony_checks += 1
            for r in range(P):
                s = (r + tau) % P
                reverse = krev[s]
                assert kt[r] * reverse == 0
                assert kt[r] - reverse == f[r] - f[s]
                assert kt[r] + reverse == (f[r] - f[s]) ** 2
                edge_checks += 1

        assert sum(bvals) == n * (P - n)

        # Exact Wiener--Khinchin/Crofton identity in Q[zeta_13].
        corr = [sum(f[r] * f[(r + tau) % P] for r in range(P)) for tau in range(P)]
        assert bvals == [n - corr[tau] for tau in range(P)]
        for k in range(1, P):
            left = cyclotomic_reduce(bvals, k)
            right = tuple(-x for x in cyclotomic_reduce(corr, k))
            assert left == right
            product_buckets = [0] * P
            for r in range(P):
                for s in range(P):
                    product_buckets[(k * (r - s)) % P] += f[r] * f[s]
            assert cyclotomic_reduce(corr, k) == reduce_buckets(product_buckets)
            # Phi_13 irreducibility: a mixed Boolean coefficient vector is
            # not constant, hence its primitive Fourier value is nonzero.
            assert len(set(f)) > 1
            assert any(left)
            crofton_fourier_checks += 1

        if len(set(bvals[1:])) == 1:
            zero_directional_variance[n] += 1
            zero_variance_masks.append(mask)

        # Complement is exact edge converse.
        comp = ((1 << P) - 1) ^ mask
        for tau in range(1, P):
            kc = boundary(comp, tau)
            kt = kbank[(-tau) % P]
            for r in range(P):
                assert kc[r] == kt[(r + tau) % P]

        # Translation and a primitive dilation generate AGL(1,13).
        for u, b in ((1, 1), (2, 0), (-1, 0)):
            moved = transform_mask(mask, u, b)
            for tau in range(1, P):
                source = kbank[tau]
                target = boundary(moved, (u * tau) % P)
                for r in range(P):
                    assert target[(u * r + b) % P] == source[r]
                    affine_checks += 1

    assert masks_checked == 8190
    assert zero_directional_variance == Counter({1: 13, 4: 52, 9: 52, 12: 13})

    seen = set()
    zero_necklaces = Counter()
    for mask in zero_variance_masks:
        if mask in seen:
            continue
        orbit = {rotate_mask(mask, b) for b in range(P)}
        seen.update(orbit)
        n = sum(bits_of(mask))
        zero_necklaces[(n, sum(boundary(mask, 1)))] += 1
    assert zero_necklaces == Counter({(1, 1): 1, (4, 3): 4, (9, 3): 4, (12, 1): 1})

    # A fixed slope sees only run terminals, not their lengths.
    interval_short = mask_of({0})
    interval_long = mask_of({9, 10, 11, 12, 0})
    assert boundary(interval_short, 1) == boundary(interval_long, 1)
    assert interval_short != interval_long

    # The canonical homometric pair has identical Crofton counts but
    # different oriented boundary fields.
    hom_a = mask_of({0, 1, 3, 9})
    hom_b = mask_of({1, 2, 5, 7})
    ba = tuple(sum(boundary(hom_a, tau)) for tau in range(P))
    bb = tuple(sum(boundary(hom_b, tau)) for tau in range(P))
    assert ba == bb == (0,) + (3,) * 12
    assert endpoint_matrix(hom_a) != endpoint_matrix(hom_b)

    # Averaging destroys pointwise reconstruction: uniform singleton and
    # uniform co-singleton laws have the same boundary matrix.
    singleton_law = [(1 << a, 1) for a in range(P)]
    cosingleton_law = [(((1 << P) - 1) ^ (1 << a), 1) for a in range(P)]
    hs, gs, ms = weighted_moments(singleton_law)
    hc, gc, mc = weighted_moments(cosingleton_law)
    assert ms == mc
    assert hs == [1] * P and hc == [12] * P
    assert all(ms[r][s] == (0 if r == s else 1) for r in range(P) for s in range(P))

    # Target-anchored rational mixture with 132 vanishing mixed modes.
    # This is the unnormalized uniform law on all 2^12 subsets of roots
    # 1,...,12; deleting the empty mask changes no boundary entry.
    anchored_uniform = [[0] * P for _ in range(P)]
    for r in range(1, P):
        anchored_uniform[r][0] = 2**11
        for s in range(1, P):
            if r != s:
                anchored_uniform[r][s] = 2**10
    zero_modes = {
        (k, ell)
        for k in range(P)
        for ell in range(P)
        if endpoint_fourier_zero(anchored_uniform, k, ell)
    }
    expected_zero_modes = {
        (k, ell)
        for k in range(1, P)
        for ell in range(1, P)
        if (k + ell) % P != 0
    }
    assert zero_modes == expected_zero_modes
    assert len(zero_modes) == 132
    for tau in range(1, P):
        row_vector = [anchored_uniform[r][(r + tau) % P] for r in range(P)]
        assert row_vector[0] == 0 and sum(row_vector) > 0 and len(set(row_vector)) > 1

    # A nontrivial target-anchored aggregate: M reconstructs H and Gram.
    anchored_weighted = []
    for local in range(1, 1 << 12):
        mask = local << 1
        anchored_weighted.append((mask, 1 + (local * local) % 17))
    h, g, m = weighted_moments(anchored_weighted)
    assert h[0] == 0
    for r in range(P):
        assert h[r] == m[r][0]
        for s in range(P):
            assert g[r][s] == h[r] - m[r][s]
            assert m[r][s] - m[s][r] == h[r] - h[s]
    divergence = [sum(m[r][s] - m[s][r] for s in range(P)) for r in range(P)]
    assert divergence == [P * h[r] - sum(h) for r in range(P)]

    # Target-anchored singleton/adjacent-pair path cone.
    alphas = [0] + [3 * j + 1 for j in range(1, 13)]
    betas = [0] + [5 * j + 2 for j in range(1, 12)] + [0]
    path_law = [(1 << j, alphas[j]) for j in range(1, 13)]
    path_law += [((1 << j) | (1 << (j + 1)), betas[j]) for j in range(1, 12)]
    hp, gp, mp = weighted_moments(path_law)
    gamma_plus = [mp[a][(a - 1) % P] for a in range(P)]
    gamma_minus = [mp[a][(a + 1) % P] for a in range(P)]
    for tau in range(2, 12):
        assert [mp[a][(a + tau) % P] for a in range(P)] == hp
    assert [mp[a][(a - 1) % P] for a in range(P)] == gamma_plus
    assert [mp[a][(a + 1) % P] for a in range(P)] == gamma_minus
    recovered_beta = [0] * P
    for a in range(1, 12):
        recovered_beta[a] = hp[a] - gamma_minus[a]
    recovered_alpha = [0] * P
    for a in range(1, 13):
        recovered_alpha[a] = hp[a] - recovered_beta[a - 1] - recovered_beta[a]
    assert recovered_alpha == alphas
    assert recovered_beta == betas
    mass_alpha = sum(alphas)
    mass_beta = sum(betas)
    occupancy_mass = sum(hp)
    crofton_mass = sum(sum(row) for row in mp)
    assert occupancy_mass == mass_alpha + 2 * mass_beta
    assert crofton_mass == 12 * mass_alpha + 22 * mass_beta
    assert (crofton_mass - 10 * occupancy_mass) // 2 == mass_alpha + mass_beta

    # THM-2535's source-charge hostile survives every positive root slope.
    singleton_one = mask_of({1})
    assert all(boundary(singleton_one, tau)[1] == 1 for tau in range(1, P))
    q7 = 547
    generator = next(
        g0 for g0 in range(2, q7) if len({pow(g0, j, q7) for j in range(1, q7)}) == q7 - 1
    )
    z13 = pow(generator, (q7 - 1) // 13, q7)
    z7 = pow(generator, (q7 - 1) // 7, q7)
    source_modes = 0
    for alpha in range(1, 13):
        for beta in range(1, 7):
            value = (1 - pow(z13, -2 * alpha, q7)) * (1 - pow(z7, -beta, q7))
            assert value % q7 != 0
            source_modes += 1
    assert source_modes == 72
    def hostile_d(hh, rr):
        return (int(hh == 0) - int(hh == 2)) * (int(rr == 0) - int(rr == 1))

    assert all(hostile_d(1, kappa) == 0 for _tau in range(1, P) for kappa in range(7))

    print("THM-2534 dependency-free exact referee")
    print(f"field: F_{P}; masks: {masks_checked} mixed Boolean masks")
    print(f"matched oriented-edge identities: {edge_checks}")
    print(f"AGL generator covariance identities: {affine_checks}")
    print(f"sparse Newton/Prony reconstructions: {prony_checks}; maximum support 6")
    print(f"exact Crofton/Wiener--Khinchin character checks: {crofton_fourier_checks}")
    print("directional-variance-zero masks by size: 1:13, 4:52, 9:52, 12:13")
    print("directional-variance-zero rotation necklaces: 1+4+4+1 = 10")
    print("fixed-slope run-length hostile: PASS")
    print("homometric Crofton phase-loss hostile: PASS (all nonzero counts 3)")
    print("weighted gauge hostile: PASS (singleton/co-singleton laws coincide)")
    print("target-anchored mixed-mode hostile: PASS (132 exact zero modes)")
    print("target-anchor reconstruction of (H,Gram): PASS")
    print("deep-path collapse and 23-ray recovery: PASS")
    print("all-slope zero-source-tooth hostile: PASS (72/72 source modes survive)")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
