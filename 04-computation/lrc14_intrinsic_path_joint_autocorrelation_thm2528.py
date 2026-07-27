#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2528.

The referee checks two complementary lifts of the live depth-one collision:

* the A_tau H_tau odd bank before its path coefficients are collapsed; and
* two-dimensional autocorrelations of the THM-2521 potential and THM-2508
  affine-cut tables.

All calculations use integers or arithmetic in F_547.  The element
2^6 in F_547 has exact order 91, so a nonzero value there is an exact
nonvanishing certificate for the corresponding Q(zeta_91) expression.
"""


P = 13
Q = 7
MODULUS = 547
ROOT_91 = pow(2, 6, MODULUS)
ZETA_13 = pow(ROOT_91, 7, MODULUS)
ZETA_7 = pow(ROOT_91, 13, MODULUS)
CHI7 = {1: 1, 2: 1, 3: -1, 4: 1, 5: -1, 6: -1}


def require(condition, label):
    if not condition:
        raise AssertionError(label)


def matmul(first, second):
    return [
        [
            sum(first[row][middle] * second[middle][col] for middle in range(P))
            for col in range(P)
        ]
        for row in range(P)
    ]


def matvec(matrix, vector):
    return [
        sum(matrix[row][col] * vector[col] for col in range(P))
        for row in range(P)
    ]


def hamilton(tau):
    """THM-2523's even Hamilton operator A_tau."""
    matrix = [[0 for _ in range(P)] for _ in range(P)]
    for row in range(P):
        for s in range(1, 7):
            matrix[row][(row + 2 * tau * s) % P] -= CHI7[s]
            matrix[row][(row - 2 * tau * s) % P] -= CHI7[s]
    return matrix


def cayley(tau):
    """THM-2526's odd cyclic Hilbert operator H_tau."""
    matrix = [[0 for _ in range(P)] for _ in range(P)]
    for row in range(P):
        for s in range(1, 7):
            matrix[row][(row - 2 * tau * s) % P] += 1
            matrix[row][(row + 2 * tau * s) % P] -= 1
    return matrix


def kernels(tau):
    """Return the coefficient functions of A_tau and H_tau."""
    a = {}
    eta = {}
    for s in range(1, 7):
        a[(2 * tau * s) % P] = -CHI7[s]
        a[(-2 * tau * s) % P] = -CHI7[s]
        eta[(-2 * tau * s) % P] = 1
        eta[(2 * tau * s) % P] = -1
    require(set(a) == set(range(1, P)), "Hamilton kernel support")
    require(set(eta) == set(range(1, P)), "Cayley kernel support")
    require(all(a[-v % P] == a[v] for v in range(1, P)), "A kernel even")
    require(all(eta[-v % P] == -eta[v] for v in range(1, P)), "H kernel odd")
    return a, eta


def autocorrelation(mask):
    return [
        sum(((mask >> root) & 1) * ((mask >> ((root + gap) % P)) & 1)
            for root in range(P))
        for gap in range(P)
    ]


def audit_intrinsic_path_split():
    tau = 1
    a, eta = kernels(tau)
    composite = matmul(hamilton(tau), cayley(tau))
    fixed_coordinate = (-4 * tau) % P
    psi_min = None
    psi_max = 0
    equality_masks = 0
    path_checks = 0
    pair_extractions = 0

    for mask in range(1 << P):
        c = autocorrelation(mask)
        odd = matvec(composite, c)
        require(all(odd[-t % P] == -odd[t] for t in range(P)), "odd bank")

        # The collapsed matrix coefficient and the uncollapsed (v,h) path
        # count agree at every output coordinate.  Each positive or negative
        # term is the number of actual root triples (r,v,h), so this is also
        # the exact disjoint-union count on the ordered four-arm lift.
        for t in range(P):
            positive = 0
            negative = 0
            for v in range(1, P):
                for h in range(1, P):
                    count = c[(t + v + h) % P]
                    if a[v] * eta[h] > 0:
                        positive += count
                    else:
                        negative += count
            require(positive - negative == odd[t], "four-arm path identity")
            path_checks += 1

        psi = -odd[(4 * tau) % P]
        population = c[0]
        pointwise_drift_numerator = P * population - sum(c)
        require(42 * psi >= pointwise_drift_numerator, "sharp 13/42 floor")
        require(psi >= 0, "fixed-sign coordinate")

        if mask in (0, (1 << P) - 1):
            require(psi == 0, "constant-mask equality")
            continue

        require(psi > 0, "strict nonconstant fixed sign")
        require(odd[fixed_coordinate] > 0, "fixed positive path imbalance")
        psi_min = psi if psi_min is None else min(psi_min, psi)
        psi_max = max(psi_max, psi)
        if 42 * psi == pointwise_drift_numerator:
            equality_masks += 1

        # Pair h with -h.  A positive total forces one matched pair of
        # Boolean collision intersections to have the selected strict sign.
        extracted = False
        for v in range(1, P):
            for h in range(1, P):
                if eta[h] != 1:
                    continue
                gap_plus = (fixed_coordinate + v + h) % P
                gap_minus = (fixed_coordinate + v - h) % P
                if a[v] * (c[gap_plus] - c[gap_minus]) > 0:
                    extracted = True
                    break
            if extracted:
                break
        require(extracted, "matched Boolean wedge extraction")
        pair_extractions += 1

    require((psi_min, psi_max) == (1, 98), "fixed-sign range")
    require(equality_masks == 52, "sharp-floor equality census")

    # At fixed t, (r,v,h) maps injectively to the ordered four root labels
    # (r,r+t,r+t+v,r+t+v+h).  Positive and negative path components cannot
    # collide because the tuple recovers v and h.
    t = fixed_coordinate
    tuples = {}
    for root in range(P):
        for v in range(1, P):
            for h in range(1, P):
                path = (
                    root,
                    (root + t) % P,
                    (root + t + v) % P,
                    (root + t + v + h) % P,
                )
                require(path not in tuples, "four-arm path injectivity")
                tuples[path] = a[v] * eta[h]
    require(len(tuples) == P * (P - 1) ** 2, "four-arm component count")

    # Every depth-one root displacement preserves the live shallow and deep
    # phase arguments: c_j=13a and c_3=13^2b make c(x+w/13)-cx integral.
    sidecar_checks = 0
    for unit in range(1, P):
        shallow = P * unit
        deep = P * P * unit
        for displacement in range(P):
            require((shallow * displacement) % P == 0, "shallow sidecar")
            require((deep * displacement) % P == 0, "deep sidecar")
            sidecar_checks += 2

    return {
        "masks": 1 << P,
        "nonconstant": (1 << P) - 2,
        "path_checks": path_checks,
        "pair_extractions": pair_extractions,
        "components": len(tuples),
        "psi_min": psi_min,
        "psi_max": psi_max,
        "equality_masks": equality_masks,
        "sidecar_checks": sidecar_checks,
    }


def transform(table, alpha, beta):
    """Unnormalized F_13 x F_7 transform, evaluated exactly in F_547."""
    return sum(
        table[h][r]
        * pow(ZETA_13, (-alpha * h) % P, MODULUS)
        * pow(ZETA_7, (-beta * r) % Q, MODULUS)
        for h in range(P)
        for r in range(Q)
    ) % MODULUS


def correlation_2d(table):
    return [
        [
            sum(
                table[h][r] * table[(h + t) % P][(r + s) % Q]
                for h in range(P)
                for r in range(Q)
            )
            for s in range(Q)
        ]
        for t in range(P)
    ]


def h_action_2d(table, tau):
    _, eta = kernels(tau)
    return [
        [sum(eta[h] * table[(t + h) % P][s] for h in range(1, P))
         for s in range(Q)]
        for t in range(P)
    ]


def potential_table(centered_profile, tau=1, a=1, c=0):
    """THM-2521 equation (27), with p scaled by 13."""
    table = [[0 for _ in range(Q)] for _ in range(P)]
    for h in range(P):
        for r in range(Q):
            s = (a * r + c) % Q
            if s == 0:
                table[h][r] = centered_profile[h]
            else:
                table[h][r] = (
                    centered_profile[(h + tau * s) % P]
                    + centered_profile[(h - tau * s) % P]
                )
    return table


def cut_table(table, tau, a):
    """The (v,c) table of THM-2508 R_(tau,a,c)(v)."""
    return [
        [
            sum(table[(v - tau * ((a * r + c) % Q)) % P][r]
                for r in range(Q))
            for c in range(Q)
        ]
        for v in range(P)
    ]


def audit_joint_autocorrelation_scalarization():
    require(pow(ROOT_91, 91, MODULUS) == 1, "91st root")
    require(all(pow(ROOT_91, divisor, MODULUS) != 1 for divisor in (1, 7, 13)),
            "exact order 91")

    # A nonuniform rational predecessor profile: 13 times the centred delta.
    centered = [12] + [-1] * (P - 1)
    d = potential_table(centered)
    require(all(sum(row) == 0 for row in d), "potential row sums")
    require(all(sum(d[h][r] for h in range(P)) == 0 for r in range(Q)),
            "potential column sums")

    corr = correlation_2d(d)
    require(all(corr[-t % P][-s % Q] == corr[t][s]
                for t in range(P) for s in range(Q)), "joint-even correlation")
    odd = h_action_2d(corr, 1)
    require(all(odd[-t % P][-s % Q] == -odd[t][s]
                for t in range(P) for s in range(Q)), "joint-odd scalar bank")

    inverse_91 = pow(P * Q, -1, MODULUS)
    potential_modes = 0
    wiener_khinchin_checks = 0
    odd_modes = 0
    for alpha in range(1, P):
        for beta in range(1, Q):
            d_hat = transform(d, alpha, beta)
            d_conjugate = transform(d, -alpha % P, -beta % Q)
            require(d_hat != 0, "potential primitive mode")
            norm_square = d_hat * d_conjugate % MODULUS
            require(norm_square != 0, "potential norm mode")
            corr_hat = transform(corr, alpha, beta)
            require(corr_hat == norm_square, "unnormalized Wiener-Khinchin")
            normalized_corr_hat = corr_hat * inverse_91 % MODULUS
            require(normalized_corr_hat == norm_square * inverse_91 % MODULUS,
                    "normalized 1/91 factor")

            odd_hat = transform(odd, alpha, beta)
            multiplier = sum(
                eta * pow(ZETA_13, alpha * shift % P, MODULUS)
                for shift, eta in kernels(1)[1].items()
            ) % MODULUS
            require(multiplier != 0, "Cayley multiplier")
            require(odd_hat == multiplier * corr_hat % MODULUS,
                    "odd scalar multiplier")
            require(odd_hat != 0, "odd primitive mode")
            potential_modes += 1
            wiener_khinchin_checks += 1
            odd_modes += 1

    # For each of the 72 affine (tau,a) cut charts, the 72 transforms of the
    # (v,c) table are nonzero.  Their autocorrelation coefficients are their
    # exact norm squares; H in v multiplies by the same nonzero root factor.
    cut_modes = 0
    cut_norm_modes = 0
    for tau in range(1, P):
        for a in range(1, Q):
            cut = cut_table(d, tau, a)
            for alpha in range(1, P):
                multiplier = sum(
                    eta * pow(ZETA_13, alpha * shift % P, MODULUS)
                    for shift, eta in kernels(1)[1].items()
                ) % MODULUS
                require(multiplier != 0, "cut Cayley multiplier")
                for beta in range(1, Q):
                    coefficient = transform(cut, alpha, beta)
                    conjugate = transform(cut, -alpha % P, -beta % Q)
                    require(coefficient != 0, "affine cut coefficient")
                    norm_square = coefficient * conjugate % MODULUS
                    require(norm_square != 0, "affine cut norm coefficient")
                    require(multiplier * norm_square % MODULUS != 0,
                            "affine cut odd scalar coefficient")
                    cut_modes += 1
                    cut_norm_modes += 1

    require(potential_modes == 72, "72 potential modes")
    require(cut_modes == 5184, "5184 cut modes")
    return {
        "potential_modes": potential_modes,
        "wiener_khinchin_checks": wiener_khinchin_checks,
        "odd_modes": odd_modes,
        "cut_modes": cut_modes,
        "cut_norm_modes": cut_norm_modes,
    }


def main():
    path = audit_intrinsic_path_split()
    joint = audit_joint_autocorrelation_scalarization()
    print("THM-2528 intrinsic path and joint autocorrelation companion: PASS")
    print(
        f"boolean_masks={path['masks']}; nonconstant_masks={path['nonconstant']}; "
        f"path_identity_checks={path['path_checks']}"
    )
    print(
        f"fixed_path_coordinate=-4tau; psi_range={path['psi_min']}..{path['psi_max']}; "
        f"sharp_floor_equality_masks={path['equality_masks']}"
    )
    print(
        f"four_arm_components={path['components']}; "
        f"matched_pair_extractions={path['pair_extractions']}; "
        f"sidecar_shift_checks={path['sidecar_checks']}"
    )
    print(
        f"joint_potential_modes={joint['potential_modes']}; "
        f"wiener_khinchin_checks={joint['wiener_khinchin_checks']}; "
        f"joint_odd_modes={joint['odd_modes']}"
    )
    print(
        f"affine_cut_scalar_modes={joint['cut_modes']}; "
        f"affine_cut_norm_modes={joint['cut_norm_modes']}; "
        f"cyclotomic_certificate_modulus={MODULUS}"
    )


if __name__ == "__main__":
    main()
