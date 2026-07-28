#!/usr/bin/env python3
"""Exact companion for THM-2820.

The script checks three distinct structures without identifying them:

1. pointwise Boolean function algebras, whose idempotent tangent spaces
   and relative derivations vanish;
2. modular convolution group rings F_p[C_p], where the orbit norm is the
   top augmentation power (g-1)^(p-1); and
3. THM-2806's integral central allocation profile at p=13.

It uses explicit exception gates, no assertions, no floating point, and no
external packages.
"""

from itertools import product
from math import comb


PRIMES = (2, 3, 5, 7, 11, 13)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def add_poly(left, right, p):
    return tuple((left[i] + right[i]) % p for i in range(p))


def mul_poly(left, right, p):
    out = [0] * p
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            if i + j < p:
                out[i + j] = (out[i + j] + a * b) % p
    return tuple(out)


def pow_poly(base, exponent, p):
    out = (1,) + (0,) * (p - 1)
    power = base
    n = exponent
    while n:
        if n & 1:
            out = mul_poly(out, power, p)
        power = mul_poly(power, power, p)
        n //= 2
    return out


def compose_poly(outer, inner, p):
    out = (0,) * p
    power = (1,) + (0,) * (p - 1)
    for coefficient in outer:
        if coefficient:
            out = add_poly(
                out,
                tuple(coefficient * value % p for value in power),
                p,
            )
        power = mul_poly(power, inner, p)
    return out


def eta_basis_group_ring_checks(p):
    one = (1,) + (0,) * (p - 1)
    eta = (0, 1) + (0,) * (p - 2)
    generator = add_poly(one, eta, p)
    norm = (0,) * p
    for exponent in range(p):
        norm = add_poly(norm, pow_poly(generator, exponent, p), p)
    top = (0,) * (p - 1) + (1,)
    require(norm == top, f"norm-top identity failed at p={p}")
    require(pow_poly(eta, p, p) == (0,) * p, f"eta^p failed at p={p}")

    fixed_norms = 0
    cotangent_multipliers = []
    for unit in range(1, p):
        image_eta = add_poly(pow_poly(generator, unit, p), tuple(-x % p for x in one), p)
        require(image_eta[0] == 0, f"augmentation failed at p={p}, u={unit}")
        require(
            image_eta[1] == unit,
            f"cotangent multiplier failed at p={p}, u={unit}",
        )
        image_norm = compose_poly(norm, image_eta, p)
        require(image_norm == norm, f"norm invariance failed at p={p}, u={unit}")
        fixed_norms += 1
        cotangent_multipliers.append(image_eta[1])

    fixed_cotangent = [
        value
        for value in range(p)
        if all(unit * value % p == value for unit in range(1, p))
    ]
    if p > 2:
        require(
            fixed_cotangent == [0],
            f"unexpected nonzero canonical jet at p={p}",
        )

    return fixed_norms, tuple(cotangent_multipliers), tuple(fixed_cotangent)


def pointwise_idempotent_tangent_checks(p, size=4):
    # Idempotents in F_p^size are exactly the Boolean vectors.  In a
    # commutative square-zero extension, e+epsilon*u is idempotent iff
    # (2e-1)u=0 coordinatewise.  Every multiplier 2e-1 is +/-1.
    idempotents = tuple(product((0, 1), repeat=size))
    zero_kernel_count = 0
    for idempotent in idempotents:
        tangent_count = 0
        for tangent in product(range(p), repeat=size):
            equation = tuple(
                ((2 * e - 1) * u) % p
                for e, u in zip(idempotent, tangent)
            )
            if equation == (0,) * size:
                tangent_count += 1
        require(
            tangent_count == 1,
            f"idempotent tangent space nonzero at p={p}, e={idempotent}",
        )
        zero_kernel_count += 1
    return len(idempotents), zero_kernel_count


def mul_bivariate(left, right, p):
    out = {}
    for (i, j), a in left.items():
        for (k, ell), b in right.items():
            if i + k < p and j + ell < p:
                key = (i + k, j + ell)
                out[key] = (out.get(key, 0) + a * b) % p
    return {key: value for key, value in out.items() if value}


def add_bivariate(left, right, p):
    out = dict(left)
    for key, value in right.items():
        out[key] = (out.get(key, 0) + value) % p
        if out[key] == 0:
            del out[key]
    return out


def scale_bivariate(poly, scalar, p):
    return {
        key: scalar * value % p
        for key, value in poly.items()
        if scalar * value % p
    }


def translate_bivariate(poly, a, b, p):
    translate = {
        (i, j): comb(a, i) * comb(b, j) % p
        for i in range(a + 1)
        for j in range(b + 1)
    }
    return mul_bivariate(poly, translate, p)


def homogeneous_piece(poly, degree):
    return {
        key: value
        for key, value in poly.items()
        if sum(key) == degree and value
    }


def linear_times_piece(piece, a, b, p):
    linear = {}
    if a:
        linear[(1, 0)] = a
    if b:
        linear[(0, 1)] = b
    return mul_bivariate(linear, piece, p)


def vector_rank_mod_p(vectors, keys, p):
    matrix = [[vector.get(key, 0) % p for key in keys] for vector in vectors]
    rank = 0
    column = 0
    while rank < len(matrix) and column < len(keys):
        pivot = next(
            (row for row in range(rank, len(matrix)) if matrix[row][column]),
            None,
        )
        if pivot is None:
            column += 1
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], -1, p)
        matrix[rank] = [inverse * value % p for value in matrix[rank]]
        for row in range(len(matrix)):
            if row == rank or matrix[row][column] == 0:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                (left - factor * right) % p
                for left, right in zip(matrix[row], matrix[rank])
            ]
        rank += 1
        column += 1
    return rank


def first_normal_translation_checks(polynomials, p):
    ranks = {}
    identity_checks = 0
    for name, poly in polynomials.items():
        degree = min(sum(key) for key in poly)
        leading = homogeneous_piece(poly, degree)
        basis_images = []
        for a, b in ((1, 0), (0, 1)):
            difference = add_bivariate(
                translate_bivariate(poly, a, b, p),
                scale_bivariate(poly, -1, p),
                p,
            )
            basis_images.append(homogeneous_piece(difference, degree + 1))
        keys = tuple(
            (i, degree + 1 - i)
            for i in range(p)
            if 0 <= degree + 1 - i < p
        )
        ranks[name] = vector_rank_mod_p(basis_images, keys, p)

        for a in range(p):
            for b in range(p):
                difference = add_bivariate(
                    translate_bivariate(poly, a, b, p),
                    scale_bivariate(poly, -1, p),
                    p,
                )
                actual = homogeneous_piece(difference, degree + 1)
                expected = linear_times_piece(leading, a, b, p)
                require(
                    actual == expected,
                    f"first-normal translation identity drift for {name}, {(a, b)}",
                )
                identity_checks += 1
    require(
        ranks == {"B": 0, "P": 1, "Q": 1, "H": 2, "Omega": 2},
        f"first-normal rank classification drift: {ranks}",
    )
    return tuple(ranks.items()), identity_checks


def carrier_gauge_checks(p):
    # THM-2763 has (ell,a,b)~(ell+sW,a+s,b-s).  After choosing
    # z.W=1, only z.ell is needed to verify the invariant
    # kappa_z=(a-z.ell,b+z.ell).
    invariant_checks = 0
    for ell_dot in range(p):
        for a in range(p):
            for b in range(p):
                kappa = ((a - ell_dot) % p, (b + ell_dot) % p)
                for s in range(p):
                    gauged = (
                        (a + s - (ell_dot + s)) % p,
                        (b - s + ell_dot + s) % p,
                    )
                    require(gauged == kappa, "carrier gauge invariant drift")
                    invariant_checks += 1

    effective_nonzero = 0
    family_checks = 0
    for l_dot in range(p):
        for q1 in range(p):
            for q2 in range(p):
                effective = ((q1 - l_dot) % p, (q2 + l_dot) % p)
                if effective != (0, 0):
                    effective_nonzero += 1
                for t in range(p):
                    observed = (
                        (t * q1 - t * l_dot) % p,
                        (t * q2 + t * l_dot) % p,
                    )
                    expected = (
                        t * effective[0] % p,
                        t * effective[1] % p,
                    )
                    require(observed == expected, "effective transverse vector drift")
                    family_checks += 1

    pure_gauge = ((1 - 1) % p, ((-1) + 1) % p)
    require(pure_gauge == (0, 0), "marked pure gauge unexpectedly visible")
    forgotten_visible = sum(
        (q1 + q2) % p != 0
        for q1 in range(p)
        for q2 in range(p)
    )
    require(forgotten_visible == p * (p - 1), "forgotten-address test drift")
    return (
        invariant_checks,
        family_checks,
        effective_nonzero,
        p**3,
        forgotten_visible,
        p**2,
    )


def integral_inverse_push_checks(p):
    # In Z[C_p^2], Omega=(N-1)x(N-1) and the positive lift of its
    # characteristic-p inverse is Theta=(N+1)x(N+1).  Convolve after a
    # translate, then push along lambda(i,j)=i+j.
    omega = {
        (i, j): int(i != 0 and j != 0)
        for i in range(p)
        for j in range(p)
    }
    theta = {
        (i, j): (1 + int(i == 0)) * (1 + int(j == 0))
        for i in range(p)
        for j in range(p)
    }
    baseline = p * (p * p - 2)
    checks = 0
    for a in range(p):
        for b in range(p):
            translated = {
                ((i + a) % p, (j + b) % p): value
                for (i, j), value in omega.items()
            }
            convolution = {(i, j): 0 for i in range(p) for j in range(p)}
            for (i, j), left in theta.items():
                for (k, ell), right in translated.items():
                    convolution[((i + k) % p, (j + ell) % p)] += left * right
            pushed = [0] * p
            for (i, j), value in convolution.items():
                pushed[(i + j) % p] += value
            expected = [baseline] * p
            expected[(a + b) % p] += 1
            require(
                pushed == expected,
                f"integral inverse/push profile drift at {(a, b)}",
            )
            checks += p
    require(set(theta.values()) == {1, 2, 4}, "positive inverse weights drift")
    return checks, baseline, tuple(sorted(set(theta.values())))


def allocation_square_p13():
    p = 13
    # Under the coefficient-vector identification with F_p[C_p x C_p],
    # constant functions become N and delta_0 becomes 1.
    allocation_terms = {
        "B": ((p - 1, p - 1),),
        "P": ((0, p - 1),),
        "Q": ((p - 1, 0),),
        "H": ((0, 0),),
    }
    augmentation_degrees = tuple(
        sum(allocation_terms[key][0]) for key in ("B", "P", "Q", "H")
    )
    require(
        augmentation_degrees == (24, 12, 12, 0),
        "tensor norm filtration drift",
    )

    omega = {
        (p - 1, p - 1): 1,
        (0, p - 1): -1 % p,
        (p - 1, 0): -1 % p,
        (0, 0): 1,
    }
    require(
        set(omega) == {(12, 12), (0, 12), (12, 0), (0, 0)},
        "allocation Mobius support drift",
    )
    require(
        not any(i + j == 1 for i, j in omega),
        "top-norm allocation unexpectedly acquired a first jet",
    )
    omega_inverse = {
        (0, 0): 1,
        (p - 1, 0): 1,
        (0, p - 1): 1,
        (p - 1, p - 1): 1,
    }
    require(
        mul_bivariate(omega, omega_inverse, p) == {(0, 0): 1},
        "raw face unit inverse drift",
    )
    first_normal_ranks, first_normal_checks = first_normal_translation_checks(
        {
            "B": {(p - 1, p - 1): 1},
            "P": {(0, p - 1): 1},
            "Q": {(p - 1, 0): 1},
            "H": {(0, 0): 1},
            "Omega": omega,
        },
        p,
    )

    # Although N alone has zero augmentation and no normalized first
    # coordinate, the rooted Mobius face has J_00=1.  THM-2201's triangular
    # translation law therefore gives an affine barycenter.  Check the full
    # two-axis table directly rather than importing the formula.
    barycenters = {}
    for a in range(p):
        for b in range(p):
            shifted = translate_bivariate(omega, a, b, p)
            j00 = shifted.get((0, 0), 0)
            j10 = shifted.get((1, 0), 0)
            j01 = shifted.get((0, 1), 0)
            require(j00 == 1, "translated raw face lost augmentation")
            barycenter = (j10 * pow(j00, -1, p) % p, j01 * pow(j00, -1, p) % p)
            require(barycenter == (a, b), "normalized Hasse barycenter drift")
            barycenters[(a, b)] = barycenter
    require(len(set(barycenters.values())) == p * p, "barycenter lost translations")

    gauge_difference_checks = 0
    gauge_origin = barycenters[(1, p - 1)]
    for a in range(p):
        for b in range(p):
            gauged = ((a + 1) % p, (b - 1) % p)
            difference = (
                (barycenters[(a, b)][0] - barycenters[(0, 0)][0]) % p,
                (barycenters[(a, b)][1] - barycenters[(0, 0)][1]) % p,
            )
            gauged_difference = (
                (barycenters[gauged][0] - gauge_origin[0]) % p,
                (barycenters[gauged][1] - gauge_origin[1]) % p,
            )
            require(
                difference == gauged_difference,
                "relative barycenter failed marked-gauge invariance",
            )
            gauge_difference_checks += 1

    central = (p * p, p, p, 1)
    hadamard = (
        sum(central),
        central[0] + central[1] - central[2] - central[3],
        central[0] - central[1] + central[2] - central[3],
        central[0] - central[1] - central[2] + central[3],
    )
    require(central == (169, 13, 13, 1), "central profile drift")
    require(hadamard == (196, 168, 168, 144), "Hadamard profile drift")
    require(hadamard[3] % p == 1, "D3 residue drift")

    valuations = []
    for value in central:
        exponent = 0
        while value % p == 0:
            value //= p
            exponent += 1
        valuations.append(exponent)
    require(tuple(valuations) == (2, 1, 1, 0), "p-adic filtration drift")

    # The tempting interpolation e+epsilon*(1-e) on the 13-point delta
    # has a nonzero idempotence defect -epsilon*(1-e) on twelve atoms.
    idempotent = (1,) + (0,) * (p - 1)
    tangent = tuple(1 - value for value in idempotent)
    first_order_defect = tuple(
        ((2 * e - 1) * u) % p for e, u in zip(idempotent, tangent)
    )
    require(
        sum(value != 0 for value in first_order_defect) == 12,
        "non-idempotent interpolation hostile drift",
    )

    return (
        augmentation_degrees,
        central,
        tuple(valuations),
        hadamard,
        sum(value != 0 for value in first_order_defect),
        len(barycenters),
        gauge_difference_checks,
        first_normal_ranks,
        first_normal_checks,
        1,
    )


def main():
    tangent_rows = []
    norm_rows = []
    for p in PRIMES:
        idempotents, zero_kernels = pointwise_idempotent_tangent_checks(p)
        fixed_norms, multipliers, fixed_cotangent = eta_basis_group_ring_checks(p)
        tangent_rows.append((p, idempotents, zero_kernels))
        norm_rows.append((p, fixed_norms, multipliers, fixed_cotangent))

    (
        augmentation_degrees,
        central,
        valuations,
        hadamard,
        interpolation_defect,
        barycenter_count,
        gauge_difference_checks,
        first_normal_ranks,
        first_normal_checks,
        omega_unit_checks,
    ) = allocation_square_p13()
    (
        carrier_invariant_checks,
        carrier_family_checks,
        effective_nonzero,
        effective_total,
        forgotten_visible,
        forgotten_total,
    ) = carrier_gauge_checks(13)
    inverse_push_checks, inverse_push_baseline, theta_weights = (
        integral_inverse_push_checks(13)
    )

    p13_norm = next(row for row in norm_rows if row[0] == 13)
    require(p13_norm[1] == 12, "p13 automorphism count drift")
    require(p13_norm[2] == tuple(range(1, 13)), "p13 jet action drift")
    require(p13_norm[3] == (0,), "p13 fixed cotangent drift")

    print("THM-2820 BOOLEAN IDEMPOTENT RIGIDITY AND NORM-TOP JET NO-GO")
    print(f"pointwise_primes={PRIMES}; tangent_rows={tuple(tangent_rows)}")
    print("relative_derivations=zero; idempotent_dual_number_tangents=zero")
    print("group_ring_identity=N=(g-1)^(p-1); checked_primes=(2,3,5,7,11,13)")
    print(
        "p13_automorphisms=12; norm_top_fixed=12; "
        "cotangent_multipliers=(1,2,3,4,5,6,7,8,9,10,11,12)"
    )
    print("p13_common_fixed_cotangent=(0,); canonical_nonzero_jet=no")
    print(
        f"allocation_augmentation_degrees={augmentation_degrees}; "
        "mobius_degrees=((12,12),(0,12),(12,0),(0,0))"
    )
    print(
        f"central={central}; p13_valuations={valuations}; "
        f"hadamard={hadamard}; D3_over_H_mod13=1"
    )
    print(
        f"nonidempotent_interpolation_defect_atoms={interpolation_defect}/13"
    )
    print(
        f"normalized_raw_face_barycenters={barycenter_count}/169; "
        f"marked_gauge_relative_checks={gauge_difference_checks}; "
        "b(g^a h^b Omega)=(a,b)"
    )
    print(
        f"first_normal_ranks={first_normal_ranks}; "
        f"first_normal_translation_checks={first_normal_checks}"
    )
    print(
        f"omega_unit_inverse_checks={omega_unit_checks}; "
        "omega_translation_orbit=regular_169"
    )
    print(
        f"carrier_gauge_invariant_checks={carrier_invariant_checks}; "
        f"effective_family_checks={carrier_family_checks}; "
        f"effective_transverse_nonzero={effective_nonzero}/{effective_total}"
    )
    print(
        "marked_(L,q)=(W,(1,-1))_is_pure_gauge=yes; "
        f"forgotten_address_visible={forgotten_visible}/{forgotten_total}; "
        "test=q1+q2"
    )
    print(
        f"integral_inverse_push_checks={inverse_push_checks}; "
        f"pushed_profile={inverse_push_baseline}*N+delta_(q1+q2); "
        f"theta_weights={theta_weights}"
    )
    print(
        "scope=norm/Rees data alone retain no nonzero jet; a rooted translated "
        "aggregate has an exact affine Hasse coordinate, but the raw mixed "
        "face is joint-absent and a lawful transverse physical translation "
        "remains required"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
