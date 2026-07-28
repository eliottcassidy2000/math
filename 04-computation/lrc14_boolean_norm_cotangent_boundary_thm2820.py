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

    # Although N alone has zero augmentation and no normalized first
    # coordinate, the rooted Mobius face has J_00=1.  THM-2201's triangular
    # translation law therefore gives an affine barycenter.  Check the full
    # two-axis table directly rather than importing the formula.
    barycenters = {}
    for a in range(p):
        for b in range(p):
            translate = {
                (i, j): comb(a, i) * comb(b, j) % p
                for i in range(a + 1)
                for j in range(b + 1)
            }
            shifted = mul_bivariate(omega, translate, p)
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
    ) = allocation_square_p13()

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
        "scope=norm/Rees data alone retain no nonzero jet; a rooted translated "
        "raw face has an exact affine Hasse coordinate, but a lawful common "
        "off-sheet physical translation remains required"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
