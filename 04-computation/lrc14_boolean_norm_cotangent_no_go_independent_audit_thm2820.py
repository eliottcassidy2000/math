#!/usr/bin/env python3
"""Independent hostile audit for THM-2820.

This audit deliberately does not import the primary companion.  The cyclic
group algebra is represented in the group-element basis, rather than in the
truncated augmentation-coordinate basis used by the primary computation.
"""

import ast
from math import comb
from pathlib import Path


def require(condition, message):
    if not bool(condition):
        raise RuntimeError(message)


def contains_assert(path):
    tree = ast.parse(Path(path).read_text(encoding="utf-8"))
    return any(isinstance(node, ast.Assert) for node in ast.walk(tree))


def prime(n):
    if n < 2:
        return False
    return all(n % d for d in range(2, int(n**0.5) + 1))


def cyclic_product(left, right, p):
    """Product in F_p[C_p], in the basis 1,g,...,g^(p-1)."""
    answer = [0] * p
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            answer[(i + j) % p] = (answer[(i + j) % p] + a * b) % p
    return tuple(answer)


def cyclic_power(base, exponent, p):
    answer = (1,) + (0,) * (p - 1)
    factor = base
    power = exponent
    while power:
        if power & 1:
            answer = cyclic_product(answer, factor, p)
        factor = cyclic_product(factor, factor, p)
        power //= 2
    return answer


def multiplier_image(vector, u, p):
    """The group automorphism g -> g^u, in the group-element basis."""
    answer = [0] * p
    for exponent, coefficient in enumerate(vector):
        answer[(u * exponent) % p] = (
            answer[(u * exponent) % p] + coefficient
        ) % p
    return tuple(answer)


def augmentation(vector, p):
    return sum(vector) % p


def cotangent_coordinate(vector, p):
    """Coefficient of h=g-1 modulo h^2 on the augmentation ideal."""
    require(augmentation(vector, p) == 0, "jet requested off augmentation ideal")
    return sum(j * coefficient for j, coefficient in enumerate(vector)) % p


def check_ring_independence():
    """Hostile finite-ring controls, including characteristic two.

    For every idempotent e, (2e-1)^2=1.  Consequently multiplication by
    2e-1 has zero kernel, which is simultaneously the derivation-basis and
    dual-number-lift gate in the proof.
    """
    idempotent_rows = 0
    tangent_rows = 0
    for modulus in range(2, 25):
        idempotents = [
            e for e in range(modulus) if (e * e - e) % modulus == 0
        ]
        for e in idempotents:
            involution = (2 * e - 1) % modulus
            require(
                involution * involution % modulus == 1 % modulus,
                f"failed involution in Z/{modulus} at e={e}",
            )
            kernel = [
                m
                for m in range(modulus)
                if involution * m % modulus == 0
            ]
            require(
                kernel == [0],
                f"nonzero dual-number tangent in Z/{modulus} at e={e}",
            )
            tangent_rows += modulus
            idempotent_rows += 1

    # Coordinate idempotents form a K-basis of K^X.  Test the actual
    # coordinatewise multiplier on nonzero hostile candidate images in
    # modules M=(Z/nZ)^X, including n=2.
    coordinate_rows = 0
    for modulus in range(2, 13):
        for size in range(1, 4):
            for coordinate in range(size):
                signs = [
                    (1 if j == coordinate else -1) % modulus
                    for j in range(size)
                ]
                for encoded in range(modulus**size):
                    value = encoded
                    candidate = []
                    for _ in range(size):
                        candidate.append(value % modulus)
                        value //= modulus
                    killed = all(
                        signs[j] * candidate[j] % modulus == 0
                        for j in range(size)
                    )
                    require(
                        not killed or all(entry == 0 for entry in candidate),
                        "nonzero coordinate-idempotent derivation image",
                    )
                    coordinate_rows += 1
    return idempotent_rows, tangent_rows, coordinate_rows


def check_group_rings():
    """Audit N=(g-1)^(p-1) and the top/cotangent character mismatch."""
    prime_rows = 0
    automorphism_rows = 0
    jet_basis_rows = 0
    hom_rows = 0
    p2_fixed_nonzero = False
    for p in [n for n in range(2, 38) if prime(n)]:
        one = (1,) + (0,) * (p - 1)
        g = (0, 1) + (0,) * (p - 2)
        h = tuple((g[j] - one[j]) % p for j in range(p))
        norm = (1,) * p
        zero = (0,) * p

        require(cyclic_power(h, p, p) == zero, f"h^p != 0 at p={p}")
        require(
            cyclic_power(h, p - 1, p) == norm,
            f"N != h^(p-1) at p={p}",
        )
        require(
            cyclic_product(h, norm, p) == zero,
            f"hN != 0 at p={p}",
        )

        for u in range(1, p):
            require(
                multiplier_image(norm, u, p) == norm,
                f"norm not fixed at p={p}, u={u}",
            )
            image_h = multiplier_image(h, u, p)
            require(
                cotangent_coordinate(image_h, p) == u,
                f"wrong cotangent character at p={p}, u={u}",
            )

            # The elements g^j-1 span the augmentation ideal.  Checking
            # their jet coordinates audits phi_u(v)=u v on I/I^2 without
            # choosing the primary script's h-adic representation.
            for j in range(1, p):
                basis = list(zero)
                basis[j] = 1
                basis[0] = -1 % p
                basis = tuple(basis)
                require(
                    cotangent_coordinate(multiplier_image(basis, u, p), p)
                    == u * cotangent_coordinate(basis, p) % p,
                    f"jet covariance at p={p}, u={u}, j={j}",
                )
                jet_basis_rows += 1
            automorphism_rows += 1

        fixed_scalars = [
            c
            for c in range(p)
            if all(u * c % p == c for u in range(1, p))
        ]
        if p == 2:
            require(fixed_scalars == [0, 1], "p=2 boundary disappeared")
            p2_fixed_nonzero = True
        else:
            require(
                fixed_scalars == [0],
                f"nonzero equivariant top-to-jet scalar at p={p}",
            )
        hom_rows += p * max(1, p - 1)
        prime_rows += 1
    require(p2_fixed_nonzero, "odd-prime hypothesis was not tested sharply")
    return prime_rows, automorphism_rows, jet_basis_rows, hom_rows


def poly2_product(left, right, p):
    """Product in F_p[e,h]/(e^p,h^p), stored sparsely."""
    answer = {}
    for (i, j), a in left.items():
        for (k, ell), b in right.items():
            if i + k >= p or j + ell >= p:
                continue
            key = (i + k, j + ell)
            answer[key] = (answer.get(key, 0) + a * b) % p
    return {key: value for key, value in answer.items() if value}


def poly2_power(base, exponent, p):
    answer = {(0, 0): 1}
    factor = base
    power = exponent
    while power:
        if power & 1:
            answer = poly2_product(answer, factor, p)
        factor = poly2_product(factor, factor, p)
        power //= 2
    return answer


def image_axis(axis, u, p):
    """(1+axis)^u-1 as an exact truncated bivariate polynomial."""
    return {
        ((degree, 0) if axis == 0 else (0, degree)): comb(u, degree) % p
        for degree in range(1, u + 1)
        if comb(u, degree) % p
    }


def check_tensor_thirteen():
    p = 13
    top = {(p - 1, p - 1): 1}
    axis_rows = 0
    for u in range(1, p):
        image_e = image_axis(0, u, p)
        for v in range(1, p):
            image_h = image_axis(1, v, p)
            image_top = poly2_product(
                poly2_power(image_e, p - 1, p),
                poly2_power(image_h, p - 1, p),
                p,
            )
            require(image_top == top, f"tensor top moved at ({u},{v})")
            require(image_e.get((1, 0), 0) == u, "first axis character")
            require(image_h.get((0, 1), 0) == v, "second axis character")
            mixed = poly2_product(image_e, image_h, p)
            require(
                mixed.get((1, 1), 0) == u * v % p,
                "mixed cotangent character",
            )
            axis_rows += 1

    fixed_axes = []
    for a in range(p):
        for b in range(p):
            if all(
                (u * a % p, v * b % p) == (a, b)
                for u in range(1, p)
                for v in range(1, p)
            ):
                fixed_axes.append((a, b))
    require(fixed_axes == [(0, 0)], "nonzero invariant cotangent axis")

    fixed_mixed = [
        c
        for c in range(p)
        if all(
            u * v * c % p == c
            for u in range(1, p)
            for v in range(1, p)
        )
    ]
    require(fixed_mixed == [0], "nonzero invariant mixed first interaction")
    return axis_rows, len(fixed_axes), len(fixed_mixed)


def main():
    require(not contains_assert(Path(__file__)), "truth-bearing assert node")
    idempotents, tangents, coordinate_rows = check_ring_independence()
    primes, automorphisms, jet_rows, hom_rows = check_group_rings()
    tensor_rows, fixed_axes, fixed_mixed = check_tensor_thirteen()
    print("THM-2820 independent Boolean/norm hostile audit")
    print(
        f"finite_ring_idempotents={idempotents} tangent_candidates={tangents} "
        f"coordinate_module_rows={coordinate_rows}: arbitrary-characteristic "
        "involution gate PASS"
    )
    print(
        f"cyclic_basis_primes={primes} automorphisms={automorphisms} "
        f"augmentation_spanning_rows={jet_rows} hom_rows={hom_rows}: "
        "norm/top and cotangent-character split PASS"
    )
    print(
        f"F13_tensor_pairs={tensor_rows} fixed_axis_vectors={fixed_axes} "
        f"fixed_mixed_scalars={fixed_mixed}: independent-character no-go PASS"
    )
    print("p=2 boundary PASS: odd-p Hom no-go is not overclaimed")
    print(
        "PASS scope=algebraic equivariance only; THM-2806 D3 is not "
        "identified with a group norm; no LRC row is excluded"
    )


if __name__ == "__main__":
    main()
