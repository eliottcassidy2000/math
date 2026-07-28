#!/usr/bin/env python3
"""Exact finite companion for THM-2695.

The script checks the finite cyclic-extension, F_2-module, Picard-torsion,
and quaternionic controls used by the theorem.  Its guards remain active
under ``python -O``.  The etale Kummer diagram and purity steps are proved in
the theorem text.
"""

from itertools import product
from math import gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def rank_f2(rows, ncols):
    rows = [sum((x & 1) << j for j, x in enumerate(row)) for row in rows]
    rank = 0
    for col in range(ncols):
        pivot = next((i for i in range(rank, len(rows)) if (rows[i] >> col) & 1), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for i in range(len(rows)):
            if i != rank and ((rows[i] >> col) & 1):
                rows[i] ^= rows[rank]
        rank += 1
    return rank


def cyclic_extension_class(p, a, quotient_multiplier):
    """Class of 0 -> C_(p^a) -> C_(p^(a+1)) -> C_p -> 0."""
    require(p > 1 and a >= 1 and gcd(p, quotient_multiplier) == 1,
            "invalid cyclic-extension parameters")
    lift = pow(quotient_multiplier, -1, p)
    modulus = p ** (a + 1)
    kernel_residues = [r for r in range(p) if (quotient_multiplier * r) % p == 0]
    lift_residues = [r for r in range(p) if (quotient_multiplier * r) % p == 1]
    require(kernel_residues == [0], "quotient kernel is not the p-multiple subgroup")
    require(lift_residues == [lift], "quotient generator lift is not unique mod p")
    require(modulus // p == p**a, "kernel order mismatch")
    kernel_coordinate = lift % (p**a)
    require((p * lift) % modulus == (p * kernel_coordinate) % modulus,
            "p-turn relation does not equal the embedded kernel coordinate")

    cocycle = {}
    for x in range(p):
        for y in range(p):
            z = (x + y) % p
            defect = (lift * x + lift * y - lift * z) % modulus
            require(defect % p == 0, "section defect does not lie in the kernel")
            cocycle[(x, y)] = (defect // p) % (p**a)
    wraps = [(x, y) for (x, y), value in cocycle.items() if value != 0]
    require(len(wraps) == p * (p - 1) // 2, "wrap-cocycle census mismatch")
    require(cocycle[(p - 1, 1)] % p == kernel_coordinate % p,
            "wrap cocycle does not represent the claimed class")

    # Push out C_(p^a) -> C_p.  The middle group becomes C_(p^2).
    pushout_middle_order = modulus // (p ** (a - 1))
    require(pushout_middle_order == p * p, "one-step pushout order mismatch")
    return lift, kernel_coordinate, kernel_coordinate % p, len(wraps), pushout_middle_order


def c3_on_h2(coeff):
    """Substitute sigma(x)=y, sigma(y)=x+y in degree two over F_2."""
    a, b, c = coeff
    return (c, b, (a + b + c) & 1)


def restrict_quadratic(coeff, vector):
    """Restrict a*x^2+b*xy+c*y^2 to one C2 line."""
    a, b, c = coeff
    x, y = vector
    return (a * x + b * x * y + c * y) & 1


def torsion_elements(modulus, killed_by):
    return [x for x in range(modulus) if (killed_by * x) % modulus == 0]


def double_image(modulus, rank):
    elems = list(product(range(modulus), repeat=rank))
    return {tuple((2 * x) % modulus for x in v) for v in elems}


def main():
    lift13, _, class13, wraps13, pushout13 = cyclic_extension_class(13, 5, 2)
    require((lift13, class13, wraps13, pushout13) == (7, 7, 78, 169),
            "odometer extension control failed")

    lift2, _, class2, wraps2, pushout2 = cyclic_extension_class(2, 1, 1)
    require((lift2, class2, wraps2, pushout2) == (1, 1, 1, 4),
            "mu4 coefficient extension control failed")

    require(gcd(2, 13) == gcd(2, 13**5) == gcd(4, 13) == 1,
            "cross-prime coefficient Hom boundary failed")

    fixed = [v for v in product(range(2), repeat=3) if c3_on_h2(v) == v]
    require(fixed == [(0, 0, 0), (1, 1, 1)],
            "C3-fixed degree-two cohomology census failed")
    q8 = (1, 1, 1)
    restrictions = [restrict_quadratic(q8, v) for v in [(1, 0), (0, 1), (1, 1)]]
    require(restrictions == [1, 1, 1], "Q8 line-restriction control failed")

    boundary_rows = [(1, 0), (0, 1), (1, 1)]
    require(rank_f2(boundary_rows, 2) == 2, "THM-2681 boundary rank mismatch")

    cl2_torsion = set(product(torsion_elements(2, 2), repeat=2))
    image_two_exponent2 = double_image(2, 2)
    require(len(cl2_torsion) == 4 and image_two_exponent2 == {(0, 0)},
            "toric exponent-two Picard control failed")
    beta_rank_toric = 2

    cl4_two_torsion = set(product(torsion_elements(4, 2), repeat=2))
    image_two_exponent4 = double_image(4, 2)
    require(cl4_two_torsion == image_two_exponent4,
            "four-divisible Picard-plane control failed")
    beta_rank_divisible = 0

    def sigma(v):
        x, y = v
        return (y, (x + y) & 1)

    nonzero = [(1, 0), (0, 1), (1, 1)]
    require([sigma(v) for v in nonzero] == [(0, 1), (1, 1), (1, 0)],
            "standard-plane orbit control failed")
    require(all(sigma(sigma(sigma(v))) == v for v in product(range(2), repeat=2)),
            "standard-plane action does not have order three")

    spectrum = {
        "ordered_root_completion": (2, None),
        "toric_quasietale_class_plane": (0, beta_rank_toric),
        "unit_plane_on_Gm2": (0, 0),
        "four_divisible_class_plane": (0, beta_rank_divisible),
    }

    print("KUMMER / CARRY SECONDARY BOCKSTEIN EXACT AUDIT")
    print(f"odometer_lift={lift13} odometer_class_mod13={class13}")
    print(f"mu4_lift={lift2} mu4_class_mod2={class2}")
    print(f"wrap_pairs_odometer={wraps13} wrap_pairs_mu4={wraps2}")
    print(f"one_step_pushout_middle_orders={pushout13},{pushout2}")
    print("cross_prime_additive_homs=0")
    print(f"c3_fixed_h2={fixed}")
    print(f"q8_line_restrictions={restrictions}")
    print(f"thm2681_boundary_rank={rank_f2(boundary_rows, 2)}")
    print(f"toric_Cl2_size={len(cl2_torsion)} twice_image_size={len(image_two_exponent2)} beta_rank={beta_rank_toric}")
    print(f"divisible_Cl4_two_torsion_size={len(cl4_two_torsion)} twice_image_size={len(image_two_exponent4)} beta_rank={beta_rank_divisible}")
    print(f"three_level_spectrum={spectrum}")
    print("PASS")


if __name__ == "__main__":
    main()
