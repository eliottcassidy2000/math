#!/usr/bin/env python3
"""THM-3202 exact moving-jet audit for C3[A^join r,A^join r,A^join r]."""

from functools import lru_cache
from hashlib import sha256
from math import comb, factorial


def require(condition, payload):
    if not condition:
        raise AssertionError(payload)


def finite_difference_profile(q_value, seed_order, r):
    """F_{A^join r}(c)=Delta^c(Q_A(j)^r)|_{j=0}."""
    total_order = seed_order * r
    values = [q_value(j) ** r for j in range(total_order + 1)]
    profile = [0] * (total_order + 2)
    row = values
    for c in range(1, total_order + 1):
        row = [row[j + 1] - row[j] for j in range(len(row) - 1)]
        profile[c] = row[0]
    return tuple(profile)


def c3_hamilton_from_factor_profile(profile):
    total_order = len(profile) - 2
    return 3 * sum(
        profile[c] ** 3
        + profile[c + 1] * profile[c] ** 2
        + profile[c + 1] ** 2 * profile[c]
        for c in range(1, total_order + 1)
    )


@lru_cache(maxsize=None)
def numerator_power(d):
    """Coefficients of (x+y+z+xy+xz+yz+3xyz)^d."""
    base = {
        (1, 0, 0): 1,
        (0, 1, 0): 1,
        (0, 0, 1): 1,
        (1, 1, 0): 1,
        (1, 0, 1): 1,
        (0, 1, 1): 1,
        (1, 1, 1): 3,
    }
    polynomial = {(0, 0, 0): 1}
    for _ in range(d):
        product = {}
        for left, left_value in polynomial.items():
            for right, right_value in base.items():
                exponent = tuple(left[i] + right[i] for i in range(3))
                product[exponent] = (
                    product.get(exponent, 0) + left_value * right_value
                )
        polynomial = product
    return polynomial


def c3_walk_family_coefficient(d, c1, c2, c3):
    """Return [x^c1 y^c2 z^c3] W_C3(x,y,z)^d."""
    polynomial = numerator_power(d)
    answer = 0
    for q in range(min(c1, c2, c3) + 1):
        residual = (c1 - q, c2 - q, c3 - q)
        coefficient = polynomial.get(residual, 0)
        if coefficient:
            answer += coefficient * comb(d + q - 1, q)
    return answer


def c3_output_ordered_profile_coordinate(factor_profile, d):
    """F_{C3[S,S,S]}(d) from the THM-3121 quotient-walk kernel."""
    total_order = len(factor_profile) - 2
    answer = 0
    for c1 in range(1, total_order + 1):
        for c2 in range(1, total_order + 1):
            for c3 in range(1, total_order + 1):
                weight = c3_walk_family_coefficient(d, c1, c2, c3)
                if weight:
                    answer += (
                        weight
                        * factor_profile[c1]
                        * factor_profile[c2]
                        * factor_profile[c3]
                    )
    return answer


def digits_hash(value):
    decimal = str(value)
    return len(decimal), sha256(decimal.encode()).hexdigest()


def q_k1(j):
    return j


def q_c3(j):
    return j**3 + 2 * j


def audit_seed(name, seed_order, q_value, levels):
    rows = []
    coordinate_checks = 0
    for r in levels:
        total_order = seed_order * r
        profile = finite_difference_profile(q_value, seed_order, r)
        require(
            profile[total_order] == factorial(total_order),
            (name, r, "singleton endpoint", profile[total_order]),
        )
        hamilton = c3_hamilton_from_factor_profile(profile)
        require(
            hamilton == c3_output_ordered_profile_coordinate(profile, 1),
            (name, r, "C3 d=1 kernel"),
        )
        fixed = []
        for d in range(1, min(3, 3 * total_order) + 1):
            ordered_coordinate = c3_output_ordered_profile_coordinate(profile, d)
            if total_order >= d:
                lower = (
                    3**d
                    * comb(total_order - 1, d - 1)
                    * factorial(total_order) ** 3
                )
                require(
                    ordered_coordinate >= lower,
                    (name, r, d, ordered_coordinate, lower),
                )
            require(
                ordered_coordinate % factorial(d) == 0,
                (name, r, d, "unordered integrality"),
            )
            fixed.append(ordered_coordinate // factorial(d))
            coordinate_checks += 1
        rows.append((r, total_order, hamilton, tuple(fixed)))
    return rows, coordinate_checks


def main():
    k1_rows, k1_checks = audit_seed("K1", 1, q_k1, range(1, 9))
    c3_rows, c3_checks = audit_seed("C3", 3, q_c3, range(1, 6))

    # K1 gives the ordered Stirling/Fubini layer exactly.
    for r, total_order, _, _ in k1_rows:
        profile = finite_difference_profile(q_k1, 1, r)
        for c in range(1, total_order + 1):
            require(
                profile[c]
                == sum(
                    (-1) ** (c - j) * comb(c, j) * j**r
                    for j in range(c + 1)
                ),
                ("K1 Stirling boundary", r, c),
            )

    # Same-H scalar hostile from THM-3134.
    profile_40 = (0, 15, 78, 198, 240, 120, 0)
    profile_76 = (0, 15, 90, 210, 240, 120, 0)
    hostile_40 = c3_hamilton_from_factor_profile(profile_40)
    hostile_76 = c3_hamilton_from_factor_profile(profile_76)
    require(hostile_40 != hostile_76, (hostile_40, hostile_76))

    print("C3 REPEATED-JOIN MOVING JET -- THM-3202 EXACT COMPANION")
    print(
        "status=PROVED+VERIFIED-EXACT;"
        "no_canonical_runtime_imports=true"
    )
    print("identity=F_Br(c)=sum_j(-1)^(c-j)binom(c,j)Q_A(j)^r")
    print("identity=F_Ur(d)=sum_c [x^c]W_C3^d product_i F_Br(c_i)")
    print("H_special=3sum_c(F_c^3+F_(c+1)F_c^2+F_(c+1)^2F_c)")
    print(
        "superexponential_floor=F_Ur(d)>=3^d*binom(ar-1,d-1)*"
        "((ar)!)^3;fixed_d>=1"
    )
    print(
        "verdict=no_fixed_finite_polynomial_exponential;"
        "no_eventual_constant_coefficient_recurrence;"
        "no_rational_ordinary_GF"
    )
    print(f"K1_levels=1..8;fixed_d_checks={k1_checks}")
    for r, order, hamilton, fixed in k1_rows[:5]:
        print(f"K1_r={r};block_order={order};H={hamilton};p_d_1_to_3={fixed}")
    k1_last = k1_rows[-1]
    print(
        f"K1_r={k1_last[0]};block_order={k1_last[1]};"
        f"H_digits_sha256={digits_hash(k1_last[2])}"
    )
    print(f"C3_levels=1..5;fixed_d_checks={c3_checks}")
    for r, order, hamilton, fixed in c3_rows[:3]:
        print(f"C3_r={r};block_order={order};H={hamilton};p_d_1_to_3={fixed}")
    for r, order, hamilton, _ in c3_rows[3:]:
        print(
            f"C3_r={r};block_order={order};"
            f"H_digits_sha256={digits_hash(hamilton)}"
        )
    print(
        f"same_H_hostile=mask40_vs_mask76;H_factor=15;"
        f"C3_equal_lift_H=({hostile_40},{hostile_76});scalar_H_REFUTED"
    )


if __name__ == "__main__":
    main()
