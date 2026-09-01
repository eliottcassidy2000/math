#!/usr/bin/env python3
"""Dependency-free Q(sqrt(265)) audit for THM-4304."""

from __future__ import annotations

from fractions import Fraction as F


Q2 = tuple[F, F]


def add(x: Q2, y: Q2) -> Q2:
    return (x[0] + y[0], x[1] + y[1])


def mul(x: Q2, y: Q2) -> Q2:
    return (x[0] * y[0] + 265 * x[1] * y[1],
            x[0] * y[1] + x[1] * y[0])


def scale(x: Q2, c: F | int) -> Q2:
    c = F(c)
    return (c * x[0], c * x[1])


def power(x: Q2, n: int) -> Q2:
    out = (F(1), F(0))
    for _ in range(n):
        out = mul(out, x)
    return out


def require(condition: bool, label: str) -> None:
    if not condition:
        raise AssertionError(label)


def main() -> None:
    delta = F(2048, 45)
    u_base = F(-315392, 3645)
    u_sqrt = F(217088, 18225)
    discriminant_rows = []
    for epsilon in (1, -1):
        u = (u_base, epsilon * u_sqrt)
        root = (F(-15, 256), F(-3 * epsilon, 256))
        value = add(add(mul(u, power(root, 3)), scale(power(root, 2), delta)),
                    add(scale(root, F(8, 3)), (F(-1, 2), F(0))))
        derivative = add(add(scale(mul(u, power(root, 2)), 3),
                             scale(root, 2 * delta)), (F(8, 3), F(0)))
        second = add(scale(mul(u, root), 6), (2 * delta, F(0)))
        require(value == (0, 0), f"root value epsilon={epsilon}")
        require(derivative == (0, 0), f"root derivative epsilon={epsilon}")
        require(second != (0, 0), f"not triple epsilon={epsilon}")
        numerator = add(add(scale(power(u, 2), 820125), scale(u, 141926400)),
                        (F(-24696061952), F(0)))
        require(numerator == (0, 0), f"U root epsilon={epsilon}")
        discriminant_rows.append((epsilon, u, root))

    numerator_u1 = F(820125 + 141926400 - 24696061952)
    discriminant_u1 = -numerator_u1 / F(121500)
    require(discriminant_u1 == F(24553315427, 121500), "U=1 control")

    forced_r = F(9, 16)
    forced_u = F(1, 2) / forced_r**3
    require(forced_u == F(2048, 729), "triple U")
    require(-3 * forced_u * forced_r == F(-128, 27), "triple y2")
    require(F(-128, 27) != delta, "triple forbidden")

    s_value, beta_value, gamma_value = 1, 3, 16
    tau = s_value + beta_value
    weights = (3 * gamma_value, 2 * gamma_value + 4 * tau,
               gamma_value + 8 * tau, 12 * tau,
               gamma_value + 12 * beta_value)
    require(weights == (48, 48, 48, 48, 52), "hostile weight ledger")
    regime_b = [(k, F(12 - k, 2), F(12 + k, 2)) for k in (1, 2, 3)]
    require(all(gamma > 4 for _, gamma, _ in regime_b), "B cubic absent")
    require(all(required not in {1, 2, 3, 4, 5, 6, 8, 10}
                for _, _, required in regime_b), "B missing q-linear row")

    print("THM4304_REPEATED_FIRST_FACE_INDEPENDENT_V1")
    print("FIELD Q(sqrt265); DEGREE repeated_factor_linear; carrier_rational_graph")
    for epsilon, u, root in discriminant_rows:
        print(f"epsilon={epsilon:+d} U=({u[0]})+({u[1]})sqrt265 root=({root[0]})+({root[1]})sqrt265 multiplicity=2")
    print(f"HOSTILE weights={(s_value,beta_value,gamma_value)} face_weights={weights}")
    print(f"CONTROL U=1 discriminant={discriminant_u1}; triple=forbidden")
    print("VERDICT PASS REDUCED_CARRIERS_RATIONAL REFINEMENTS_OPEN")


if __name__ == "__main__":
    main()
