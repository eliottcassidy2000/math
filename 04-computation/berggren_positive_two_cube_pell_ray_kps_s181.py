#!/usr/bin/env python3
"""Exact audit of the positive Berggren/two-cube Pell ray (THM-3375).

All theorem identities are checked with Python integers.  The infinite-family
proof itself is algebraic: a norm-one unit preserves the Pell conic and a
mod-29 subsequence preserves the integral parameter u=h/29.
"""

from __future__ import annotations

from hashlib import sha256
from math import gcd, isqrt
from pathlib import Path


D = 621
EPSILON = (7775, 312)
SEED = (5059, 203)
MODULUS = 29
ALLOWED_R63 = (2, 4, 10, 11, 16, 24, 25, 31, 37, 38, 46, 51, 52, 58, 60)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def mul(z: tuple[int, int], w: tuple[int, int]) -> tuple[int, int]:
    a, b = z
    c, d = w
    return a * c + D * b * d, a * d + b * c


def power(z: tuple[int, int], n: int) -> tuple[int, int]:
    out = (1, 0)
    while n:
        if n & 1:
            out = mul(out, z)
        z = mul(z, z)
        n >>= 1
    return out


def norm(z: tuple[int, int]) -> int:
    return z[0] * z[0] - D * z[1] * z[1]


def family_member(j: int, epsilon5: tuple[int, int]) -> dict[str, int]:
    w, h = mul(SEED, power(epsilon5, j))
    require(h % MODULUS == 0, "h is not divisible by 29")
    u = h // MODULUS
    d = 841 * u * u + 2
    e = u * w
    a = 21866 * u**3 + 81 * u
    q = 568516 * u**4 + 2860 * u * u + 1
    require((d + e) % 2 == 0 and (d - e) % 2 == 0, "cube parity")
    x = (d + e) // 2
    y = (d - e) // 2
    require((a - 1) % 2 == 0, "depth parity")
    r = (a - 1) // 2
    return {"j": j, "w": w, "h": h, "u": u, "d": d, "e": e,
            "a": a, "q": q, "x": x, "y": y, "r": r}


def cube_sum_support(modulus: int) -> set[int]:
    cubes = {pow(x, 3, modulus) for x in range(modulus)}
    return {(x + y) % modulus for x in cubes for y in cubes}


def q_value(r: int) -> int:
    return (2 * r + 1) ** 2 + 2


def orbit_residue_period(
    modulus: int, epsilon5: tuple[int, int]
) -> tuple[int, set[int]]:
    """Period and r-residue support of the distinguished (W,u) orbit."""
    a_unit, b_unit = epsilon5
    require(b_unit % MODULUS == 0, "B/29 is not integral")
    b_div = b_unit // MODULUS
    state = (SEED[0] % modulus, (SEED[1] // MODULUS) % modulus)
    initial = state
    residues: set[int] = set()
    period = 0
    while True:
        w, u = state
        a = 21866 * u**3 + 81 * u
        residues.add(((a - 1) * pow(2, -1, modulus)) % modulus)
        state = (
            (a_unit * w + D * b_unit * MODULUS * u) % modulus,
            (b_div * w + a_unit * u) % modulus,
        )
        period += 1
        if state == initial:
            return period, residues
        require(period <= modulus * modulus, "modular orbit did not close")


def main() -> None:
    require(norm(EPSILON) == 1, "unit norm")
    require(norm(SEED) == 2692, "seed norm")
    epsilon5 = power(EPSILON, 5)
    require(epsilon5[0] % MODULUS == MODULUS - 1, "unit scalar mod 29")
    require(epsilon5[1] % MODULUS == 0, "unit radical coefficient mod 29")

    # Coefficientwise verification of the two polynomial identities in u.
    require(21866**2 == 841 * 568516, "u^6 identity")
    require(2 * 21866 * 81 == 841 * 2860 + 2 * 568516, "u^4 identity")
    require(81**2 == 841 + 2 * 2860, "u^2 identity")
    require(4 * 568516 - 841**2 == 3 * 522261, "Pell u^4 identity")
    require(4 * 2860 - 4 * 841 == 3 * 2692, "Pell u^2 identity")
    require(522261 == 621 * 29**2, "Pell discriminant factor")

    members = [family_member(j, epsilon5) for j in range(8)]
    for member in members:
        w = member["w"]
        h = member["h"]
        u = member["u"]
        d = member["d"]
        e = member["e"]
        a = member["a"]
        q = member["q"]
        x = member["x"]
        y = member["y"]
        r = member["r"]
        require(w * w - D * h * h == 2692, "Pell orbit")
        require(w * w - 522261 * u * u == 2692, "rescaled Pell orbit")
        require(a * a + 2 == d * q, "factor identity")
        require(4 * q - d * d == 3 * e * e, "square identity")
        require(x > y > 0 and x - y == e and x + y == d, "positive chamber")
        require(x**3 + y**3 == a * a + 2 == q_value(r), "cube/Berggren identity")
        require(gcd(d, gcd(2 * a, q)) == 1, "primitive discriminant form")
        require((2 * a) ** 2 - 4 * d * q == -8, "discriminant -8")
        require(w < MODULUS * MODULUS * u, "W/h<29")

    require(all(members[i]["h"] < members[i + 1]["h"] for i in range(7)),
            "strict orbit growth")
    first = members[0]
    require(
        (first["x"], first["y"], first["a"], first["r"])
        == (38312, 2899, 7500605, 3750302),
        "first family member",
    )

    support7 = cube_sum_support(7)
    support9 = cube_sum_support(9)
    require(support7 == {0, 1, 2, 5, 6}, "cube-sum support mod 7")
    require(support9 == {0, 1, 2, 7, 8}, "cube-sum support mod 9")
    allowed = tuple(
        r for r in range(63)
        if q_value(r) % 7 in support7 and q_value(r) % 9 in support9
    )
    require(allowed == ALLOWED_R63, "CRT depth classes")
    require(all(member["r"] % 63 in ALLOWED_R63 for member in members),
            "family misses local support")

    period63, orbit_r63 = orbit_residue_period(63, epsilon5)
    period5, orbit_r5 = orbit_residue_period(5, epsilon5)
    require(period63 == 24 and orbit_r63 == {38, 51, 60}, "orbit mod 63 scar")
    require(orbit_r5 == {2}, "orbit mod 5 scar")
    require(epsilon5[0] % 11 == 1 and epsilon5[1] % 11 == 0,
            "unit mod 11 scar")
    require((SEED[0] % 11, SEED[1] % 11) == (10, 5), "seed mod 11 scar")

    # Exact rational brackets used in the positivity proof.
    require(24 * SEED[1] < SEED[0] < 25 * SEED[1], "seed ratio bracket")
    require(SEED[0] < MODULUS * SEED[1], "seed ratio below 29")
    require(isqrt(D) == 24 and 24 * 24 < D < 25 * 25, "sqrt bracket")

    semantic = ";".join(
        f"{m['j']}:{m['r']}:{m['x']}:{m['y']}:{m['d']}:{m['a']}"
        for m in members
    )
    source = Path(__file__).read_bytes().replace(b"\r\n", b"\n")
    print("BERGGREN POSITIVE TWO-CUBE PELL RAY")
    print("status=PROVED+VERIFIED-EXACT")
    print(f"epsilon={EPSILON[0]}+{EPSILON[1]}*sqrt({D});norm={norm(EPSILON)}")
    print(
        f"epsilon5_scalar_digits={len(str(epsilon5[0]))};"
        f"epsilon5_radical_digits={len(str(epsilon5[1]))};"
        f"epsilon5_mod29=({epsilon5[0] % 29},{epsilon5[1] % 29})"
    )
    print(f"seed={SEED[0]}+{SEED[1]}*sqrt({D});norm={norm(SEED)};h_div29=true")
    print("pell=W^2-522261*u^2=2692;h=29*u")
    print(
        "compiler=d=841*u^2+2;e=u*W;a=21866*u^3+81*u;"
        "q=568516*u^4+2860*u^2+1"
    )
    print("identities=a^2+2=d*q;4*q-d^2=3*e^2;verified_coefficientwise=true")
    print(
        f"first=j=0;r={first['r']};x={first['x']};y={first['y']};"
        f"a={first['a']};Q={first['a']**2+2}"
    )
    print("positivity=0<e<d;ratio_proof=sqrt(621)<W/h<=5059/203<29")
    print("limit=e/d->sqrt(621)/29<1;family_growth=exponential")
    print(f"local_depth_classes_mod63={list(ALLOWED_R63)};count=15;density=5/21")
    print("ambient_local_sieve=mod7_and_mod9_are_the_complete_p_adic_cube_sum_obstruction_for_Q_r")
    print(
        f"distinguished_orbit_scar=period_mod63={period63};"
        f"r_mod63={sorted(orbit_r63)};r_mod5={sorted(orbit_r5)};"
        "epsilon5_mod11=(1,0)"
    )
    print("form_compiler=primitive_[d,2a,q]_has_discriminant_-8_and_class_number_one")
    print(f"first8_semantic_sha256={sha256(semantic.encode('ascii')).hexdigest()}")
    print(
        "boundary=no_asymptotic_for_all_positive_intersections;no_fixed_d_"
        "positive_infinite_orbit;no_LRC_FC_JC_or_AMM_consequence"
    )
    print(f"source_lf_sha256={sha256(source).hexdigest()}")


if __name__ == "__main__":
    main()
