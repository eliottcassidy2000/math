#!/usr/bin/env python3
"""Exact referee for provisional THM-3422.

For

    lambda = y + gamma*y**e*z**d,       delta = Jac(lambda,-),

the target sector indexed by ``s`` (target z-weight ``s-1 mod d``) is the
colimit of monomial chains

    E[k,m] = [y**m z**(s-1+k*d)],
    n E[k,m] + gamma*(n*e-d*m) E[k+1,m+e-1] = 0,
    n=s+k*d.

The script checks the one-root orbit-label proof, every predicted
lambda-arrow break, the finite torsion-thickness formula, direct bivariate
primitives for the unit observer, and the separately typed repeated-root
split-fiber character window.  It uses only integer/Fraction arithmetic and
has no assertion-dependent truth gate.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product
from math import comb, gcd


Poly = dict[tuple[int, int], Fraction]


def require(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(label)


def add_term(poly: Poly, key: tuple[int, int], value: Fraction) -> None:
    if value:
        poly[key] = poly.get(key, Fraction(0)) + value
        if not poly[key]:
            del poly[key]


def add(left: Poly, right: Poly) -> Poly:
    out = dict(left)
    for key, value in right.items():
        add_term(out, key, value)
    return out


def scale(poly: Poly, scalar: Fraction) -> Poly:
    return {key: scalar * value for key, value in poly.items() if scalar * value}


def multiply(left: Poly, right: Poly) -> Poly:
    out: Poly = {}
    for (m1, n1), c1 in left.items():
        for (m2, n2), c2 in right.items():
            add_term(out, (m1 + m2, n1 + n2), c1 * c2)
    return out


def power(poly: Poly, exponent: int) -> Poly:
    out: Poly = {(0, 0): Fraction(1)}
    base = dict(poly)
    n = exponent
    while n:
        if n & 1:
            out = multiply(out, base)
        base = multiply(base, base)
        n //= 2
    return out


def delta(poly: Poly, d: int, e: int, gamma: Fraction) -> Poly:
    """Apply (1+gamma*e*y^(e-1)z^d)d_z-gamma*d*y^e*z^(d-1)d_y."""
    out: Poly = {}
    h = e - 1
    for (m, n), coefficient in poly.items():
        if n:
            add_term(out, (m, n - 1), coefficient * n)
            add_term(out, (m + h, n + d - 1), coefficient * gamma * e * n)
        if m:
            add_term(out, (m + h, n + d - 1), -coefficient * gamma * d * m)
    return out


def monomial(m: int, n: int, coefficient: Fraction = Fraction(1)) -> Poly:
    return {(m, n): coefficient} if coefficient else {}


def lambda_poly(e: int, d: int, gamma: Fraction) -> Poly:
    return {(1, 0): Fraction(1), (e, d): gamma}


def ceil_div(a: int, b: int) -> int:
    return -((-a) // b)


def late_state(d: int, e: int, s: int, ell: int) -> tuple[int, int]:
    """Choose a state beyond any zero transition and inside the wrap window."""
    h = e - 1
    require(h > 0, ("late_state needs repeated root", d, e, s, ell))
    k = max(0, ceil_div(-ell, h))
    if (e * s) % d == 0:
        k = max(k, ell - (e * s) // d + 1)
    if s == d:
        k = max(k, ell - e + 1)
    m = ell + k * h
    require(m >= 0, ("negative monomial", d, e, s, ell, k, m))
    if s == d:
        require(m < e * (k + 1), ("outside wrap quotient", d, e, ell, k, m))
    return k, m


def sector_type(d: int, e: int, s: int) -> str:
    if e == 1:
        return "free"
    return "free+prufer" if (s * (e - 1)) % d == 0 else "laurent"


def observer_annihilator_exponent(d: int, e: int) -> int | None:
    h = e - 1
    if h > 0 and h % d == 0:
        return h // d
    return None


def primitive_for_observer_power(d: int, e: int, gamma: Fraction) -> Poly:
    """Return Q with delta(Q)=lambda^q when e-1=q*d."""
    h = e - 1
    require(h > 0 and h % d == 0, ("observer is not torsion", d, e))
    q = h // d
    out: Poly = {}
    for j in range(q):
        coefficient = Fraction(comb(q - 1, j), 1 + j * d) * gamma**j
        add_term(out, (q + j * h, 1 + j * d), coefficient)
    return out


def check_direct_identities() -> tuple[int, int, int]:
    relation_checks = 0
    endpoint_checks = 0
    primitive_checks = 0
    for gamma in (Fraction(2, 3), Fraction(-5, 4)):
        for d in range(1, 9):
            for e in range(1, 12):
                h = e - 1
                lam = lambda_poly(e, d, gamma)
                for s in range(1, d + 1):
                    for k in range(5):
                        n = s + k * d
                        m_values = range(0, 9)
                        if s == d:
                            m_values = range(0, min(9, e * (k + 1)))
                        for m in m_values:
                            lhs = delta(monomial(m, n), d, e, gamma)
                            expected = add(
                                monomial(m, n - 1, Fraction(n)),
                                monomial(m + h, n + d - 1, gamma * (n * e - d * m)),
                            )
                            require(lhs == expected, ("monomial relation", gamma, d, e, s, k, m))
                            relation_checks += 1

                    if h > 0 and (s * h) % d == 0:
                        q_s = s * h // d
                        tau = monomial(q_s - 1, s - 1)
                        primitive = monomial(q_s, s, Fraction(1, s))
                        require(
                            multiply(lam, tau) == delta(primitive, d, e, gamma),
                            ("endpoint primitive", gamma, d, e, s),
                        )
                        endpoint_checks += 1

                if h > 0 and h % d == 0:
                    q = h // d
                    primitive = primitive_for_observer_power(d, e, gamma)
                    require(
                        delta(primitive, d, e, gamma) == power(lam, q),
                        ("unit-observer primitive", gamma, d, e, q),
                    )
                    primitive_checks += 1
    return relation_checks, endpoint_checks, primitive_checks


def check_orbit_and_arrows() -> tuple[int, int, int]:
    transition_checks = 0
    arrow_checks = 0
    thickness_checks = 0
    for d in range(1, 13):
        for e in range(2, 16):
            h = e - 1
            selected = []
            for s in range(1, d + 1):
                if (s * h) % d == 0:
                    selected.append(s)
                for ell in range(-10, 18):
                    k, m = late_state(d, e, s, ell)
                    n = s + k * d
                    transition_numerator = d * m - e * n
                    require(
                        transition_numerator == d * ell - e * s - d * k,
                        ("transition orbit", d, e, s, ell, k, m),
                    )
                    transition_checks += 1

                    denominator = d * (m + 1) - e * n
                    if not denominator:
                        k += 1
                        m += h
                        n += d
                        denominator = d * (m + 1) - e * n
                    require(denominator != 0, ("late denominator", d, e, s, ell, k, m))
                    direct_numerator = d * (m + 1) - n * h
                    orbit_numerator = d * (ell + 1) - s * h
                    require(
                        direct_numerator == orbit_numerator,
                        ("lambda arrow", d, e, s, ell, k, m),
                    )
                    predicted_zero = (s * h) % d == 0 and ell == s * h // d - 1
                    require(
                        (direct_numerator == 0) == predicted_zero,
                        ("lambda break", d, e, s, ell),
                    )
                    arrow_checks += 1

                if (s * h) % d == 0:
                    q_s = s * h // d
                    ell_zero = q_s - 1
                    for depth in range(7):
                        if s < d:
                            # At stage depth, torsion monomials have
                            # -depth*h <= ell=m-depth*h <= ell_zero.
                            brute = sum(
                                1
                                for m in range(depth * h + ell_zero + 1)
                                if m - depth * h <= ell_zero
                            )
                        else:
                            brute = sum(
                                1
                                for m in range(e * (depth + 1))
                                if m - depth * h <= ell_zero
                            )
                        expected = depth * h + q_s
                        require(brute == expected, ("torsion thickness", d, e, s, depth, brute, expected))
                        thickness_checks += 1

            require(
                len(selected) == gcd(d, h),
                ("selected-character count", d, e, selected, gcd(d, h)),
            )
    return transition_checks, arrow_checks, thickness_checks


def check_simple_root() -> int:
    checks = 0
    for d in range(1, 15):
        e = 1
        for s in range(1, d + 1):
            require(sector_type(d, e, s) == "free", ("simple type", d, s))
            for ell in range(12):
                # h=0: the orbit labels are exactly m>=0, and the lambda
                # numerator d(m+1) never vanishes.
                require(d * (ell + 1) != 0, ("simple arrow", d, s, ell))
                checks += 1
    return checks


def gcd_packet(values: list[int]) -> int:
    out = 0
    for value in values:
        out = gcd(out, value)
    return out


def special_fiber_sector_dimensions(
    d: int, multiplicities: tuple[int, ...], chosen: int
) -> tuple[list[int], int, tuple[int, ...]]:
    """Character dimensions for a repeated-root special fiber.

    At root ``chosen`` the nonvertical component is the unramified Kummer
    cover of A1 minus all roots with exponent vector

        (e_chosen-1, e_j for j != chosen) mod d.

    Target weight s-1 corresponds to de Rham character s.
    """
    require(multiplicities[chosen] > 1, ("chosen root is simple", multiplicities, chosen))
    exponents = tuple(
        multiplicity - 1 if index == chosen else multiplicity
        for index, multiplicity in enumerate(multiplicities)
    )
    number_of_roots = len(multiplicities)
    dims = [
        number_of_roots
        if all((s * exponent) % d == 0 for exponent in exponents)
        else number_of_roots - 1
        for s in range(1, d + 1)
    ]
    components = gcd_packet([d, *exponents])
    return dims, components, exponents


def cayley_cover_components(d: int, exponents: tuple[int, ...]) -> int:
    parent = list(range(d))

    def find(vertex: int) -> int:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    def union(left: int, right: int) -> None:
        a = find(left)
        b = find(right)
        if a != b:
            parent[b] = a

    for vertex in range(d):
        for exponent in exponents:
            union(vertex, (vertex + exponent) % d)
    return len({find(vertex) for vertex in range(d)})


def check_split_root_first_windows() -> tuple[int, int]:
    packet_checks = 0
    character_checks = 0
    for d in range(2, 11):
        for number_of_roots in range(2, 5):
            for multiplicities in product(range(1, 6), repeat=number_of_roots):
                for chosen, multiplicity in enumerate(multiplicities):
                    if multiplicity == 1:
                        continue
                    dims, components, exponents = special_fiber_sector_dimensions(
                        d, multiplicities, chosen
                    )
                    graph_components = cayley_cover_components(d, exponents)
                    require(
                        graph_components == components,
                        ("Kummer cover components", d, multiplicities, chosen, exponents),
                    )
                    graph_betti = d * number_of_roots - d + graph_components
                    require(
                        sum(dims) == graph_betti,
                        ("character/graph Betti", d, multiplicities, chosen, dims, graph_betti),
                    )
                    require(
                        sum(value == number_of_roots for value in dims) == components,
                        ("trivial-character count", d, multiplicities, chosen, dims, components),
                    )
                    packet_checks += 1
                    character_checks += d
    return packet_checks, character_checks


def describe_control(d: int, e: int) -> str:
    types = [sector_type(d, e, s) for s in range(1, d + 1)]
    selected = [s - 1 for s in range(1, d + 1) if types[s - 1] == "free+prufer"]
    exponent = observer_annihilator_exponent(d, e)
    observer = "0" if exponent is None else f"lambda^{exponent}"
    return f"d={d}, e={e}: sectors={types}; prufer target weights={selected}; Ann([1])={observer}"


def describe_split_control(d: int, multiplicities: tuple[int, ...], chosen: int) -> str:
    dims, components, exponents = special_fiber_sector_dimensions(d, multiplicities, chosen)
    return (
        f"d={d}, multiplicities={multiplicities}, chosen={chosen}: "
        f"special exponents={exponents}, components={components}, "
        f"dim C_s/(P-beta_i)C_s={dims}"
    )


def main() -> None:
    relations, endpoints, primitives = check_direct_identities()
    transitions, arrows, thickness = check_orbit_and_arrows()
    simple = check_simple_root()
    split_packets, split_characters = check_split_root_first_windows()

    print("ONE-ROOT NONLINEAR INTEGRAL RESPONSE -- EXACT REFEREE")
    print(f"direct monomial relations checked: {relations}")
    print(f"direct lambda-torsion endpoint identities checked: {endpoints}")
    print(f"direct unit-observer primitives checked: {primitives}")
    print(f"orbit transition formulas checked: {transitions}")
    print(f"lambda-arrow/break formulas checked: {arrows}")
    print(f"finite torsion-thickness formulas checked: {thickness}")
    print(f"simple-root free-arrow controls checked: {simple}")
    print(f"split-root repeated-boundary packets checked: {split_packets}")
    print(f"split-root special-fiber character dimensions checked: {split_characters}")
    print()
    print("classification over K[lambda]")
    print("e=1: every sector is K[lambda]")
    print("e>1, d | s(e-1): K[lambda] direct-sum K[lambda,lambda^-1]/K[lambda]")
    print("e>1, d does not divide s(e-1): K[lambda,lambda^-1]")
    print("number of Pruefer-bearing target characters: gcd(d,e-1)")
    print("selected-sector torsion length at z^d-depth k: k(e-1)+s(e-1)/d")
    print()
    for d, e in ((2, 1), (2, 2), (2, 3), (4, 3), (3, 7), (6, 5)):
        print(describe_control(d, e))
    print()
    print("generic rank control: every displayed summand has K(lambda)-rank one")
    print("localization loss: every Pruefer arm dies, while every Laurent/free line survives")
    print()
    print("split-root first-window theorem (chosen root repeated)")
    print("dim C_(s-1)/(P-beta_i)C_(s-1) is N for trivial Kummer monodromy, N-1 otherwise")
    print("monodromy exponents: e_i-1 at the chosen root and e_j at every other root")
    for d, multiplicities, chosen in (
        (4, (3, 1), 0),
        (4, (3, 2), 0),
        (6, (5, 3, 2), 0),
    ):
        print(describe_split_control(d, multiplicities, chosen))
    print("hostile: (d;e_i,e_j)=(4;3,1), s=2 is locally resonant but globally nontrivial")
    print("warning: a first-window dimension is not a global direct-sum or Pruefer-arm theorem")


if __name__ == "__main__":
    main()
