#!/usr/bin/env python3
"""Exact companion for THM-2552: flat deck gain and coherent-face hostile."""

from __future__ import annotations

from collections import Counter
from itertools import product


P = 13
FIELD = 547
CHECKS = 0


def require(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def prime_factors(n: int) -> list[int]:
    out: list[int] = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out


def primitive_root(p: int) -> int:
    factors = prime_factors(p - 1)
    for g in range(2, p):
        if all(pow(g, (p - 1) // ell, p) != 1 for ell in factors):
            return g
    raise AssertionError("no primitive root")


def face(t: int) -> list[tuple[int, int, int]]:
    return [(x, y, (t - x - y) % P) for x in range(P) for y in range(P)]


def marginal(atoms: list[tuple[int, int, int]], axes: tuple[int, ...]) -> Counter[tuple[int, ...]]:
    return Counter(tuple(atom[i] for i in axes) for atom in atoms)


def main() -> None:
    # Every triangular closed deck path has zero telescoping gain.  Affine
    # vertex gauges exhaust a nontrivial exact family of gauge controls.
    deck_checks = 0
    gauge_checks = 0
    for r, s, u in product(range(P), repeat=3):
        gain = ((s - r) + (u - s) + (r - u)) % P
        require(gain == 0, "pure deck triangle has nonzero gain")
        deck_checks += 1
        for a in range(P):
            # phi(v)=a*v; constants cancel automatically.
            gauged = ((s - r) + a * (s - r)
                       + (u - s) + a * (u - s)
                       + (r - u) + a * (r - u)) % P
            require(gauged == 0, "vertex gauge changed closed gain")
            gauge_checks += 1

    nonzero = set(range(1, P))
    two_wall = {(x + y) % P for x in nonzero for y in nonzero}
    require(two_wall == set(range(P)), "two free walls did not fill F13")
    seven_wall = {0}
    for _ in range(7):
        seven_wall = {(x + y) % P for x in seven_wall for y in nonzero}
    require(seven_wall == set(range(P)), "seven free walls did not fill F13")

    faces = {t: face(t) for t in (0, 1)}
    require(all(len(atoms) == P * P for atoms in faces.values()),
            "coherent face has wrong size")
    marginal_checks = 0
    for axes in ((0,), (1,), (2,), (0, 1), (0, 2), (1, 2)):
        m0 = marginal(faces[0], axes)
        m1 = marginal(faces[1], axes)
        require(m0 == m1, "proper marginals distinguish curvature")
        expected = P ** (2 - len(axes))
        require(set(m0.values()) == {expected}, "proper marginal is not uniform")
        marginal_checks += len(m0)
    require(all((x + y + z) % P == 0 for x, y, z in faces[0]),
            "zero face has nonzero curvature")
    require(all((x + y + z) % P == 1 for x, y, z in faces[1]),
            "unit face has wrong curvature")

    generator = primitive_root(FIELD)
    zeta = pow(generator, (FIELD - 1) // P, FIELD)
    fourier_checks = 0
    nonzero_face_modes = 0
    for t in (0, 1):
        atoms = faces[t]
        for alpha, beta, gamma in product(range(P), repeat=3):
            total = sum(pow(zeta,
                            (alpha * x + beta * y + gamma * z) % P,
                            FIELD)
                        for x, y, z in atoms) % FIELD
            expected = (P * P * pow(zeta, (gamma * t) % P, FIELD)) % FIELD if alpha == beta == gamma else 0
            require(total == expected, "coherent-face Fourier formula failed")
            if total:
                nonzero_face_modes += 1
            fourier_checks += 1
    require(nonzero_face_modes == 2 * P, "wrong diagonal Fourier support")

    print("THM-2552 exact common-base gain/coherent-face referee")
    print(f"deck_triangle_checks={deck_checks}")
    print(f"affine_gauge_triangle_checks={gauge_checks}")
    print(f"two_wall_sumset={len(two_wall)} seven_wall_sumset={len(seven_wall)}")
    print(f"face_atoms_each={len(faces[0])} proper_marginal_cells_checked={marginal_checks}")
    print(f"face_curvatures=0,1")
    print(f"fourier_checks={fourier_checks} nonzero_diagonal_modes={nonzero_face_modes} nonconstant_three_edge_modes={nonzero_face_modes - 2}")
    print(f"explicit_require_checks={CHECKS}")
    print("ALL CHECKS PASS")


if __name__ == "__main__":
    main()
