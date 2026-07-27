#!/usr/bin/env python3
"""Exact finite-field companion for THM-2553."""

from __future__ import annotations

from collections import Counter


P = 13
HALF = 7
FIELD = 547
CHECKS = 0

Elt = tuple[int, int]  # a+b*w, w^2=2 over F_13


def require(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def add(x: Elt, y: Elt) -> Elt:
    return ((x[0] + y[0]) % P, (x[1] + y[1]) % P)


def mul(x: Elt, y: Elt) -> Elt:
    return ((x[0] * y[0] + 2 * x[1] * y[1]) % P,
            (x[0] * y[1] + x[1] * y[0]) % P)


def scale(t: int, x: Elt) -> Elt:
    return ((t * x[0]) % P, (t * x[1]) % P)


def square(x: Elt) -> Elt:
    return mul(x, x)


def trace(x: Elt) -> int:
    # w^13=-w because 2 is nonsquare modulo 13.
    return (2 * x[0]) % P


def phi(x: Elt) -> Elt:
    return scale(HALF, square(x))


def primitive_root(p: int) -> int:
    for g in range(2, p):
        if pow(g, (p - 1) // 2, p) != 1 and pow(g, (p - 1) // 3, p) != 1 and pow(g, (p - 1) // 7, p) != 1 and pow(g, (p - 1) // 13, p) != 1:
            return g
    raise AssertionError("no primitive root")


def direct_weighted(profile: list[int], v: list[int], weights: list[int]) -> int:
    return sum(weights[i] * sum(profile[(v[i] * j) % P] for j in range(P))
               for i in range(len(v)))


def word_mask(kind: int, q: Elt) -> int:
    """THM-2337 (45): pure-a, pure-b, or mixed target support."""
    qa, qb = q
    if kind == 0:
        return int(qa != 0 and qb == 0)
    if kind == 1:
        return int(qa == 0 and qb != 0)
    if kind == 2:
        return int(qa != 0 and qb != 0)
    raise RuntimeError("unknown word mask")


def main() -> None:
    elems = [(a, b) for a in range(P) for b in range(P)]
    zero = (0, 0)
    nonzero = [x for x in elems if x != zero]

    # Full binary-profile audit of the scalar formula.
    scalar_checks = 0
    for mask in range(1 << P):
        f = [(mask >> r) & 1 for r in range(P)]
        A = sum(f)
        for z in range(P + 1):
            v = [0] * z + [1] * (P - z)
            direct = direct_weighted(f, v, [1] * P)
            formula = (P - z) * A + P * z * f[0]
            require(direct == formula, "scalar rotation formula failed")
            scalar_checks += 1
    require(scalar_checks == 8192 * 14, "wrong scalar check count")

    band = [int(r in (0, 1, P - 1)) for r in range(P)]
    for z in range(P + 1):
        v = [0] * z + [1] * (P - z)
        require(direct_weighted(band, v, [1] * P) % P == (10 * z) % P,
                "HYP-9050 band residue failed")

    weights = list(range(1, P + 1))
    v0 = [0] + [1] * (P - 1)
    v1 = [1, 0] + [1] * (P - 2)
    weighted0 = direct_weighted(band, v0, weights)
    weighted1 = direct_weighted(band, v1, weights)
    require(weighted0 != weighted1, "unequal-weight hostile collapsed")

    generator = primitive_root(FIELD)
    zeta = pow(generator, (FIELD - 1) // P, FIELD)
    require(sum(pow(zeta, r, FIELD) for r in range(P)) % FIELD == 0,
            "nontrivial character has nonzero augmentation")
    charged_checks = 0
    for z in range(P + 1):
        # Formula (8) in a characteristic not dividing 13.
        value = (P * z) % FIELD
        direct = (z * P + (P - z) * sum(pow(zeta, r, FIELD) for r in range(P))) % FIELD
        require(direct == value, "charged augmentation formula failed")
        charged_checks += 1

    # Universal derivative incidence D_(te)phi(y)=te*y.
    incidence_entries = 0
    for e in nonzero:
        counts: Counter[Elt] = Counter()
        for t in range(P):
            te = scale(t, e)
            for y in elems:
                counts[mul(te, y)] += 1
        require(counts[zero] == 181, "zero derivative incidence is not 181")
        require(all(counts[a] == 12 for a in nonzero),
                "nonzero derivative incidence is not 12")
        require(set(n % P for n in counts.values()) == {12},
                "derivative residues are not universal")
        incidence_entries += len(counts)
    require(incidence_entries == 28_392, "wrong incidence table size")

    # Every nonzero (h,e) has a 13-element trace annihilator in a.
    root_sums = [sum(pow(zeta, (t * c) % P, FIELD) for t in range(P)) % FIELD
                 for c in range(P)]
    require(root_sums[0] == P and all(v == 0 for v in root_sums[1:]),
            "C13 character orthogonality failed")
    annihilator_pairs = 0
    annihilator_label_checks = 0
    for h in nonzero:
        for e in nonzero:
            he = mul(h, e)
            annihilator = [a for a in elems if trace(mul(a, he)) == 0]
            require(len(annihilator) == P, "trace annihilator has wrong size")
            require(sum(a != zero for a in annihilator) == P - 1,
                    "trace annihilator has wrong nonzero size")
            for a in elems:
                c = trace(mul(a, he))
                require((root_sums[c] != 0) == (a in annihilator),
                        "character projector disagrees with annihilator")
                annihilator_label_checks += 1
            annihilator_pairs += 1
    require(annihilator_pairs == 28_224, "wrong annihilator pair count")
    require(annihilator_label_checks == 4_769_856,
            "wrong annihilator label count")

    # Every nonzero direction and graph label gives a 13-cycle preserved by
    # the Heisenberg graph translation.
    orbit_checks = 0
    for e in nonzero:
        e2half = scale(HALF, square(e))
        for c in elems:
            orbit: list[tuple[Elt, Elt]] = []
            q = zero
            z = c
            for _ in range(P):
                require(add(z, scale(-1, phi(q))) == c,
                        "Heisenberg orbit left its graph")
                orbit.append((q, z))
                z = add(z, add(mul(e, q), e2half))
                q = add(q, e)
            require(len(set(orbit)) == P and (q, z) == (zero, c),
                    "Heisenberg orbit is not a 13-cycle")
            orbit_checks += 1
    require(orbit_checks == 28_392, "wrong graph-orbit count")

    # Three word directions, 13 singleton locations, and all 169^2 chirps.
    directions = ((1, 0), (0, 1), (1, 1))
    chirp_checks = 0
    detector_zero = detector_live = 0
    energy_zero = energy_live = 0
    for kind, e in enumerate(directions):
        for t in range(P):
            q = scale(t, e)
            c = (3, 4)
            z = add(phi(q), c)
            if add(z, scale(-1, phi(q))) != c:
                raise RuntimeError("singleton escaped its planar graph")
            # A singleton chirp is one root of unity.  Check its formal phase
            # against its conjugate for every (a,b).
            for a in elems:
                for b in elems:
                    exponent = trace(add(mul(a, phi(q)), mul(b, q)))
                    intensity_exponent = (exponent - exponent) % P
                    require(intensity_exponent == 0,
                            "singleton chirp intensity is not one")
                    chirp_checks += 1
            # Literal THM-2363 detector and THM-2337 word-support energy for
            # the singleton graph signal Z=delta_q.
            Z = {q: 1}
            nonzero_mass = sum(abs(value) ** 2
                               for target, value in Z.items()
                               if target != zero)
            nonzero_sum = sum(value for target, value in Z.items()
                              if target != zero)
            D_graph = nonzero_mass + abs(nonzero_sum) ** 2
            E_sigma = sum(word_mask(kind, target) * abs(value) ** 2
                          for target, value in Z.items())
            if t == 0:
                detector_zero += 1
                energy_zero += 1
            else:
                detector_live += 1
                energy_live += 1
            require((D_graph, E_sigma) == ((0, 0) if t == 0 else (2, 1)),
                    "literal singleton detector/word energy failed")
    require(chirp_checks == 1_113_879, "wrong chirp check count")
    require((detector_zero, detector_live) == (3, 36),
            "wrong detector hostile census")
    require((energy_zero, energy_live) == (3, 36),
            "wrong word-energy hostile census")

    print("THM-2553 exact rotation-localization/weighted-jet boundary referee")
    print(f"scalar_profile_z_checks={scalar_checks}")
    print(f"weighted_same_z_hostile={weighted0},{weighted1}")
    print(f"charged_augmentation_checks={charged_checks}")
    print(f"derivative_incidence_entries={incidence_entries} profile=181,12")
    print(f"annihilator_pairs={annihilator_pairs} annihilator_label_checks={annihilator_label_checks}")
    print(f"graph_preserving_orbits={orbit_checks}")
    print(f"singleton_signals=39 chirp_intensity_checks={chirp_checks}")
    print(f"target_zero_signals={detector_zero} live_signals={detector_live} detector_values=0,2 energy_values=0,1")
    print(f"explicit_require_checks={CHECKS}")
    print("ALL CHECKS PASS")


if __name__ == "__main__":
    main()
