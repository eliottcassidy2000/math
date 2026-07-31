#!/usr/bin/env python3
"""Exact finite controls for THM-2611's principal-C13 holotopy theorem."""

from itertools import product


P = 13
Q = 7
N = P * Q


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def iota(r: int) -> int:
    return (14 * r) % N


def delta(c: int) -> int:
    return (78 + iota(c)) % N


def c91_exact_sequence() -> tuple[int, int, int]:
    kernel = {z for z in range(N) if z % Q == 0}
    image = {iota(r) for r in range(P)}
    require(kernel == image and len(kernel) == P,
            "the C13 embedding is not the full C91->C7 kernel")

    fibres = []
    for b in range(Q):
        fibre = {z for z in range(N) if z % Q == b}
        require(len(fibre) == P, "a quotient fibre is not principal C13")
        for x in fibre:
            require({(x + iota(r)) % N for r in range(P)} == fibre,
                    "the kernel action is not transitive on a fibre")
        fibres.append(fibre)

    edge_lifts = {delta(c) for c in range(P)}
    require(edge_lifts == fibres[1] and len(edge_lifts) == P,
            "positive-clock edge lifts are not one C13 torsor")
    require(delta(0) == 78 and delta(7) == 85,
            "constant-six rail lifts changed")
    require({(P * delta(c)) % N for c in range(P)} == {13},
            "Koopman multiplication did not collapse all edge lifts")
    require({z for z in range(N) if (P * z) % N == 0} == kernel,
            "multiplication by thirteen has the wrong kernel")
    return len(fibres), len(edge_lifts), len(kernel)


def equivariant_map_counts() -> tuple[int, int]:
    fixed = []
    semilinear = []
    for kappa in range(P):
        for c in range(P):
            phi = tuple((kappa * r + c) % P for r in range(P))
            if all(phi[(r + 1) % P] == (phi[r] + 1) % P
                   for r in range(P)):
                fixed.append((kappa, c))
            if kappa != 0 and len(set(phi)) == P:
                require(all(phi[(r + 1) % P] ==
                            (phi[r] + kappa) % P for r in range(P)),
                        "an affine semilinear map lost covariance")
                semilinear.append((kappa, c))
    require(fixed == [(1, c) for c in range(P)],
            "fixed-action identification count changed")
    require(len(semilinear) == P * (P - 1),
            "semilinear affine identification count changed")

    # A K-equivariant injection of one regular orbit has a free image orbit.
    # Exhaust all translation-stable subsets of the regular 13-set: only the
    # empty set and full set exist, which is the sharp finite equality case.
    invariant_subsets = []
    for mask in range(1 << P):
        subset = {r for r in range(P) if (mask >> r) & 1}
        if all({(x + k) % P for x in subset} == subset for k in range(P)):
            invariant_subsets.append(len(subset))
    require(invariant_subsets == [0, P],
            "the regular C13 action acquired a smaller nonempty orbit")
    return len(fixed), len(semilinear)


def binary_cycle_sections() -> tuple[int, int, int]:
    pairs = 0
    closed = 0
    sections = 0
    for cs in product((0, 7), repeat=Q):
        for a in range(P):
            pairs += 1
            hol = (Q * a + sum(cs)) % P
            # Directly solve h_i-h_(i-1)=-(a+c_i) for every base origin.
            solutions = 0
            for base in range(P):
                h = [None] * Q
                previous = base
                for ell in range(Q - 1):
                    h[ell] = (previous - a - cs[ell]) % P
                    previous = h[ell]
                h[-1] = base
                if all((a + cs[ell] + h[ell] - h[(ell - 1) % Q]) % P == 0
                       for ell in range(Q)):
                    solutions += 1
            require(solutions == (P if hol == 0 else 0),
                    "cycle section count disagrees with holonomy")

            lifted = (sum(delta(c) for c in cs) + Q * iota(a)) % N
            require((lifted == 0) == (hol == 0),
                    "C91 lifted closure disagrees with C13 holonomy")
            if hol == 0:
                closed += 1
                sections += solutions
    require(pairs == (2**Q) * P, "binary cycle census changed")
    require(closed == 2**Q and sections == (2**Q) * P,
            "binary cycle equality census changed")
    return pairs, closed, sections


def thm2601_intertwiners() -> tuple[int, tuple[tuple[int, int], ...]]:
    # THM-2601's exact target values t_q and successor S(q)=t_(q+1) in
    # scalar coordinates.  S is printed here as a map on scalar values.
    target_values = (9, 3, 10, 1, 8, 0, 7, 2, 12, 11, 4, 6, 5)
    position = {value: q for q, value in enumerate(target_values)}
    require(len(position) == P, "THM-2601 scalar coordinate is not bijective")
    successor = tuple(target_values[(position[c] + 1) % P]
                      for c in range(P))
    require(successor == (7, 8, 12, 10, 6, 9, 5, 2, 0, 3, 1, 4, 11),
            "THM-2601 successor changed")

    maps = []
    references = []
    for c in range(P):
        values = [c]
        for _ in range(1, P):
            values.append(successor[values[-1]])
        phi = tuple(values)
        require(len(set(phi)) == P,
                "an equivariant root-to-target map is not bijective")
        require(all(phi[(r + 1) % P] == successor[phi[r]]
                    for r in range(P)),
                "THM-2601 intertwiner lost equivariance")
        maps.append(phi)
        references.append((c, position[c]))
    require(len(set(maps)) == P, "THM-2601 does not have thirteen splittings")
    expected = ((0, 5), (1, 3), (2, 7), (3, 1), (4, 10), (5, 12),
                (6, 11), (7, 6), (8, 4), (9, 0), (10, 2), (11, 9),
                (12, 8))
    require(tuple(references) == expected,
            "THM-2601 reference-origin pairs changed")
    return len(maps), tuple(references)


def main() -> None:
    fibres, edge_lifts, kernel = c91_exact_sequence()
    fixed, semilinear = equivariant_map_counts()
    pairs, closed, sections = binary_cycle_sections()
    maps, references = thm2601_intertwiners()

    print("THM-2611 exact principal-C13 bibundle controls")
    print(f"c91_fibres={fibres} fibre_size={kernel} positive_edge_lifts={edge_lifts}")
    print("crt_idempotents=(14,78) rail_lifts=(78,85) koopman_common_image=13")
    print(f"fixed_action_identifications={fixed} semilinear_affine_identifications={semilinear}")
    print("regular_c13_invariant_subset_sizes=(0,13) minimum_faithful_sidecar=13")
    print(f"binary_cycle_marker_pairs={pairs} closed_pairs={closed} parallel_sections={sections}")
    print(f"thm2601_equivariant_maps={maps} references={references}")
    print("verdict=PASS: zero kernel holonomy iff thirteen sections; physical bibundle remains absent")


if __name__ == "__main__":
    main()
