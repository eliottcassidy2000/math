#!/usr/bin/env python3
"""Exact companion for THM-3046.

The script checks the labelled quartic/resolvent difference identities,
the valuation-clutch map on all residue classes modulo six, the rational
projector integrality criterion, S4 equivariance, and sharp 5-adic controls.
All truth-bearing checks use explicit exceptions and survive ``python -O``.
"""

from fractions import Fraction
from itertools import permutations, product
from hashlib import sha256

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


VERTICES = range(4)
EDGES = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
EDGE_INDEX = {edge: i for i, edge in enumerate(EDGES)}
MATCHINGS = (
    ((0, 1), (2, 3)),
    ((0, 2), (1, 3)),
    ((0, 3), (1, 2)),
)
MATCHING_INDEX = {
    frozenset((tuple(sorted(e)), tuple(sorted(f)))): i
    for i, (e, f) in enumerate(MATCHINGS)
}


def matching_sums(x):
    return tuple(
        x[EDGE_INDEX[tuple(sorted(e))]] + x[EDGE_INDEX[tuple(sorted(f))]]
        for e, f in MATCHINGS
    )


def clutch(x):
    s = matching_sums(x)
    return tuple(value % 2 for value in s), sum(x) % 3


def projectors(x):
    total = sum(x)
    p0 = tuple(Fraction(total, 6) for _ in EDGES)
    p31 = [Fraction(0) for _ in EDGES]
    pplus = [Fraction(0) for _ in EDGES]
    for e, f in MATCHINGS:
        i = EDGE_INDEX[tuple(sorted(e))]
        j = EDGE_INDEX[tuple(sorted(f))]
        p31[i] = Fraction(x[i] - x[j], 2)
        p31[j] = Fraction(x[j] - x[i], 2)
        pplus[i] = pplus[j] = Fraction(x[i] + x[j], 2)
    p22 = tuple(pplus[i] - p0[i] for i in range(6))
    return p0, p22, tuple(p31)


def integral(values):
    return all(value.denominator == 1 for value in values)


def permute_edge_vector(x, perm):
    out = [None] * 6
    for edge, value in zip(EDGES, x):
        image = tuple(sorted((perm[edge[0]], perm[edge[1]])))
        out[EDGE_INDEX[image]] = value
    require(all(value is not None for value in out), "edge permutation incomplete")
    return tuple(out)


def matching_permutation(perm):
    out = []
    for e, f in MATCHINGS:
        image = frozenset(
            (
                tuple(sorted((perm[e[0]], perm[e[1]]))),
                tuple(sorted((perm[f[0]], perm[f[1]]))),
            )
        )
        out.append(MATCHING_INDEX[image])
    return tuple(out)


def vp(nonzero_integer, prime):
    require(nonzero_integer != 0, "valuation requested at zero")
    value = abs(nonzero_integer)
    order = 0
    while value % prime == 0:
        value //= prime
        order += 1
    return order


def root_valuation_vector(roots, prime):
    return tuple(vp(roots[i] - roots[j], prime) for i, j in EDGES)


def base_roots_for_t_a(t, a, prime=5):
    """Realize matching sums (t+a,t,t) with four distinct rational roots."""
    require(t >= 0 and a >= 0, "negative cluster depth")
    if t == 0 and a == 0:
        return (0, 1, 2, 3)
    if t == 0:
        return (0, prime**a, 1, 2)
    if a == 0:
        return (0, prime**t, 2 * prime**t, 1)
    return (0, prime ** (t + a), prime**t, 1)


def realize_clutch_class(kappa, tau, prime=5):
    """Return honest roots realizing a prescribed F2^3 x F3 class."""
    special = next(
        i for i in range(3) if kappa[(i + 1) % 3] == kappa[(i + 2) % 3]
    )
    common = kappa[(special + 1) % 3]
    t = common
    parity_a = (kappa[special] - common) % 2
    a = next(value for value in range(6) if value % 2 == parity_a and value % 3 == tau)
    base = base_roots_for_t_a(t, a, prime)
    for perm in permutations(VERTICES):
        roots = tuple(base[perm[i]] for i in VERTICES)
        if clutch(root_valuation_vector(roots, prime)) == (kappa, tau):
            return roots
    raise RuntimeError("failed to place special matching coordinate")


# Root-level polynomial identities.
z = sp.symbols("z0:4")
u = (
    z[0] * z[1] + z[2] * z[3],
    z[0] * z[2] + z[1] * z[3],
    z[0] * z[3] + z[1] * z[2],
)
identities = (
    u[1] - u[2] - (z[0] - z[1]) * (z[2] - z[3]),
    u[0] - u[2] - (z[0] - z[2]) * (z[1] - z[3]),
    u[0] - u[1] - (z[0] - z[3]) * (z[1] - z[2]),
)
for number, identity in enumerate(identities, 1):
    require(sp.expand(identity) == 0, f"resolvent difference identity {number}")

plucker_terms = (
    (z[0] - z[1]) * (z[2] - z[3]),
    -(z[0] - z[2]) * (z[1] - z[3]),
    (z[0] - z[3]) * (z[1] - z[2]),
)
require(sp.expand(sum(plucker_terms)) == 0, "Pluecker relation")

quartic_vandermonde = sp.prod(z[i] - z[j] for i in range(4) for j in range(i + 1, 4))
resolvent_vandermonde = (u[0] - u[1]) * (u[0] - u[2]) * (u[1] - u[2])
require(
    sp.expand(resolvent_vandermonde**2 - quartic_vandermonde**2) == 0,
    "quartic/resolvent discriminant identity",
)

# The residue universe is exactly the denominator-six obstruction universe.
fibre_counts = {}
projector_good = 0
for x in product(range(6), repeat=6):
    kappa, tau = clutch(x)
    fibre_counts[(kappa, tau)] = fibre_counts.get((kappa, tau), 0) + 1
    s = matching_sums(x)
    v_disc = 2 * sum(x)
    require(tau == sum(s) % 3, "ternary matching-sum mismatch")
    require(tau == (v_disc // 2) % 3, "half-discriminant mismatch")
    p0, p22, p31 = projectors(x)
    all_integral = integral(p0) and integral(p22) and integral(p31)
    expected = kappa == (0, 0, 0) and tau == 0
    require(all_integral == expected, "projector integrality mismatch")
    projector_good += int(all_integral)

require(len(fibre_counts) == 24, "clutch map is not onto all 24 classes")
require(set(fibre_counts.values()) == {6**6 // 24}, "clutch fibres are not uniform")
require(projector_good == 6**6 // 24, "wrong integral-projector census")

# The tropical Pluecker restriction still realizes all 24 clutch classes.
root_realizations = {}
for kappa in product(range(2), repeat=3):
    for tau in range(3):
        roots = realize_clutch_class(kappa, tau)
        x = root_valuation_vector(roots, 5)
        s = matching_sums(x)
        require(clutch(x) == (kappa, tau), "root clutch realization mismatch")
        require(s.count(min(s)) >= 2, "tropical Pluecker minimum is not repeated")
        root_realizations[(kappa, tau)] = (roots, s)
require(len(root_realizations) == 24, "not all clutch classes are root-realized")

# S4 relabelling permutes the three binary matching coordinates and fixes tau.
for perm in permutations(VERTICES):
    sigma = matching_permutation(perm)
    for basis_index in range(6):
        x = tuple(int(i == basis_index) for i in range(6))
        kappa, tau = clutch(x)
        image_kappa, image_tau = clutch(permute_edge_vector(x, perm))
        expected_kappa = [None] * 3
        for i, value in enumerate(kappa):
            expected_kappa[sigma[i]] = value
        require(image_kappa == tuple(expected_kappa), "S4 matching equivariance failure")
        require(image_tau == tau, "S4 ternary invariance failure")

# Sharp root-realizable controls over Q_5.
two_pairs_depth_1 = root_valuation_vector((0, 5, 1, 6), 5)
two_pairs_depth_2 = root_valuation_vector((0, 25, 1, 26), 5)
triple_refinement = root_valuation_vector((0, 25, 5, 1), 5)

require(matching_sums(two_pairs_depth_1) == (2, 0, 0), "depth-one control")
require(matching_sums(two_pairs_depth_2) == (4, 0, 0), "depth-two control")
require(matching_sums(triple_refinement) == (2, 1, 1), "triple control")
require(
    2 * sum(two_pairs_depth_2) == 2 * sum(triple_refinement) == 8,
    "same-discriminant control",
)
require(
    clutch(two_pairs_depth_2)[0] != clutch(triple_refinement)[0],
    "discriminant unexpectedly recovers binary clutch",
)
require(
    clutch(two_pairs_depth_1)[0] == clutch(two_pairs_depth_2)[0]
    and clutch(two_pairs_depth_1)[1] != clutch(two_pairs_depth_2)[1],
    "binary clutch unexpectedly recovers ternary clutch",
)

semantic_payload = repr(
    (
        EDGES,
        MATCHINGS,
        sorted(fibre_counts.items()),
        two_pairs_depth_1,
        two_pairs_depth_2,
        triple_refinement,
        sorted(root_realizations.items()),
    )
).encode("ascii")

print("THM-3046 QUARTIC RESOLVENT VALUATION CLUTCH")
print("root_difference_identities=3")
print("discriminant_identity=quartic_equals_resolvent")
print("tropical_pluecker_identity=PASS")
print("residue_vectors=46656")
print("clutch_classes=24")
print("uniform_fibre_size=1944")
print(f"integral_projector_classes={projector_good}")
print("root_realized_clutch_classes=24")
print("root_realization_prime=5")
print("s4_permutations=24")
print("s4_basis_equivariance_checks=144")
print(f"depth1_matching_sums={matching_sums(two_pairs_depth_1)}")
print(f"depth2_matching_sums={matching_sums(two_pairs_depth_2)}")
print(f"triple_matching_sums={matching_sums(triple_refinement)}")
print("same_discriminant_binary_hostile=PASS")
print("same_binary_ternary_hostile=PASS")
print(f"semantic_sha256={sha256(semantic_payload).hexdigest()}")
print("ALL_EXACT_CHECKS_PASS")
