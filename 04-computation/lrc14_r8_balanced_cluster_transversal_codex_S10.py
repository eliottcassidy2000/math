#!/usr/bin/env python3
"""Exact certificate for THM-786's balanced-cluster transversal bound.

For the six companions C, a balanced cluster B is any nonempty subset with

    sum_(c in B) c^(-1) = -g^(-1)  (mod 7).

A nonnegative fractional transversal lambda satisfies sum_(c in B) lambda_c
>= 1 for every balanced cluster.  Its speed weight W=sum lambda_c*c and mass
Lambda=sum lambda_c replace the crude all-companion quantities in the wall-
grid density proof.  This script certifies both sides of the exact threshold:
at g=9 the old condition sum(C)<g fails but W<g, while the smallest admissible
profile g=8 has W>g and supplies a capacity-compatible dual obstruction.

Tournament Analysis / assumption challenge
-------------------------------------------
The theorem-bearing object is the bipartite incidence hypergraph between
balanced clusters and companion speeds.  Ordering clusters by cardinality,
speed sum, and lexicographic order gives a transitive telemetry tournament
(one Hamiltonian path, no cycles), but that order forgets the hitting
predicate.  The fractional-cover weights and speed labels must remain as a
sidecar; runner vertices or a bare cluster tournament do not preserve the
density certificate.
"""

from __future__ import annotations

from fractions import Fraction as F
from hashlib import sha256


MODULUS = 7
COMPANIONS = (1, 2, 3, 4, 5, 6)
EXPECTED_DIGEST = "b87a919e44f7af0903d0ac812185c5e39b400962bae03c539210c8578b4ab778"


def inverse_mod7(value: int) -> int:
    assert value % MODULUS
    return pow(value % MODULUS, -1, MODULUS)


def subset(mask: int) -> tuple[int, ...]:
    return tuple(
        speed for index, speed in enumerate(COMPANIONS) if (mask >> index) & 1
    )


def is_balanced_cluster(g: int, mask: int) -> bool:
    return (
        inverse_mod7(g)
        + sum(inverse_mod7(speed) for speed in subset(mask))
    ) % MODULUS == 0


def is_transversal(mask: int, clusters: tuple[int, ...]) -> bool:
    return all(mask & cluster for cluster in clusters)


def analyze(g: int):
    clusters = tuple(mask for mask in range(1, 1 << len(COMPANIONS))
                     if is_balanced_cluster(g, mask))
    transversals = tuple(mask for mask in range(1 << len(COMPANIONS))
                         if is_transversal(mask, clusters))
    ranked = sorted(
        (sum(subset(mask)), len(subset(mask)), subset(mask), mask)
        for mask in transversals
    )
    return clusters, transversals, ranked


def main() -> None:
    soft_g, fastest = 9, 10
    soft_clusters, soft_transversals, soft_ranked = analyze(soft_g)
    soft_weight, soft_size, soft_set, _ = soft_ranked[0]

    # Three disjoint balanced clusters give the matching fractional dual lower
    # bound 5+1+2=8: every fractional cover must pay at least one unit across
    # each, and their cheapest speeds have weights 5, 1, and 2 respectively.
    soft_dual_clusters = ((5,), (1, 4), (2, 6))
    assert all(
        tuple_ in tuple(subset(mask) for mask in soft_clusters)
        for tuple_ in soft_dual_clusters
    )
    soft_dual_lower_bound = sum(min(cluster) for cluster in soft_dual_clusters)

    assert len(soft_clusters) == 9
    assert len(soft_transversals) == 17
    assert (soft_weight, soft_size, soft_set) == (8, 3, (1, 2, 5))
    assert soft_dual_lower_bound == soft_weight
    assert sum(COMPANIONS) == 21 >= soft_g > soft_weight

    wall_bound = (
        F(soft_size, 1) - F(soft_weight, soft_g)
        + F(2 * soft_weight, fastest)
    ) / (1 - F(soft_weight, soft_g))
    assert wall_bound == F(167, 5)

    # The smallest valid six-companion profile is genuinely high-transversal.
    hard_g = 8
    hard_clusters, hard_transversals, hard_ranked = analyze(hard_g)
    hard_weight, hard_size, hard_set, _ = hard_ranked[0]
    hard_dual_clusters = ((1, 3), (2, 4), (6,))
    assert all(
        tuple_ in tuple(subset(mask) for mask in hard_clusters)
        for tuple_ in hard_dual_clusters
    )
    hard_dual_weights = (F(1), F(2), F(6))
    hard_dual_value = sum(hard_dual_weights)
    assert len(hard_clusters) == 9
    assert len(hard_transversals) == 17
    assert (hard_weight, hard_size, hard_set) == (9, 3, (1, 2, 6))
    assert hard_dual_value == hard_weight > hard_g

    # Scale the dual packing to total mass g.  Dividing by g gives a
    # probability distribution on balanced clusters whose companion marginals
    # fit below the individual grid capacities c/g.
    scale = F(hard_g, hard_dual_value)
    hard_probabilities = tuple(weight * scale / hard_g
                               for weight in hard_dual_weights)
    assert sum(hard_probabilities) == 1
    for companion in COMPANIONS:
        marginal = sum(
            probability
            for cluster, probability in zip(hard_dual_clusters, hard_probabilities)
            if companion in cluster
        )
        assert marginal <= F(companion, hard_g)

    digest = sha256()
    for g, clusters, ranked in (
        (soft_g, soft_clusters, soft_ranked),
        (hard_g, hard_clusters, hard_ranked),
    ):
        digest.update(f"G|{g}\n".encode())
        for mask in clusters:
            digest.update(("B|" + ",".join(map(str, subset(mask))) + "\n").encode())
        for total, size, speeds, _ in ranked:
            digest.update(
                (f"T|{total}|{size}|" + ",".join(map(str, speeds)) + "\n").encode()
            )
    hexdigest = digest.hexdigest()

    print("THM-786 balanced-cluster transversal certificate")
    print("strict-extension profile")
    print(f"g={soft_g} f={fastest} companions={COMPANIONS} sum_companions={sum(COMPANIONS)}")
    print(
        "inverse_residues="
        + str(tuple(inverse_mod7(speed) for speed in (soft_g,) + COMPANIONS))
        + " order=(g,C)"
    )
    print(f"balanced_cluster_count={len(soft_clusters)}")
    print("balanced_clusters=" + str(tuple(subset(mask) for mask in soft_clusters)))
    print(f"transversal_count={len(soft_transversals)}")
    print(
        f"minimum_speed_weight={soft_weight} minimum_transversal={soft_set} "
        f"mass={soft_size}"
    )
    print(
        f"fractional_optimality_dual={soft_dual_clusters} "
        f"lower_bound={soft_dual_lower_bound}"
    )
    print(f"strict_extension=sum(C)={sum(COMPANIONS)}>=g={soft_g}>W={soft_weight}")
    print(f"for_f={fastest}: M<{wall_bound} hence integer_M<={wall_bound.numerator // wall_bound.denominator}")
    print()
    print("high-transversal residual profile")
    print(f"g={hard_g} companions={COMPANIONS} sum_companions={sum(COMPANIONS)}")
    print(f"balanced_cluster_count={len(hard_clusters)}")
    print("balanced_clusters=" + str(tuple(subset(mask) for mask in hard_clusters)))
    print(
        f"minimum_speed_weight={hard_weight} minimum_transversal={hard_set} "
        f"mass={hard_size}"
    )
    print(
        f"fractional_optimality_dual={hard_dual_clusters} "
        f"weights={tuple(int(weight) for weight in hard_dual_weights)} "
        f"value={hard_dual_value}"
    )
    print(
        "capacity_distribution=("
        + ",".join(map(str, hard_probabilities))
        + ") "
        "on_clusters=((1,3),(2,4),(6,))"
    )
    print(f"genuine_residual=W={hard_weight}>g={hard_g}")
    print(f"incidence_digest={hexdigest}")
    print()
    print("Tournament Analysis / assumption challenge:")
    print("  vertices: nine balanced companion clusters in each residue profile")
    print("  gauge: (cardinality, speed sum, lexicographic order)")
    print("  score_histogram={0:1,...,8:1} cycles=0 sccs=9 hp_count=1")
    print("  information lost by tournament: cluster-speed incidence and cover weight")
    print("  exact carrier: labelled cluster hypergraph + fractional transversal")

    if EXPECTED_DIGEST != "TO_BE_FILLED":
        assert hexdigest == EXPECTED_DIGEST
    print("PASS")


if __name__ == "__main__":
    main()
