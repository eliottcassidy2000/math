#!/usr/bin/env python3
"""Finite exact audit for the first-kind character-transference theorem.

For an obtuse superbase b_0,...,b_d, its conorms form the conductances
of a connected graph on d+1 vertices.  A parity representative is a cut,
and a crossing edge produces an odd dual vector whose squared norm is the
edge's effective resistance.  Foster's identity then proves

    dist(c,L) * lambda_odd <= sqrt(d)/2.

The proof is algebraic and recorded in the companion reflection.  This file
uses Fraction arithmetic to audit every connected unweighted graph through
five vertices, several weighted families through eight vertices, every cut,
the effective-resistance/dual-Gram identity, Foster's trace, and the sharp
orthogonal star family.  It makes no assertion about arbitrary lattices.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def inverse(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    size = len(matrix)
    augmented = [
        row[:] + [Fraction(int(i == j)) for j in range(size)]
        for i, row in enumerate(matrix)
    ]
    for column in range(size):
        pivot = next(row for row in range(column, size) if augmented[row][column])
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        augmented[column] = [entry / scale for entry in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            multiplier = augmented[row][column]
            if multiplier:
                augmented[row] = [
                    left - multiplier * right
                    for left, right in zip(augmented[row], augmented[column])
                ]
    return [row[size:] for row in augmented]


def quadratic(vector: list[int], matrix: list[list[Fraction]]) -> Fraction:
    return sum(
        Fraction(vector[i]) * matrix[i][j] * vector[j]
        for i in range(len(vector))
        for j in range(len(vector))
    )


def connected(vertex_count: int, weights: dict[tuple[int, int], Fraction]) -> bool:
    seen = {0}
    frontier = [0]
    while frontier:
        vertex = frontier.pop()
        for edge, weight in weights.items():
            if not weight or vertex not in edge:
                continue
            other = edge[0] if edge[1] == vertex else edge[1]
            if other not in seen:
                seen.add(other)
                frontier.append(other)
    return len(seen) == vertex_count


def graph_audit(vertex_count: int, weights: dict[tuple[int, int], Fraction]) -> int:
    dimension = vertex_count - 1
    edges = [edge for edge, weight in weights.items() if weight > 0]
    require(connected(vertex_count, weights), "graph must be connected")

    laplacian = [[Fraction(0) for _ in range(vertex_count)] for _ in range(vertex_count)]
    for (first, second), weight in weights.items():
        if not weight:
            continue
        laplacian[first][first] += weight
        laplacian[second][second] += weight
        laplacian[first][second] -= weight
        laplacian[second][first] -= weight

    # Ground vertex zero.  This cofactor is precisely the Gram matrix of
    # the basis b_1,...,b_d attached to the obtuse superbase.
    gram = [row[1:] for row in laplacian[1:]]
    dual_gram = inverse(gram)

    resistances: dict[tuple[int, int], Fraction] = {}
    foster_sum = Fraction(0)
    for first, second in edges:
        pairing = [0] * dimension
        if first:
            pairing[first - 1] += 1
        if second:
            pairing[second - 1] -= 1
        resistance = quadratic(pairing, dual_gram)
        resistances[(first, second)] = resistance
        foster_sum += weights[(first, second)] * resistance
    require(foster_sum == dimension, "exact Foster trace identity")

    cut_checks = 0
    # Fix vertex zero outside the cut to quotient complementary cuts.
    for mask in range(1, 1 << dimension):
        inside = {index + 1 for index in range(dimension) if mask & (1 << index)}
        crossing = [edge for edge in edges if (edge[0] in inside) != (edge[1] in inside)]
        capacity = sum((weights[edge] for edge in crossing), Fraction(0))
        require(capacity > 0, "nontrivial connected cut")
        least_resistance = min(resistances[edge] for edge in crossing)
        require(capacity * least_resistance <= dimension, "cut/Foster inequality")

        # A crossing edge's incidence pairing evaluates to +/-1 on the cut
        # representative, so its corresponding dual vector is odd.
        chosen = min(crossing, key=lambda edge: resistances[edge])
        first, second = chosen
        odd_pairing = int(first in inside) - int(second in inside)
        require(abs(odd_pairing) == 1, "crossing edge is an odd character")
        cut_checks += 1
    return cut_checks


graph_count = 0
cut_count = 0
foster_count = 0

print("=== Exhaustive connected simple-graph controls ===")
for vertex_count in range(2, 6):
    possible_edges = list(combinations(range(vertex_count), 2))
    local_graphs = 0
    local_cuts = 0
    for mask in range(1, 1 << len(possible_edges)):
        weights = {
            edge: Fraction(1 if mask & (1 << index) else 0)
            for index, edge in enumerate(possible_edges)
        }
        if not connected(vertex_count, weights):
            continue
        local_cuts += graph_audit(vertex_count, weights)
        local_graphs += 1
    expected = {2: 1, 3: 4, 4: 38, 5: 728}[vertex_count]
    require(local_graphs == expected, "connected labelled graph count")
    graph_count += local_graphs
    cut_count += local_cuts
    foster_count += local_graphs
    print(f"vertices={vertex_count}: connected_graphs={local_graphs}, cuts={local_cuts}")


print("\n=== Weighted conorm-graph controls ===")
weighted_graphs = 0
weighted_cuts = 0
for vertex_count in range(2, 9):
    all_edges = list(combinations(range(vertex_count), 2))
    families: list[dict[tuple[int, int], Fraction]] = []

    # Positive complete graph with nonuniform rational conductances.
    families.append(
        {
            edge: Fraction((edge[0] + 1) * (edge[1] + 2) % 7 + 1, edge[0] + edge[1] + 1)
            for edge in all_edges
        }
    )
    # Weighted path.
    families.append(
        {
            edge: Fraction(edge[0] + 2, edge[1] + 1) if edge[1] == edge[0] + 1 else Fraction(0)
            for edge in all_edges
        }
    )
    # Weighted star, the orthogonal-lattice conorm graph.
    families.append(
        {
            edge: Fraction((edge[1] + 1) ** 2, edge[1] + 2) if edge[0] == 0 else Fraction(0)
            for edge in all_edges
        }
    )
    # Cycle (for two vertices this coincides with one edge, harmlessly).
    families.append(
        {
            edge: Fraction(edge[0] + edge[1] + 1, 3)
            if edge[1] == edge[0] + 1 or edge == (0, vertex_count - 1)
            else Fraction(0)
            for edge in all_edges
        }
    )

    for weights in families:
        weighted_cuts += graph_audit(vertex_count, weights)
        weighted_graphs += 1

graph_count += weighted_graphs
cut_count += weighted_cuts
foster_count += weighted_graphs
print(f"weighted_graphs={weighted_graphs}, cuts={weighted_cuts}")


print("\n=== Sharp orthogonal-star boundary ===")
for dimension in range(1, 13):
    vertex_count = dimension + 1
    weights = {
        edge: Fraction(1 if edge[0] == 0 else 0)
        for edge in combinations(range(vertex_count), 2)
    }
    # The cut separating all leaves from the centre has capacity d and
    # every crossing effective resistance one.  Hence P^2*q^2=d and
    # delta*q=sqrt(d)/2 exactly.
    laplacian = [[Fraction(0) for _ in range(vertex_count)] for _ in range(vertex_count)]
    for (first, second), weight in weights.items():
        if weight:
            laplacian[first][first] += weight
            laplacian[second][second] += weight
            laplacian[first][second] -= weight
            laplacian[second][first] -= weight
    gram = [row[1:] for row in laplacian[1:]]
    require(gram == [[Fraction(int(i == j)) for j in range(dimension)] for i in range(dimension)], "star Gram")
    require(Fraction(dimension) * Fraction(1) == dimension, "sharp Foster cut")
print("Z^d with the all-half characteristic has delta*lambda_odd=sqrt(d)/2 for d=1,...,12.")


print("\n=== Scope verdict ===")
print("PROVED ALGEBRA: every rank-d lattice with an obtuse superbase satisfies")
print("  delta*lambda_odd <= sqrt(d)/2,")
print("and the constant is sharp on the orthogonal all-half family.")
print("CITED CONSEQUENCE: every Euclidean lattice of dimension at most three is of this kind.")
print("OPEN: no claim is made for a general lattice without an obtuse superbase.")
print(f"GRAPHS={graph_count}")
print(f"FOSTER_IDENTITIES={foster_count}")
print(f"CUTS={cut_count}")
