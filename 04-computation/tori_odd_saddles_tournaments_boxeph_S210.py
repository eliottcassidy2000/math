#!/usr/bin/env python3
"""Exact controls for corrected HYP-8835.

The original S210 script promoted several analogies to one antisymmetry
theorem.  This version checks the independent facts and prints the scope
failures: odd tournament-game support uses a mod-2 tournament coordinate;
pure optimal play means Condorcet winner, not transitivity; the RPS invariant
comes from a positive kernel vector; and the cutting bagel is a solid torus,
not its boundary T^2.
"""

from fractions import Fraction
from itertools import combinations, product
from math import comb, pi, sin


def section(title: str) -> None:
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def matrix_rank(matrix: list[list[Fraction | int]]) -> int:
    work = [[Fraction(x) for x in row] for row in matrix]
    rows = len(work)
    cols = len(work[0]) if rows else 0
    rank = 0
    for col in range(cols):
        pivot = next((i for i in range(rank, rows) if work[i][col]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        scale = work[rank][col]
        work[rank] = [x / scale for x in work[rank]]
        for i in range(rows):
            if i != rank and work[i][col]:
                factor = work[i][col]
                work[i] = [work[i][j] - factor * work[rank][j] for j in range(cols)]
        rank += 1
    return rank


def rank_mod2(matrix: list[list[int]]) -> int:
    work = [[x & 1 for x in row] for row in matrix]
    rows = len(work)
    cols = len(work[0]) if rows else 0
    rank = 0
    for col in range(cols):
        pivot = next((i for i in range(rank, rows) if work[i][col]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        for i in range(rows):
            if i != rank and work[i][col]:
                work[i] = [a ^ b for a, b in zip(work[i], work[rank])]
        rank += 1
    return rank


def nullspace(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    work = [row[:] for row in matrix]
    rows = len(work)
    cols = len(work[0]) if rows else 0
    pivots: list[int] = []
    rank = 0
    for col in range(cols):
        pivot = next((i for i in range(rank, rows) if work[i][col]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        scale = work[rank][col]
        work[rank] = [x / scale for x in work[rank]]
        for i in range(rows):
            if i != rank and work[i][col]:
                factor = work[i][col]
                work[i] = [work[i][j] - factor * work[rank][j] for j in range(cols)]
        pivots.append(col)
        rank += 1
    basis: list[list[Fraction]] = []
    for free in (col for col in range(cols) if col not in pivots):
        vector = [Fraction(0)] * cols
        vector[free] = Fraction(1)
        for row, pivot in enumerate(pivots):
            vector[pivot] = -work[row][free]
        basis.append(vector)
    return basis


def all_tournaments(n: int):
    edges = list(combinations(range(n), 2))
    for bits in product((0, 1), repeat=len(edges)):
        adjacency = [[0] * n for _ in range(n)]
        for (i, j), bit in zip(edges, bits):
            adjacency[i][j] = bit
            adjacency[j][i] = 1 - bit
        yield adjacency


def payoff(adjacency: list[list[int]]) -> list[list[int]]:
    n = len(adjacency)
    return [
        [adjacency[i][j] - adjacency[j][i] for j in range(n)]
        for i in range(n)
    ]


def optimal_supports(adjacency: list[list[int]]) -> list[tuple[int, ...]]:
    """Enumerate strict supports satisfying Mp<=0 and Mp=0 on support."""
    n = len(adjacency)
    matrix = [[Fraction(x) for x in row] for row in payoff(adjacency)]
    found: list[tuple[int, ...]] = []
    for size in range(1, n + 1):
        for support in combinations(range(n), size):
            principal = [[matrix[i][j] for j in support] for i in support]
            kernel = nullspace(principal)
            if len(kernel) != 1:
                continue
            vector = kernel[0]
            if not (all(x > 0 for x in vector) or all(x < 0 for x in vector)):
                continue
            vector = [abs(x) for x in vector]
            total = sum(vector)
            strategy = [Fraction(0)] * n
            for local, vertex in enumerate(support):
                strategy[vertex] = vector[local] / total
            image = [
                sum(matrix[i][j] * strategy[j] for j in range(n))
                for i in range(n)
            ]
            if all(x <= 0 for x in image):
                found.append(support)
        if found:
            break
    return found


section("P1  Even skew rank is not by itself the odd-support theorem")
for size in (2, 4, 6):
    mod2_principal = [[int(i != j) for j in range(size)] for i in range(size)]
    assert rank_mod2(mod2_principal) == size
    print(f"even support size {size}: tournament payoff block mod 2 has full rank {size}")
print("An optimal support block must be singular; the mod-2 tournament form excludes even size.")

for n in (3, 4, 5):
    support_sizes: set[int] = set()
    count = 0
    for adjacency in all_tournaments(n):
        count += 1
        matrix = payoff(adjacency)
        assert matrix_rank(matrix) % 2 == 0
        supports = optimal_supports(adjacency)
        assert len(supports) == 1
        support_sizes.add(len(supports[0]))
        assert len(supports[0]) % 2 == 1
    print(f"n={n}: all {count} tournaments; unique-support sizes={sorted(support_sizes)}")


section("P2  Pure optimum means Condorcet winner, not transitivity")
# Vertex 0 beats everyone; vertices 1,2,3 form 1->2->3->1.
witness = [
    [0, 1, 1, 1],
    [0, 0, 1, 0],
    [0, 0, 0, 1],
    [0, 1, 0, 0],
]
matrix = payoff(witness)
assert optimal_supports(witness) == [(0,)]
assert witness[1][2] == witness[2][3] == witness[3][1] == 1
print("4-vertex witness: 0 is Condorcet winner; 1->2->3->1 is a directed cycle.")
print("The tournament is intransitive but its unique optimal support is {0}.")
print("For interior replicator states, x_0' = x_0(1-x_0), so x_0 tends to 1.")


section("P3  The correct dynamical hinge is a positive kernel vector")
rps = [[0, 1, -1], [-1, 0, 1], [1, -1, 0]]
p = [Fraction(1, 3)] * 3
assert [sum(rps[i][j] * p[j] for j in range(3)) for i in range(3)] == [0, 0, 0]
assert [sum(rps[i][j] for i in range(3)) for j in range(3)] == [0, 0, 0]
print("RPS has p=(1/3,1/3,1/3) in ker(M).")
print("Hence d/dt sum_i p_i log(x_i)=p^T Mx=0; equivalently x0*x1*x2 is conserved.")
print("Positive regular levels inside Delta^2 are periodic circles S^1, not a surface T^2.")

# A skew matrix need not have zero column sums.
generic_skew = [[0, 1, 1], [-1, 0, 1], [-1, -1, 0]]
column_sums = [sum(generic_skew[i][j] for i in range(3)) for j in range(3)]
assert column_sums != [0, 0, 0]
print(f"hostile skew control column sums={column_sums}; antisymmetry alone does not conserve the product.")


section("P4  Closed torus, solid torus, and cutting sequence are different objects")
for n in range(1, 6):
    betti = [comb(n, k) for k in range(n + 1)]
    chi = sum((-1) ** k * betti[k] for k in range(n + 1))
    assert chi == 0
    print(f"T^{n}: Betti={betti}, chi={chi}")
print("solid torus D^2 x S^1: Betti=[1,1], chi=0")
print("boundary T^2: Betti=[1,2,1], chi=0; a perfect Morse function has counts 1,2,1")
print("No computation here identifies bagel(n)-cake(n)=T_n-1 with either Euler characteristic.")


section("P5  Three different sign actions")
for n in (1, 2, 3):
    fixed = list(product((Fraction(0), Fraction(1, 2)), repeat=n))
    assert len(fixed) == 2**n
    print(f"torus inversion on T^{n}: {len(fixed)} fixed 2-torsion points")

delta = 0.2
for k in (1, 2, 3):
    ghat = lambda h: -sin(2 * pi * h * delta) / (pi * h)
    assert abs(ghat(k) - ghat(-k)) < 1e-12
print("LRC far-set Fourier weights are even, not odd.")


def vandermonde(values: list[int]) -> int:
    result = 1
    for i, j in combinations(range(len(values)), 2):
        result *= values[j] - values[i]
    return result


values = [1, 2, 4, 7]
assert vandermonde([2, 1, 4, 7]) == -vandermonde(values)
print("Vandermonde changes sign under a coordinate transposition.")
print("Matrix complement M->-M, torus inversion, and permutation parity are distinct Z/2 actions.")


section("SUMMARY")
print(
    "Exact survivors: skew rank parity; tournament-specific odd optimal support; torus Betti data;\n"
    "the RPS positive-kernel invariant; torus-inversion fixed points; even LRC weights; and\n"
    "alternating Vandermonde.  The one-hinge synthesis, pure<=>transitive claim, bagel Euler\n"
    "bridge, and intransitivity=>toroidal recurrence claim are refuted."
)
