#!/usr/bin/env python3
"""Exact certificate for phase-aware extension beyond matched residues.

The script has two purposes.

1. It gives a minimal one-coordinate, one-unit counterexample to carrying the
   canonical reflected high-pair floor to arbitrary residue packets.
2. It verifies a strict positive noncanonical packet through the phase-aware
   L1/Hunter stability lemma stated in the matching scratch note.

All interval arithmetic is exact Fraction arithmetic and this file is
self-contained.
"""
from fractions import Fraction as Q
from itertools import combinations


DMAX = Q(186_636_088_362, 11_773_143_757_375)
TARGET = DMAX / 5


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def clip(left, right):
    left = max(Q(), left)
    right = min(Q(1), right)
    return (left, right) if left < right else None


def arcs(L, z, j):
    """The exact open radius-1/14 clause on cell t=(j+u)/L."""
    require(L > 0 and z > 0 and 0 <= j < L, (L, z, j))
    lo = (z * j) // L - 1
    hi = (z * (j + 1)) // L + 1
    rows = []
    for n in range(lo, hi + 1):
        piece = clip(Q(L * (14 * n - 1), 14 * z) - j,
                     Q(L * (14 * n + 1), 14 * z) - j)
        if piece is not None:
            rows.append(piece)
    return tuple(sorted(set(rows)))


def merge(rows):
    answer = []
    for left, right in sorted(rows):
        if answer and left <= answer[-1][1]:
            if right > answer[-1][1]:
                answer[-1] = (answer[-1][0], right)
        else:
            answer.append((left, right))
    return tuple(answer)


def mass(rows):
    return sum((right - left for left, right in merge(rows)), Q())


def overlap(first, second):
    first, second = merge(first), merge(second)
    i = j = 0
    answer = Q()
    while i < len(first) and j < len(second):
        answer += max(Q(), min(first[i][1], second[j][1])
                      - max(first[i][0], second[j][0]))
        if first[i][1] < second[j][1]:
            i += 1
        else:
            j += 1
    return answer


def symmetric_difference(first, second):
    return mass(first) + mass(second) - 2 * overlap(first, second)


def maximum_spanning_tree(weights, edges):
    parent = list(range(6))

    def root(vertex):
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    tree = []
    credit = Q()
    for edge in sorted(edges, key=lambda item: (weights[item], item), reverse=True):
        left, right = map(root, edge)
        if left == right:
            continue
        parent[left] = right
        tree.append(edge)
        credit += weights[edge]
        if len(tree) == 5:
            return tuple(tree), credit
    raise RuntimeError("disconnected edge bank")


def hostile_pair():
    # Upper-median cell for the canonical hostile body (1,2,3,4,6,12).
    L, j = 168, 90
    z = 3 * L - 12
    w = 5 * L - 4
    perturbed = z - 1
    first = arcs(L, z, j)
    second = arcs(L, w, j)
    moved = arcs(L, perturbed, j)
    canonical_overlap = overlap(first, second)
    moved_overlap = overlap(moved, second)
    require(canonical_overlap == Q(6, 209) > TARGET, canonical_overlap)
    require(moved_overlap == 0, moved_overlap)
    require(first == ((Q(5, 41), Q(7, 41)),
                      (Q(19, 41), Q(21, 41)),
                      (Q(33, 41), Q(35, 41))), first)
    require(moved == ((Q(0), Q(6, 491)),
                      (Q(150, 491), Q(174, 491)),
                      (Q(318, 491), Q(342, 491)),
                      (Q(486, 491), Q(1))), moved)
    # h=-1 and m=-1 give Delta=h*j-m*L=78, a 13/28 phase jump.
    h, m = -1, -1
    delta = h * j - m * L
    require(delta == 78 and Q(abs(delta), L) == Q(13, 28), delta)
    return {
        "L": L,
        "j": j,
        "canonical": (z, w),
        "perturbed": (perturbed, w),
        "canonical_overlap": canonical_overlap,
        "perturbed_overlap": moved_overlap,
        "phase_jump": Q(abs(delta), L),
        "canonical_single": mass(first),
        "perturbed_single": mass(moved),
        "second_single": mass(second),
    }


def positive_tree():
    # A 649-bank upper-median body with an exact cell-phase stabilizer:
    # 28*j = 15*L.  The audit restricts to cross edges between the two
    # repeated-level classes, hence to a K_3,3 certificate.
    body = (1, 2, 3, 8, 9, 14)
    levels = (3, 3, 3, 5, 5, 5)
    L, j = 7056, 3780
    drifts = tuple(q * L - e for q, e in zip(levels, body))
    clauses = tuple(arcs(L, z, j) for z in drifts)
    edges = tuple((a, b) for a in range(3) for b in range(3, 6))
    weights = {edge: overlap(clauses[edge[0]], clauses[edge[1]])
               for edge in edges}
    tree, credit = maximum_spanning_tree(weights, edges)
    debt = sum((mass(clause) for clause in clauses), Q()) - Q(6, 7)
    margin = credit - debt
    expected_tree = ((1, 4), (1, 5), (0, 5), (2, 3), (0, 3))
    require(tree == expected_tree, tree)
    require(margin == Q(71_440_713_312_252_560_278,
                        527_394_888_495_258_135_905) > 0, margin)

    # Perturb only z_0 by h=28.  Its cell phase is unchanged because
    # 28*j=15*L, while its slope changes by 1/252.
    vertex, h, phase_integer = 0, 28, 15
    require(h * j == phase_integer * L, (h, j, L))
    moved_drifts = list(drifts)
    moved_drifts[vertex] += h
    moved_clauses = list(clauses)
    moved_clauses[vertex] = arcs(L, moved_drifts[vertex], j)
    moved_clauses = tuple(moved_clauses)
    sigma = symmetric_difference(clauses[vertex], moved_clauses[vertex])

    # General phase-aware bound.  Here Delta=0 and eta=h/L=1/252.
    delta = h * j - phase_integer * L
    eta = Q(abs(delta) + abs(h), L)
    sigma_bound = min(Q(1), 4 * eta * (1 + Q(L, drifts[vertex])))
    require(eta == Q(1, 252), eta)
    require(sigma == Q(67_424, 16_616_095), sigma)
    require(sigma_bound == Q(28_223, 1_333_521), sigma_bound)
    require(sigma <= sigma_bound, (sigma, sigma_bound))

    degree = sum(vertex in edge for edge in tree)
    require(degree == 2, degree)
    exact_stability_floor = margin - (degree + 1) * sigma
    coarse_stability_floor = margin - (degree + 1) * sigma_bound
    require(exact_stability_floor == Q(
        2_041_646_828_324_165_638_790,
        16_560_199_498_751_105_467_417) > 0, exact_stability_floor)
    require(coarse_stability_floor == Q(
        113_864_784_228_062_404_699,
        1_582_184_665_485_774_407_715) > 0, coarse_stability_floor)

    moved_credit = sum(
        (overlap(moved_clauses[a], moved_clauses[b]) for a, b in tree), Q()
    )
    moved_debt = sum((mass(clause) for clause in moved_clauses), Q()) - Q(6, 7)
    moved_margin = moved_credit - moved_debt
    require(moved_margin == Q(
        1_585_980_049_263_141_094,
        11_735_389_638_648_206_265) > exact_stability_floor,
        moved_margin)
    return {
        "body": body,
        "L": L,
        "j": j,
        "levels": levels,
        "tree": tree,
        "margin": margin,
        "perturbed_vertex": vertex,
        "h": h,
        "delta": delta,
        "eta": eta,
        "sigma": sigma,
        "sigma_bound": sigma_bound,
        "exact_stability_floor": exact_stability_floor,
        "coarse_stability_floor": coarse_stability_floor,
        "actual_moved_margin": moved_margin,
    }


def main():
    hostile = hostile_pair()
    positive = positive_tree()
    print("LRC14 REFLECTED-RESIDUE PHASE/STABILITY AUDIT")
    print("target", TARGET)
    print("hostile", hostile)
    print("positive", positive)
    print("map=B_(z+h)=phi_h(B_z), phi_h(u)=(z*u-Delta)/(z+h), Delta=h*j-m*L")
    print("L1_bound=sigma_i<=min(1,4*eta_i*(1+L/z_i)), eta_i=(abs(Delta_i)+abs(h_i))/L")
    print("tree_bound=M_prime>=M-sum_i(deg_T(i)+1)*sigma_i")
    print("conclusion=arbitrary-residue pair floor REFUTED; phase-aware frozen-tree chamber PROVED")
    print("status=PASS")


if __name__ == "__main__":
    main()
