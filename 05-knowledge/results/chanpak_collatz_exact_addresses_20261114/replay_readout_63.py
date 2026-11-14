"""Literal finite checks for the order-eight H63 readout core.

This script verifies the example, not the low-count classification, which is
inherited from THM-344 and the audited all-order low-count inventory.
"""
from itertools import permutations

ROWS = (238, 220, 216, 16, 225, 206, 8, 72)
P = (0, 1, 2, 3, 4, 5, 7, 6)
Q = (0, 2, 4, 5, 1, 7, 6, 3)
ORDER = (3, 6, 7, 2, 5, 0, 4, 1)
TRIANGLES = ((0, 1, 4), (0, 2, 4), (1, 4, 5), (2, 4, 5), (3, 4, 7))


def arcs(p):
    return set(zip(p, p[1:]))


def check(condition, description):
    if not condition:
        raise ValueError(description)


def main():
    n = len(ROWS)
    all_arcs = {(i, j) for i in range(n) for j in range(n) if ROWS[i] >> j & 1}
    check(not any((i, i) in all_arcs for i in range(n)), "loop")
    check(all(((i, j) in all_arcs) != ((j, i) in all_arcs)
              for i in range(n) for j in range(i + 1, n)), "not a tournament")
    triangle_arcs = [arcs(t + (t[0],)) for t in TRIANGLES]
    check(all(t <= all_arcs for t in triangle_arcs), "not a directed triangle")
    hamiltonian = [w for w in permutations(range(n)) if arcs(w) <= all_arcs]
    check(len(hamiltonian) == 63, "wrong Hamiltonian count")
    check(all(any(not arcs(w) & t for t in triangle_arcs) for w in hamiltonian),
          "one Hamiltonian path meets all five triangles")
    check(P in hamiltonian and Q in hamiltonian, "not two Hamiltonian paths")
    pos = {v: i for i, v in enumerate(ORDER)}
    check(all(pos[i] > pos[j] for i, j in all_arcs - arcs(P) - arcs(Q)),
          "remaining arcs not all descending")
    print("63 Hamiltonian paths; each misses one of five triangles; two listed paths suffice.")


if __name__ == "__main__":
    main()
