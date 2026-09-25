"""Finite exact controls for the four-core tournament substitution note."""
from itertools import combinations, product
from math import comb


def arc(u, v):
    assert u != v
    if u == 9:
        return True
    if v == 9:
        return False
    bu, iu = divmod(u, 3)
    bv, iv = divmod(v, 3)
    return (iv - iu) % 3 == 1 if bu == bv else (bv - bu) % 3 == 1


def hamiltonian_path(vertices, edge):
    path = []
    for v in vertices:
        position = next((i for i, w in enumerate(path) if edge(v, w)), len(path))
        path.insert(position, v)
    assert len(set(path)) == len(vertices)
    assert all(edge(u, v) for u, v in zip(path, path[1:]))
    return path


def main():
    for A in range(1, 26):
        vertices = [(i, j) for i in range(3) for j in range(A)] + [(3, 0)]
        counts = [0, 0, 0, 0]
        for (i, j), (k, ell) in combinations(vertices, 2):
            if i == k:
                counts[0] += 1
            elif j and ell:
                counts[1] += 1
            elif j == ell == 0:
                counts[3] += 1
            else:
                counts[2] += 1
        expected = [3 * comb(A, 2), 3 * (A - 1)**2, 9 * (A - 1), 6]
        assert counts == expected and sum(counts) == comb(3 * A + 1, 2)
        original = expected[0] + expected[1] + 9 * A + 6
        assert original - sum(counts) == 9
    print("Edge partition: all A=1..25, unfiltered; direct pair classification agrees with formula.")
    print("Original count exceeds actual count by9; A=1 gives15 vs6, A=3 gives54 vs45.")

    vertices = list(range(10))
    assert all(arc(u, v) != arc(v, u) for u, v in combinations(vertices, 2))
    for reps in product(range(3), repeat=3):
        marked = [3 * i + reps[i] for i in range(3)] + [9]
        for i, j in combinations(range(4), 2):
            expected = False if j == 3 else (j - i) % 3 == 1
            assert arc(marked[i], marked[j]) == expected
    print("Positive control: Q[C3,C3,C3,1], Q=source over C3; all27 representative cores equalQ.")

    pair_modules = [(u, v) for u, v in combinations(vertices, 2)
                    if all(arc(w, u) == arc(w, v) for w in vertices if w not in (u, v))]
    assert pair_modules == []
    print("Halving hostile: this same10-vertex substitution has0 two-vertex modules among45 pairs.")
    forward = hamiltonian_path(vertices, arc)
    converse = lambda u, v: arc(v, u)
    backward = hamiltonian_path(vertices, converse)
    assert all(not converse(9, v) for v in range(9))
    print(f"Hamiltonian positive path: {forward}")
    print(f"Rooted-route hostile: converse has Hamiltonian path{backward}, but root9 has no outgoing arc.")
    print("FINITE-EXACT controls only; no universal Collatz certificate is claimed.")


if __name__ == "__main__":
    main()
