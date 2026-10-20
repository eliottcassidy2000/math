"""Exact seven-port junction tomography, with independent tree controls.

Universe: all 546 ternary words of length seven with all three multiplicities
odd, quotiented by the free global colour-permutation action (91 orbits).
No numerical rank or graph-isomorphism package is used.
"""
from collections import Counter
from itertools import combinations, product
import json
from pathlib import Path


def require(test, message):
    if not test:
        raise AssertionError(message)


def canonical_word(word):
    mapping = {}
    return tuple(mapping.setdefault(c, len(mapping)) for c in word)


def matchings(items):
    if not items:
        yield ()
        return
    a = items[0]
    for i in range(1, len(items)):
        for rest in matchings(items[1:i] + items[i + 1:]):
            yield ((a, items[i]),) + rest


def rank_mod(rows, prime=101):
    a = [list(row) for row in rows]
    r = 0
    pivot_columns = []
    for col in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][col] % prime),
                     None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inv = pow(a[r][col] % prime, -1, prime)
        a[r] = [(x * inv) % prime for x in a[r]]
        for i in range(r + 1, len(a)):
            scale = a[i][col] % prime
            if scale:
                a[i] = [(x - scale * y) % prime
                        for x, y in zip(a[i], a[r])]
        pivot_columns.append(col)
        r += 1
    return r, pivot_columns


PORT_VERTEX = (0, 0, 1, 2, 3, 4, 4)


def enumerate_path_tree():
    """Enumerate internal colours, independently of the cut-XOR criterion."""
    counts = Counter()
    for internal in product(range(3), repeat=4):
        if any(internal[i] == internal[i + 1] for i in range(3)):
            continue
        middle = [3 - internal[i] - internal[i + 1] for i in range(3)]
        left = [c for c in range(3) if c != internal[0]]
        right = [c for c in range(3) if c != internal[-1]]
        for l in (left, left[::-1]):
            for r in (right, right[::-1]):
                counts[tuple(l + middle + r)] += 1
    return counts


def xor_tree_count(word):
    # Colours 0,1,2 correspond to the nonzero elements 1,2,3 of F_2^2.
    if any(word.count(c) % 2 != 1 for c in range(3)):
        return 0
    for cut in range(4):
        flow = 0
        for port, vertex in enumerate(PORT_VERTEX):
            if vertex <= cut:
                flow ^= word[port] + 1
        if flow == 0:
            return 0
    return 1


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def main():
    physical = [word for word in product(range(3), repeat=7)
                if all(word.count(c) % 2 for c in range(3))]
    words = sorted({canonical_word(word) for word in physical})
    require(len(physical) == 546 and len(words) == 91, "word universe")
    require(set(Counter(canonical_word(w) for w in physical).values()) == {6},
            "colour action must be free")
    kernels = []
    for triple in combinations(range(7), 3):
        rest = tuple(v for v in range(7) if v not in triple)
        for pairs in matchings(rest):
            kernels.append([int(len({word[v] for v in triple}) == 3
                                and all(word[a] == word[b] for a, b in pairs))
                            for word in words])
    require(len(kernels) == 105, "junction universe")
    require({sum(row) for row in kernels} == {9}, "J support")
    defects = [[(1 if sorted(Counter(w).values()) == [1, 1, 5] else -1)
                if w.count(w[i]) == 1 else 0 for w in words]
               for i in range(7)]
    require(all(dot(row, w) == 0 for row in kernels for w in defects),
            "every J annihilates every singleton defect")
    require(rank_mod(defects)[0] == 7, "independent defects")
    rank, pivots = rank_mod(kernels)
    require(rank == 84, "rank over F_101")
    # Seven explicit rational null vectors bound Q-rank by 84, and the
    # nonzero mod-101 minor bounds it below by 84: exact rational rank.
    enumerated = enumerate_path_tree()
    require(len(enumerated) == 96 and set(enumerated.values()) == {1},
            "P5 tree extension counts")
    require(all(enumerated[w] == xor_tree_count(w)
                for w in product(range(3), repeat=7)),
            "independent cut-XOR vs internal-edge enumeration")
    trees = []
    for i in range(7):
        def read(word):
            w = list(word)
            w[i], w[3] = w[3], w[i]
            return enumerated[tuple(w)]
        trees.append([read(w) for w in words])
    pairings = [[dot(tree, w) for w in defects] for tree in trees]
    require(pairings == [[-4 * int(i == j) for j in range(7)]
                         for i in range(7)], "tree-defect diagonal pairing")
    full_rank = rank_mod(kernels + trees)[0]
    require(full_rank == 91, "complete tomography")
    receipt = {
        "status": "FINITE-EXACT; exact rational rank follows from integer bounds",
        "physical_parity_words": len(physical),
        "colour_orbits": len(words),
        "junction_kernels": len(kernels),
        "junction_word_support": 54,
        "junction_rank_Q": rank,
        "nullity_Q": 91 - rank,
        "modular_prime": 101,
        "modular_pivot_columns": pivots,
        "column_sums_by_pattern": {
            str(pattern): sorted({sum(row[k] for row in kernels)
                                  for k, w in enumerate(words)
                                  if tuple(sorted(Counter(w).values())) == pattern})
            for pattern in ((1, 1, 5), (1, 3, 3))
        },
        "tree_internal_enumeration_equals_cut_XOR_on_all_2187_words": True,
        "tree_support": len(enumerated),
        "tree_defect_pairings": pairings,
        "augmented_rank_Q": full_rank,
        "augmented_physical_closures": len(kernels) + len(trees),
    }
    target = Path(__file__).with_name("junctions.json")
    target.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
