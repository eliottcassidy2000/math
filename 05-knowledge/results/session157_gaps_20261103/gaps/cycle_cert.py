#!/usr/bin/env python3
"""Direct finite certificate for path counts of strong tournaments of order 3--6.

Every strong tournament has a directed Hamiltonian cycle. Fix that cycle,
orient each remaining chord in all possible ways, then test permutations.
"""
from collections import Counter
from itertools import permutations

EXPECTED = {
  3: {3: 1},
  4: {5: 4},
  5: {9: 10, 11: 5, 13: 10, 15: 7},
  6: {15: 12, 17: 12, 19: 12, 23: 36, 25: 24, 27: 12,
      29: 66, 31: 36, 33: 72, 37: 126, 41: 24, 43: 60, 45: 20},
}

def spectrum(n):
    cycle = {(i, (i + 1) % n) for i in range(n)}
    chords = [(i, j) for i in range(n) for j in range(i + 1, n)
              if (i, j) not in cycle and (j, i) not in cycle]
    counts = Counter()
    for mask in range(1 << len(chords)):
        edges = cycle | {(j, i) if (mask >> k) & 1 else (i, j)
                         for k, (i, j) in enumerate(chords)}
        h = sum(all((p[i], p[i + 1]) in edges for i in range(n - 1))
                for p in permutations(range(n)))
        counts[h] += 1
    if counts != EXPECTED[n]: raise AssertionError((n,counts))
    print(f'n={n}: {sorted(counts.items())}')
    return counts

if __name__ == '__main__':
    for n in range(3, 7):
        counts = spectrum(n)
        if 7 in counts or 21 in counts: raise AssertionError((n,counts))
    print('All assertions passed.')
