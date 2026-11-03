"""Finite 16-case lemma: two complementary directed triangle pairs imply C>=10."""
from itertools import combinations, permutations

edges = list(combinations(range(6), 2))
# Directed triangles: (0,1,2), (3,4,5), (0,1,3), (2,4,5).
forced_arcs = [(0,1), (1,2), (2,0), (3,4), (4,5), (5,3),
               (1,3), (3,0), (2,4), (5,2)]
forced = {(min(a,b), max(a,b)): a < b for a,b in forced_arcs}
free = [(0,4), (0,5), (1,4), (1,5), (2,3)]
if len(forced) != 10 or set(forced) | set(free) != set(edges):
    raise RuntimeError('bad partition of arcs')
expected = [
    [(6,4), (6,6), (7,7), (6,8)],
    [(7,7), (6,8), (8,8), (6,6)],
    [(8,6), (8,8), (8,6), (7,7)],
    [(8,6), (7,7), (8,6), (6,4)],
]
for mask in range(32):
    orientation = dict(forced)
    orientation.update((e, bool(mask & (1 << i))) for i,e in enumerate(free))
    def arc(a, b):
        return orientation[min(a,b), max(a,b)] == (a < b)
    cycles = []
    for k in (3,5):
        for subset in combinations(range(6), k):
            for perm in permutations(subset[1:]):
                cycle = (subset[0],) + perm
                if all(arc(cycle[i], cycle[(i+1) % k]) for i in range(k)):
                    cycles.append(cycle)
    result = tuple(sum(len(cycle) == k for cycle in cycles) for k in (3,5))
    row, col = mask & 3, (mask >> 2) & 3
    if result != expected[row][col]:
        raise RuntimeError(f'mismatch for mask {mask}: {result} != {expected[row][col]}')
    if sum(result) < 10:
        raise RuntimeError(f'counterexample for mask {mask}')
print('PASS: all 32 orientations match the 4x4 table and satisfy C >= 10.')
