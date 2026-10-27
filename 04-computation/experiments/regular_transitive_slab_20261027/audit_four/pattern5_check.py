"""pattern5_check.py -- independent check of the 'First data' paragraph appended to the sigma-four reflection (section 5,
item 1): which of the twelve 5-vertex type counts are functions of (sorted scores, tr A^4, tr A^5), at orders 5..8,
and the pattern-free class counts per type.  Written independently of pattern5.py: the 5-vertex types are identified
by orbit enumeration (the S_5-orbit of each labelled 5-tournament, built by BFS over transpositions), and labelled by
(sorted scores, c3, c4, c5, H, |Aut|) so that the result can be compared with any numbering.
Usage: python3 pattern5_check.py nmax
"""
import sys, subprocess, itertools
import numpy as np
nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 8
P5 = [(i, j) for i in range(5) for j in range(i + 1, 5)]
IDX = {p: k for k, p in enumerate(P5)}
def adj_from_bits(b):
    A = np.zeros((5, 5), dtype=int)
    for k, (i, j) in enumerate(P5):
        if b >> k & 1: A[i, j] = 1
        else: A[j, i] = 1
    return A
def bits_from_adj(A):
    return sum(1 << k for k, (i, j) in enumerate(P5) if A[i, j])
def relabel(b, perm):
    A = adj_from_bits(b); B = np.zeros_like(A)
    for i in range(5):
        for j in range(5):
            B[perm[i], perm[j]] = A[i, j]
    return bits_from_adj(B)
# orbits of S_5 on the 1024 labelled 5-tournaments, by closure under the four adjacent transpositions
orbit_of = [-1] * 1024; orbits = []
gens = [tuple(list(range(t)) + [t + 1, t] + list(range(t + 2, 5))) for t in range(4)]
for b in range(1024):
    if orbit_of[b] >= 0: continue
    k = len(orbits); stack = [b]; orbit_of[b] = k; members = [b]
    while stack:
        x = stack.pop()
        for g in gens:
            y = relabel(x, g)
            if orbit_of[y] < 0: orbit_of[y] = k; members.append(y); stack.append(y)
    orbits.append(members)
assert len(orbits) == 12, len(orbits)
def invariants(A):
    n = len(A)
    sc = tuple(sorted(int(x) for x in A.sum(1)))
    c3 = int(np.trace(np.linalg.matrix_power(A, 3))) // 3
    c4 = int(np.trace(np.linalg.matrix_power(A, 4))) // 4
    c5 = int(np.trace(np.linalg.matrix_power(A, 5))) // 5
    H = sum(1 for p in itertools.permutations(range(n)) if all(A[p[i], p[i + 1]] for i in range(n - 1)))
    return sc, c3, c4, c5, H
labels = []
for k, mem in enumerate(orbits):
    A = adj_from_bits(mem[0]); sc, c3, c4, c5, H = invariants(A)
    labels.append(f"scores{sc} c3={c3} c4={c4} c5={c5} H={H} |Aut|={120 // len(mem)}")
order = sorted(range(12), key=lambda k: labels[k])
print("the twelve 5-vertex types:")
for k in order: print(f"   T{k:02d}: {labels[k]}")
def classes(n):
    out = subprocess.run(["gentourng", "-q", str(n)], capture_output=True, text=True).stdout.split()
    for line in out:
        A = np.zeros((n, n), dtype=int); k = 0
        for i in range(n):
            for j in range(i + 1, n):
                if line[k] == '1': A[i, j] = 1
                else: A[j, i] = 1
                k += 1
        yield A
for n in range(5, nmax + 1):
    groups = {}; groups_full = {}; free = [0] * 12; total = 0
    subsets = list(itertools.combinations(range(n), 5))
    for A in classes(n):
        total += 1
        cnt = [0] * 12
        for q in subsets:
            b = 0
            for k, (i, j) in enumerate(P5):
                if A[q[i], q[j]]: b |= 1 << k
            cnt[orbit_of[b]] += 1
        key = (tuple(sorted(int(x) for x in A.sum(1))), int(np.trace(np.linalg.matrix_power(A, 4))), int(np.trace(np.linalg.matrix_power(A, 5))))
        groups.setdefault(key, []).append(tuple(cnt))
        # stronger statistic: sorted scores together with ALL trace moments tr A^k, k <= n (the full spectrum, by Newton)
        keyfull = (key[0],) + tuple(int(np.trace(np.linalg.matrix_power(A, k))) for k in range(3, n + 1))
        groups_full.setdefault(keyfull, []).append(tuple(cnt))
        for t in range(12):
            if cnt[t] == 0: free[t] += 1
    determined = [t for t in range(12) if all(len({c[t] for c in v}) == 1 for v in groups.values())]
    print(f"n={n}: classes {total}; statistic groups {len(groups)}; types whose count is a function of (scores, trA^4, trA^5): {[f'T{t:02d}' for t in determined]}")
    detfull = [t for t in range(12) if all(len({c[t] for c in v}) == 1 for v in groups_full.values())]
    print(f"      with (scores, all of tr A^3..tr A^n): groups {len(groups_full)}; determined types: {[f'T{t:02d}' for t in detfull]}")
    print(f"      pattern-free classes per type: " + ", ".join(f"T{t:02d}:{free[t]}" for t in order) + f"; min {min(free)}, max {max(free)}")
