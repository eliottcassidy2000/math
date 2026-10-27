"""The twelve 5-vertex tournaments inside a tournament (data for the reflection's next step 1, after THM-4823):
which of the 12 type counts are functions of (sorted scores, tr A^4, tr A^5)?  For every class of orders 5..nmax
(gentourng), count the 5-subsets of each isomorphism type, group the classes by the statistic, and report per type
whether its count is constant on every group.  Also reports the pattern-free class counts per type.
Usage: python3 pattern5.py nmax"""
import sys, subprocess
from itertools import combinations, permutations
from collections import defaultdict
nmax = int(sys.argv[1])
# canonical type of a 5-tournament from its 10 arc bits (lexicographically least relabelled bit string)
pairs5 = [(i, j) for i in range(5) for j in range(i + 1, 5)]
def canon5(adj):
    best = None
    for p in permutations(range(5)):
        key = tuple(1 if adj[p[i]][p[j]] else 0 for i, j in pairs5)
        if best is None or key < best: best = key
    return best
types = {}
for bits in range(1 << 10):
    A = [[0] * 5 for _ in range(5)]
    for k, (i, j) in enumerate(pairs5):
        if bits >> k & 1: A[i][j] = 1
        else: A[j][i] = 1
    c = canon5(A)
    if c not in types: types[c] = len(types)
code = {}
for bits in range(1 << 10):
    A = [[0] * 5 for _ in range(5)]
    for k, (i, j) in enumerate(pairs5):
        if bits >> k & 1: A[i][j] = 1
        else: A[j][i] = 1
    code[bits] = types[canon5(A)]
T = len(types); assert T == 12, T
def type_scores(t):
    c = [k for k, v in types.items() if v == t][0]
    A = [[0] * 5 for _ in range(5)]
    for k, (i, j) in enumerate(pairs5):
        if c[k]: A[i][j] = 1
        else: A[j][i] = 1
    return tuple(sorted(sum(r) for r in A))
def classes(n):
    for line in subprocess.run(["gentourng", "-q", str(n)], capture_output=True, text=True).stdout.split():
        A = [[0] * n for _ in range(n)]; k = 0
        for i in range(n):
            for j in range(i + 1, n):
                if line[k] == '1': A[i][j] = 1
                else: A[j][i] = 1
                k += 1
        yield A
def mat(X, Y):
    n = len(X); return [[sum(X[i][k] * Y[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
print("types (index: sorted scores):", {t: type_scores(t) for t in range(T)})
for n in range(5, nmax + 1):
    groups = defaultdict(list); free = [0] * T; total = 0
    for A in classes(n):
        total += 1
        cnt = [0] * T
        for q in combinations(range(n), 5):
            bits = 0
            for k, (i, j) in enumerate(pairs5):
                if A[q[i]][q[j]]: bits |= 1 << k
            cnt[code[bits]] += 1
        for t in range(T):
            if cnt[t] == 0: free[t] += 1
        A2 = mat(A, A); A4 = mat(A2, A2); A5 = mat(A4, A)
        stat = (tuple(sorted(sum(r) for r in A)), sum(A4[i][i] for i in range(n)), sum(A5[i][i] for i in range(n)))
        groups[stat].append(cnt)
    determined = [all(len({g[t] for g in grp}) == 1 for grp in groups.values()) for t in range(T)]
    print(f"n={n}: classes {total}; statistic groups {len(groups)}; type counts determined by (scores, tr A^4, tr A^5): "
          f"{[t for t in range(T) if determined[t]]}; not determined: {[t for t in range(T) if not determined[t]]}; pattern-free classes per type: {free}")
