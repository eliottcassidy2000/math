"""The four 4-vertex tournaments inside a tournament: which counts the scores determine, and the one trace.
Types: TT4 (transitive), W4 (a source over a 3-cycle), L4 (a 3-cycle over a sink), S4 (the strong 4-tournament).
With P = sum C(s_v, 3) (4-sets with a source) and Q = sum C(n-1-s_v, 3) (with a sink), N = C(n,4):
    #TT4 + #W4 = P,  #TT4 + #L4 = Q,  #TT4 + #W4 + #L4 + #S4 = N,
and S4 is the only 4-tournament with a directed 4-cycle (exactly one), while an oriented graph has no degenerate
closed 4-walk, so #S4 = tr(A^4)/4.  Checked on every class of orders 4..n (gentourng), together with:
  - the number of score sequences at which #S4 is not constant (so one matrix product is needed);
  - pattern-free classes for each of the four patterns (S4-free = transitive stacks of K1 and C3 blocks, since a strong
    tournament on >= 4 vertices has a 4-cycle; TT4-free only up to order 7 (every 8-tournament contains TT4)).
Usage: python3 pattern4.py nmax"""
import sys, subprocess
from itertools import combinations
from collections import defaultdict
nmax = int(sys.argv[1])
def classes(n):
    out = subprocess.run(["gentourng", "-q", str(n)], capture_output=True, text=True).stdout.split()
    for line in out:
        A = [[0]*n for _ in range(n)]; k = 0
        for i in range(n):
            for j in range(i+1, n):
                if line[k] == '1': A[i][j] = 1
                else: A[j][i] = 1
                k += 1
        yield A
def type4(A, q):
    sc = sorted(sum(A[u][v] for v in q if v != u) for u in q)
    if sc == [0,1,2,3]: return 'TT4'
    if sc == [1,1,1,3]: return 'W4'
    if sc == [0,2,2,2]: return 'L4'
    return 'S4'   # [1,1,2,2]
def matmul(X, Y):
    n = len(X); return [[sum(X[i][k]*Y[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
for n in range(4, nmax + 1):
    from math import comb
    free = defaultdict(int); byscores = defaultdict(set); ok = True; total = 0
    for A in classes(n):
        total += 1
        cnt = defaultdict(int)
        for q in combinations(range(n), 4): cnt[type4(A, q)] += 1
        s = [sum(r) for r in A]
        P = sum(comb(x, 3) for x in s); Q = sum(comb(n-1-x, 3) for x in s)
        A2 = matmul(A, A); A4 = matmul(A2, A2); tr4 = sum(A4[i][i] for i in range(n))
        if not (cnt['TT4'] + cnt['W4'] == P and cnt['TT4'] + cnt['L4'] == Q and sum(cnt.values()) == comb(n, 4) and 4 * cnt['S4'] == tr4):
            ok = False
        byscores[tuple(sorted(s))].add(cnt['S4'])
        for t in ('TT4', 'W4', 'L4', 'S4'):
            if cnt[t] == 0: free[t] += 1
    nonconst = sum(1 for v in byscores.values() if len(v) > 1)
    print(f"n={n}: classes {total}; identities hold: {ok}; score sequences {len(byscores)}, with #S4 not determined by the scores: {nonconst}; "
          f"pattern-free classes: TT4 {free['TT4']}, W4 {free['W4']}, L4 {free['L4']}, S4 {free['S4']}")
