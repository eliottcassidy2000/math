"""Rooted digit ratio B(T,t)/H(T), B(T,t) = sum over paths P starting at t of H(T - V(P)) (wave 56; near-join H = h + 2^k B).
7h is a near-join count iff some rooted realizer of h has B = 3h (k = 1 is forced: 6h/2^k must be an integer).
Scan all classes of order <= 7."""
import sys
sys.path.insert(0, '.')
from w45_hk import H_of, sub
from collections import defaultdict
def allH(n, out):
    """H of every induced subtournament (by mask) and path-start counts dp[S][v] (paths covering S starting at v)."""
    N = 1 << n
    st = [[0] * n for _ in range(N)]   # paths covering S that START at v
    for v in range(n): st[1 << v][v] = 1
    # build by prepending: path starting at v covering S = v -> (path starting at w covering S - v) with v -> w
    for S in range(1, N):
        for v in range(n):
            if not S >> v & 1 or S == 1 << v: continue
            R = S & ~(1 << v); tot = 0; o = out[v] & R
            while o:
                w = (o & -o).bit_length() - 1; o &= o - 1; tot += st[R][w]
            st[S][v] = tot
    Hs = [1] + [sum(st[S]) for S in range(1, N)]
    return Hs, st
def B_all(n, out):
    Hs, st = allH(n, out); full = (1 << n) - 1
    B = [0] * n
    for S in range(1, 1 << n):
        rest = full & ~S
        for v in range(n):
            if st[S][v]: B[v] += st[S][v] * Hs[rest]
    return Hs[full], B
pairs = defaultdict(set); ratio3 = defaultdict(list)
for line in open('w57_reps.txt'):
    t = list(map(int, line.split())); n = t[0]; out = t[1:]
    h, B = B_all(n, out)
    for v in range(n):
        pairs[h].add(B[v])
        if B[v] == 3 * h: ratio3[h].append((n, out, v))
print("B values available at small counts (orders <= 7):")
for h in (1, 3, 5, 9, 11, 13, 15):
    print(f"  h={h}: B in {sorted(pairs[h])[:30]}{' ...' if len(pairs[h]) > 30 else ''}")
print("counts h with a ratio-3 root (B = 3h) at orders <= 7:", sorted(ratio3)[:60])
print("odd h <= 45 realized at orders <= 7 WITHOUT a ratio-3 root:", [h for h in sorted(pairs) if h <= 45 and h not in ratio3])
