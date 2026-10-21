"""Exact gluings of two interval chains (session opus-2026-10-21-S50).  Note: 05-knowledge/results/
h_spectrum_multiplicities_20261021.md; THM-4745, HYP-9356.

For chains A = IC(s), B = IC(t) - all ordered pairs with parts in {2,3,4,5}, at most 6 vertices each and
(|A|-1)(|B|-1) <= 16, plus the four pairs of C_3 or R_4 with the optimal 7-vertex chain IC(3,3) - every choice of
the shared vertex v (position i in A, j in B) and every orientation of the (|A|-1)(|B|-1) cross arcs between A - v and
B - v, compute H(T) for T = A u B by subset DP and record the exact gluings H(T) = H(A) H(B).
Checks:
  (X1) Busch's merge inequality H(T) >= H(A) H(B) holds for every gluing;
  (X2) the one-way lemma: in every exact gluing all cross arcs point the same way;
  (X3) every exact gluing is the end-to-start concatenation: with cross arcs A' -> B', v is the last spine vertex of A
       and the first of B (any vertex when the chain is the 3-cycle, whose rotations are automorphisms), and
       symmetrically for B' -> A';
  (X4) conversely every pair (A, B) has its two concatenations among the exact gluings.
Prints ALL CHECKS PASSED.  Runtime about 10 minutes (numba, parallel over chunks of cross patterns).  The range covers
every pair of extremal chain orders that enters the induction of THM-4745 (4) up to order 10: (3,3), (3,4), (4,4),
(3,6), (4,5), (4,6), (3,7), (4,7)."""
import sys, os, time
import numpy as np
import numba as nb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from h_multiplicities import ham_one                                     # noqa: E402

FAIL = []
T0 = time.time()


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m, flush=True)
    if not c:
        FAIL.append(m)


def chain(c):
    n = 1 + sum(c)
    adj = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            adj[i] |= 1 << j
    l = 0
    for p in c:
        r = l + p
        adj[l] &= ~(1 << r); adj[r] |= 1 << l
        l = r
    return adj


def comps(total, allowed):
    if total == 0:
        yield ()
        return
    for p in allowed:
        if p <= total:
            for rest in comps(total - p, allowed):
                yield (p,) + rest


@nb.njit(cache=True)
def ham_buf(n, adj, dp):
    N = 1 << n
    for S in range(N):
        for v in range(n):
            dp[S, v] = 0
    for v in range(n):
        dp[1 << v, v] = 1
    for S in range(1, N):
        for v in range(n):
            c = dp[S, v]
            if c == 0:
                continue
            for u in range(n):
                if (S >> u) & 1:
                    continue
                if (adj[v] >> u) & 1:
                    dp[S | (1 << u), u] += c
    tot = 0
    for v in range(n):
        tot += dp[N - 1, v]
    return tot


@nb.njit(parallel=True, cache=True)
def glue_all(adjA, adjB, i, j, out):
    a = adjA.shape[0]; b = adjB.shape[0]
    n = a + b - 1
    Aidx = np.empty(a, dtype=np.int64); Bidx = np.empty(b, dtype=np.int64)
    k = 0
    for x in range(a):
        if x == i:
            Aidx[x] = a - 1
        else:
            Aidx[x] = k; k += 1
    k = a
    for y in range(b):
        if y == j:
            Bidx[y] = a - 1
        else:
            Bidx[y] = k; k += 1
    base = np.zeros(n, dtype=np.int64)
    for x in range(a):
        for x2 in range(a):
            if (adjA[x] >> x2) & 1:
                base[Aidx[x]] |= 1 << Aidx[x2]
    for y in range(b):
        for y2 in range(b):
            if (adjB[y] >> y2) & 1:
                base[Bidx[y]] |= 1 << Bidx[y2]
    ncross = (a - 1) * (b - 1)
    total = 1 << ncross
    nch = 96
    for ch in nb.prange(nch):
        dp = np.zeros((1 << n, n), dtype=np.int64)
        adj = np.zeros(n, dtype=np.int64)
        for mask in range(ch, total, nch):
            for t in range(n):
                adj[t] = base[t]
            for x in range(a - 1):
                for y in range(b - 1):
                    if (mask >> (x * (b - 1) + y)) & 1:
                        adj[x] |= 1 << (a + y)
                    else:
                        adj[a + y] |= 1 << x
            out[mask] = ham_buf(n, adj, dp)


def H(adj):
    return int(ham_one(len(adj), np.array(adj, dtype=np.int64)))


chains = [c for t in range(2, 6) for c in comps(t, (2, 3, 4, 5))]
extra = [((2,), (3, 3)), ((3, 3), (2,)), ((3,), (3, 3)), ((3, 3), (3,))]
below = 0; mixed = 0; nonconcat = []; exact_total = 0; pairs = 0; missing = []
todo = [(cA, cB) for cA in chains for cB in chains if (len(chain(cA)) - 1) * (len(chain(cB)) - 1) <= 16] + extra
for cA, cB in todo:
    A = chain(cA); B = chain(cB)
    a, b = len(A), len(B)
    pairs += 1
    hA, hB = H(A), H(B)
    full = (1 << ((a - 1) * (b - 1))) - 1
    seen_ab = seen_ba = False
    for i in range(a):
        for j in range(b):
            out = np.zeros(full + 1, dtype=np.int64)
            glue_all(np.array(A, dtype=np.int64), np.array(B, dtype=np.int64), i, j, out)
            below += int(np.sum(out < hA * hB))
            for m in np.where(out == hA * hB)[0]:
                exact_total += 1
                if m == full:          # every cross arc A' -> B'
                    ok = (i == a - 1 or cA == (2,)) and (j == 0 or cB == (2,))
                    seen_ab |= (i == a - 1 and j == 0)
                elif m == 0:           # every cross arc B' -> A'
                    ok = (j == b - 1 or cB == (2,)) and (i == 0 or cA == (2,))
                    seen_ba |= (j == b - 1 and i == 0)
                else:
                    mixed += 1
                    ok = False
                if not ok:
                    nonconcat.append((cA, cB, i, j, int(m)))
    if not (seen_ab and seen_ba):
        missing.append((cA, cB))
    print(f'     A = {cA}, B = {cB}: done [{time.time() - T0:.0f}s]', flush=True)
print(f'     {pairs} ordered chain pairs, {exact_total} exact gluings')
check(below == 0, "(X1) Busch's merge inequality H(T) >= H(A)H(B) for every gluing tested")
check(mixed == 0, '(X2) one-way lemma: every exact gluing has all cross arcs in one direction')
check(not nonconcat, '(X3) every exact gluing is an end-to-start concatenation (up to the rotations of a 3-cycle factor)'
      + ('' if not nonconcat else f'; exceptions {nonconcat[:5]}'))
check(not missing, '(X4) both concatenations of every pair are exact (consistent with the product law)')

print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f' ({time.time() - T0:.0f}s)')
