"""audit_side_checks.py -- S10 auditor side checks for THM-4868 / HYP-9460 (pure arithmetic and small combinatorics).

1. HYP-9460 P7: recompute the recurrence rows N(n, r) to n = 18 from N(4,.) = (1, 0), N(5,.) = (1, 1, 1); compare with
   the rows printed in HYP-9460 (orders 14-18) and with the note's section-4 table (orders 5-13); row sums vs F(2n-6);
   the generating function t^2 x^4 (1 + (t^2 - t) x) / (1 - (1 + 2t) x + t x^2) expanded as a power series;
   the fitted columns / top diagonals claimed to follow (n <= 40).
2. The (F) formulas F(2n-10)+1, F(2n-9)-1, F(2n-8)+1 sum to F(2n-6); a + 2b = F(2n-5); HYP-9460 P5 values.
3. Proof of (B): index bookkeeping of the even-n second cycle for every even n in [6, 200]; and every round
   tournament (N+(u) = {u+1, ..., u+k(u)} along a Hamiltonian cycle) with all k(u) >= 2, n = 5..13: the constructed
   walk is a Hamiltonian cycle different from C, and a vertex with k >= 3 exists at even n.
Run: python3 -I audit_side_checks.py
"""
import itertools, sys, time

t0 = time.time()
def F(k):
    if k < 0:
        return (-1) ** (k + 1) * F(-k)
    a, b = 0, 1
    for _ in range(k):
        a, b = b, a + b
    return a

ok_all = True
def report(name, ok):
    global ok_all
    ok_all &= ok
    print(f"{name}: {'OK' if ok else 'FAIL'}")

# ---------------- 1. P7 rows ----------------
N = {4: {2: 1, 3: 0}, 5: {2: 1, 3: 1, 4: 1}}
for n in range(6, 41):
    N[n] = {r: N[n - 1].get(r, 0) + 2 * N[n - 1].get(r - 1, 0) - N[n - 2].get(r - 1, 0) for r in range(2, n)}
hyp_rows = {
    14: [1, 10, 55, 210, 606, 1372, 2471, 3530, 3920, 3232, 1792, 512],
    15: [1, 11, 66, 275, 870, 2178, 4396, 7183, 9420, 9680, 7424, 3840, 1024],
    16: [1, 12, 78, 352, 1210, 3312, 7380, 13504, 20256, 24600, 23552, 16896, 8192, 2048],
    17: [1, 13, 91, 442, 1639, 4862, 11826, 23868, 40081, 55692, 63072, 56576, 38144, 17408, 4096],
    18: [1, 14, 105, 546, 2171, 6930, 18238, 40140, 74313, 115598, 149856, 159168, 134400, 85504, 36864, 8192],
}
note_rows = {
    5: [1, 1, 1], 6: [1, 2, 3, 2], 7: [1, 3, 6, 7, 4], 8: [1, 4, 10, 16, 16, 8], 9: [1, 5, 15, 30, 41, 36, 16],
    10: [1, 6, 21, 50, 85, 102, 80, 32], 11: [1, 7, 28, 77, 155, 231, 248, 176, 64],
    12: [1, 8, 36, 112, 259, 456, 608, 592, 384, 128], 13: [1, 9, 45, 156, 406, 819, 1289, 1560, 1392, 832, 256],
}
report("P7 recurrence reproduces the note's section-4 rows (5-13)",
       all([N[n][r] for r in range(2, n)] == row for n, row in note_rows.items()))
report("P7 recurrence reproduces HYP-9460's predicted rows (14-18)",
       all([N[n][r] for r in range(2, n)] == row for n, row in hyp_rows.items()))
report("row sums = F(2n-6), n = 4..40", all(sum(N[n].values()) == F(2 * n - 6) for n in range(4, 41)))
s16 = [2048, 8192, 16896, 23552, 24600, 20256, 13504, 7380, 3312, 1210, 352, 78, 12, 1]
report("order-16 recurrence row reversed = S74's published T(16,k) row", [N[16][17 - k] for k in range(2, 16)] == s16)
# generating function: series of t^2 x^4 (1 + (t^2 - t) x) / (1 - (1 + 2t) x + t x^2), polynomials in t as dicts
def padd(a, b, c=1):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + c * v
    return {k: v for k, v in out.items() if v}
def pmul(a, b):
    out = {}
    for i, u in a.items():
        for j, v in b.items():
            out[i + j] = out.get(i + j, 0) + u * v
    return {k: v for k, v in out.items() if v}
num = {4: {2: 1}, 5: {4: 1, 3: -1}}          # coefficient of x^n in the numerator
G = {}
for n in range(0, 41):                       # G_n = num_n + (1 + 2t) G_{n-1} - t G_{n-2}
    g = dict(num.get(n, {}))
    if n - 1 in G:
        g = padd(g, pmul({0: 1, 1: 2}, G[n - 1]))
    if n - 2 in G:
        g = padd(g, pmul({1: 1}, G[n - 2]), -1)
    G[n] = g
report("generating function coefficients = N(n, r), n = 4..40",
       all(G[n] == {r: v for r, v in N[n].items() if v} for n in range(4, 41)) and all(not G[n] for n in range(0, 4)))
from math import comb
report("columns T(n,2)=2^(n-5) (n>=5), T(n,3)=n 2^(n-7) (n>=6), T(n,4)=(n^2+n-8) 2^(n-10) (n>=7), n <= 40",
       all(N[n][n - 1] * 2 ** 0 == 2 ** (n - 5) for n in range(5, 41))
       and all(N[n][n - 2] * 2 ** 7 == n * 2 ** n for n in range(6, 41))
       and all(N[n][n - 3] * 2 ** 10 == (n * n + n - 8) * 2 ** n for n in range(7, 41)))
report("top diagonals N(n,2)=1, N(n,3)=n-4, N(n,4)=C(n-3,2), n <= 40",
       all(N[n][2] == 1 for n in range(4, 41)) and all(N[n][3] == n - 4 for n in range(5, 41))
       and all(N[n][4] == comb(n - 3, 2) for n in range(6, 41)))

# ---------------- 2. (F) arithmetic ----------------
report("(F) table sums to F(2n-6) and a + 2b = F(2n-5), n = 5..60",
       all((F(2*n-10)+1) + 2*(F(2*n-9)-1) + (F(2*n-8)+1) == F(2*n-6) and F(2*n-8) + 2*F(2*n-7) == F(2*n-5) for n in range(5, 61)))
p5 = {14: (2585, 4180, 6766), 15: (6766, 10945, 17712), 16: (17712, 28656, 46369), 17: (46369, 75024, 121394), 18: (121394, 196417, 317812)}
report("HYP-9460 P5 values = formulas", all(p5[n] == (F(2*n-10)+1, F(2*n-9)-1, F(2*n-8)+1) for n in p5))
p1 = [17711, 46368, 121393, 317811, 832040]
report("HYP-9460 P1 values = F(2n-6), n = 14..18", p1 == [F(2 * n - 6) for n in range(14, 19)])

# ---------------- 3. (B) second-cycle bookkeeping ----------------
def even_walk(n):  # indices relative to u_0 = w+++ ; w = u_(n-3)
    return list(range(0, n - 1, 2)) + [n - 1] + list(range(1, n - 2, 2))
ok = True
for n in range(6, 201, 2):
    W = even_walk(n)
    steps = [(W[(i + 1) % n] - W[i]) % n for i in range(n)]
    if sorted(W) != list(range(n)) or steps.count(2) != n - 2 or steps.count(1) != 1 or steps.count(3) != 1:
        ok = False
    if (W[n - 1], W[0]) != (n - 3, 0) or steps[W.index(n - 2)] != 1:
        ok = False
report("even-n walk: permutation of Z_n with steps +2 (n-2 times), +1 at u_(n-2), +3 at u_(n-3) -> u_0, n = 6..200", ok)

def round_tournaments(n):
    """all k in [1, n-2]^n such that u -> u+d (d <= k(u)) defines a tournament (backtracking)."""
    k = [0] * n
    def consistent(u, v):
        d = (v - u) % n
        return (d <= k[u]) != ((n - d) <= k[v])
    def rec(i):
        if i == n:
            yield tuple(k)
            return
        for val in range(1, n - 1):
            k[i] = val
            if all(consistent(j, i) for j in range(i)):
                yield from rec(i + 1)
    yield from rec(0)

def hc_count(adj, n):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    dp[1][0] = 1
    for m in range(1, 1 << n, 2):
        for v in range(n):
            c = dp[m][v]
            if c:
                for w in range(n):
                    if adj[v] >> w & 1 and not m >> w & 1:
                        dp[m | 1 << w][w] += c
    return sum(dp[full][v] for v in range(1, n) if adj[v] & 1)

for n in range(5, 14):
    cnt = bad = 0; minhc = None
    for k in round_tournaments(n):
        if min(k) < 2:
            continue
        cnt += 1
        adj = [sum(1 << ((u + d) % n) for d in range(1, k[u] + 1)) for u in range(n)]
        if n % 2:
            walk = [(2 * i) % n for i in range(n)]
        else:
            ws = [u for u in range(n) if k[u] >= 3]
            if not ws:
                bad += 1; continue
            w = ws[0]; u0 = (w + 3) % n
            walk = [(u0 + i) % n for i in even_walk(n)]
        good = sorted(walk) == list(range(n)) and all(adj[walk[i]] >> walk[(i + 1) % n] & 1 for i in range(n))
        good &= any((walk[(i + 1) % n] - walk[i]) % n != 1 for i in range(n))
        if not good:
            bad += 1
        if n <= 11:
            h = hc_count(adj, n)
            minhc = h if minhc is None else min(minhc, h)
    print(f"   round tournaments n={n} with all k(u) >= 2: {cnt}; walk not a second Hamiltonian cycle: {bad}"
          + (f"; min hc = {minhc}" if minhc is not None else ""))
    ok_all &= (bad == 0)

print("ALL OK" if ok_all else "SOME CHECK FAILED", f"({time.time() - t0:.1f} s)")
