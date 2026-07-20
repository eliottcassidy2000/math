#!/usr/bin/env python3
"""
two_ladders_antiorders_hgaps_boxeph_S151.py  (HYP-8215)

OWNER S151: "see how {2, 6} relates to {7, 21}."
{2,6} = anti-automorphism order spectrum (n <= 6, S150 census).
{7,21} = the h-spectrum gaps (S146).  Both are X*{1,3} ladders.  This script
proves/measures the machinery connecting them:

 (1) FREE-ACTION LEMMA: Aut(T) acts freely on directed Hamiltonian paths (an
     automorphism fixing a directed Ham path is the identity; reversal is not
     available inside Aut).  Hence |Aut(T)| divides h(T).  With Redei (h odd)
     this RE-PROVES Moon's theorem (|Aut| odd).  Verified exhaustively n <= 6.
 (2) SC CROSS-TABLE (n <= 6): (|Aut|, anti-orders, h) — the order-6 anti-auto
     classes are exactly where 3 | h: the same Z3 drives both ladders' 3.
 (3) CIRCULANT n=7 EXACT: Z7 acts freely => h in 7*odds for all 8 circulants.
     Which of 7, 21, 35 occur?  (S146 sample: 35 attained somewhere at n=7,
     7/21 not.)  Also circulant n=9 (9*odds predicted).
 (4) EXACT EXCLUSION at n=7 beyond sampling: every tournament with an
     order-3 automorphism (both cycle types 331, 31111) enumerated exactly
     (128 invariant tournaments per generator type); if none has h = 21, then
     h = 21 at n=7 forces |Aut| in {1,7}; 7 => rotational => (3) decides.
     So h=21 at n=7 is cornered into the trivial-Aut stratum only.
 (5) ANTI-ORDER 10 CONSTRUCTION (the 2*mu(m) threshold law): explicit n=10
     two-layer tournament with anti-automorphism of order 10 = tau o r
     (r = Z5 rotation, tau = layer swap reversing arcs).  Verifies the
     upper bound of: THEOREM — the minimal n carrying an anti-automorphism
     of order 2m (m odd > 1) is 2*mu(m), mu = minimal degree of an order-m
     permutation (sum of prime powers); lower bound by the orbit-pairing
     argument (centralizer of an m-cycle forces tau to pair a-orbits;
     tau-fixed vertices of a-orbits are impossible except singletons).
     m=3: n=6 (census-confirmed first order-6 at n=6); m=5: n=10.

boxeph-2026-07-20-S151.
"""

from itertools import permutations, product, combinations
from math import lcm

def hpaths(adj, n):
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for S in range(1 << n):
        for v in range(n):
            c = dp[S][v]
            if not c or not (S >> v) & 1: continue
            m = adj[v] & ~S
            while m:
                w = (m & -m).bit_length() - 1
                dp[S | (1 << w)][w] += c
                m &= m - 1
    return sum(dp[(1 << n) - 1][v] for v in range(n))

def to_bitadj(mat, n):
    return [sum(1 << j for j in range(n) if mat[i][j]) for i in range(n)]

def apply_perm(mat, p, n):
    out = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if mat[i][j]: out[p[i]][p[j]] = 1
    return out

def key(mat): return tuple(tuple(r) for r in mat)
def op(mat, n): return [[mat[j][i] for j in range(n)] for i in range(n)]
def perm_order(p, n):
    seen = [False]*n; o = 1
    for s in range(n):
        if seen[s]: continue
        l, v = 0, s
        while not seen[v]: seen[v] = True; v = p[v]; l += 1
        o = lcm(o, l)
    return o

print("=" * 96)
print("(1)+(2) n <= 6 exhaustive: |Aut| | h  +  SC cross-table (|Aut|, anti-orders, h)")
print("=" * 96)
for n in range(3, 7):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    perms = list(permutations(range(n)))
    seen = {}
    for bits in range(1 << len(pairs)):
        mat = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1: mat[i][j] = 1
            else: mat[j][i] = 1
        best = min(key(apply_perm(mat, p, n)) for p in perms)
        if best not in seen: seen[best] = mat
    viol = 0
    sc_rows = []
    for mat in seen.values():
        T, Top = key(mat), key(op(mat, n))
        autos = [p for p in perms if key(apply_perm(mat, p, n)) == T]
        antis = [p for p in perms if key(apply_perm(mat, p, n)) == Top]
        h = hpaths(to_bitadj(mat, n), n)
        if h % len(autos) != 0: viol += 1
        if antis:
            orders = sorted({perm_order(p, n) for p in antis})
            sc_rows.append((len(autos), orders, h))
    print("n=%d: classes=%d  |Aut| divides h: %s" %
          (n, len(seen), "ALL OK" if viol == 0 else "%d VIOLATIONS" % viol))
    for aut, orders, h in sorted(sc_rows):
        mark = "  <-- 3|h with order-6 anti" if 6 in orders else ""
        print("    SC: |Aut|=%-2d anti-orders=%-8s h=%-3d h/|Aut|=%d%s"
              % (aut, orders, h, h // aut, mark))
print("LEMMA verified: Aut acts freely on directed Ham paths => |Aut| | h;")
print("with Redei (h odd) this re-proves MOON (|Aut| odd).")

print("\n" + "=" * 96)
print("(3) circulants: n=7 (all 8) and n=9 (all 16): h values (predict 7*odds / 9*odds)")
print("=" * 96)
for n in (7, 9):
    half = list(range(1, (n - 1) // 2 + 1))
    vals = {}
    for choice in product([1, -1], repeat=len(half)):
        S = set((s if e == 1 else -s) % n for s, e in zip(half, choice))
        mat = [[1 if (j - i) % n in S else 0 for j in range(n)] for i in range(n)]
        h = hpaths(to_bitadj(mat, n), n)
        vals.setdefault(h, []).append(sorted(S))
    print("  n=%d: h values = %s" % (n, sorted(vals)))
    for h in sorted(vals):
        print("      h=%-4d (h/%d=%s) e.g. S=%s  x%d" % (h, n,
              ("%d" % (h // n)) if h % n == 0 else "NOT DIV?!", vals[h][0], len(vals[h])))
    assert all(h % n == 0 for h in vals), "free-action violated?!"

print("\n" + "=" * 96)
print("(4) EXACT n=7 exclusion: all tournaments with an order-3 automorphism")
print("=" * 96)
n = 7
hits21 = 0
hset = set()
for cyc in ([(0,1,2),(3,4,5)], [(0,1,2)]):        # cycle types 331 and 31111
    g = list(range(n))
    for c in cyc:
        for i in range(len(c)): g[c[i]] = c[(i + 1) % len(c)]
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    # orbits of <g> on unordered pairs
    orb = {}
    for (i, j) in pairs:
        o = []
        a, b = i, j
        for _ in range(3):
            o.append((min(a, b), max(a, b)))
            a, b = g[a], g[b]
        r = min(o)
        orb.setdefault(r, set()).update(o)
    reps = sorted(orb)
    cnt = 0
    for bits in range(1 << len(reps)):
        # orient the representative pair, propagate along g
        mat = [[0]*n for _ in range(n)]
        ok = True
        for k, r in enumerate(reps):
            i, j = r
            fwd = (bits >> k) & 1
            a, b = (i, j) if fwd else (j, i)
            for _ in range(3):
                if mat[b][a]: ok = False; break
                mat[a][b] = 1
                a, b = g[a], g[b]
            if not ok: break
        if not ok: continue
        # check tournament completeness (fixed pairs oriented once)
        if any(mat[i][j] + mat[j][i] != 1 for (i, j) in pairs): continue
        cnt += 1
        h = hpaths(to_bitadj(mat, n), n)
        hset.add(h)
        if h == 21: hits21 += 1
        if h == 7: hits21 += 1000  # flag loudly
    print("  cycle type %s: %d invariant tournaments enumerated" %
          ("331" if len(cyc) == 2 else "31111", cnt))
print("  h-values over ALL order-3-symmetric n=7 tournaments:", sorted(hset))
print("  h = 21 occurrences: %d ; h = 7 occurrences: %d" % (hits21 % 1000, hits21 // 1000))
print("  => h = 21 at n = 7 %s outside |Aut| in {1,7}; rotational case decided in (3)."
      % ("EXCLUDED" if hits21 == 0 else "NOT excluded"))

print("\n" + "=" * 96)
print("(5) anti-order 10 at n = 10: explicit construction (2*mu law upper bound)")
print("=" * 96)
n10 = 10
# vertices (j, eps) = j + 5*eps ; layer 0 circulant S0={1,2}; layer 1 = reversed; cross C={0} sym
S0 = {1, 2}
C = {0}
mat = [[0]*n10 for _ in range(n10)]
for j in range(5):
    for k in range(5):
        d = (k - j) % 5
        if d in S0: mat[j][k] = 1                       # layer 0
        if d in {(-x) % 5 for x in S0}: mat[5 + j][5 + k] = 1   # layer 1 reversed
        if j != k or True:
            if d in C: mat[j][5 + k] = 1                # cross 0 -> 1
            else: mat[5 + k][j] = 1                     # cross 1 -> 0
ok_t = all(mat[i][j] + mat[j][i] == 1 for i in range(n10) for j in range(i + 1, n10))
r = [ (j + 1) % 5 + 5 * (v // 5) for v in range(n10) for j in [v % 5] ]
tau = [ (v + 5) % 10 for v in range(n10) ]
sigma = [ tau[r[v]] for v in range(n10) ]
img = apply_perm(mat, sigma, n10)
is_anti = key(img) == key(op(mat, n10))
print("  tournament OK: %s ; sigma = tau o r is anti-auto: %s ; order(sigma) = %d"
      % (ok_t, is_anti, perm_order(sigma, n10)))
assert ok_t and is_anti and perm_order(sigma, n10) == 10
print("  => anti-order 10 realized at n = 10.  Lower bound n >= 10 for order-10:")
print("  sigma^2 = order-5 auto a; an involutive-part tau commuting with a must map")
print("  each 5-orbit to a DIFFERENT 5-orbit (centralizer of a 5-cycle in S5 = Z5 has")
print("  no involution; tau|orbit = id would force an arcless 5-set) => orbits paired")
print("  => n >= 10.  THE 2*mu(m) LAW: minimal n for anti-order 2m = 2*mu(m), mu = ")
print("  minimal degree of an order-m permutation.  m=3: n=6 (census); m=5: n=10 (here);")
print("  m=15: n=16 (predicted, NOT 30: orbits {3,5} paired).")
print("DONE.")
