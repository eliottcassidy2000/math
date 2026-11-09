"""AUDIT2 for THM-5042: rebuild the proof's witnesses from scratch and verify them with independent code.

Base: my own exhaustive census (a2_census.c) for n <= 7 (witness masks of odd graphs on exactly n vertices).
Induction (THM-5042 text), for n >= 8:
  k <= 2^(n-1) - 1           : pad a witness from n - 1
  k = 2^(n-1) + x:
    x >= 2^(n-2) + 1         : Lemma A, K_1 v (K_1 v T), T threshold on n - 2 vertices with d(T) = x - 2^(n-2)
    x <= 2^(n-2) - 1, x not in {1,5,37} : Lemma B, K_1 v (odd witness of x, padded to n - 1)
    x in {1, 5, 37}          : Lemma C, the THM-5042 seeds at m = 8, 9, then SR at the tracked leaf
Every new value at orders 8..NMAX is built and checked by a2_check.exe: exact order, d (subset DP) and an odd
automorphism found by an independent backtracking search.  A random subset is re-checked with networkx VF2 and a
direct subset-by-subset domination count.  The Lemma C chain is checked to order 16, including that the tracked leaf
q is fixed by an odd automorphism (search restricted to the stabiliser of q) and that there is no universal vertex.
usage: python a2_thm5042.py NMAX
"""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a2_lib import *

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 11
EXC = (1, 5, 37)

def census_witnesses(n):
    pairs = [(u, v) for u in range(n) for v in range(u + 1, n)]
    W = {}
    for line in open(os.path.join(HERE, 'a2_census_n%d.out' % n)):
        f = line.split()
        if f and f[0] == 'V':
            k = int(f[1]); w = int(f[12])
            if w >= 0:
                W[k] = from_edges(n, [pairs[e] for e in range(len(pairs)) if w >> e & 1])
    return W

BASE = {n: census_witnesses(n) for n in range(1, 8)}
SEEDS = {(8, 1): (3, '0102030406071516172757'), (8, 5): (3, '01020304051617262767'),
         (8, 37): (5, '010203041526273637464767'), (9, 1): (4, '0102030405071618273868'),
         (9, 5): (6, '0102030405071627384858'), (9, 37): (6, '010203040507081617182838485878')}

def seed(m, c):
    if m in (8, 9):
        q, s = SEEDS[(m, c)]
        return from_hex_pairs(m, s), q
    g, q = seed(m - 2, c)
    return SR(g, q)

memo = {}
def build(n, k):
    """odd graph on exactly n vertices with d = k; also returns the lemma used"""
    if (n, k) in memo: return memo[(n, k)]
    assert k % 2 == 1 and k <= 2**n - 1 and k not in EXC
    if k <= 2**(n - 1) - 1:
        g, how = build(n - 1, k); res = (pad(g, n), 'pad')
    elif n <= 7:
        res = (BASE[n][k], 'census')
    else:
        x = k - 2**(n - 1)
        if x >= 2**(n - 2) + 1:
            res = (cone(cone(threshold(x - 2**(n - 2), n - 2))), 'A')
        elif x not in EXC:
            g, _ = build(n - 2, x)        # x <= 2^(n-2) - 1: a witness on n - 2 vertices exists
            res = (cone(pad(g, n - 1)), 'B')
        else:
            res = (seed(n, x)[0], 'C')
    memo[(n, k)] = res
    return res

ok = True
# 1. seeds
items = []
for (m, c), (q, s) in sorted(SEEDS.items()):
    g = from_hex_pairs(m, s)
    items.append((g, q))
res = run_check(items)
for ((m, c), (q, s)), (g, _), r in zip(sorted(SEEDS.items()), items, res):
    db = d_brute(g)
    nxodd, cnt = nx_sign_odd(g, fix=q)
    deg_q = bin(g[1][q]).count('1')
    good = (r[0] == db == 2**(m - 1) + c and r[1] == 'ODD' and r[2] == 'QODD' and r[3] == 0 and deg_q == 1 and nxodd)
    ok &= good
    print("seed m=%d d=%d (2^%d+%d): brute d=%d, C %s/%s, universal=%d, deg(q=%d)=%d, networkx odd aut fixing q: %s -> %s"
          % (m, r[0], m - 1, c, db, r[1], r[2], r[3], q, deg_q, nxodd, 'OK' if good else 'BAD'))

# 2. Lemma C chain to order 16
items = []; labels = []
for c in EXC:
    for m in range(8, 17):
        g, q = seed(m, c); items.append((g, q)); labels.append((m, c))
res = run_check(items)
badC = 0
for (m, c), (g, q), r in zip(labels, items, res):
    good = r[0] == 2**(m - 1) + c and r[1] == 'ODD' and r[2] == 'QODD' and r[3] == 0 and bin(g[1][q]).count('1') == 1 and g[0] == m
    badC += not good
    if not good: print("Lemma C chain FAIL", m, c, r)
print("Lemma C chain m = 8..16, c in {1,5,37}: %d graphs, %d failures" % (len(items), badC))
ok &= badC == 0

# 3. every new value at each order 8..NMAX
rng = random.Random(2026)
for n in range(8, NMAX + 1):
    ks = [k for k in range(2**(n - 1) + 1, 2**n, 2)]
    items = []; hows = []
    for k in ks:
        g, how = build(n, k); items.append((g, -1)); hows.append(how)
    res = run_check(items)
    fails = [(k, r) for k, (g, _), r in zip(ks, items, res) if not (g[0] == n and r[0] == k and r[1] == 'ODD')]
    cnt = {h: hows.count(h) for h in set(hows)}
    # independent re-check of a random subset (networkx VF2 + direct domination count)
    sub = rng.sample(range(len(ks)), min(12, len(ks))) if n <= 10 else []
    nxbad = 0
    for i in sub:
        g = items[i][0]
        if d_brute(g) != ks[i] or not nx_sign_odd(g, limit=2 * 10**6)[0]: nxbad += 1
    ok &= not fails and nxbad == 0
    print("n=%d: %d new values built (%s); failures %d; networkx/direct re-check of %d: %d bad" % (n, len(ks), cnt, len(fails), len(sub), nxbad))
    for k, r in fails[:10]: print("   FAIL", k, r)
print("ALL:", "PASS" if ok else "FAIL")
