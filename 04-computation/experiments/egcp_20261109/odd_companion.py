"""THM-5042 (odd-graph companion of THM-5040): checks and witnesses.

(1) From the exhaustive census output (dom511.out, `dom_even 511 16`, parity column O/E), the values of odd graphs on
    <= n vertices for n <= 9 are exactly the odd numbers <= 2^n - 1 minus {1, 5, 37}.  (Run with the census file.)
(2) Seeds: odd graphs on 8 and 9 vertices with d = 2^(n-1) + c, c in {1, 5, 37}, a leaf q, no universal vertex, and an
    automorphism fixing q that reverses an odd number of edges (listed below, re-verified here).
(3) Witnesses for the inductive step, built as in the proof, brute-forced for d and oddness at n = 10, 11.

usage: python odd_companion.py [dom511.out]
"""
import sys, random
from egcp_lib import (from_hex, dcount, automorphisms, edge_sign, degs, pad, cone, SR, union, no_universal, has_leaf)

SEEDS = {  # (n, d): (leaf q, hex edges) - from the census
    (8, 129): (3, '0102030406071516172757'),
    (8, 133): (3, '01020304051617262767'),
    (8, 165): (5, '010203041526273637464767'),
    (9, 257): (4, '0102030405071618273868'),
    (9, 261): (6, '0102030405071627384858'),
    (9, 293): (6, '010203040507081617182838485878'),
}

def is_odd(g):
    return any(edge_sign(g, p) == -1 for p in automorphisms(g))

def odd_fixing(g, q):
    return any(p[q] == q and edge_sign(g, p) == -1 for p in automorphisms(g))

def threshold(x, m):
    """a graph on exactly m vertices with d = x (x odd, x <= 2^m - 1): Lemma B1 threshold graphs (parity irrelevant)"""
    if m == 1:
        assert x == 1; return (1, [0])
    if x <= 2**(m - 1) - 1:
        return pad(threshold(x, m - 1), m)
    return cone_any(threshold(x - 2**(m - 1), m - 1))

def cone_any(g):
    n, a = g
    adj = a[:] + [(1 << n) - 1]
    for v in range(n): adj[v] |= 1 << n
    return (n + 1, adj)

def build_odd(n, k, base):
    """odd graph on exactly n vertices with d = k (k odd, k <= 2^n - 1, k not in {1, 5, 37})"""
    if k <= 2**(n - 1) - 1:
        return pad(build_odd(n - 1, k, base), n)
    if n <= 9:
        return pad(base[k], n)
    x = k - 2**(n - 1)
    if x >= 2**(n - 2) + 1:                      # cone over a cone: two adjacent universal vertices
        inner = threshold(x - 2**(n - 2), n - 2)
        return cone_any(cone_any(inner))
    if x not in (1, 5, 37):                      # padded odd graph, coned
        return cone(pad(build_odd(n - 2, x, base), n - 1))
    g = seed(n, x)                               # SR chain from the seeds
    return g

def seed(n, c):
    if n <= 9:
        q, e = SEEDS[(n, 2**(n - 1) + c)]
        return from_hex(n, e)
    g = seed(n - 2, c)
    q = next(v for v in range(g[0]) if degs(g)[v] == 1 and odd_fixing(g, v))
    return SR(g, q)

if __name__ == "__main__":
    ok = True
    for (n, d), (q, e) in SEEDS.items():
        g = from_hex(n, e)
        good = dcount(g) == d and has_leaf(g) and no_universal(g) and degs(g)[q] == 1 and odd_fixing(g, q)
        ok &= good
        print("seed n=%d d=%d leaf=%d: %s" % (n, d, q, "OK" if good else "BAD"))
    base = {}
    if len(sys.argv) <= 1:
        print("census check skipped (pass dom511.out to run it); only the seeds were re-verified")
    if len(sys.argv) > 1:
        from collections import defaultdict
        connE = defaultdict(set); connO = defaultdict(set); conn_graph = {}
        for line in open(sys.argv[1]):
            f = line.split(); n = int(f[0]); d = int(f[2]); par = f[3]
            if n > 9: continue
            (connE if par == 'E' else connO)[n].add(d)
            if par == 'O' and d not in conn_graph: conn_graph[d] = from_hex(n, f[5] if len(f) > 5 else '')
        tab = [{False: set(), True: set()} for _ in range(10)]; tab[0][False] = {1}
        for s in range(1, 10):
            for size in range(1, s + 1):
                for flag in (False, True):
                    for x in tab[s - size][flag]:
                        for y in connE[size]: tab[s][flag].add(x * y)
                        for y in connO[size]: tab[s][True].add(x * y)
        for n in range(2, 10):
            Eo = set().union(*[tab[m][True] for m in range(1, n + 1)])
            target = set(range(1, 2**n, 2)) - {1, 5, 37}
            r = Eo == target; ok &= r
            print("odd-graph values on <= %d vertices == odd <= 2^n - 1 minus {1,5,37}: %s" % (n, r))
        # base witnesses for build_odd: connected odd realisers, or K2 (x) cone products
        for k in range(3, 512, 2):
            if k in (5, 37): continue
            if k in conn_graph: base[k] = conn_graph[k]
            elif k % 3 == 0 and (k // 3) in conn_graph: base[k] = union(conn_graph[k // 3], (2, [2, 1]))
        rng = random.Random(42)
        for n in (10, 11):
            new = [k for k in range(2**(n - 1) + 1, 2**n, 2)]
            sample = sorted(set(rng.sample(new, 25)) | {2**(n - 1) + c for c in (1, 3, 5, 7, 37)} | {2**n - 1, 2**(n - 1) + 2**(n - 2) + 1, 2**(n - 1) + 2**(n - 2) - 1})
            for k in sample:
                try:
                    g = build_odd(n, k, base)
                except KeyError as ex:
                    print("missing base for", ex); ok = False; continue
                good = g[0] == n and dcount(g) == k and is_odd(g)
                ok &= good
                if not good: print("FAIL n=%d k=%d" % (n, k))
            print("n=%d: %d odd witnesses checked" % (n, len(sample)))
    print("ALL:", "PASS" if ok else "FAIL")
