"""Extract the census base used by the S89 proof from the output of dom_even (THM-4996, Lemma B3).

Input : even511.txt = lines of `dom_even 511 16` with parity E (awk '$4=="E"'), format "n m d E I2 edges(hex pairs)".
Output: base_witnesses.txt  - W lines: a minimal-order even realiser of every value in E(9);
                              L lines: leaf graphs without universal vertex on exactly m in {8, 9} vertices
                                       for the top band T(m) and the offsets 2^(m-1) + {1,3,5,7,11,13}
        base_sets.json      - E(n), E0(n) (no universal vertex), C(n) (complement connected), L(n) for n <= 9
"""
import json, sys
from collections import defaultdict
src = sys.argv[1] if len(sys.argv) > 1 else 'even511.txt'

def parse(line):
    f = line.split(); n = int(f[0]); d = int(f[2]); e = f[5] if len(f) > 5 else ''
    adj = [0] * n
    for i in range(0, len(e), 2):
        u, v = int(e[i], 16), int(e[i + 1], 16); adj[u] |= 1 << v; adj[v] |= 1 << u
    return n, d, adj
def edges_str(n, adj):
    return ''.join('%x%x' % (u, v) for u in range(n) for v in range(u + 1, n) if adj[u] >> v & 1)
def comp_connected(n, adj):
    full = (1 << n) - 1; cadj = [(full ^ adj[v]) & ~(1 << v) for v in range(n)]
    seen = 1; st = [0]
    while st:
        v = st.pop(); nb = cadj[v] & ~seen
        while nb:
            w = (nb & -nb).bit_length() - 1; nb &= nb - 1; seen |= 1 << w; st.append(w)
    return seen == full
def union(gs):
    n = 0; adj = []
    for (m, a) in gs:
        adj += [x << n for x in a]; n += m
    return n, adj

conn = defaultdict(dict)        # value -> {n: adj} first realiser per order
conn_leaf = defaultdict(dict)   # value -> {n: adj} with a leaf
conn_leaf_nouniv = defaultdict(dict)
sets_conn = defaultdict(set); sets_conn0 = defaultdict(set); sets_connC = defaultdict(set); sets_connL = defaultdict(set)
leafany = defaultdict(set)
for line in open(src):
    n, d, adj = parse(line)
    if n > 9: continue
    deg = [bin(a).count('1') for a in adj]
    conn[d].setdefault(n, adj)
    sets_conn[n].add(d)
    univ = n >= 2 and max(deg) == n - 1
    if n == 1 or not univ: sets_conn0[n].add(d)
    if n == 1 or comp_connected(n, adj): sets_connC[n].add(d)
    if 1 in deg:
        leafany[n].add(d); conn_leaf[d].setdefault(n, adj)
        if not univ: sets_connL[n].add(d); conn_leaf_nouniv[d].setdefault(n, adj)

# value sets for n <= 9 (unions of components)
NM = 9
tab = [defaultdict(set) for _ in range(NM + 1)]   # tab[s][(ncomp>=2?, hasleafcomp?)] values
tab[0][(0, False)] = {1}
for s in range(1, NM + 1):
    for size in range(1, s + 1):
        for (c, hl), vals in list(tab[s - size].items()):
            for y in sets_conn[size]:
                for x in vals:
                    tab[s][(min(c + 1, 2), hl or (y in leafany[size]))].add(x * y)
E = {}; E0 = {}; C = {}; L = {}
for n in range(1, NM + 1):
    exact = set().union(*tab[n].values())
    multi = set().union(*[v for (c, hl), v in tab[n].items() if c >= 2])
    multileaf = set().union(*[v for (c, hl), v in tab[n].items() if c >= 2 and hl] or [set()])
    E[n] = set().union(*[set().union(*tab[m].values()) for m in range(1, n + 1)])
    prev = E[n - 1] if n > 1 else set()
    E0[n] = prev | multi | (sets_conn0[n] if n > 1 else set())
    C[n] = prev | multi | (sets_connC[n] if n > 1 else set())
    L[n] = multileaf | sets_connL[n]
json.dump({k: {str(n): sorted(v) for n, v in s.items()} for k, s in (('E', E), ('E0', E0), ('C', C), ('L', L))},
          open('base_sets.json', 'w'))

out = open('base_witnesses.txt', 'w')
out.write('# W n d edges : minimal-order even realiser of d (n <= 9), hex vertex pairs\n')
out.write('# L m d edges : even leaf graph without universal vertex on exactly m vertices\n')
minconn = {v: min(conn[v]) for v in conn}
for v in sorted(E[9]):
    best = None
    if v in minconn: best = (minconn[v], [(minconn[v], conn[v][minconn[v]])])
    for v1 in minconn:
        if v1 > 1 and v % v1 == 0:
            v2 = v // v1
            if v2 in minconn and v2 > 1:
                cand = (minconn[v1] + minconn[v2], [(minconn[v1], conn[v1][minconn[v1]]), (minconn[v2], conn[v2][minconn[v2]])])
                if best is None or cand[0] < best[0]: best = cand
            for v2 in minconn:
                if v2 > 1 and (v // v1) % v2 == 0 and (v // v1 // v2) in minconn and v // v1 // v2 > 1:
                    v3 = v // v1 // v2
                    cand = (minconn[v1] + minconn[v2] + minconn[v3],
                            [(minconn[v1], conn[v1][minconn[v1]]), (minconn[v2], conn[v2][minconn[v2]]), (minconn[v3], conn[v3][minconn[v3]])])
                    if best is None or cand[0] < best[0]: best = cand
    if v == 1: best = (1, [(1, [0])])
    n, adj = union(best[1])
    out.write('W %d %d %s\n' % (n, v, edges_str(n, adj)))
need = {8: {2**7 - j for j in (1, 3, 5, 9)} | {2**7 + c for c in (1, 3, 5, 7, 11, 13)},
        9: {2**8 - j for j in (1, 3, 5, 7)} | {2**8 + c for c in (1, 3, 5, 7, 11, 13)}}
for m, vals in need.items():
    for v in sorted(vals):
        g = None
        if m in conn_leaf_nouniv.get(v, {}): g = (m, conn_leaf_nouniv[v][m])
        elif (m - 1) in conn_leaf.get(v, {}): g = union([(m - 1, conn_leaf[v][m - 1]), (1, [0])])
        if g is None: raise SystemExit("no leaf base for m=%d v=%d" % (m, v))
        out.write('L %d %d %s\n' % (g[0], v, edges_str(*g)))
out.close()
print({n: len(E[n]) for n in E}, 'leaf base written')
