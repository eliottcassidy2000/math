"""Auditor's check of the proof of HYP-9512 (D): at odd order n = 2m + 1 >= 3 the self-converse classes split exactly
evenly between H = 1 and H = 3 (mod 4).

Model: vertices 0 (the fixed point v0) and pairs (a_i, b_i) = (2i - 1, 2i), i = 1..m; sigma swaps a_i <-> b_i.
X = tournaments with sigma as an anti-automorphism (u -> w  iff  sigma(w) -> sigma(u)); |X| = 2^(m(m+1)).
P(x) = number of 'pair-distinct' directed paths v0 -> x_1 -> ... -> x_k = x (k >= 1, at most one vertex per pair).
O = {i : P(a_i) + P(b_i) odd};  Phi = reverse the internal arc of every pair in O.
Claims checked: (i) S = sum_x P(x) is odd; (ii) H(T) = 1 + 2F (mod 4) with F = sum_i P(head_i), head_i the tail of the
internal arc of pair i (i.e. the x with x -> sigma(x)); (iii) H(Phi T) = H(T) + 2 (mod 4); (iv) Phi is an involution.
n = 3, 5, 7: all of X; n = 9, 11: random samples.  usage: python3 aud_sc_mod4.py [samples]
"""
import sys, random, itertools

def build(m, bits):
    """bits: dict with keys ('v', i), ('p', i), ('x', i, j), ('y', i, j) for i < j; returns out-neighbour sets."""
    n = 2 * m + 1
    adj = [[False] * n for _ in range(n)]
    def arc(u, w):
        adj[u][w] = True
    sig = lambda u: 0 if u == 0 else (u + 1 if u % 2 == 1 else u - 1)
    def setpair(u, w):           # u -> w and its sigma-image sigma(w) -> sigma(u)
        arc(u, w); arc(sig(w), sig(u))
    for i in range(1, m + 1):
        a, b = 2 * i - 1, 2 * i
        if bits[('v', i)]: setpair(0, a)
        else: setpair(a, 0)
        if bits[('p', i)]: arc(a, b)
        else: arc(b, a)
        for j in range(i + 1, m + 1):
            c, d = 2 * j - 1, 2 * j
            if bits[('x', i, j)]: setpair(a, c)
            else: setpair(c, a)
            if bits[('y', i, j)]: setpair(a, d)
            else: setpair(d, a)
    for u in range(n):
        for w in range(n):
            if u != w: assert adj[u][w] != adj[w][u]
    return adj

def hpaths(adj):
    n = len(adj); F = 1 << n
    f = [[0] * n for _ in range(F)]
    for v in range(n): f[1 << v][v] = 1
    for S in range(1, F):
        for v in range(n):
            x = f[S][v]
            if x:
                for u in range(n):
                    if not (S >> u) & 1 and adj[v][u]: f[S | 1 << u][u] += x
    return sum(f[F - 1])

def pcounts(adj, m):
    """P(x) for every x != 0: pair-distinct paths from 0 ending at x (k >= 1)."""
    pair = lambda u: (u + 1) // 2
    from collections import defaultdict
    cur = {}
    for x in range(1, 2 * m + 1):
        if adj[0][x]: cur[(1 << pair(x), x)] = 1
    P = defaultdict(int)
    while cur:
        nxt = defaultdict(int)
        for (mask, x), c in cur.items():
            P[x] += c
            for y in range(1, 2 * m + 1):
                if not (mask >> pair(y)) & 1 and adj[x][y]:
                    nxt[(mask | 1 << pair(y), y)] += c
        cur = nxt
    return P

def check(m, bits):
    adj = build(m, bits)
    H = hpaths(adj); P = pcounts(adj, m)
    S = sum(P.values())
    F = 0
    for i in range(1, m + 1):
        a, b = 2 * i - 1, 2 * i
        F += P[a] if adj[a][b] else P[b]
    O = [i for i in range(1, m + 1) if (P[2 * i - 1] + P[2 * i]) % 2]
    bits2 = dict(bits)
    for i in O: bits2[('p', i)] ^= 1
    adj2 = build(m, bits2); H2 = hpaths(adj2); P2 = pcounts(adj2, m)
    O2 = [i for i in range(1, m + 1) if (P2[2 * i - 1] + P2[2 * i]) % 2]
    return (S % 2 == 1, H % 4 == (1 + 2 * F) % 4, (H2 - H) % 4 == 2, O2 == O, H % 4)

def keys(m):
    ks = [('v', i) for i in range(1, m + 1)] + [('p', i) for i in range(1, m + 1)]
    for i in range(1, m + 1):
        for j in range(i + 1, m + 1):
            ks += [('x', i, j), ('y', i, j)]
    return ks

samples = int(sys.argv[1]) if len(sys.argv) > 1 else 300
random.seed(20261107)
for m in (1, 2, 3, 4, 5):
    ks = keys(m); tot = 0; good = [0, 0, 0, 0]; res = {1: 0, 3: 0}
    if m <= 3:
        it = (dict(zip(ks, v)) for v in itertools.product((0, 1), repeat=len(ks)))
    else:
        it = ({k: random.randint(0, 1) for k in ks} for _ in range(samples if m == 4 else max(30, samples // 10)))
    for bits in it:
        r = check(m, bits); tot += 1
        for t in range(4): good[t] += r[t]
        res[r[4]] += 1
    print('n=%d (%s, %d tournaments in X): S odd %d, H=1+2F mod 4 %d, Phi flips H mod 4 %d, Phi involution %d; labelled H mod 4: %s'
          % (2 * m + 1, 'all of X' if m <= 3 else 'random sample', tot, good[0], good[1], good[2], good[3], res))
