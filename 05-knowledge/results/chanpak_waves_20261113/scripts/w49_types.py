"""Converse types from strong components.
A tournament = word of strong components (ordinal sum); H multiplies; converse = reverse the word and converse each letter.
So the classes with count h at order n = words of strong letters (singleton letter '1' of size 1) with product of H's = h
and total size n; self-converse classes = words fixed by the converse.  Strong classes up to order 7 give the complete
strong spectrum for h <= 44 (strong minima 3, 5, 9, 15, 25, 45 at orders 3..8, THM-1370), hence every type of every
h <= 44 at every order.  Verified against brute force at n = 5, 6, 7."""
import itertools, sys, pickle
sys.path.insert(0, '.')
from w45_hk import H_of
def classes(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]; m = len(pairs)
    idx = {pr: t for t, pr in enumerate(pairs)}; perms = list(itertools.permutations(range(n)))
    seen = bytearray(1 << m); out = []
    def relabel(bits, p):
        b = 0
        for t, (i, j) in enumerate(pairs):
            x = bits >> t & 1; a, c = p[i], p[j]
            if a < c: b |= x << idx[(a, c)]
            else: b |= (1 - x) << idx[(c, a)]
        return b
    for bits in range(1 << m):
        if seen[bits]: continue
        orb = set(relabel(bits, p) for p in perms)
        for o in orb: seen[o] = 1
        canon = min(orb); conv = min(relabel((1 << m) - 1 ^ bits, p) for p in perms)
        adj = [0] * n
        for t, (i, j) in enumerate(pairs):
            if bits >> t & 1: adj[i] |= 1 << j
            else: adj[j] |= 1 << i
        # strong?
        def reach(v, fwd=True):
            r = 1 << v; st = [v]
            while st:
                u = st.pop()
                nb = adj[u] if fwd else sum(1 << w for w in range(n) if adj[w] >> u & 1)
                t2 = nb & ~r
                while t2:
                    w = (t2 & -t2).bit_length() - 1; t2 &= t2 - 1; r |= 1 << w; st.append(w)
            return r
        strong = n == 1 or (reach(0) == (1 << n) - 1 and reach(0, False) == (1 << n) - 1)
        out.append(dict(n=n, canon=canon, conv=conv, H=H_of(n, adj), strong=strong))
    return out
allc = {n: classes(n) for n in range(1, 8)}
strong = [c for n in allc for c in allc[n] if c['strong']]
# letters: (H, size, id, conv_id)
key = {(c['n'], c['canon']): i for i, c in enumerate(strong)}
letters = [(c['H'], c['n'], i, key[(c['n'], c['conv'])]) for i, c in enumerate(strong)]
by_h = {}
for (h, sz, i, ci) in letters:
    by_h.setdefault(h, []).append((sz, i, ci))
print("strong spectrum (h <= 45): h: [(order, self-converse?)]")
for h in sorted(by_h):
    if h <= 45: print("  ", h, [(sz, i == ci) for sz, i, ci in sorted(by_h[h])])
def words(h, n):
    """all words of non-singleton strong letters with product h, plus singleton gaps, total size n"""
    res = []
    def fact(rem, seq):
        if rem == 1:
            yield list(seq); return
        for (hh, sz, i, ci) in letters:
            if hh > 1 and rem % hh == 0:
                seq.append((hh, sz, i, ci)); yield from fact(rem // hh, seq); seq.pop()
    for seq in fact(h, []):
        s = sum(x[1] for x in seq); free = n - s
        if free < 0: continue
        m = len(seq)
        for gaps in itertools.combinations_with_replacement(range(m + 1), free):
            g = [0] * (m + 1)
            for x in gaps: g[x] += 1
            res.append((tuple(x[2] for x in seq), tuple(g)))
    return set(res)
conv_of = {i: ci for (_, _, i, ci) in letters}
def typ(h, n):
    W = words(h, n); s = 0
    for (ls, g) in W:
        cw = (tuple(conv_of[x] for x in reversed(ls)), tuple(reversed(g)))
        if cw == (ls, g): s += 1
    return s, (len(W) - s) // 2
# verification against brute force at n = 5, 6, 7
bad = 0
for n in (5, 6, 7):
    brute = {}
    for c in allc[n]:
        brute.setdefault(c['H'], [0, 0])
        if c['canon'] == c['conv']: brute[c['H']][0] += 1
        else: brute[c['H']][1] += 1
    for h, (s, t2) in brute.items():
        if typ(h, n) != (s, t2 // 2): bad += 1
print("word formula vs brute force at n = 5, 6, 7: mismatches", bad)
# 1+2+2 type for h <= 44 at all orders (n <= 12 suffices: classes grow past 5)
plat = {}
for h in range(1, 45, 2):
    for n in range(1, 13):
        if typ(h, n) == (1, 2): plat.setdefault(h, []).append(n)
print("h <= 44 with type 1 + 2 + 2, and the orders where it occurs:", plat)
print("never 1+2+2 (h <= 44, h realized):", [h for h in range(1, 45, 2) if h not in plat and h not in (7, 21)])
targets = {'(2,2) four-dimensional polytopes': (2, 2), '(1,1) plane tilings / >=5D polytopes / 4D honeycombs': (1, 1),
           '(1,0) 3D and >=5D honeycombs': (1, 0)}
for name, t in targets.items():
    hits = {}
    for h in range(1, 45, 2):
        for n in range(1, 13):
            if typ(h, n) == t: hits.setdefault(h, []).append(n)
    print(name, "->", hits)
# H = 3 realizes every (1, k) at odd orders and every (0, k) at even orders
print("H=3 types n=3..12:", [typ(3, n) for n in range(3, 13)])
# smoothness: largest prime p with 7p <= 20405 (census bound, wave 46)
def isprime(x):
    return x > 1 and all(x % d for d in range(2, int(x ** .5) + 1))
pmax = max(p for p in range(3, 20405 // 7 + 1) if isprime(p))
print("largest prime p with 7p <= 20405:", pmax, "-> every odd", pmax, "-smooth number except 7, 21 is a count")
