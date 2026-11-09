"""audit D2 (THM-5056): brute-force tests of the lemmas used by Part II, on graphs of every parity.
  L4.1   d(PH(G,u,l)) = d(G) + (2^r - 1) 2^(n+l-r), r = |N[u]|, for ALL graphs G (odd ones included), all u, l = 1..3
  K2J    d(K_2 v G) = 3 * 2^|G| + d(G)  (used in case 1, stated without proof in Part II)
  CONE   d(K_1 v G) = 2^|G| + d(G) for |G| >= 1
  O1(b)  K_1 v H is odd for every odd H  (exhaustive over all labelled graphs H on <= 5 vertices, random on 6-7)
  O1(c)  K_2 v G is odd for every G
  O1(d)  PH(G,u,l) keeps adjacent twins x, y when u is not in {x, y}, hence is odd
Parity is decided by a complete backtracking enumeration of Aut and the reversed-edge count of every element."""
import random, itertools
rng = random.Random(77)

def dcount(G):
    n = len(G); cl = [G[v] | (1 << v) for v in range(n)]; full = (1 << n) - 1; c = 0
    for S in range(1 << n):
        dom = 0
        for v in range(n):
            if (S >> v) & 1: dom |= cl[v]
        c += dom == full
    return c
def auts(G):
    n = len(G); dg = [bin(a).count('1') for a in G]; img = [-1] * n; used = [False] * n; out = []
    def bt(v):
        if v == n: out.append(img[:]); return
        for c in range(n):
            if used[c] or dg[c] != dg[v]: continue
            if any(((G[v] >> j) & 1) != ((G[c] >> img[j]) & 1) for j in range(v)): continue
            img[v] = c; used[c] = True; bt(v + 1); used[c] = False; img[v] = -1
    bt(0); return out
def is_odd(G):
    n = len(G)
    for g in auts(G):
        if sum(1 for a in range(n) for b in range(a + 1, n) if (G[a] >> b) & 1 and g[a] > g[b]) % 2: return True
    return False
def cone(G): n = len(G); return [a | (1 << n) for a in G] + [(1 << n) - 1]
def k2join(G):
    n = len(G); allG = (1 << n) - 1
    return [x | (3 << n) for x in G] + [allG | (1 << (n + 1)), allG | (1 << n)]
def ph(G, u, l):
    n = len(G); H = list(G); hm = 0
    for v in range(n):
        if v != u: H[v] |= 1 << n; hm |= 1 << v
    for i in range(l): hm |= 1 << (n + 1 + i)
    return H + [hm] + [1 << n] * l
def rand_graph(n, p):
    G = [0] * n
    for a in range(n):
        for b in range(a + 1, n):
            if rng.random() < p: G[a] |= 1 << b; G[b] |= 1 << a
    return G
def all_graphs(n):
    pairs = list(itertools.combinations(range(n), 2))
    for mask in range(1 << len(pairs)):
        G = [0] * n
        for i, (a, b) in enumerate(pairs):
            if (mask >> i) & 1: G[a] |= 1 << b; G[b] |= 1 << a
        yield G
def twins_pairs(G):
    n = len(G)
    return [(x, y) for x in range(n) for y in range(x + 1, n) if (G[x] >> y) & 1 and (G[x] | 1 << x) == (G[y] | 1 << y)]

stats = dict(L41=0, K2J=0, CONE=0, O1b=0, O1b_univ=0, O1c=0, O1d=0)
fail = []
# L4.1, K2J, CONE on random graphs of all parities (n = 1..8)
for trial in range(400):
    n = rng.randint(1, 8); G = rand_graph(n, rng.choice((0.2, 0.4, 0.6, 0.8))); dG = dcount(G)
    for u in range(n):
        r = bin(G[u] | 1 << u).count('1')
        for l in (1, 2, 3):
            if len(G) + 1 + l > 12: continue
            stats['L41'] += 1
            if dcount(ph(G, u, l)) != dG + (2 ** r - 1) * 2 ** (n + l - r): fail.append(('L41', G, u, l))
    stats['K2J'] += 1
    if dcount(k2join(G)) != 3 * 2 ** n + dG: fail.append(('K2J', G))
    stats['CONE'] += 1
    if dcount(cone(G)) != 2 ** n + dG: fail.append(('CONE', G))
# O1(b): exhaustive over labelled H on <= 5 vertices
for n in range(1, 6):
    for H in all_graphs(n):
        if not is_odd(H): continue
        stats['O1b'] += 1
        if any(bin(a).count('1') == n - 1 for a in H): stats['O1b_univ'] += 1
        if not is_odd(cone(H)): fail.append(('O1b', H))
# O1(b) random on 6-7 vertices, O1(c) random
for trial in range(300):
    n = rng.randint(6, 7); H = rand_graph(n, rng.choice((0.3, 0.5, 0.7, 0.9)))
    if is_odd(H):
        stats['O1b'] += 1
        if any(bin(a).count('1') == n - 1 for a in H): stats['O1b_univ'] += 1
        if not is_odd(cone(H)): fail.append(('O1b', H))
for trial in range(150):
    G = rand_graph(rng.randint(1, 6), rng.random()); stats['O1c'] += 1
    if not is_odd(k2join(G)): fail.append(('O1c', G))
# O1(d): random graphs with an adjacent twin pair; PH with u outside the pair keeps the pair and is odd
for trial in range(300):
    n = rng.randint(3, 7); G = rand_graph(n, rng.random())
    tp = twins_pairs(G)
    if not tp: continue
    x, y = rng.choice(tp)
    for u in range(n):
        if u in (x, y): continue
        l = rng.randint(1, 2); H = ph(G, u, l); stats['O1d'] += 1
        if (x, y) not in twins_pairs(H) or not is_odd(H): fail.append(('O1d', G, x, y, u, l))
print('tests:', stats)
print('failures:', len(fail), fail[:3])
