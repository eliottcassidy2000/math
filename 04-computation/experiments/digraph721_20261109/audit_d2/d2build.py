"""audit D2 (THM-5056), item 6: an independent re-implementation of the THM-5056 construction.

Written for the audit (own graph operations, own seeds, own threshold graphs).  Differences from oddbuild.py that
make this an adversarial test of the proof rather than a re-run of it:
  * the base below 512 is taken from d2base_witnesses.txt (realisers found and verified by this audit, random among
    the census candidates), not from oddbuild.py's first-found census graph;
  * offset chains (Lemma O2) make RANDOM valid choices: a random order of the +2 (R) and +3 (PH(3,2)) moves, a random
    excluded pendant (outside the twin pair, neighbour not universal) for R, a random excluded vertex with |N[w]| = 3
    (outside the twin pair) for PH(3,2), a random admissible seed and a random seed vertex w;
  * after every move the odd-good invariants of Lemma O2 are asserted (twins x, y; >= 2 pendants outside {x, y}
    with non-universal neighbours; a vertex outside {x, y} with |N[w]| = 3), and the value formula of Lemma 4.1 is
    asserted symbolically (d tracked as an integer, checked later by brute force in d2check).
Output: wcheck-format lines (set O) for d2check.
usage: python3 d2build.py MODE OUT [seed]
  MODE = all:LO:HI       every odd k in [LO, HI)
         rand:LO:HI:N    N random odd k in [LO, HI)
         special:LMAX    2^(L-1)+c, 2^(L-1)+2^(L-3)+c, 3*2^(L-2)+c (c in {1,5,37}), 10 <= L <= LMAX, plus random
                         deep cone chains ending in an offset, and offset graphs at every order in Gamma_odd(c)
"""
import sys, random

def bl(k): return k.bit_length()
def deg(G, v): return bin(G[v]).count('1')
def pad(G, N): return G + [0] * (N - len(G))
def cone(G):
    n = len(G); return [a | (1 << n) for a in G] + [(1 << n) - 1]
def union(G, H):
    n = len(G); return G + [a << n for a in H]
def k2join(G):                       # K_2 v G: two new vertices a, b adjacent to each other and to all of G
    n = len(G); allG = (1 << n) - 1; a, b = n, n + 1
    return [x | (1 << a) | (1 << b) for x in G] + [allG | (1 << b), allG | (1 << a)]
def ph(G, u, l):                     # hub adjacent to every vertex of G except u, plus l leaves on the hub
    n = len(G); h = n; H = list(G); hm = 0
    for v in range(n):
        if v != u: H[v] |= 1 << h; hm |= 1 << v
    for i in range(l): hm |= 1 << (n + 1 + i)
    return H + [hm] + [1 << h] * l
def threshold(m, N):
    if N <= 1: assert m == 1; return [0] * N
    if m >= 2 ** (N - 1): return cone(threshold(m - 2 ** (N - 1), N - 1))
    return threshold(m, N - 1) + [0]
def edges_of(pairs, n):
    G = [0] * n
    for u, v in pairs: G[u] |= 1 << v; G[v] |= 1 << u
    return G
def adj_twins(G, x, y):
    return (G[x] >> y) & 1 and (G[x] | (1 << x)) == (G[y] | (1 << y))
def universal(G, v): return deg(G, v) == len(G) - 1
def pendants_ok(G, tw):
    return [p for p in range(len(G)) if p not in tw and deg(G, p) == 1
            and not universal(G, (G[p] & -G[p]).bit_length() - 1)]
def deg2(G, tw): return [w for w in range(len(G)) if w not in tw and deg(G, w) == 2]
def assert_odd_good(G, tw):
    assert adj_twins(G, *tw), 'twins lost'
    assert len(pendants_ok(G, tw)) >= 2, 'fewer than two admissible pendants'
    assert len(deg2(G, tw)) >= 1, 'no admissible degree-2 vertex'

# seeds: (graph, twin pair, admissible w's, d)
def seeds(c):
    S = []
    if c == 5:
        G = edges_of([(0, 1), (0, 2), (0, 3), (1, 2), (1, 3)], 4); S.append((G, (0, 1), 13))          # K_4 - e
    if c == 1:
        G = edges_of([(0, 1), (2, 3), (3, 4), (4, 5), (5, 2)], 6); S.append((G, (0, 1), 33))          # K_2 u C_4
    if c == 37:
        pairs = [(0, 1)] + [(a, v) for a in (0, 1) for v in range(2, 7)] + [(2, 3), (2, 4)]
        S.append((edges_of(pairs, 7), (0, 1), 101))                                                   # K_2 v (P_3 u E_2)
        S.append((union(edges_of([(0, 1)], 2), threshold(55, 6)), (0, 1), 165))                       # K_2 u T_55
    return S

def offset_graph(c, n, rng):
    """an odd-good graph on n vertices with d = 2^(n-1) + c, by a random admissible seed and random move order."""
    opts = []
    for G0, tw, d0 in seeds(c):
        n0 = len(G0); assert d0 == 2 ** (n0 - 1) + c
        r = n - (n0 + 3)
        if r >= 0 and r != 1: opts.append((G0, tw, d0, r))
    if not opts: return None
    G0, tw, d0, r = rng.choice(opts)
    ws = deg2(G0, tw); assert ws, 'seed has no admissible w'
    w = rng.choice(ws); n0 = len(G0)
    G = ph(G0, w, 2); d = d0 + 7 * 2 ** (n0 - 1)
    assert_odd_good(G, tw)
    # random composition of r into parts 2 and 3
    moves = []
    while r > 0:
        choices = [s for s in (2, 3) if r - s >= 0 and r - s != 1]
        s = rng.choice(choices); moves.append(s); r -= s
    for s in moves:
        m = len(G)
        if s == 2:
            p = rng.choice(pendants_ok(G, tw)); G = ph(G, p, 1); d += 3 * 2 ** (m - 1)
        else:
            w = rng.choice(deg2(G, tw)); G = ph(G, w, 2); d += 7 * 2 ** (m - 1)
        assert_odd_good(G, tw)
        assert d == 2 ** (len(G) - 1) + c
    assert len(G) == n
    return G

BASE = {}
def load_base(path):
    for line in open(path):
        p = line.split(); n = int(p[2]); k = int(p[3]); G = [0] * n
        for e in line.split('|')[1].split():
            u, v = map(int, e.split('-')); G[u] |= 1 << v; G[v] |= 1 << u
        BASE[k] = G

def build(k, rng):
    L = bl(k)
    if k < 512: return BASE[k]
    if k >= 3 * 2 ** (L - 2):
        t = k - 3 * 2 ** (L - 2); assert t % 2 == 1 and 1 <= t < 2 ** (L - 2)
        return k2join(pad(threshold(t, bl(t)), L - 2))
    x = k - 2 ** (L - 1); assert x % 2 == 1 and x < 2 ** (L - 2)
    if x not in (1, 5, 37):
        R = build(x, rng); assert len(R) == bl(x) <= L - 2
        return cone(pad(R, L - 1))
    if (x, L) == (1, 10):
        return union(edges_of([(0, 1)], 2), threshold(171, 8))
    G = offset_graph(x, L, rng); assert G is not None, (k, x, L)
    return G

def emit(out, tag, G, k):
    n = len(G)
    out.write(f'{tag} O {n} {k} | ' + ' '.join(f'{u}-{v}' for u in range(n) for v in range(u + 1, n) if (G[u] >> v) & 1) + '\n')

def main():
    mode = sys.argv[1]; out = open(sys.argv[2], 'w'); seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    rng = random.Random(seed); load_base('d2base_witnesses.txt')
    p = mode.split(':'); cnt = 0
    if p[0] in ('all', 'rand'):
        lo, hi = int(p[1]), int(p[2]); ks = [k for k in range(lo | 1, hi, 2) if k not in (1, 5, 37)]
        if p[0] == 'rand': ks = sorted(rng.sample(ks, int(p[3])))
        for k in ks:
            G = build(k, rng); assert len(G) == bl(k); emit(out, f'k{k}', G, k); cnt += 1
    elif p[0] == 'special':
        LM = int(p[1]); ks = set()
        for L in range(10, LM + 1):
            for c in (1, 5, 37):
                ks.add(2 ** (L - 1) + c); ks.add(3 * 2 ** (L - 2) + c)
                if L >= 13: ks.add(2 ** (L - 1) + 2 ** (L - 3) + c)
        # random deep cone chains: leading 10 at every level, ending in an offset value 2^(L0-1) + c, L0 >= 10
        for _ in range(60):
            c = rng.choice((1, 5, 37)); L = rng.randint(10, 12); k = 2 ** (L - 1) + c
            while True:
                L2 = bl(k) + rng.randint(2, 3)
                if L2 > LM: break
                k = 2 ** (L2 - 1) + k
            ks.add(k)
        for k in sorted(ks):
            G = build(k, rng); assert len(G) == bl(k); emit(out, f's{k}', G, k); cnt += 1
        # offset graphs at every order of Gamma_odd(c) up to LM, three random constructions each
        for c in (1, 5, 37):
            for n in range(4, LM + 1):
                for t in range(3):
                    G = offset_graph(c, n, rng)
                    if G is None: break
                    emit(out, f'off{c}_n{n}_t{t}', G, 2 ** (n - 1) + c); cnt += 1
    out.close(); print('d2build', mode, 'seed', seed, 'graphs', cnt, '(all invariants asserted; orders = L(k))')

if __name__ == '__main__':
    main()
