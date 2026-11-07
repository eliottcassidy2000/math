"""Targeted check of THM-4892 (D): order joins of strong tournaments, Q = Pk2/global swap vs Q^+ = prod Q(C_i)."""
import itertools, random, sys
sys.path.insert(0, '.')
from sheaf_check import random_tournament, ham_paths, strong_components, odd_cycles, packings

def join(blocks):
    n = sum(len(b) for b in blocks); adj = [[False] * n for _ in range(n)]; off = 0; ranges = []
    for b in blocks:
        k = len(b)
        for i in range(k):
            for j in range(k):
                adj[off + i][off + j] = b[i][j]
        ranges.append(range(off, off + k)); off += k
    for x in range(len(ranges)):
        for y in range(x + 1, len(ranges)):
            for i in ranges[x]:
                for j in ranges[y]:
                    adj[i][j] = True
    return adj

rng = random.Random(11)
def strong_random(k):
    while True:
        t = random_tournament(k, rng)
        if len(strong_components(t, range(k))) == 1:
            return t
ok = bad = 0; seen = set()
for _ in range(80):
    k = rng.choice([2, 2, 3]); sizes = [rng.choice([3, 3, 4]) for _ in range(k)] if k == 2 else [3, 3, rng.choice([3, 4])]; parts = [strong_random(s) for s in sizes]
    adj = join(parts); n = len(adj); V = list(range(n))
    comps = strong_components(adj, V)
    cyc = odd_cycles(adj, n); pk = packings(cyc, V)
    H = ham_paths(adj, V)
    q_presheaf = sum(1 if not p else 2 ** (len(p) - 1) for p in pk)          # classes mod global swap
    q_sheaf = 1
    for c in comps:
        q_sheaf *= (ham_paths(adj, c) + 1) // 2
    # fibre sizes: for each packing, (#global classes) / (#local images) must be 2^(m-1)
    fib_ok = True
    for p in pk:
        if not p: continue
        comp_of = [next(k for k, C in enumerate(comps) if set(cyc[i]) <= C) for i in p]
        m = len(set(comp_of))
        glob = 2 ** (len(p) - 1)
        loc = 1
        for k in set(comp_of):
            loc *= 2 ** (comp_of.count(k) - 1)
        fib_ok &= (glob == loc * 2 ** (m - 1))
    qs = [(ham_paths(adj, c) - 1) // 2 for c in comps]
    if len(qs) == 2:
        a, b = qs
        law_ok = (q_presheaf - 1 == a + b + 2 * a * b) and (q_sheaf - 1 == a + b + a * b)
    else:
        law_ok = True
    if q_presheaf == (H + 1) // 2 and fib_ok and law_ok and q_presheaf > q_sheaf: ok += 1
    else: bad += 1; print('BAD', qs, q_presheaf, q_sheaf, H)
print('order joins of 2-3 strong pieces: |Q| = (H+1)/2, fibres 2^(m-1), join laws, Q not separated:', ok, 'ok,', bad, 'bad')
