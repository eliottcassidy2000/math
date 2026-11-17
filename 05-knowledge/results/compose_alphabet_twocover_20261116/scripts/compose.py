#!/usr/bin/env python3
"""compose.py: flexible cores composed with the one-port alphabet.  From profiles8.npz (classes of orders 1-8 with h, e, B,
M): (1) nontrivial balanced flips (labelg) and, for every root, the change of the L3 profile (H, g, B, e) across the
flip; (2) the census of rooted L3 profiles: permanent +2 pairs of kind K1 (same (g, e, B), H' = H + 2) and K2 (same
(H, g, B), e' = e + 2), their centres; (3) profile richness per fibre vs neutral edges; (4) the tribonacci check of
iterated N."""
import numpy as np, subprocess
from collections import defaultdict, Counter
exec(open('nearjoin_check.py').read().split('# (1)')[0])
P = np.load('profiles8.npz')
def d6_encode(n, M):
    bits = ''.join(str(int(M[i][j])) for i in range(n) for j in range(n)); bits += '0' * ((-len(bits)) % 6)
    return '&' + chr(n + 63) + ''.join(chr(int(bits[i:i+6], 2) + 63) for i in range(0, len(bits), 6))
def canon(lines):
    out = subprocess.run(['labelg', '-q'], input='\n'.join(lines) + '\n', capture_output=True, text=True).stdout.split(); assert len(out) == len(lines); return out
def rev(T): n = len(T); return [[T[j][i] for j in range(n)] for i in range(n)]
def rooted(T):
    """(H, g_r, e_r, B_r) for every root r."""
    n = len(T); N = 1 << n; full = N - 1; dp, Hm = ham_all(T); dpr, _ = ham_all(rev(T))
    out = []
    for r in range(n):
        B = sum(dpr[m][r] * Hm[full & ~m] for m in range(1, N) if (m >> r) & 1)
        out.append((Hm[full], Hm[full & ~(1 << r)], dp[full][r], B))
    return out
# (1) flips
flip_stats = Counter(); endpoint_B0 = 0; endpoint_total = 0; examples = []
for n in range(3, 9):
    A, h, M = P[f'A{n}'], P[f'h{n}'], P[f'M{n}']; k = len(h)
    cn = canon([d6_encode(n, A[i]) for i in range(k)]); index = {s: i for i, s in enumerate(cn)}
    flips = []; owner = []
    for i in range(k):
        for a in range(n):
            for b in range(n):
                if A[i][a][b] and M[i][a][b] == M[i][b][a]:
                    T2 = A[i].copy(); T2[a][b] = 0; T2[b][a] = 1; flips.append(d6_encode(n, T2)); owner.append((i, a, b, T2))
    cf = canon(flips) if flips else []
    for s, (i, a, b, T2) in zip(cf, owner):
        if index[s] == i: continue
        p1 = rooted(A[i].tolist()); p2 = rooted(T2.tolist()); assert p1[0][0] == p2[0][0]
        for r in range(n):
            H, g, e, B = p1[r]; H2, g2, e2, B2 = p2[r]
            if r in (a, b):
                endpoint_total += 1; endpoint_B0 += (B2 == B)
            else:
                key = ('dg0' if g2 == g else 'dg1', 'dB0' if B2 == B else 'dB1', 'de0' if e2 == e else 'de1'); flip_stats[key] += 1
                if g2 == g and B2 == B and e2 != e and len(examples) < 6: examples.append((n, i, a, b, r, (H, g, B), e, e2))
tot_nonend = sum(flip_stats.values()); stored = flip_stats[('dg0', 'dB0', 'de1')]; same_l2 = stored + flip_stats[('dg0', 'dB0', 'de0')]
print(f"(1) nontrivial balanced flips, roots off the flipped arc: {tot_nonend} (flip, root) pairs; with Delta g = 0 and Delta B = 0: {same_l2}, of which Delta e != 0 (stored choice): {stored} ({stored/max(1,same_l2):.3f}); full split {dict(flip_stats)}")
print(f"    endpoints of the flipped arc: {endpoint_total} (flip, root) pairs, Delta B = 0 in {endpoint_B0}; examples of stored choices (n, class, a, b, root, (H, g, B), e, e'): {examples}")
# (2) rooted profile census and permanent pairs
prof = defaultdict(list)   # (H, g, e, B) -> list of (n, i, r)
L = defaultdict(set)
for n in range(1, 9):
    A, h = P[f'A{n}'], P[f'h{n}']
    for i in range(len(h)):
        for r, (H, g, e, B) in enumerate(rooted(A[i].tolist())):
            prof[(int(H), int(g), int(e), int(B))].append((n, i, r)); L[int(H)].add((int(g), int(e), int(B)))
keys = set(prof)
K1 = set(); K2 = set()
for (H, g, e, B) in keys:
    if (H + 2, g, e, B) in keys: K1.add(H + 1)
    if (H, g, e + 2, B) in keys: K2.add(H + B - e - 1)
LIM = 1000; evens = list(range(4, LIM + 1, 2))
c1 = sum(1 for c in evens if c in K1); c2 = sum(1 for c in evens if c in K2); cu = sum(1 for c in evens if c in K1 or c in K2)
print(f"(2) rooted L3 profiles (orders <= 8): {len(keys)} distinct; permanent pairs K1 (same (g, e, B), H' = H + 2): {len(K1)} centres, K2 (same (H, g, B), e' = e + 2, permanent after one N): {len(K2)} centres")
print(f"    even centres in [4, {LIM}] covered: K1 {c1}, K2 {c2}, union {cu} of {len(evens)}; missed by the union: {[c for c in evens if c not in K1 and c not in K2][:40]}")
print(f"    K1 centres <= 60: {sorted(c for c in K1 if c <= 60)}; K2 centres <= 60: {sorted(c for c in K2 if c <= 60)}")
# (3) richness vs neutral edges (recompute edges per fibre from (1)'s data is heavy; reuse neutral.out numbers)
edges = {}
for line in open('neutral_fibres.txt'):
    hh, ne = map(int, line.split()); edges[hh] = ne
rows = [(hh, len(L[hh]), edges[hh]) for hh in sorted(L) if hh <= 125 and hh in edges]
from scipy.stats import spearmanr
print(f"(3) profile richness |L3 profiles of F_h| vs nontrivial neutral edges over {len(rows)} fibres: Spearman {spearmanr([r[1] for r in rows], [r[2] for r in rows]).correlation:.3f}; richness of 11, 13, 37, 47, 67: {[len(L[x]) for x in (11, 13, 37, 47, 67)]}")
# (4) tribonacci check of iterated N
bad = tot = 0
for n in range(3, 7):
    A = P[f'A{n}']
    for i in range(len(A)):
        T = A[i].tolist()
        for r in range(n):
            seq = [ham_all(T)[1][-1]]; U = [row[:] for row in T]; root = r
            for step in range(5):
                m = len(U); U = [row[:] + [0 if j == root else 1] for j, row in enumerate(U)] + [[1 if j == root else 0 for j in range(m)] + [0]]
                root = m  # the new vertex becomes the root (the tower)
                seq.append(ham_all(U)[1][-1])
            for kk in range(3, len(seq)):
                tot += 1
                if seq[kk] != seq[kk-1] + seq[kk-2] + seq[kk-3]: bad += 1
print(f"(4) iterated N at the newest vertex (the tower): {tot} tribonacci checks, {bad} failures")
