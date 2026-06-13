#!/usr/bin/env python3
"""
n19_alpha12.py -- alpha_1, alpha_2 for the n=19 maximizer.
alpha_3 is too slow (|nonzero cc|^2 loop), skip for now.
"""

import time

def circulant(n, S):
    adj = [0]*n
    for v in range(n):
        for k in S:
            adj[v] |= 1 << ((v+k) % n)
    return adj

def cycle_cc(adj, n):
    full = (1 << n)-1
    cc = [0] * (1 << n)
    for s in range(n):
        s_bit = 1 << s
        hi_mask = full & ~(s_bit - 1)
        queue = {(s_bit, s): 1}
        while queue:
            nq = {}
            for (mask, v), cnt in queue.items():
                L = bin(mask).count('1')
                if L >= 3 and L % 2 == 1:
                    if (adj[v] >> s) & 1:
                        cc[mask] += cnt
                cands = adj[v] & ~mask & hi_mask
                while cands:
                    ub = cands & -cands; u = ub.bit_length()-1
                    key = (mask|ub, u)
                    nq[key] = nq.get(key, 0) + cnt
                    cands ^= ub
            queue = nq
    return cc

def sos(f, n):
    f = list(f)
    for i in range(n):
        for mask in range(1 << n):
            if (mask >> i) & 1:
                f[mask] += f[mask ^ (1 << i)]
    return f

n = 19
S = list(range(1, (n+1)//2))   # cyclic interval {1,...,9}
adj = circulant(n, S)
full = (1 << n) - 1

print(f"n={n}, S={S}")

t0 = time.time()
cc = cycle_cc(adj, n)
print(f"cycle_cc: {time.time()-t0:.1f}s, |nonzero|={sum(1 for c in cc if c):,}")

t0 = time.time()
f = sos(cc, n)
print(f"SOS: {time.time()-t0:.1f}s")

alpha1 = sum(cc)
print(f"α₁ = {alpha1:,}")

# cycle length distribution
lc = {}
for mask, cnt in enumerate(cc):
    if cnt:
        L = bin(mask).count('1')
        lc[L] = lc.get(L, 0) + cnt
print(f"Cycle lengths: {dict(sorted(lc.items()))}")

t0 = time.time()
alpha2 = sum(cc[m] * f[(~m) & full] for m in range(1 << n)) // 2
print(f"α₂ = {alpha2:,}  ({time.time()-t0:.1f}s)")

H_ref = 1_184_212_824_763
ratio = alpha1 / (2 * alpha2)
print(f"\nH (ref) = {H_ref:,}")
print(f"α₁/(2α₂) = {ratio:.5f}")
print(f"α₃/α₂ ratio not computed (too slow with naive loop)")

# Also compare Paley alpha_1, alpha_2
print("\n--- Also Paley T_19 ---")
paley_S = sorted({i*i % n for i in range(1, n)} - {0})
print(f"Paley S = {paley_S}")
adj_p = circulant(n, paley_S)
t0 = time.time()
cc_p = cycle_cc(adj_p, n)
f_p = sos(cc_p, n)
a1_p = sum(cc_p)
a2_p = sum(cc_p[m] * f_p[(~m) & full] for m in range(1 << n)) // 2
print(f"Paley: α₁={a1_p:,}, α₂={a2_p:,}, ratio={a1_p/(2*a2_p):.5f}  ({time.time()-t0:.1f}s)")
lc_p = {}
for mask, cnt in enumerate(cc_p):
    if cnt:
        L = bin(mask).count('1')
        lc_p[L] = lc_p.get(L, 0) + cnt
print(f"Paley cycle lengths: {dict(sorted(lc_p.items()))}")

print(f"\nDifference breakdown:")
print(f"  α₁ (cyclic - Paley): {alpha1 - a1_p:+,}")
print(f"  α₂ (cyclic - Paley): {alpha2 - a2_p:+,}")
print(f"  2α₁ diff: {2*(alpha1-a1_p):+,}")
print(f"  4α₂ diff: {4*(alpha2-a2_p):+,}")
H_pal_ref = 1_172_695_746_915
print(f"  H diff: {H_ref - H_pal_ref:+,}")
