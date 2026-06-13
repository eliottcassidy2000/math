#!/usr/bin/env python3
"""
n17_alpha2.py -- Corrected alpha_4 including [3,3,5,5] quadruples.

Missing from previous: [3,3,5,5] at n=17 (3+3+5+5=16 ≤ 17).
At n=15 this wasn't needed (16 > 15). At n=17 it is.

Strategy:
  alpha4 = [3,3,3,3] + [3,3,3,L≥5] + [3,3,5,5]
  [3,3,3,L≥5]: covered by raw_all - raw_3s (triple of 3-cycles × non-3-cycle in complement)
  [3,3,5,5]:   NEW — iterate pairs of 3-cycles, pairs of 5-cycles in complement
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

n = 17
S = [1,2,3,4,5,6,7,8]
adj = circulant(n, S)
full = (1 << n) - 1

# Load from previous run (already computed, just need cc)
print("Computing cc...")
t0 = time.time()
cc = cycle_cc(adj, n)
f = sos(cc, n)
print(f"Done: {time.time()-t0:.1f}s")

alpha1 = sum(cc)
alpha2 = sum(cc[m] * f[(~m) & full] for m in range(1 << n)) // 2
print(f"α₁={alpha1:,}  α₂={alpha2:,}")

# Pre-compute cc3, cc5, their SOS and lists
cc3 = [cc[m] if bin(m).count('1') == 3 else 0 for m in range(1 << n)]
cc5 = [cc[m] if bin(m).count('1') == 5 else 0 for m in range(1 << n)]
f3 = sos(cc3, n)
f5 = sos(cc5, n)
tc3 = [(m, cc3[m]) for m in range(1 << n) if cc3[m]]
tc5 = [(m, cc5[m]) for m in range(1 << n) if cc5[m]]
print(f"|tc3|={len(tc3)}, |tc5|={len(tc5)}")

# ── alpha_4: [3,3,3,3] and [3,3,3,L≥5] ───────────────────────────────────────
t0 = time.time()
raw_3s = 0
raw_all = 0
for i, (m1, c1) in enumerate(tc3):
    comp1 = (~m1) & full
    for j, (m2, c2) in enumerate(tc3):
        if j <= i: continue
        if m1 & m2: continue
        m12 = m1 | m2
        comp12 = comp1 & (~m2) & full
        for k, (m3, c3) in enumerate(tc3):
            if k <= j: continue
            if m12 & m3: continue
            comp123 = comp12 & (~m3) & full
            raw_3s  += c1 * c2 * c3 * f3[comp123]
            raw_all += c1 * c2 * c3 * f[comp123]

alpha4_3333 = raw_3s // 4
alpha4_3_35plus = raw_all - raw_3s   # [3,3,3,L≥5] with 4th cycle of odd length ≥ 5
print(f"[3,3,3,3]={alpha4_3333:,}, [3,3,3,L≥5]={alpha4_3_35plus:,}  ({time.time()-t0:.1f}s)")

# ── alpha_4: [3,3,5,5] ────────────────────────────────────────────────────────
# Iterate unordered pairs of 3-cycles × pairs of 5-cycles in complement.
# Each [3,3,5,5] quadruple:
#   pair of 3-cycles counted once (unordered by i<j)
#   5-cycle pair: iterate m3 over tc5 ∩ comp12, lookup f5[comp12 & ~m3]
#   Each quadruple counted twice (once with each 5-cycle as "m3"), divide by 2.
t0 = time.time()
raw_3355 = 0
pair_count = 0
for i, (m1, c1) in enumerate(tc3):
    for j in range(i+1, len(tc3)):
        m2, c2 = tc3[j]
        if m1 & m2: continue
        pair_count += 1
        m12 = m1 | m2
        comp12 = (~m12) & full
        for m3, c3 in tc5:
            if m3 & m12: continue   # 5-cycle must be in complement
            raw_3355 += c1 * c2 * c3 * f5[comp12 & ~m3]

alpha4_3355 = raw_3355 // 2
print(f"[3,3,5,5]={alpha4_3355:,}  ({time.time()-t0:.1f}s, {pair_count} 3-cycle pairs)")

alpha4 = alpha4_3333 + alpha4_3_35plus + alpha4_3355
print(f"α₄ = {alpha4:,}")
print(f"  Expected: 45,997,104")
print(f"  Match: {alpha4 == 45_997_104}")

# ── Verify ───────────────────────────────────────────────────────────────────
alpha3 = 458_011_858  # from previous run (372s computation)
alpha5 = 1_800_368

H_ocf = 1 + 2*alpha1 + 4*alpha2 + 8*alpha3 + 16*alpha4 + 32*alpha5
H_ref = 13_689_269_499
print(f"\nH (OCF) = {H_ocf:,}")
print(f"H (ref) = {H_ref:,}")
print(f"Match: {H_ocf == H_ref}")
