#!/usr/bin/env python3
"""
n17_alpha.py -- Full alpha decomposition for n=17 maximizer.
S = {1,2,3,4,5,6,7,8} (consecutive-half tournament), H = 13,689,269,499
"""

import time
from math import factorial

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

print(f"n={n}, S={S}")
print(f"kmax = {n//3}")

# ── Step 1: cycle cc ─────────────────────────────────────────────
t0 = time.time()
cc = cycle_cc(adj, n)
print(f"cycle_cc: {time.time()-t0:.2f}s")

# ── Step 2: SOS ──────────────────────────────────────────────────
t0 = time.time()
f = sos(cc, n)
print(f"SOS:      {time.time()-t0:.2f}s")

# ── alpha_1 ──────────────────────────────────────────────────────
alpha1 = sum(cc)
print(f"\nα₁ = {alpha1:,}")

# cycle length distribution
len_counts = {}
for mask, cnt in enumerate(cc):
    if cnt:
        L = bin(mask).count('1')
        len_counts[L] = len_counts.get(L, 0) + cnt
print(f"Cycle lengths: {dict(sorted(len_counts.items()))}")

# ── alpha_2 ──────────────────────────────────────────────────────
t0 = time.time()
alpha2 = sum(cc[m] * f[(~m) & full] for m in range(1 << n)) // 2
print(f"α₂ = {alpha2:,}   ({time.time()-t0:.2f}s)")

# ── alpha_3 ──────────────────────────────────────────────────────
t0 = time.time()
nonzero = [(m, cc[m]) for m in range(1 << n) if cc[m]]
print(f"|nonzero cc| = {len(nonzero):,}")

alpha3 = 0
for i, (m1, c1) in enumerate(nonzero):
    comp1 = (~m1) & full
    for m2, c2 in nonzero:
        if m2 > m1 and not (m1 & m2):
            alpha3 += c1 * c2 * f[comp1 & (~m2) & full]
alpha3 //= 3
print(f"α₃ = {alpha3:,}   ({time.time()-t0:.2f}s)")

# ── alpha_4 ─────────────────────────────────────────────────────
# At n=17: kmax=5. For k=4: cycles could be [3,3,3,3] (12≤17), [3,3,3,5] (14≤17), [3,3,5,5] (16≤17), [3,3,3,7] (16≤17), [3,5,9] etc.
# Use the split approach: f3 for only 3-cycles, f_all = f for all odd cycles.
# But for n=17 with 4-tuples, we need to be careful about which cycle sizes sum to ≤ 17.
# Use general approach: iterate triples of cycles, look up f[complement] and f3[complement].
t0 = time.time()
cc3 = [cc[m] if bin(m).count('1') == 3 else 0 for m in range(1 << n)]
f3 = sos(cc3, n)
tc3 = [(m, cc3[m]) for m in range(1 << n) if cc3[m]]
print(f"|tc3| = {len(tc3):,}")

# alpha_4: sum over ORDERED 4-cycles / 4!  but account for different size combinations
# For each ordered quadruple of disjoint odd cycles:
# - [3,3,3,3]: counted 4! = 24 times in ordered sum, divide by 4 for unordered
# - [3,3,3,L≠3]: counted 3! = 6 for the 3-cycle order × 1 for the L-cycle = 6 times? No...
#
# Actually: α_4 = (1/4) * sum_{disjoint pairs (m1,m2,m3)} c1*c2*c3 * f[complement]
# where f[complement] = (sum of all odd cycles disjoint from m1,m2,m3)
# But that gives 4-tuples where the 4th cycle can be any odd cycle.
# The correct formula: iterate over 4-subsets of disjoint cycles (unordered).
# Faster: raw = sum_{ordered triples} c1*c2*c3 * f[comp] / 3 for the triple-part,
# then the 4th element contributes f[comp] which SOS already includes all sizes.
# Actually the standard SOS formula is:
#   α_4 = (1/4!) * sum_{ordered 4-tuples of disjoint cycles} product_counts
# But we can compute via:
#   sum_{disjoint ordered triples m1<m2<m3} c1*c2*c3 * f[~(m1|m2|m3)&full] gives
#   sum over unordered triples times the number of 4th cycles in complement = sum_3 * (disjoint cycles in comp)
# Hmm, this is getting complex. Let's just do:
#   raw_all = sum over ordered triples (using nonzero) × f[comp]
#   raw_3s  = sum over ordered triples (using nonzero) × f3[comp]
# Then 4-tuples of [3,3,3,3] = raw_3s / (3! * ... wait)
#
# Better: use the direct nested loop approach (expensive but correct)
# alpha_4 = (1/24) * sum_{all ordered 4-tuples of disjoint odd cycles} prod c_i
# Too slow. Use partial SOS:
# alpha_4 = (1/4) * sum_{ordered triples of disjoint cycles, m1<m2<m3} c1*c2*c3 * f[comp]
# NO — the "1/4" doesn't work simply when cycles have different sizes.
#
# Correct clean approach:
#   I_4 = sum_{k>=4} alpha_k * x^k contributes to I(Omega, x) at x^4
# Via inclusion-exclusion on SOS:
#   2*alpha_1 = sum cc = I_1
#   2*alpha_2 = (sum_m cc[m] * f[~m]) / 2 wait...
#
# Actually the SOS formula: f[T] = sum_{S subset T} cc[S] = number of ordered sequences of
# disjoint cycles with union ⊆ T (where empty = 1)... no that's not right.
# f[T] = I(Omega restricted to T, 1) = sum_{k>=0} alpha_k(T) where alpha_k(T) = k-packings in T
#
# The generating approach for alpha_4:
#   Treat the 3-cycle SOS as f3. Then:
#   g4 = sum_{m} cc[m] * f3[~m & full]  -> counts ordered pairs (L-cycle, 3-cycle-tuple) disjoint
#   But we need FOUR DISJOINT cycles.
#
# Let's use the layer-by-layer approach:
#   h1[T] = sum_{m ⊆ T} cc[m]  = f[T]   (all odd cycles in T)
#   h2[T] = sum_{disjoint m1,m2 ⊆ T} cc[m1]*cc[m2] / 2! = alpha_2 restricted to T
#   h3[T] = ... = alpha_3 restricted to T
# Then alpha_4 = sum_{m: cc[m]>0} cc[m] * h3[~m & full] / 4
# This requires h3 as an array, which is what f3-after-SOS gives IF h3 were computed recursively.
#
# For now: use the 3-cycle-only approach, getting alpha_4 for all-3-cycle quadruples,
# plus correction for mixed quadruples.

raw_3s = 0
for i, (m1, c1) in enumerate(tc3):
    comp1 = (~m1) & full
    for j, (m2, c2) in enumerate(tc3):
        if j <= i: continue
        if m1 & m2: continue
        m12 = m1 | m2
        comp12 = comp1 & (~m2) & full
        for m3, c3 in tc3[j+1:]:
            if not (m12 & m3):
                raw_3s += c1 * c2 * c3 * f3[comp12 & (~m3) & full]

alpha4_3333 = raw_3s // 4   # ordered triples of 3-cycles × f3 / 4 = all-3-cycle quads (unordered)

# For the mixed contribution (at least one cycle of length ≥ 5):
# raw_all = same triples but using full f
raw_all = 0
for i, (m1, c1) in enumerate(tc3):
    comp1 = (~m1) & full
    for j, (m2, c2) in enumerate(tc3):
        if j <= i: continue
        if m1 & m2: continue
        m12 = m1 | m2
        comp12 = comp1 & (~m2) & full
        for m3, c3 in tc3[j+1:]:
            if not (m12 & m3):
                raw_all += c1 * c2 * c3 * f[comp12 & (~m3) & full]

# raw_all counts: triples of 3-cycles × all cycles in complement
# = (4 * alpha4_3333) + (quadruples with exactly one non-3-cycle, counted once)
# = 4 * alpha4_3333 + alpha4_mixed
alpha4_mixed = raw_all - 4 * alpha4_3333
alpha4 = alpha4_3333 + alpha4_mixed
print(f"α₄ = {alpha4:,}  (3333: {alpha4_3333:,}, mixed: {alpha4_mixed:,})  ({time.time()-t0:.2f}s)")

# ── alpha_5 ─────────────────────────────────────────────────────
# For n=17, kmax=5, so we need 5 disjoint odd cycles covering all 17 vertices.
# But 5*3=15 < 17 and 3*4+5=17, so perfect factorizations include [3,3,3,3,5] only.
# Pure [3,3,3,3,3] = 15 < 17 vertices — NOT a perfect factorization!
# So alpha_5 counts [3,3,3,3,5] partitions (15 vertices + one 5-cycle? = 20 > 17... wait)
# Actually alpha_5 = 5-packings of disjoint odd cycles, NOT necessarily covering all vertices.
# kmax = floor(17/3) = 5, so 5 disjoint cycles using 5*3=15 vertices minimum, ≤ 17.
# Possible: 5 3-cycles (15 verts), 4 3-cycles + 1 5-cycle (17 verts = perfect!).
#
# Use 0/1 knapsack for all-3-cycle 5-tuples:
t0 = time.time()
dp = [0] * (1 << n); dp[0] = 1
for m, c in tc3:
    for mask in range((1 << n)-1, -1, -1):
        if dp[mask] and not (mask & m):
            dp[mask|m] += dp[mask] * c

# alpha5_33333 = sum over 15-vertex masks of dp[mask] (5 disjoint 3-cycles)
# alpha5_33335 = sum over all-17-vertex masks (perfect [3,3,3,3,5] coverings) computed below
alpha5_33333 = sum(dp[mask] for mask in range(1<<n) if bin(mask).count('1') == 15 and dp[mask]) // 5
# Wait — knapsack counts ORDERED 5-tuples? No — it builds up by adding cycles one at a time
# but each mask gets multiple cycle choices, so it's counting ordered... actually no:
# dp[0]=1. After processing each 3-cycle (m,c), we allow adding it or not.
# This is the 0/1 knapsack: each MASK can be added AT MOST ONCE. But cc[m] counts multiple
# cycles with the same vertex set... and we have tc3 with distinct masks.
# So dp[mask] after all = sum over subsets of {cycle masks} that are disjoint and union = mask,
# weighted by product of cc values. This counts UNORDERED sets of cycle masks.
# For 5 disjoint 3-cycles covering 15 vertices: dp[15-vertex mask] already = unordered count.
# alpha5_33333 = sum of dp[15-vertex masks] (no division needed)

alpha5_33333 = sum(dp[mask] for mask in range(1 << n) if bin(mask).count('1') == 15)

# For alpha5_33335 ([3,3,3,3,5]): use dp from above (4 disjoint 3-cycles covering 12 verts),
# then count 5-cycles in the remaining 5 vertices.
cc5 = [cc[m] if bin(m).count('1') == 5 else 0 for m in range(1 << n)]
# dp_4 = 4-3-cycle knapsack
dp4 = [0] * (1 << n); dp4[0] = 1
for m, c in tc3:
    for mask in range((1 << n)-1, -1, -1):
        if dp4[mask] and not (mask & m):
            dp4[mask|m] += dp4[mask] * c

alpha5_33335 = 0
for mask12 in range(1 << n):
    if bin(mask12).count('1') != 12: continue
    if not dp4[mask12]: continue
    comp = (~mask12) & full
    if bin(comp).count('1') != 5: continue
    alpha5_33335 += dp4[mask12] * cc5[comp]
# dp4[mask12] already counts 4-cycle unordered sets; each [3,3,3,3,5] is counted 4 times
# (once for each choice of which 4 3-cycles form the "base" — but wait, the 5-cycle is fixed
# as the complement). Actually alpha5_33335 as written = unordered 4-tuples of 3-cycles × cc5[comp]
# This counts each [3,3,3,3,5] exactly once: the 4 3-cycles form an unordered set (from dp4),
# and the 5-cycle is uniquely determined by the 5-vertex complement.

alpha5 = alpha5_33333 + alpha5_33335
print(f"α₅ = {alpha5:,}  (33333: {alpha5_33333:,}, 33335: {alpha5_33335:,})  ({time.time()-t0:.2f}s)")

# ── Verify H ─────────────────────────────────────────────────────
alphas = {1: alpha1, 2: alpha2, 3: alpha3, 4: alpha4, 5: alpha5}
H_ocf = 1 + sum(2**k * v for k, v in alphas.items())
H_ref = 13_689_269_499
print(f"\nH (OCF) = {H_ocf:,}")
print(f"H (ref) = {H_ref:,}")
print(f"Match: {H_ocf == H_ref}")

# ── Mechanism table ──────────────────────────────────────────────
print(f"\nTerm breakdown:")
print(f"  2^0 * α₀ = {1:>20,}")
for k, v in sorted(alphas.items()):
    print(f"  2^{k} * α_{k} = {2**k * v:>20,}  (α_{k} = {v:,})")
print(f"  Total    = {H_ocf:>20,}")

ratio = alpha1 / (2 * alpha2)
print(f"\nα₁/(2α₂) = {ratio:.4f}  (< 1: α₂-dominant)")
print(f"α₃/(α₂/2) = {alpha3 / (alpha2/2):.4f}")
