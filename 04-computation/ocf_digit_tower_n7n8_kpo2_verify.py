#!/usr/bin/env python3
"""
ocf_digit_tower_n7n8_kpo2_verify.py  -- ADVERSARIAL VERIFIER (kind-pasteur-2026-06-10-S2)

Independent verification of THM-466 claims at n=7 (FULL labeled census, 2^21)
and n=8 (random sample; REDUCED to 20,000 from the worker's 200,000 for the
15-minute verifier budget -- stated explicitly).

FRESH METHODS (different from worker's Held-Karp DP):
  * directed odd cycles: per odd-size subset, all cyclic orders with min fixed
    first (DIRECTED cycles; a subset can host several -- repo mistake guard);
  * H by inclusion-exclusion walk counting (exact int64), cross-checked by pure
    permutation brute force;
  * n=7 v2(H-1) full distribution derived from the alpha census via the digit
    identity H-1 = 2*(alpha1 + 2*alpha2) (the identity itself is verified at
    n=7 on 20,000 sampled codes + 2 transitive codes with IE-computed H, and on
    ALL tournaments n<=6 in the companion script);
  * GLOBAL analytic cross-check at n=7: Sum(1+2*alpha1+4*alpha2) over all 2^21
    codes must equal Sum H = 7! * 2^(21-6) = 165,150,720 (walk-count total),
    tying my alpha census to the known labeled Hamiltonian-path total.

Encoding (mine): bit p(u,v) for u<v lexicographic = 1 iff u->v. T^op = bitwise
complement => reversal invariance <=> palindromic alpha arrays.
"""
import sys, math, itertools, time
import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
FAIL = 0
T0 = time.time()

def chk(cond, msg):
    global FAIL
    if cond:
        print("  PASS  " + msg)
    else:
        FAIL += 1
        print("  *** FAIL *** " + msg)

def pair_index(n):
    idx, p = {}, 0
    for u in range(n):
        for v in range(u + 1, n):
            idx[(u, v)] = p; p += 1
    return idx, p

def adj_from_code(code, n, idx):
    A = [[0] * n for _ in range(n)]
    for (u, v), p in idx.items():
        if (code >> p) & 1: A[u][v] = 1
        else: A[v][u] = 1
    return A

def ham_paths_perm(A, n):
    cnt = 0
    for perm in itertools.permutations(range(n)):
        if all(A[perm[i]][perm[i + 1]] for i in range(n - 1)):
            cnt += 1
    return cnt

def batched_power(M, e):
    result, base = None, M
    while e:
        if e & 1:
            result = base if result is None else np.matmul(result, base)
        e >>= 1
        if e: base = np.matmul(base, base)
    return result

def H_IE_batch(codes, n, idx):
    N = len(codes)
    A = np.zeros((N, n, n), dtype=np.int64)
    for (u, v), p in idx.items():
        b = (codes >> p) & 1
        A[:, u, v] = b
        A[:, v, u] = 1 - b
    H = np.zeros(N, dtype=np.int64)
    for r in range(1, n + 1):
        sign = (-1) ** (n - r)
        for S in itertools.combinations(range(n), r):
            P = batched_power(A[:, S, :][:, :, S], n - 1)
            H += sign * P.sum(axis=(1, 2))
    return H

def v2_of(x):
    lo = (x & -x).astype(np.float64)
    return np.round(np.log2(lo)).astype(np.int64)

# =====================================================================
n = 7
idx7, m7 = pair_index(n)
N = 2 ** m7
print("=" * 78)
print(f"PART A: n=7 FULL labeled census ({N} tournaments), vectorized alpha census")
print("=" * 78)
codes = np.arange(N, dtype=np.int64)
B = {}
for (u, v), p in idx7.items():
    b = ((codes >> p) & 1).astype(np.uint8)
    B[(u, v)] = b
    B[(v, u)] = (1 - b).astype(np.uint8)

# c3 with per-triple indicators (needed for alpha2)
cyc3 = {}
c3 = np.zeros(N, dtype=np.int32)
for t in itertools.combinations(range(7), 3):
    i, j, k = t
    cyc3[t] = (B[(i, j)] & B[(j, k)] & B[(k, i)]) | (B[(i, k)] & B[(k, j)] & B[(j, i)])
    c3 += cyc3[t]
print(f"  [t={time.time()-T0:.0f}s] c3 census done")

# c5: 21 subsets x 24 directed cyclic orders (totals only; no (3,5) pairs at n=7)
c5 = np.zeros(N, dtype=np.int32)
for sub in itertools.combinations(range(7), 5):
    first = sub[0]
    for perm in itertools.permutations(sub[1:]):
        seq = (first,) + perm
        ind = B[(seq[0], seq[1])].copy()
        for i in range(1, 5):
            ind &= B[(seq[i], seq[(i + 1) % 5])]
        c5 += ind
print(f"  [t={time.time()-T0:.0f}s] c5 census done")

# c7: 720 directed cyclic orders on all 7 vertices
c7 = np.zeros(N, dtype=np.int32)
sub = tuple(range(7))
for perm in itertools.permutations(sub[1:]):
    seq = (0,) + perm
    ind = B[(seq[0], seq[1])].copy()
    for i in range(1, 7):
        ind &= B[(seq[i], seq[(i + 1) % 7])]
    c7 += ind
print(f"  [t={time.time()-T0:.0f}s] c7 census done")

alpha1 = c3 + c5 + c7
# alpha2: only (3,3) pairs possible at n=7 (3+5=8>7)
alpha2 = np.zeros(N, dtype=np.int32)
trips = list(cyc3.keys())
npairs33 = 0
for a in range(len(trips)):
    for bidx in range(a + 1, len(trips)):
        if not (set(trips[a]) & set(trips[bidx])):
            alpha2 += (cyc3[trips[a]] & cyc3[trips[bidx]]).astype(np.int32)
            npairs33 += 1
print(f"  [t={time.time()-T0:.0f}s] alpha2 census done ({npairs33} disjoint triple pairs; expect 70)")
chk(npairs33 == 70, "exactly 70 disjoint (3,3) subset pairs at n=7")

# aggregate identities (analytic, pre-registered)
chk(int(c3.sum()) == 35 * 2 ** 19, f"Sum c3 == C(7,3)*2^19 == {35*2**19}")
chk(int(c5.sum()) == 21 * 24 * 2 ** 16, f"Sum c5 == C(7,5)*24*2^16 == {21*24*2**16}")
chk(int(c7.sum()) == 720 * 2 ** 14, f"Sum c7 == 720*2^14 == {720*2**14}")
chk(int(alpha2.sum()) == 70 * 2 ** 17, f"Sum alpha2 == 70*4*2^15 == 70*2^17 == {70*2**17}")
tot = N + 2 * int(alpha1.sum()) + 4 * int(alpha2.sum())
chk(tot == math.factorial(7) * 2 ** 15,
    f"GLOBAL: Sum(1+2*a1+4*a2) == 7!*2^15 == {math.factorial(7)*2**15} (known Sum H), got {tot}")

# reversal invariance, full census
chk(bool(np.array_equal(alpha1, alpha1[::-1])) and bool(np.array_equal(alpha2, alpha2[::-1])),
    "FULL-census reversal invariance at n=7: alpha1 and alpha2 are palindromes under T -> T^op")

# transitive count
ntrans = int(((alpha1 == 0) & (alpha2 == 0)).sum())
chk(ntrans == math.factorial(7), f"#(alpha-vector == 0) == 7! == {math.factorial(7)} (transitive)")

# v2(H-1) distribution via identity (identity verified by sampled IE-H below + full n<=6)
t_ = alpha1.astype(np.int64) + 2 * alpha2.astype(np.int64)
pos = t_[t_ > 0]
v2H = 1 + v2_of(pos)
distr = {int(v): int(c) for v, c in zip(*np.unique(v2H, return_counts=True))}
print("  v2(H-1) distribution n=7 (full census): " +
      ", ".join(f"v2={k}: {v}" for k, v in sorted(distr.items())) + f", inf: {ntrans}")
p_odd = int((alpha1 % 2 == 1).sum())
print(f"  P(alpha1 odd) = {p_odd}/{N} = {p_odd/N:.6f}")
chk(distr.get(1, 0) == p_odd, "#(v2(H-1)==1) == #(alpha1 odd)")
c57odd = int(((c5 + c7) % 2 == 1).sum())
print(f"  #((c5+c7) odd) = {c57odd}/{N} = {c57odd/N:.6f}")
c5odd = int((c5 % 2 == 1).sum()); c7odd = int((c7 % 2 == 1).sum())
print(f"  #(c5 odd) = {c5odd}/{N};  #(c7 odd) = {c7odd}/{N}")

# Kendall + score classes
scores = np.zeros((N, 7), dtype=np.int64)
for u in range(7):
    for v in range(7):
        if u != v: scores[:, u] += B[(u, v)]
kend = math.comb(7, 3) - (scores * (scores - 1) // 2).sum(axis=1)
chk(bool((kend == c3).all()), "Kendall: c3 == C(7,3) - sum C(s_i,2) on all 2^21")
skey = np.sort(scores, axis=1)
key = np.zeros(N, dtype=np.int64)
for i in range(7):
    key = key * 8 + skey[:, i]
uniq, inv = np.unique(key, return_inverse=True)
cnt = np.bincount(inv).astype(np.float64)
def noncon(par):
    odd = np.bincount(inv, weights=par.astype(np.float64))
    return int(((odd > 0.5) & (odd < cnt - 0.5)).sum())
nc_a1 = noncon((alpha1 % 2)); nc_c3 = noncon((c3 % 2)); nc_c5 = noncon((c5 % 2)); nc_c7 = noncon((c7 % 2))
print(f"  score classes at n=7: {len(uniq)} (expect 59 = A000571(7))")
print(f"  non-constant parity classes: alpha1: {nc_a1}/59, c5: {nc_c5}/59, c7: {nc_c7}/59, c3: {nc_c3}/59")
chk(len(uniq) == 59, "59 distinct sorted score sequences at n=7")
chk(nc_c3 == 0, "c3 parity constant on every score class")

# free memory before IE
del B, cyc3
print(f"  [t={time.time()-T0:.0f}s] PART A done")

print("=" * 78)
print("PART B: n=7 identity spot-verification: H by inclusion-exclusion on 20,002 codes")
print("=" * 78)
rng = np.random.default_rng(314159)
sample = np.unique(np.concatenate([rng.integers(0, N, size=20000),
                                   np.array([0, N - 1], dtype=np.int64)]))
H_s = H_IE_batch(sample, 7, idx7)
rhs = 1 + 2 * alpha1[sample].astype(np.int64) + 4 * alpha2[sample].astype(np.int64)
bad = int((H_s != rhs).sum())
chk(bad == 0, f"H == 1+2*alpha1+4*alpha2 on all {len(sample)} sampled n=7 codes "
              f"(incl. both transitive extremes; failures={bad})")
chk(bool((H_s % 2 == 1).all()) and int(H_s.min()) == 1, "sampled H odd, min == 1 (transitive included)")
# brute-force permutation cross-check of the IE engine
ok = True
for code in sample[rng.integers(0, len(sample), size=8)]:
    A = adj_from_code(int(code), 7, idx7)
    if ham_paths_perm(A, 7) != int(H_s[np.searchsorted(sample, code)]):
        ok = False
chk(ok, "8 random sampled codes: IE walk-count H == pure 7!-permutation brute force")
print(f"  [t={time.time()-T0:.0f}s] PART B done")
del alpha1, alpha2, c3, c5, c7, t_, scores, kend, skey, key, inv

print("=" * 78)
print("PART C: n=8, 20,000 random tournaments (REDUCED from worker's 200,000; seed 271828)")
print("=" * 78)
n8 = 8
idx8, m8 = pair_index(n8)
rng8 = np.random.default_rng(271828)
codes8 = rng8.integers(0, 2 ** m8, size=20000, dtype=np.int64)
N8 = len(codes8)
B8 = {}
for (u, v), p in idx8.items():
    b = ((codes8 >> p) & 1).astype(np.uint8)
    B8[(u, v)] = b
    B8[(v, u)] = (1 - b).astype(np.uint8)
cyc3_8 = {}
c3_8 = np.zeros(N8, dtype=np.int32)
for t in itertools.combinations(range(8), 3):
    i, j, k = t
    cyc3_8[t] = (B8[(i, j)] & B8[(j, k)] & B8[(k, i)]) | (B8[(i, k)] & B8[(k, j)] & B8[(j, i)])
    c3_8 += cyc3_8[t]
cyc5cnt = {}
c5_8 = np.zeros(N8, dtype=np.int32)
for sub in itertools.combinations(range(8), 5):
    acc = np.zeros(N8, dtype=np.uint8)
    first = sub[0]
    for perm in itertools.permutations(sub[1:]):
        seq = (first,) + perm
        ind = B8[(seq[0], seq[1])].copy()
        for i in range(1, 5):
            ind &= B8[(seq[i], seq[(i + 1) % 5])]
        acc += ind
    cyc5cnt[sub] = acc
    c5_8 += acc
c7_8 = np.zeros(N8, dtype=np.int32)
for sub in itertools.combinations(range(8), 7):
    first = sub[0]
    for perm in itertools.permutations(sub[1:]):
        seq = (first,) + perm
        ind = B8[(seq[0], seq[1])].copy()
        for i in range(1, 7):
            ind &= B8[(seq[i], seq[(i + 1) % 7])]
        c7_8 += ind
alpha1_8 = c3_8 + c5_8 + c7_8
# alpha2 = (3,3) pairs + (3,5) pairs
alpha2_33 = np.zeros(N8, dtype=np.int32)
trips8 = list(cyc3_8.keys())
np33 = 0
for a in range(len(trips8)):
    for bb in range(a + 1, len(trips8)):
        if not (set(trips8[a]) & set(trips8[bb])):
            alpha2_33 += (cyc3_8[trips8[a]] & cyc3_8[trips8[bb]]).astype(np.int32)
            np33 += 1
alpha2_35 = np.zeros(N8, dtype=np.int32)
allv = set(range(8))
for t in trips8:
    comp = tuple(sorted(allv - set(t)))
    alpha2_35 += cyc3_8[t].astype(np.int32) * cyc5cnt[comp].astype(np.int32)
alpha2_8 = alpha2_33 + alpha2_35
chk(np33 == 280, f"280 disjoint (3,3) subset pairs at n=8 (got {np33}); 56 (3,5) pairs by complement")
print(f"  [t={time.time()-T0:.0f}s] n=8 cycle censuses done; computing IE H for 20,000 codes...")
H8 = H_IE_batch(codes8, 8, idx8)
rhs8 = 1 + 2 * alpha1_8.astype(np.int64) + 4 * alpha2_8.astype(np.int64)
bad8 = int((H8 != rhs8).sum())
chk(bad8 == 0, f"H == 1+2*alpha1+4*alpha2 on all {N8} sampled n=8 tournaments (failures={bad8})")
chk(bool((H8 % 2 == 1).all()), "all sampled n=8 H odd (Redei)")
act35 = int((alpha2_35 > 0).sum())
print(f"  (3,5)-pair channel ACTIVE (alpha2_35>0) in {act35}/{N8} = {act35/N8:.4f} "
      f"(worker: 174495/200000 = 0.8725)")
c57odd8 = int(((c5_8 + c7_8) % 2 == 1).sum())
print(f"  (c5+c7) odd in {c57odd8}/{N8} = {c57odd8/N8:.4f} (worker: 98802/200000 = 0.4940)")
# brute-force permutation cross-check at n=8
ok8 = True
for j in rng8.integers(0, N8, size=3):
    A = adj_from_code(int(codes8[j]), 8, idx8)
    if ham_paths_perm(A, 8) != int(H8[j]):
        ok8 = False
chk(ok8, "3 random n=8 codes: IE H == pure 8!-permutation brute force")
# reversal invariance spot check at n=8 (complement codes)
sub_idx = rng8.integers(0, N8, size=200)
op_codes = (2 ** m8 - 1) - codes8[sub_idx]
both = np.concatenate([codes8[sub_idx], op_codes])
Bb = {}
for (u, v), p in idx8.items():
    b = ((both >> p) & 1).astype(np.uint8)
    Bb[(u, v)] = b; Bb[(v, u)] = (1 - b).astype(np.uint8)
def quick_alpha1(Bx, NN):
    a = np.zeros(NN, dtype=np.int32)
    for t in itertools.combinations(range(8), 3):
        i, j, k = t
        a += (Bx[(i, j)] & Bx[(j, k)] & Bx[(k, i)]) | (Bx[(i, k)] & Bx[(k, j)] & Bx[(j, i)])
    for s in (5, 7):
        for sub in itertools.combinations(range(8), s):
            first = sub[0]
            for perm in itertools.permutations(sub[1:]):
                seq = (first,) + perm
                ind = Bx[(seq[0], seq[1])].copy()
                for i in range(1, s):
                    ind &= Bx[(seq[i], seq[(i + 1) % s])]
                a += ind
    return a
a1b = quick_alpha1(Bb, len(both))
chk(bool((a1b[:200] == a1b[200:]).all()), "alpha1(T) == alpha1(T^op) on 200 random n=8 pairs")

# score-class snapshot at n=8 (sample-based; counts are lower bounds only)
sc8 = np.zeros((N8, 8), dtype=np.int64)
for u in range(8):
    for v in range(8):
        if u != v: sc8[:, u] += B8[(u, v)]
sk8 = np.sort(sc8, axis=1)
k8 = np.zeros(N8, dtype=np.int64)
for i in range(8):
    k8 = k8 * 9 + sk8[:, i]
u8, inv8 = np.unique(k8, return_inverse=True)
cnt8 = np.bincount(inv8).astype(np.float64)
odd8 = np.bincount(inv8, weights=(alpha1_8 % 2).astype(np.float64))
nc8 = int(((odd8 > 0.5) & (odd8 < cnt8 - 0.5)).sum())
print(f"  n=8 sample: {len(u8)} score classes seen, alpha1-parity non-constant in >= {nc8} "
      f"(LOWER BOUND from 20k sample; worker: 114/166 from 200k)")
print(f"  [t={time.time()-T0:.0f}s] PART C done")

print("=" * 78)
print(f"SUMMARY (n7n8): {'ALL CHECKS PASSED' if FAIL == 0 else str(FAIL) + ' FAILURES'}")
print("=" * 78)
sys.exit(0 if FAIL == 0 else 1)
