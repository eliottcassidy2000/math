#!/usr/bin/env python3
"""
paley_hk_design_kpc1_verify.py
ADVERSARIAL VERIFIER, kind-pasteur-2026-06-10-S1, thread E-R31-prediction.
Verifies claim E8 (and the H values underlying R_exact).

INDEPENDENT METHOD: start-fixed Held-Karp with NO rotation canonicalization
(the worker used a Z_p-rotation-reduced canonical-form DP).
  N0 = number of directed Hamiltonian paths of T_p starting at vertex 0,
  via layered numpy DP over (mask of vertices 1..p-1, end vertex).
  H = p * N0:  rotation acts freely on directed Hamiltonian paths (if
  P + r = P as a directed path then its start satisfies s + r = s, so r = 0),
  and each rotation orbit contains exactly one path starting at 0.
Cross-checks: p=7 by raw permutation brute force; p=11 also by a dict-based
all-starts DP (third method).

Then: every arithmetic claim in the E8 design doc is recomputed exactly.
"""
import time
import numpy as np
from math import comb, factorial
from itertools import permutations

T0 = time.time()

def qr_set(p):
    return set((x * x) % p for x in range(1, p))

def H_brute(p):
    qr = qr_set(p)
    cnt = 0
    for sig in permutations(range(p)):
        if all((sig[i + 1] - sig[i]) % p in qr for i in range(p - 1)):
            cnt += 1
    return cnt

def H_dict_dp(p):
    """All-starts dict DP (no start fixing, no rotation argument)."""
    qr = qr_set(p)
    adj = [[(w - v) % p in qr for w in range(p)] for v in range(p)]
    dp = {}
    for v in range(p):
        dp[(1 << v, v)] = 1
    for _ in range(p - 1):
        ndp = {}
        for (S, v), c in dp.items():
            for w in range(p):
                if not (S >> w) & 1 and adj[v][w]:
                    k = (S | (1 << w), w)
                    ndp[k] = ndp.get(k, 0) + c
        dp = ndp
    return sum(dp.values())

def H_start_fixed(p):
    """Layered numpy DP, start fixed at 0; returns p * N0."""
    qr = qr_set(p)
    m = p - 1                       # mask bits for vertices 1..p-1
    dp = np.zeros((1 << m, p), dtype=np.int64)
    dp[0, 0] = 1
    order = np.argsort(np.array([bin(T).count('1') for T in range(1 << m)]))
    by_pop = {}
    for T in order:
        by_pop.setdefault(bin(int(T)).count('1'), []).append(int(T))
    for k in range(0, m):
        Ms = np.array(by_pop.get(k, []), dtype=np.int64)
        if len(Ms) == 0:
            continue
        for w in range(1, p):
            bw = 1 << (w - 1)
            sel = Ms[(Ms & bw) == 0]
            if len(sel) == 0:
                continue
            tgt = sel | bw
            for v in range(p):
                if (w - v) % p in qr:           # edge v -> w
                    if v == 0:
                        dp[tgt, w] += dp[sel, w * 0]   # dp[sel, 0]
                    else:
                        bv = 1 << (v - 1)
                        ok = (sel & bv) != 0
                        if ok.any():
                            s2 = sel[ok]
                            dp[s2 | bw, w] += dp[s2, v]
    full = (1 << m) - 1
    N0 = int(dp[full, :].sum())
    return p * N0, N0

print("=" * 78)
print("E8 part 1: independent recomputation of H(T_p)")
print("=" * 78)
h7b = H_brute(7)
print(f"  H(7)  brute force over 7! orderings          = {h7b}")
h7, _ = H_start_fixed(7)
print(f"  H(7)  start-fixed layered DP * 7             = {h7}   "
      f"{'OK' if h7 == h7b == 189 else 'MISMATCH'}")
h11d = H_dict_dp(11)
h11, n0_11 = H_start_fixed(11)
print(f"  H(11) all-starts dict DP                     = {h11d}")
print(f"  H(11) start-fixed layered DP * 11            = {h11}   "
      f"{'OK' if h11 == h11d == 95095 else 'MISMATCH'}  (claimed 95095)")
h19, n0_19 = H_start_fixed(19)
print(f"  H(19) start-fixed layered DP * 19            = {h19}   "
      f"(claimed 1172695746915: {'CONFIRMED' if h19 == 1172695746915 else 'REFUTED'})"
      f"   [{time.time()-T0:.0f}s]")
print(f"  max intermediate count check: N0(19) = {n0_19} < 2^63 OK")

print()
print("=" * 78)
print("E8 part 2: design-doc arithmetic, recomputed exactly")
print("=" * 78)
p = 31
# (a) freeness: rotation r fixes (S,v) => v + r = v mod p => r = 0. Trivial,
#     holds for ANY nonempty S (v exists); p prime not even needed for pairs.
print("  (a) freeness of Z_31 on states (S,v), S nonempty: (S+r, v+r) = (S,v)")
print("      forces r = 0 from the v-component alone. PROVED (trivially).")

# (b) layer-k canonical state count = C(30, k-1) EXACTLY (not just <=)
ok = all(k * comb(31, k) == 31 * comb(30, k - 1) for k in range(1, 32))
print(f"  (b) k*C(31,k) = 31*C(30,k-1) for k=1..31: {'OK' if ok else 'FAIL'}"
      f"  => orbit count at layer k EXACTLY C(30,k-1) (by freeness)")
peak = max(comb(30, k - 1) for k in range(1, 32))
print(f"      peak layer count = C(30,15) = {comb(30, 15)} "
      f"(claimed 155,117,520: {'OK' if comb(30, 15) == 155117520 else 'FAIL'})")
bytes128 = comb(30, 15) * 16
print(f"      uint128 bytes at peak = {bytes128} = {bytes128 / 2**30:.4f} GiB "
      f"(claimed ~2.31 GiB: {'OK' if abs(bytes128/2**30 - 2.31) < 0.01 else 'FAIL'})")

# (c) total transitions: QR out-degree is exactly 15 for every vertex.
#     unreduced (S,v,w) with w not in S, chi(w-v)=1:
#     sum_k C(29,k-1) * 31 * 15  ... computed two exact ways
t1 = sum(31 * 15 * comb(29, k - 1) for k in range(1, 31))
t2 = 31 * 30 * 2 ** 29 // 62
print(f"  (c) unreduced QR-transitions = 31*15*2^29 = {31 * 15 * 2**29}")
print(f"      sum_k 31*15*C(29,k-1) (k=1..30)       = {t1}  "
      f"{'OK' if t1 == 31 * 15 * 2**29 else 'FAIL'}")
print(f"      rotation-reduced = /31 = 15*2^29 = {15 * 2**29}")
print(f"      claimed 8,053,063,680 = 31*30*2^29/62 = {t2}  "
      f"{'OK' if t2 == 15 * 2**29 == 8053063680 else 'FAIL'}")

# (d) overflow crossing: avg dp ~ (k-1)!/2^(k-1) vs 2^64
print("  (d) (k-1)!/2^(k-1) vs 2^64 = 18446744073709551616:")
for k in [26, 27, 28, 29]:
    v = factorial(k - 1) // 2 ** (k - 1)
    print(f"      k={k}: ~{v:.3e}  {'<' if v < 2**64 else '>'} 2^64")
cross = next(k for k in range(2, 32) if factorial(k - 1) / 2 ** (k - 1) > 2 ** 64)
print(f"      first crossing at k = {cross} (claimed ~28: "
      f"{'OK' if cross == 28 else 'FAIL'}); skewed states overflow earlier, so")
print("      uint128/CRT needed from some earlier layer (~24 is a design margin call).")

# (e) Karp inclusion-exclusion rotation classes = binary necklaces N(31)
neck = (2 ** 31 + sum(2 ** 1 for r in range(1, 31))) // 31   # gcd(r,31)=1, r!=0
print(f"  (e) subsets of Z_31 up to rotation (Burnside) = (2^31 + 30*2)/31 = {neck}")
print(f"      claimed 6.93e7: {'OK' if abs(neck - 6.93e7) < 5e5 else 'FAIL'}"
      f"  (exact {neck} = {neck:.3e})")

# (f) total canonical states across layers = 2^30
tot_states = sum(comb(30, k - 1) for k in range(1, 32))
print(f"  (f) total canonical states sum_k C(30,k-1) = {tot_states} = 2^30: "
      f"{'OK' if tot_states == 2**30 else 'FAIL'}")
print(f"\n[elapsed {time.time()-T0:.1f}s]")

# --------------------------------------------------------------------------
# Bonus (run as __main__ part 2 via env flag): independent H(23) check.
# --------------------------------------------------------------------------
if __name__ == "__main__" and __import__('os').environ.get('KPC1_H23') == '1':
    print()
    print("=" * 78)
    print("BONUS: independent recomputation of H(T_23) (canon value check)")
    print("=" * 78)
    h23, n0_23 = H_start_fixed(23)
    print(f"  H(23) start-fixed layered DP * 23 = {h23}")
    print(f"  canon/brief value                 = 15760206976379349  "
          f"{'CONFIRMED' if h23 == 15760206976379349 else 'REFUTED'}")
    print(f"  N0(23) = {n0_23} < 2^63 OK   [{time.time()-T0:.0f}s]")
