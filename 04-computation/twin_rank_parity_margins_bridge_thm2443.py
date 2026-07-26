#!/usr/bin/env python3
r"""
twin_rank_parity_margins_bridge_thm2443.py
(kind-pasteur-2026-07-26-S131; companion to THM-2443)

Typing (MISTAKE-268): C_all = A014574 (twin centers), C_6 = C_all \ {4},
K = C_6/6 = A002822 (twin ranks: 6k-1 and 6k+1 both prime).

Exact companion for the parent-fibre parity lemma, the dyadic
representation margins, and the boundary-crossing (Seymour) bridge.

Universe: all twin ranks k with center 6k <= LIMIT = 10^8
(the same census universe as THM-2422 eq (36)).

Computes, all exact:
  1. R(k) = #{(a,b) in K^2 : a+b = k}   (ordered parent count)
     by real-FFT convolution of the K indicator, rounded to integers
     (counts << 2^52, rounding exact), CROSS-CHECKED against a direct
     double-loop brute force for every k in K, k <= 10^4.
  2. Weak law R(k) > 0 for all k in K, k >= 2, and distinct-parent law
     R(k) - [k even and k/2 in K] > 0 for all k >= 3
     (independent-path reproduction of THM-2422 (36): convolution vs scan).
  3. Parity lemma verification: R(k) odd  <=>  k even and k/2 in K.
  4. Dyadic margins: min/mean R over K-terms in [2^j, 2^{j+1}).
  5. Doubling nodes: k in K with 2k in K.
  6. Seymour bridge range: the largest V0 such that the boundary-crossing
     statement (OEIS A002822, J. Seymour comment 2026-05-20) is
     FINITE-EXACT verified for all 1 <= V <= V0 from this census.
  7. Gap partner table: counts of consecutive-rank gaps g with the
     Hardy-Littlewood-type normalization
        w(g) = prod_{p>=5, p|g} (p-2)/(p-4)
             * prod_{p>=5, p | 9g^2-1} (p-3)/(p-4)
     (primes up to 199 exactly; the tail is O(1/p) per prime).

Controls:
  * raw centers 4 and 6 have no two-earlier-center representation, and
    4 is summand-inert mod 6 (positive controls for the startup hole);
  * center 12 (k=2) has exactly the doubled representation 6+6 (R(2)=1);
  * FFT vs brute force must agree exactly on all K-terms <= 10^4.

Reproduction: python 04-computation/twin_rank_parity_margins_bridge_thm2443.py
"""
import numpy as np

LIMIT = 100_000_000

def fail(msg):
    raise SystemExit("CHECK FAILED: " + msg)

# ---- sieve and K ----
sieve = np.ones(LIMIT + 3, dtype=bool)
sieve[:2] = False
for p in range(2, int((LIMIT + 2) ** 0.5) + 1):
    if sieve[p]:
        sieve[p * p:: p] = False
mid = (np.where(sieve[:-2] & sieve[2:])[0] + 1).astype(np.int64)
if mid[0] != 4 or mid[1] != 6:
    fail("startup centers")
K = (mid[mid >= 6] // 6)
if (6 * K != mid[mid >= 6]).any():
    fail("center divisibility typing (MISTAKE-268)")
M = int(K[-1])
print(f"universe: centers <= {LIMIT}; |K| = {len(K)}; max k = {M}")

Kset = np.zeros(M + 2, dtype=bool)
Kset[K] = True

# ---- raw-unit startup controls ----
Cset = set(mid[mid <= 1000].tolist())
if any(a in Cset and (v - a) in Cset and a <= v - a for v in (4, 6) for a in range(2, v // 2 + 1)):
    fail("4/6 unexpectedly representable")
if any((c - 4) % 6 == 0 for c in Cset if c >= 12):
    fail("4 not summand-inert")
print("controls: centers 4 and 6 rep-free; 4 summand-inert mod 6: PASS")

# ---- FFT ordered parent counts ----
size = 1
while size < 2 * M + 2:
    size *= 2
f = np.zeros(size, dtype=np.float64)
f[K] = 1.0
F = np.fft.rfft(f)
conv = np.fft.irfft(F * F, size)[: M + 1]
R = np.rint(conv).astype(np.int64)
del f, F, conv

# ---- brute-force cross-check on k <= 10^4 ----
smallK = [int(k) for k in K if k <= 10_000]
smallset = set(smallK)
for k in smallK:
    r = sum(1 for a in smallK if a < k and (k - a) in smallset)
    if r != int(R[k]):
        fail(f"FFT/brute mismatch at k={k}: {R[k]} vs {r}")
print(f"FFT vs brute-force cross-check on {len(smallK)} ranks <= 10^4: PASS")

targets = K[K >= 2]
Rt = R[targets]
if int(R[2]) != 1:
    fail("R(2) != 1 (center 12 must have exactly 6+6)")
weak_fail = targets[Rt == 0]
halfok = ((targets % 2 == 0) & Kset[np.where(targets % 2 == 0, targets // 2, 1)])
dist_fail = targets[(Rt - halfok.astype(np.int64) == 0) & (targets >= 3)]
print(f"weak law failures (k>=2): {weak_fail.tolist()}")
print(f"distinct-parent law failures (k>=3): {dist_fail.tolist()}")
if len(weak_fail) or len(dist_fail):
    fail("parent law failure")

# ---- parity lemma ----
viol = int(np.count_nonzero((Rt % 2 == 1) != halfok))
print(f"parity lemma (R odd <=> k even and k/2 in K) violations: {viol}")
if viol:
    fail("parity lemma")

# ---- dyadic margins ----
print("dyadic_margins: [lo,hi) min_R argmin mean_R")
lo = 2
while lo <= M:
    hi = min(2 * lo, M + 1)
    sel = (targets >= lo) & (targets < hi)
    if sel.any():
        w = Rt[sel]
        am = int(targets[sel][int(np.argmin(w))])
        print(f"  [{lo},{hi}): {int(w.min())} {am} {w.mean():.1f}")
    lo = hi

# ---- doubling nodes ----
dbl = [int(x) for x in K if 2 * x <= M and Kset[2 * x]]
print(f"doubling_nodes (k and 2k in K): {len(dbl)}; first 25: {dbl[:25]}")

# ---- Seymour bridge range ----
# For 2 <= V <= M-1: w = least K-term > V has distinct parents (a,b), a<b<w,
# and b <= V by minimality.  V = 1: (1,1,2).  Verified range:
V0 = M - 1
print(f"seymour_boundary_crossing verified for all 1 <= V <= {V0}")

# ---- gap partner table ----
gaps = np.diff(K)
from collections import Counter
gh = Counter(gaps.tolist())
PR = [p for p in range(5, 200) if all(p % q for q in range(2, int(p ** .5) + 1))]

def w_of(g):
    w = 1.0
    for p in PR:
        if g % p == 0:
            w *= (p - 2) / (p - 4)
        elif (9 * g * g - 1) % p == 0:
            w *= (p - 3) / (p - 4)
    return w

print("gap_table: g count w count/w")
for g in sorted(gh)[:40]:
    c = gh[g]
    w = w_of(g)
    print(f"  {g} {c} {w:.4f} {c / w:.1f}")

print("ALL CHECKS PASSED")
