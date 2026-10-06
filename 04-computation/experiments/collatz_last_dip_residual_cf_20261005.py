#!/usr/bin/env python3
"""
collatz_last_dip_residual_cf_20261005.py -- residual lengths of the last-dip lemma for large X via the
continued fraction of log2 3.

After the lemma is known for every source m* <= X (sweep, or CITED: Terras's coefficient stopping time
conjecture t(n) = tau(n) verified for 2 <= n <= 2.8e19 by Rozier--Terracol 2025, Corollary 5.2; a last-dip
violation at (m*, n) is a CST failure at m*), a violation needs m* > X, hence n > X, and by (R1)
      1 - {l log2 3} < log2( (1 + 2/(3n-3)) e^{(l-1)/(3n)} ) < (l + 2)/(3 n ln 2)  (for n >= X, l << n).
For l < sqrt(3 X ln2 / 2) this forces ||l log2 3|| < 1/(2l), so l is a convergent denominator; beyond that
also intermediate fractions can qualify.  This script enumerates convergents and intermediate fractions of
log2 3 (100-digit arithmetic), keeps the UPPER approximants (p/q > log2 3, i.e. 2^p > 3^q), and tests (R1)
with n = X for X in {2^28, 2.8e19}.  For survivors it reports delta_l = 2^A/3^l - 1 and the upper end
N_l of the n-window from (R2): 2 + delta_l n <= n((1 + 1/(3n))^l - 1).

Reproduce: python3 collatz_last_dip_residual_cf_20261005.py
"""
import math
from decimal import Decimal, getcontext
getcontext().prec = 110

THETA = Decimal(3).ln() / Decimal(2).ln()
LN2 = Decimal(2).ln()

def convergents(x, nterms):
    cf = []; y = x
    h0, h1 = 1, 0; k0, k1 = 0, 1   # h = numerators, k = denominators (standard recursion)
    conv = []
    for _ in range(nterms):
        a = int(y); cf.append(a)
        h0, h1 = a * h0 + h1, h0
        k0, k1 = a * k0 + k1, k0
        conv.append((h0, k0, a))
        frac = y - a
        if frac == 0: break
        y = 1 / frac
    return cf, conv

cf, conv = convergents(THETA, 60)
print("log2 3 =", str(THETA)[:40], "...")
print("continued fraction:", cf[:40])

def upper_excess(p, q):
    # p - q log2 3 (> 0 for upper approximants)
    return Decimal(p) - Decimal(q) * THETA

def r1_bound(l, n):
    # log2((1 + 2/(3n-3)) e^{(l-1)/(3n)})
    n = Decimal(n); l = Decimal(l)
    return ((1 + 2 / (3 * n - 3)).ln() + (l - 1) / (3 * n)) / LN2

def N_l(l, delta):
    def ok(n):
        n = Decimal(n)
        e = Decimal(l) * (1 + 1 / (3 * n)).ln()
        if e > 200:
            return True                      # right side astronomically large
        return 2 + delta * n <= n * (e.exp() - 1)
    lo = None; n = 3
    while n < Decimal(10) ** 60:
        if ok(n): lo = n; break
        n *= 2
    if lo is None: return None
    hi = lo
    while ok(hi * 2): hi *= 2
    lo, hi = hi, hi * 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if ok(mid): lo = mid
        else: hi = mid
    return lo

# candidate approximants: convergents and intermediate fractions (p_{k-1} + j p_k)/(q_{k-1} + j q_k)
cands = set()
for i in range(1, len(conv) - 1):
    p0, q0, _ = conv[i - 1]; p1, q1, a_next = conv[i]
    cands.add((p1, q1))
    _, _, a2 = conv[i + 1]
    for j in range(1, a2 + 1):
        cands.add((p0 + j * p1, q0 + j * q1))
cands = sorted(c for c in cands if c[1] >= 6 and c[1] <= 10 ** 14)
uppers = [(p, q, upper_excess(p, q)) for p, q in cands]
uppers = [(p, q, e) for p, q, e in uppers if e > 0]
print(f"\nupper approximants p/q of log2 3 with 6 <= q <= 1e14: {len(uppers)} (convergents and intermediate fractions)")

for X in [2 ** 28, Decimal("2.8e19")]:
    print("\n" + "=" * 78)
    print(f"X = {X}: lengths l (= q) surviving (R1) at n = X, with window [X, N_l]")
    print("=" * 78)
    survivors = []
    for p, q, e in uppers:
        l = q; A = p
        if e < r1_bound(l, X):
            delta = Decimal(2) ** e - 1
            Nl = N_l(l, delta)
            if Nl is not None and Nl >= X:
                survivors.append((l, A, delta, Nl))
    print(f"  survivors: {len(survivors)}")
    for l, A, delta, Nl in survivors[:25]:
        print(f"   l = {l:>14d}  A = {A:>14d}  delta_l = {float(delta):.3e}  n in [{X}, {Nl}]  log2 N_l = {math.log2(int(Nl)):.1f}")
    if survivors:
        lmin = min(s[0] for s in survivors)
        print(f"  smallest residual length: l = {lmin}  (an excursion of {lmin} odd steps above n, returning to n exactly)")
print("\nCAVEAT (2026-10-06): the Beatty gate 1 - {l log2 3} < (l+2)/(3 X ln 2) is linear in l, so sums and small multiples of")
print("qualifying lengths qualify as well: the residual set is the Bohr set B(X) = {l : 1 - {l log2 3} < (l+2)/(3 X ln 2)},")
print("of which this list shows only the convergents and intermediate fractions.  The SMALLEST residual length is still the")
print("first upper convergent denominator past sqrt(3 X ln 2 / 2), because below that bound the gate forces a convergent.")
print("The brute-force script (collatz_last_dip_residual_20261005.py) enumerates B(X) exactly for l <= 1e5.")
