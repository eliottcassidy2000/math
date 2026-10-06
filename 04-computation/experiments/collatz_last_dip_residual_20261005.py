#!/usr/bin/env python3
"""
collatz_last_dip_residual_20261005.py -- the explicit residual of the last-dip lemma after a sweep to X.

PROVED (note section 4): a violating excursion (word w of length l, total valuation A, source m*, target n,
gap g = n - m* >= 2) satisfies
   (R1)  2^A/3^l = (m*/n) prod_{j=1}^{l} (1 + 1/(3 n_{j-1}))  with n_0 = m*, n_{j-1} > n for j >= 2,
         hence 1 < 2^A/3^l < (1 + 1/(3m*)) (1 + 1/(3n))^(l-1) <= (1 + 2/(3n-3)) e^{(l-1)/(3n)}   [m* > (2n-1)/3],
         so A = ceil(l log2 3) and  delta_l := 2^A/3^l - 1 < (1 + 2/(3n-3)) e^{(l-1)/(3n)} - 1 =: eps(l, n);
   (R2)  g = c_w/3^l - delta_l n >= 2 and c_w/3^l <= n((1 + 1/(3n))^l - 1), hence n <= N_l := the largest n with
         2 + delta_l n <= n((1+1/(3n))^l - 1).
After an exhaustive sweep below X (no violation with n < X), a violation needs n >= X, so by (R1)
   1 - {l log2 3} < log2(1 + eps(l, X)),
which (for l << X) is about (l+2)/(3 X ln 2): only l whose multiples l log2 3 land within that distance
below an integer -- the denominators of upper (semi)convergents of log2 3.  This script lists, for each X in
{2^20, 2^28, 2^32, 2^36}, the lengths l <= LMAX surviving (R1) and their windows [X, N_l] from (R2).

Reproduce: python3 collatz_last_dip_residual_20261005.py [LMAX]
"""
import sys, math
from fractions import Fraction
from decimal import Decimal, getcontext
getcontext().prec = 60

LMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 2_000_000
L23 = Decimal(3).ln() / Decimal(2).ln()      # log2 3 to 60 digits

def frac_part(l):
    x = l * L23
    return x - int(x)

def delta_of(l):
    # 2^ceil(l log2 3)/3^l - 1 = 2^(1 - {l log2 3}) - 1
    f = frac_part(l)
    return (Decimal(2) ** (1 - f)) - 1

def eps(l, n):
    n = Decimal(n)
    return (1 + 2 / (3 * n - 3)) * ((l - 1) / (3 * n)).exp() - 1

def N_l(l, delta):
    # largest n with 2 + delta*n <= n((1+1/(3n))^l - 1); the right side -> l/3 as n -> inf; solve by bisection
    def ok(n):
        n = Decimal(n)
        e = Decimal(l) * (1 + 1 / (3 * n)).ln()
        if e > 200:
            return True
        return 2 + delta * n <= n * (e.exp() - 1)
    lo, hi = 3, 10 ** 30
    if not ok(lo):
        # find any ok n by scanning powers of two
        found = None
        n = 4
        while n < hi:
            if ok(n):
                found = n; break
            n *= 2
        if found is None:
            return None
        lo = found
    # bisection for the largest ok n (ok(n) is true on an interval; we want its right end)
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if ok(mid):
            lo = mid
        else:
            hi = mid
    return lo

print("Continued fraction of log2 3 (convergents p/q; upper = p/q > log2 3):")
# compute convergents from the decimal expansion
x = L23
cf = []
a = int(x); cf.append(a); y = x - a
for _ in range(16):
    y = 1 / y; a = int(y); cf.append(a); y = y - a
h0, h1 = 1, cf[0]; k0, k1 = 0, 1
conv = [(h1, k1)]
for a in cf[1:]:
    h0, h1 = h1, a * h1 + h0; k0, k1 = k1, a * k1 + k0
    conv.append((h1, k1))
for p, q in conv[:12]:
    side = "upper" if Decimal(p) / Decimal(q) > L23 else "lower"
    print(f"   {p}/{q}  ({side}; 1-{{q log2 3}} = {float(1 - frac_part(q)) if side == 'upper' else float(frac_part(q)):.3e})")

print()
for X in [2 ** 20, 2 ** 28, 2 ** 32, 2 ** 36]:
    survivors = []
    for l in range(6, LMAX + 1):
        f = frac_part(l)
        gap = 1 - f                      # distance below the next integer
        # quick filter: need gap < log2(1 + eps(l, X)) ~ (l+2)/(3 X ln2)
        quick = Decimal(l + 2) / (Decimal(3) * Decimal(X) * Decimal(2).ln()) * Decimal("1.0001") + Decimal("1e-40")
        if gap >= quick:
            continue
        e = eps(l, X)
        if gap < (1 + e).ln() / Decimal(2).ln():
            d = delta_of(l)
            Nl = N_l(l, d)
            if Nl is not None and Nl >= X:
                survivors.append((l, d, Nl))
    print(f"X = 2^{int(math.log2(X))}: lengths l <= {LMAX} surviving (R1) with window [X, N_l] nonempty: {len(survivors)}")
    for l, d, Nl in survivors[:20]:
        print(f"   l = {l:8d}  A = {math.ceil(l * float(L23)):8d}  delta_l = {float(d):.3e}  n in [{X}, {Nl}]  (log2 N_l = {math.log2(Nl):.2f})")
    if len(survivors) > 20:
        print(f"   ... {len(survivors) - 20} more")
print("\nThe lemma's residual after a sweep to X is the union of these (l, [X, N_l]) windows: an excursion of l odd")
print("steps above n, returning to n exactly, with A = ceil(l log2 3) halvings and a non-rising word.")
