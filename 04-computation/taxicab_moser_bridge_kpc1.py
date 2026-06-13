#!/usr/bin/env python3
"""
taxicab_moser_bridge_kpc1.py
THREAD A (kind-pasteur-2026-06-10-S1): HYP-2367 / THM-463
The taxicab--Moser Eisenstein bridge: two-cube representations are a divisor
property on the split axis.

EXACT integer arithmetic throughout the verification logic. numpy int64 is
used ONLY for the part-(c) meet-in-the-middle enumeration (max value
2*10^12 << 2^63, no overflow possible); every collision found there is
re-derived and re-checked in pure Python ints before being reported.

(a) BIJECTION reps <-> good divisors, all n in [1, 10^6]:
    unordered {x,y} subset Z with x^3+y^3 = n   <-->   divisors d > 0 of n
    with 3d | d^3 - n and Delta = (4n - d^3)/(3d) a perfect square.
    Completeness of the |x|,|y| <= 1000 box for n <= 10^6 is PROVED:
    writing a rep with max-coordinate x >= 1 and d = x+y >= 1,
    n = d*(3x^2 - 3dx + d^2) has d-derivative 3(x-d)^2 >= 0, so
    n >= [value at d=1] = 3x^2 - 3x + 1, hence x <= 577.
    Also verifies the classical parity condition sqrt(Delta) == d (mod 2)
    is AUTOMATIC (e^2 = d^2 - 4s == d^2 mod 4 forces matching parity).

(b) SPLIT LEMMA on all primitive pairs 1 <= y <= x <= 2000:
    q = x^2 - xy + y^2 (= Eisenstein norm N(x + y*omega)) is odd;
    every prime p != 3 dividing q satisfies p == 1 (mod 3);
    v_3(q) <= 1, with equality iff 3 | x + y.

(c) PRIMITIVE positive 2-representation numbers (both reps gcd(x,y) = 1)
    up to 10^12, meet-in-the-middle on x <= 10^4. First 25 listed with
    factorizations. Empirical inert-prime test: theory (THM-463 Cor C)
    says an inert prime p == 2 (mod 3) dividing n must divide
    gcd(d_1, d_2) of the primitive divisors d_i = x_i + y_i; we check
    every doubly-primitive n in range for gcd(d_i, d_j) > 1.

(d) The 1729 slice: full good-divisor scan of 1729 (expect good = {13,19}),
    cofactors 133 = 7*19 and 91 = 7*13 (split products), B(1729) = 8,
    r_E(1729) = 48, THM-434 rosette = 12 + 48 = 60 units, record-B check
    (1729 = smallest t with B = 8; no t <= 1729 has B = 7), taxicab
    minimality (no n < 1729 with two positive reps), Ramanujan near-miss
    9^3 + 10^3 = 12^3 + 1.
"""

import sys, time
from math import isqrt, gcd
from array import array

T0 = time.time()


def icbrt_floor(m):
    """floor(cbrt(m)) for any integer m, exact."""
    if m >= 0:
        r = int(round(m ** (1.0 / 3.0)))  # float seed only; exact adjust below
        while r * r * r > m:
            r -= 1
        while (r + 1) ** 3 <= m:
            r += 1
        return r
    # m < 0: floor(cbrt(m)) = -ceil(cbrt(-m))
    return -icbrt_ceil(-m)


def icbrt_ceil(m):
    if m >= 0:
        f = icbrt_floor(m)
        return f if f * f * f == m else f + 1
    return -icbrt_floor(-m)


def banner(s):
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)


# ============================================================================
banner("(a) BIJECTION  reps <-> good divisors,  n in [1, 10^6]")
# ============================================================================
ta = time.time()
N_A = 10 ** 6
XB = 1000          # search box (PROVEN sufficient bound is 577; margin to 1000)
PROVEN_BOUND = 577

reps = bytearray(N_A + 1)
max_coord_seen = 0
n_pairs_counted = 0

# Any rep of n >= 1 has max(x,y) >= 1 (else n <= 0).  WLOG enumerate x >= y,
# x from 1..XB, y in the exact range making 1 <= n <= N_A and y >= -XB.
for x in range(1, XB + 1):
    x3 = x * x * x
    ylo = max(-XB, icbrt_ceil(1 - x3))        # n >= 1
    yhi = min(x, icbrt_floor(N_A - x3)) if x3 <= N_A else min(x, icbrt_floor(N_A - x3))
    if x3 > N_A + XB ** 3:
        break
    for y in range(ylo, yhi + 1):
        n = x3 + y * y * y
        if 1 <= n <= N_A:
            reps[n] += 1
            n_pairs_counted += 1
            mc = max(abs(x), abs(y))
            if mc > max_coord_seen:
                max_coord_seen = mc

print(f"reps side: {n_pairs_counted} unordered integer pairs (x>=y) with x^3+y^3 in [1,10^6]")
print(f"max |coordinate| over all reps found: {max_coord_seen} "
      f"(proven bound {PROVEN_BOUND}; box {XB} -> box completeness "
      f"{'CONFIRMED' if max_coord_seen <= PROVEN_BOUND else 'VIOLATED'})")

# Independent good-divisor count.  d | n, write k = n/d.  The condition
# 3d | d^3 - n  <=>  k == d^2 (mod 3)  (since n = dk and (d^3-n)/d = d^2-k...
# precisely: 3d | d(d^2-k) <=> 3 | d^2 - k when gcd... CAREFUL:
# 3d | d^3 - n = d(d^2 - k)  <=>  3 | (d^2 - k)  -- d cancels exactly.)
# Then Delta = (4n - d^3)/(3d) = (4k - d^2)/3, integral by the same congruence,
# and Delta >= 0 <=> k >= ceil(d^2/4).  Only d^3 <= 4n can contribute.
gd = bytearray(N_A + 1)
parity_violations = 0
n_good_divisors = 0
d = 1
while d * d * d <= 4 * N_A:
    d2 = d * d
    kmin = (d2 + 3) // 4                  # Delta >= 0
    if kmin < 1:
        kmin = 1
    kmax = N_A // d
    r = d2 % 3
    k = kmin + ((r - kmin) % 3)           # smallest k >= kmin with k == d^2 (mod 3)
    while k <= kmax:
        Delta = (4 * k - d2) // 3
        e = isqrt(Delta)
        if e * e == Delta:
            gd[d * k] += 1
            n_good_divisors += 1
            if (e - d) % 2 != 0:
                parity_violations += 1
        k += 3
    d += 1

print(f"divisor side: {n_good_divisors} good divisors over all n in [1,10^6]")

mismatch = 0
first_mismatches = []
for n in range(1, N_A + 1):
    if reps[n] != gd[n]:
        mismatch += 1
        if len(first_mismatches) < 10:
            first_mismatches.append((n, reps[n], gd[n]))
n_with_rep = sum(1 for n in range(1, N_A + 1) if reps[n] > 0)
mx = max(reps[1:])
argmx = [n for n in range(1, N_A + 1) if reps[n] == mx][:5]
print(f"#n in [1,10^6] that are sums of two integer cubes: {n_with_rep}")
print(f"max #reps = {mx}, first attaining n: {argmx}")
print(f"parity condition sqrt(Delta) == d (mod 2): {parity_violations} violations "
      f"(0 expected -- the condition is automatic)")
print(f"BIJECTION MISMATCHES: {mismatch}" + (f"  first: {first_mismatches}" if mismatch else ""))
print(f"(a) {'PASS' if mismatch == 0 and parity_violations == 0 and max_coord_seen <= PROVEN_BOUND else 'FAIL'}"
      f"   [{time.time()-ta:.1f}s]")

# ============================================================================
banner("(b) SPLIT LEMMA  on all primitive pairs 1 <= y <= x <= 2000")
# ============================================================================
tb = time.time()
LIM_B = 2000
QMAX = LIM_B * LIM_B          # q = x^2-xy+y^2 <= x^2 for 0 <= y <= x

# smallest-prime-factor sieve to QMAX (pure int)
spf = array('i', range(QMAX + 1))
i = 2
while i * i <= QMAX:
    if spf[i] == i:
        for j in range(i * i, QMAX + 1, i):
            if spf[j] == j:
                spf[j] = i
    i += 1

prim_pairs = 0
viol_odd = viol_split = viol_v3 = 0
for x in range(1, LIM_B + 1):
    xx = x * x
    for y in range(1, x + 1):
        if gcd(x, y) != 1:
            continue
        prim_pairs += 1
        q = xx - x * y + y * y
        if q % 2 == 0:
            viol_odd += 1
        v3 = 0
        qq = q
        while qq % 3 == 0:
            qq //= 3
            v3 += 1
        if not (v3 <= 1 and (v3 == 1) == ((x + y) % 3 == 0)):
            viol_v3 += 1
        ok = True
        while qq > 1:
            p = spf[qq]
            if p % 3 != 1:
                ok = False
                break
            while qq % p == 0:
                qq //= p
        if not ok:
            viol_split += 1

print(f"primitive pairs checked: {prim_pairs}")
print(f"violations -- q even: {viol_odd}; prime p != 3 of q with p !== 1 (mod 3): {viol_split}; "
      f"v3 law: {viol_v3}")
print(f"(b) {'PASS' if viol_odd == viol_split == viol_v3 == 0 else 'FAIL'}   [{time.time()-tb:.1f}s]")

# primes <= 10^6 for part (c) factorization, harvested from the spf sieve
PRIMES_1M = [p for p in range(2, 10 ** 6 + 1) if spf[p] == p]
print(f"(harvested {len(PRIMES_1M)} primes <= 10^6 for part (c) trial division)")


def factorize(n):
    """exact factorization for n <= 10^12 via primes <= 10^6."""
    fac = []
    m = n
    for p in PRIMES_1M:
        if p * p > m:
            break
        if m % p == 0:
            e = 0
            while m % p == 0:
                m //= p
                e += 1
            fac.append((p, e))
    if m > 1:
        fac.append((m, 1))   # remaining cofactor is prime (n <= 10^12)
    return fac


def fac_str(fac):
    def cls(p):
        if p == 3:
            return "ram"
        return "split" if p % 3 == 1 else "INERT"
    return " * ".join((f"{p}^{e}" if e > 1 else f"{p}") + f"[{cls(p)}]" for p, e in fac)


# ============================================================================
banner("(c) PRIMITIVE positive 2-rep numbers <= 10^12  (meet-in-the-middle)")
# ============================================================================
tc = time.time()
import numpy as np

N_C = 10 ** 12
X_C = 10 ** 4
cb = np.arange(0, X_C + 1, dtype=np.int64) ** 3      # exact: 10^12 < 2^63
assert int(cb[X_C]) == X_C ** 3                       # overflow guard

chunks = []
for x in range(1, X_C + 1):
    x3 = x * x * x
    rem = N_C - x3
    if rem < 1:
        break
    ymax = min(x, icbrt_floor(rem))
    if ymax >= 1:
        chunks.append(x3 + cb[1:ymax + 1])
S = np.concatenate(chunks)
del chunks
total_pairs = S.size
S.sort()
dupmask = S[1:] == S[:-1]
dup_vals = np.unique(S[1:][dupmask])
del S, dupmask
print(f"positive pairs y <= x with x^3+y^3 <= 10^12: {total_pairs}")
print(f"values with >= 2 positive representations: {dup_vals.size}")

cube_set = set(i * i * i for i in range(1, X_C + 1))


def positive_reps(n):
    """all positive reps x >= y >= 1 of n, exact pure-int."""
    out = []
    xlo = icbrt_ceil((n + 1) // 2)
    if 2 * xlo ** 3 < n:   # guard: ensure x^3 >= n - x^3 region fully covered
        xlo = max(1, xlo - 1)
    xhi = icbrt_floor(n - 1)
    for x in range(max(1, xlo - 1), xhi + 1):
        r = n - x * x * x
        if r < 1 or r > x * x * x:
            continue
        if r in cube_set:
            y = icbrt_floor(r)
            if y * y * y == r and 1 <= y <= x:
                out.append((x, y))
    return out


two_rep_all = 0
doubly_primitive = []     # (n, [(x,y) primitive reps])
shared_d = []             # (n, gcd-of-primitive-d's > 1)
for v in dup_vals:
    n = int(v)
    rr = positive_reps(n)
    assert len(rr) >= 2, (n, rr)   # re-derivation must confirm the collision
    two_rep_all += 1
    prr = [(x, y) for (x, y) in rr if gcd(x, y) == 1]
    if len(prr) >= 2:
        doubly_primitive.append((n, prr))
        ds = [x + y for (x, y) in prr]
        g = 0
        for dv in ds:
            g = gcd(g, dv)
        if g > 1:
            shared_d.append((n, prr, ds, g))

doubly_primitive.sort()
print(f"values with >= 2 positive reps (re-derived exactly): {two_rep_all}")
print(f"DOUBLY-PRIMITIVE (>= 2 reps with gcd(x,y) = 1): {len(doubly_primitive)}")
print()
print("first 25 doubly-primitive taxicab numbers with factorizations:")
print(f"{'n':>15}  reps (x,y) | d = x+y | cofactor m = n/d (all m split-axis by THM-463.B)")
inert_in_first25 = 0
for n, prr in doubly_primitive[:25]:
    fac = factorize(n)
    has_inert = any(p % 3 == 2 for p, e in fac)
    if has_inert:
        inert_in_first25 += 1
    rep_desc = "; ".join(f"{x}^3+{y}^3, d={x+y}, m={n // (x + y)}" for (x, y) in prr)
    print(f"{n:>15}  = {fac_str(fac)}")
    print(f"{'':>15}  {rep_desc}{'   <<< CONTAINS INERT PRIME' if has_inert else ''}")
    # exact split-axis check of each cofactor (verifies THM-463.B at 10^12 scale)
    for (x, y) in prr:
        m = n // (x + y)
        mf = factorize(m)
        v3m = sum(e for p, e in mf if p == 3)
        bad = [p for p, e in mf if p % 3 == 2] + (["v3>1"] if v3m > 1 else [])
        if bad:
            print(f"{'':>15}  *** SPLIT-LEMMA VIOLATION in cofactor {m}: {bad}")

print()
print(f"inert-prime test over ALL {len(doubly_primitive)} doubly-primitive n <= 10^12:")
print(f"  n with gcd of primitive d's > 1: {len(shared_d)}")
if shared_d:
    for n, prr, ds, g in shared_d[:20]:
        gf = factorize(g)
        inert = [p for p, e in gf if p % 3 == 2]
        nf = factorize(n)
        print(f"  n = {n} = {fac_str(nf)}  d's = {ds}  gcd = {g} = {fac_str(gf)}"
              f"  inert in gcd: {inert if inert else 'NONE'}")
    n_inert_total = sum(1 for n, prr, ds, g in shared_d
                        if any(p % 3 == 2 for p, e in factorize(g)))
else:
    n_inert_total = 0
# Theory check: an inert prime can divide a doubly-primitive n ONLY via gcd(d_i).
# So if no gcd > 1 contains an inert prime, NO doubly-primitive n in range has one.
print(f"  => doubly-primitive n containing an inert prime (p == 2 mod 3): "
      f"{n_inert_total} of {len(doubly_primitive)}")
print(f"  (theory THM-463.C: inert p | n forces p^v_p(n) | every primitive d, "
      f"hence p | gcd of the d's)")
print(f"(c) done   [{time.time()-tc:.1f}s]")

# ============================================================================
banner("(d) THE 1729 SLICE")
# ============================================================================
td = time.time()
n = 1729
divs = [d for d in range(1, n + 1) if n % d == 0]
print(f"1729 = {fac_str(factorize(1729))}   divisors: {divs}")
print()
print(f"{'d':>5} {'3d | d^3-n?':>12} {'s=xy':>8} {'Delta':>9} {'square?':>8} {'rep {x,y}':>14}")
good = []
for d in divs:
    if (d ** 3 - n) % (3 * d) != 0:
        print(f"{d:>5} {'no':>12} {'-':>8} {'-':>9} {'-':>8} {'-':>14}")
        continue
    s = (d ** 3 - n) // (3 * d)
    Delta = d * d - 4 * s
    if Delta < 0:
        print(f"{d:>5} {'yes':>12} {s:>8} {Delta:>9} {'(<0) no':>8} {'-':>14}")
        continue
    e = isqrt(Delta)
    if e * e == Delta:
        x, y = (d + e) // 2, (d - e) // 2
        assert x ** 3 + y ** 3 == n
        good.append(d)
        print(f"{d:>5} {'yes':>12} {s:>8} {Delta:>9} {'YES':>8} {str({x, y}):>14}")
    else:
        print(f"{d:>5} {'yes':>12} {s:>8} {Delta:>9} {'no':>8} {'-':>14}")
print(f"\ngood divisors of 1729: {good}  (expected [13, 19])")
ok_d1 = (good == [13, 19])

m13, m19 = 1729 // 13, 1729 // 19
print(f"cofactors: 1729/13 = {m13} = {fac_str(factorize(m13))};  "
      f"1729/19 = {m19} = {fac_str(factorize(m19))}  (both split products)")
print(f"gcd(d1, d2) = gcd(13, 19) = {gcd(13, 19)}  -> no inert prime can divide 1729 "
      f"(and indeed 1729 = 7*13*19 is completely split)")

# r_E(1729): representations of 1729 by THM-434's norm form x^2 + xy + y^2
rE = 0
R = isqrt(4 * 1729 // 3) + 2
for a in range(-R, R + 1):
    for b in range(-R, R + 1):
        if a * a + a * b + b * b == 1729:
            rE += 1
# B(1729) via the divisor character chi_{-3}
chi = lambda d: 0 if d % 3 == 0 else (1 if d % 3 == 1 else -1)
B1729 = sum(chi(d) for d in divs)
units = 12 + rE
print(f"r_E(1729) = {rE} (expected 48);  B(1729) = sum chi_-3(d) = {B1729} (expected 8);  "
      f"r_E = 6B: {rE == 6 * B1729}")
print(f"THM-434 rosette: #units(L_1729) = 12 + r_E = {units} (expected 60)")
ok_d2 = (rE == 48 and B1729 == 8 and units == 60)

# record-B check: 1729 is the smallest t with B(t) = 8; no t <= 1729 has B = 7
TT = 1729
Bt = [0] * (TT + 1)
for d in range(1, TT + 1):
    c = chi(d)
    if c:
        for t in range(d, TT + 1, d):
            Bt[t] += c
rec = []
best = -1
for t in range(1, TT + 1):
    if Bt[t] > best:
        best = Bt[t]
        rec.append((t, Bt[t]))
print(f"B-records up to 1729: {rec}")
has_B7 = any(Bt[t] == 7 for t in range(1, TT + 1))
ok_d3 = (rec[-1] == (1729, 8) and not has_B7)
print(f"1729 is the FIRST t with B = 8 (60-unit rosette): {rec[-1] == (1729, 8)};  "
      f"any t <= 1729 with B = 7: {has_B7} (expected False; B = 7 needs p^6, min 7^6 = {7**6})")

# taxicab minimality: no n < 1729 with two positive reps
pos2 = {}
for x in range(1, 13):
    for y in range(1, x + 1):
        v = x ** 3 + y ** 3
        if v < 1729:
            pos2[v] = pos2.get(v, 0) + 1
min2 = [v for v, c in sorted(pos2.items()) if c >= 2]
print(f"n < 1729 with two positive reps: {min2} (expected []) -> Ta(2) = 1729")
ok_d4 = (min2 == [])

print(f"reps of 1729: 1^3 + 12^3 = {1 + 12**3};  9^3 + 10^3 = {9**3 + 10**3};  "
      f"both primitive: {gcd(1,12) == 1 and gcd(9,10) == 1}")
print(f"Ramanujan near-miss: 9^3 + 10^3 - 12^3 = {9**3 + 10**3 - 12**3} "
      f"(= 1: first taxicab number is one above Klein's 1728 = j(i))")
ok_d = ok_d1 and ok_d2 and ok_d3 and ok_d4
print(f"(d) {'PASS' if ok_d else 'FAIL'}   [{time.time()-td:.1f}s]")

banner(f"TOTAL [{time.time()-T0:.1f}s]")
