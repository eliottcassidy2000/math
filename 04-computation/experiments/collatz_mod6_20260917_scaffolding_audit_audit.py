#!/usr/bin/env python3
"""Adversarial audit of lane scaffolding_audit (session collatz-mod6-20260917, mac-mini; audit written 2026-09-21).

Independent recomputation of every number quoted in
05-knowledge/results/collatz_mod6_20260917_scaffolding_audit.md, written from scratch
(different code paths from the lane script wherever possible: sympy symbolic identities instead of
hand-rolled polynomial dicts; cyclotomic-value criterion for primitive primes instead of
factorisation+order; permutation-matching census of tournaments instead of a bitmask DP;
independent peak census with the box/ladder analysis that the lane note lacked).

Every check is an explicit `raise` (survives python3 -O).  Runtime ~40 s, RAM << 1 GB.

Reproduce:
  python3 04-computation/experiments/collatz_mod6_20260917_scaffolding_audit_audit.py \
      > 05-knowledge/results/collatz_mod6_20260917_scaffolding_audit_audit.out
"""
import itertools
import math
from collections import Counter, defaultdict
from fractions import Fraction

import sympy
from sympy import Rational, binomial, expand, factorial, symbols

FAILS = []


def check(cond, msg):
    if not cond:
        FAILS.append(msg)
        raise AssertionError("AUDIT CHECK FAILED: " + msg)


def hdr(t):
    print()
    print("=" * 78)
    print(t)
    print("=" * 78)


def T(n):
    return n * (n + 1) // 2


A, B = symbols("A B")
Ts = lambda x: x * (x + 1) / 2  # symbolic triangular


print("scaffolding_audit_audit: independent recomputation for lane scaffolding_audit")
print("session collatz-mod6-20260917 (mac-mini), audit written 2026-09-21")

# ----------------------------------------------------------------------------- A1
hdr("A1. Claim 1 (core identity, repair, flows): symbolic re-derivation with sympy")
N = (A + 1) * (B + 1) - 1
lhs = expand(Ts(N) - Ts(A * B - 1))
shear = A * (B**2 - 1) + B * (A**2 - 1)
rhs = expand(Ts(A) + Ts(B) + shear)
d = expand(lhs - rhs)
print("  LHS - RHS(paste) =", d)
check(d == expand(2 * A * B + A + B), "difference must be 2AB+A+B")
check(expand(lhs - (Ts(A) + Ts(B) + A * B * (A + B + 2))) == 0, "true shear AB(A+B+2)")
check(expand(A * B * (A + B + 2) - (2 * A * Ts(B) + 2 * B * Ts(A))) == 0, "AB(A+B+2)=2A T(B)+2B T(A)")
check(expand(lhs - (A + B + 1) * (2 * A * B + A + B) / 2) == 0, "closed form (A+B+1)(2AB+A+B)/2")
# repair
check(expand(Ts(N - 1) - Ts(A * B) - rhs) == 0, "repaired identity T((A+1)(B+1)-2)-T(AB) = T(A)+T(B)+Shear")
# ALL rational solutions of the shift system (not just shifts in [-3,3])
s1, s2, t1, t2 = symbols("s1 s2 t1 t2")
gen = expand(Ts(N + s1) - Ts(A * B + s2) - Ts(A + t1) - Ts(B + t2) - shear)
eqs = sympy.Poly(gen, A, B).coeffs()
sols = sympy.solve(eqs, [s1, s2, t1, t2], dict=True)
sols_t = sorted(tuple(int(s[v]) for v in (s1, s2, t1, t2)) for s in sols)
print("  all rational solutions (s1,s2,t1,t2) of the shifted form, as polynomial identity:", sols_t)
check(sols_t == [(-3, -2, -2, -2), (-1, 0, 0, 0)], "exactly the two shift solutions, over Q not just [-3,3]")
# Lemma 2.2: no quadratic T(n)=a n^2+b n+c makes the STATED form an identity
a_, b_, c_ = symbols("a b c")
Tq = lambda x: a_ * x**2 + b_ * x + c_
gq = expand(Tq(N) - Tq(A * B - 1) - Tq(A) - Tq(B) - shear)
Pq = sympy.Poly(gq, A, B)
cA2B = Pq.coeff_monomial(A**2 * B)
cAB = Pq.coeff_monomial(A * B)
print("  quadratic reading: coeff of A^2B =", cA2B, "; coeff of AB =", cAB)
check(cA2B == 2 * a_ - 1 and cAB == 4 * a_, "a=1/2 forced by A^2B, then AB coefficient 2 != 0")
check(sympy.solve(Pq.coeffs(), [a_, b_, c_]) == [], "no quadratic T rescues the stated form")
# flows
g_h = A + B**2 + 2 * A * B + B
g_v = B + A**2 + 2 * A * B + A
S = A * B * (A + B) + binomial(A, 2) + binomial(B, 2)
S = expand(sympy.expand_func(S))
check(expand(S.subs(A, A + 1) - S - g_h) == 0 and expand(S.subs(B, B + 1) - S - g_v) == 0, "S solves both flows")
check(expand(g_h.subs(B, B + 1) - g_h - (g_v.subs(A, A + 1) - g_v)) == 0, "compatibility")
check(expand(S - rhs) == 0, "S = paste's right side")
check(expand(lhs.subs(A, A + 1) - lhs - g_h - (2 * B + 1)) == 0, "paste's LHS has horizontal flow g_h + 2B+1")
print("  CONFIRMED: stated identity false (A=B=1: %d vs %d); repair; flows' unique solution = paste RHS;"
      % (T(3) - T(0), T(1) + T(1)))
print("             LHS horizontal flow = g_h + (2B+1).")

# ----------------------------------------------------------------------------- A2
hdr("A2. Claim 2 (tetrahedral / pentatope / d-simplex laws): symbolic")
Te = lambda x: binomial(x + 2, 3)
Pt = lambda x: binomial(x + 3, 4)
ef = lambda e: expand(sympy.expand_func(e))
tet = ef(Te(A + B + 1) - Te(A) - Te(B) - (A + 1) * Ts(B) - (B + 1) * Ts(A) - (A + 1) * (B + 1))
pen = ef(Pt(A + B + 1) - Pt(A) - Pt(B) - (A + 1) * Te(B) - (B + 1) * Te(A) - Ts(A + 1) * Ts(B + 1))
two = ef(Ts(A + B + 1) - Ts(A) - Ts(B) - (A + 1) * (B + 1))
check(tet == 0 and pen == 0 and two == 0, "tetrahedral, pentatope, 2D laws are polynomial identities")
print("  tetrahedral, pentatope and 2D laws: polynomial identities in Q[A,B] (residuals 0, 0, 0)")
for dd in range(1, 9):
    for aa in range(0, dd + 1):
        bb = dd - aa
        e = ef(binomial(A + B + dd, dd) - sum(binomial(A + aa, i) * binomial(B + bb, dd - i) for i in range(dd + 1)))
        check(e == 0, "Vandermonde d=%d a=%d" % (dd, aa))
    e = ef(binomial(A + B + dd, dd) - sum(binomial(A + i, i) * binomial(B + dd - i - 1, dd - i) for i in range(dd + 1)))
    check(e == 0, "simplex convolution d=%d" % dd)
print("  Chu-Vandermonde, all splits, d<=8: polynomial identities; simplex convolution d<=8: identities.")
print("  boundary: binomial(B-1,0)=1 at B=0 is the S_0(0)=1 convention (sympy: %s)" % binomial(-1, 0))
check(binomial(-1, 0) == 1 and binomial(-1, 1) == -1, "note the lane's binom() returns 0 for C(-1,1); sympy gives -1")
# does the lane's convention matter?  S_j(0)=C(j-1,j) for j>=1: sympy value vs lane value
vals = [(j, int(binomial(j - 1, j))) for j in range(1, 6)]
print("  C(j-1,j), j=1..5 (sympy):", vals, " -- all 0, so the lane's binom() agrees on every used value")
check(all(v == 0 for _, v in vals), "C(j-1,j)=0 for j>=1")

# ----------------------------------------------------------------------------- A3
hdr("A3. Claim 3 (g-operator)")
check(all(B_**B_ != math.factorial(B_) for B_ in range(2, 10)) and 1**1 == math.factorial(1), "B^B=B! only at B=1")
# non-uniqueness without symmetry: f(A,B)=A!B!*c^(A-B-1) for A>B (free boundary f(k,0)=k! c^(k-1), k>=2)
# also solves f(A,B)=AB f(A-1,B-1) for A,B>=1 and f(1,B)=B! for all B>=0 (f(1,0)=1 is forced).
def f_alt(A_, B_, c=5):
    return math.factorial(A_) * math.factorial(B_) * (c ** (A_ - B_ - 1) if A_ > B_ else 1)
ok = all(f_alt(A_, B_) == A_ * B_ * f_alt(A_ - 1, B_ - 1) for A_ in range(1, 8) for B_ in range(1, 8))
ok = ok and all(f_alt(1, B_) == math.factorial(B_) for B_ in range(0, 8))
check(ok, "asymmetric alternative solution exists")
print("  CONFIRMED: A!B! is the unique SYMMETRIC solution; without symmetry f(A,B)=A!B!*5^(A-B-1) (A>B) also")
print("             satisfies f=AB f(A-1,B-1) and f(1,B)=B!, so 'symmetric' in Lemma 4.1 is load-bearing.")
check(all(f_alt(A_, A_) == math.factorial(A_) ** 2 for A_ in range(0, 8)), "AgA=(A!)^2 for the alternative too")
print("             (AgA=(A!)^2 holds for every solution, since f(A,A)=(A!)^2 f(0,0) and f(0,0)=1 from f(1,1)=1.)")

# ----------------------------------------------------------------------------- A4
hdr("A4. Claim 4 (tournament census by permutation matching, independent of the DP)")
try:
    import numpy as np
except ImportError:  # pragma: no cover
    np = None
for n in range(2, 7):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    idx = {p: k for k, p in enumerate(pairs)}
    m = len(pairs)
    # bit k = 1 means edge i->j for pair (i,j) with i<j
    perms = list(itertools.permutations(range(n)))
    masks, vals = [], []
    for pi in perms:
        mk = vl = 0
        for u, v in zip(pi, pi[1:]):
            if u < v:
                k = idx[(u, v)]
                mk |= 1 << k
                vl |= 1 << k
            else:
                k = idx[(v, u)]
                mk |= 1 << k
        masks.append(mk)
        vals.append(vl)
    if np is not None:
        allT = np.arange(1 << m, dtype=np.int64)
        h = np.zeros(1 << m, dtype=np.int64)
        for mk, vl in zip(masks, vals):
            h += ((allT & mk) == vl)
        hs = h.tolist()
    else:
        hs = [sum(1 for mk, vl in zip(masks, vals) if (t & mk) == vl) for t in range(1 << m)]
    tot = sum(hs)
    seen = sorted(set(hs))
    print("  n=%d: #T=%d sum_h=%d (n!2^T(n-2)=%d) max=%d seen=%s" % (n, len(hs), tot, math.factorial(n) * 2 ** T(n - 2), max(hs), seen))
    check(tot == math.factorial(n) * 2 ** T(n - 2), "double counting n=%d" % n)
    check(all(x % 2 == 1 for x in hs), "Redei n=%d" % n)
    check(7 not in seen and 21 not in seen, "holes")
    if n == 6:
        check(seen == [1, 3, 5, 9, 11, 13, 15, 17, 19, 23, 25, 27, 29, 31, 33, 37, 41, 43, 45], "n=6 spectrum")
        missing = [x for x in range(1, 46, 2) if x not in seen]
        print("  n=6: odd values <=45 NOT seen: %s (7 and 21 are THM-1745's holes; 35, 39 are merely absent at n=6)" % missing)
        check(missing == [7, 21, 35, 39], "n=6 missing odds")
    if n == 5:
        check(seen == [1, 3, 5, 9, 11, 13, 15], "n=5 spectrum")
for n in (3, 4, 5, 6, 8, 10):
    ex = math.log2(math.factorial(n)) - (n - 1)
    print("  n=%2d excess bits %.3f" % (n, ex))
check(abs(math.log2(6) - 2 - 0.585) < 1e-3 and abs(math.log2(math.factorial(10)) - 9 - 12.791) < 1e-3, "excess values")

# ----------------------------------------------------------------------------- A5
hdr("A5. Bang/Zsigmondy via cyclotomic values (no factorisation of 2^n-1 or 4^n-1)")


def phi_val(n, a):
    """Phi_n(a) as an exact integer, via Phi_n(a) = prod_{d|n} (a^d-1)^{mu(n/d)}."""
    num = den = 1
    for dd in sympy.divisors(n):
        mu = sympy.mobius(n // dd)
        if mu == 1:
            num *= a**dd - 1
        elif mu == -1:
            den *= a**dd - 1
    check(num % den == 0, "cyclotomic value integrality")
    return num // den


def has_primitive_prime(a, n):
    """a^n-1 has a primitive prime divisor iff Phi_n(a) stripped of primes dividing n exceeds 1 (n>=2 standard;
    n=1: primitive prime iff a-1>1)."""
    if n == 1:
        return a - 1 > 1
    v = phi_val(n, a)
    for p in sympy.primefactors(n):
        while v % p == 0:
            v //= p
    return v > 1


exc2 = [n for n in range(1, 41) if not has_primitive_prime(2, n)]
exc4 = [n for n in range(1, 41) if not has_primitive_prime(4, n)]
print("  2^n-1 without primitive prime, n<=40:", exc2)
print("  4^n-1 without primitive prime, n<=40:", exc4)
check(exc2 == [1, 6] and exc4 == [], "Bang exceptions [1,6]; none for base 4")
check(phi_val(6, 2) == 3 and phi_val(5, 4) == 341, "Phi_6(2)=3, Phi_5(4)=341")
print("  Phi_6(2)=3 (divides 6: the Catalan exception); Phi_5(4)=341=11*31, both primes primitive for 4^5-1")
# multiplicative orders by hand
def mord(a, p):
    if a % p == 0:
        raise ValueError("order undefined when p | a")
    k, x = 1, a % p
    while x != 1:
        x = x * a % p
        k += 1
    return k
check(mord(4, 11) == 5 and mord(4, 31) == 5 and mord(2, 9) == 6, "orders")
check(all(mord(2, p) != 6 for p in sympy.primerange(3, 100)), "no odd prime <100 of order 6 to base 2 (63=3^2*7)")
# Cipolla
for p in [5, 7, 11, 13, 17, 19, 23]:
    Np = (4**p - 1) // 3
    check(pow(2, Np - 1, Np) == 1 and not sympy.isprime(Np), "Cipolla p=%d" % p)
    check((Np - 1) % (2 * p) == 0, "2p | N_p-1")
check(pow(2, 20, 21) == 4, "p=3 fails: 2^20 = 4 mod 21")
check(11 * 31 == 341 and 43 * 127 == 5461 and 23 * 89 * 683 == 1398101 and 2731 * 8191 == 22369621, "factorisations")
check(131071 * 43691 == 5726623061 and 174763 * 524287 == 91625968981 and 47 * 178481 * 2796203 == 23456248059221, "factorisations 2")
print("  Cipolla p=5..23 confirmed; 2p | N_p-1; p=3 fails (2^20 = 4 mod 21); all seven factorisations re-multiplied.")
print("  INHERITANCE: 05-knowledge/results/collatz_mod6_20260917_wild_typing.md section 3 already proves Cipolla,")
print("    the exact criterion 6j | 4^j-4, Zsigmondy-no-exception for 4^n-1, and REFUTES '341 stalls primitive primes'.")
R_orb = [0]
for _ in range(6):
    R_orb.append(4 * R_orb[-1] + 1)
check(R_orb == [0, 1, 5, 21, 85, 341, 1365], "R-orbit")

# ----------------------------------------------------------------------------- A6
hdr("A6. Rational 3-cycle")
f = lambda x: x * x - Fraction(29, 16)
x = Fraction(-7, 4)
o = [x, f(x), f(f(x)), f(f(f(x)))]
check(o == [Fraction(-7, 4), Fraction(5, 4), Fraction(-1, 4), Fraction(-7, 4)], "3-cycle")
print("  orbit:", [str(v) for v in o], " -- confirmed; AP with difference 3/2")

# ----------------------------------------------------------------------------- A7
hdr("A7. Peaks census (independent implementation) and the ladder 'explanation'")


def C(n):
    return n // 2 if n % 2 == 0 else (3 * n + 1) // 2


def peak_shortcut(n):
    p = n
    while n != 1:
        n = C(n)
        if n > p:
            p = n
    return p


def peak_full(n):
    p = n
    while n != 1:
        n = n // 2 if n % 2 == 0 else 3 * n + 1
        if n > p:
            p = n
    return p


LIM = 10**5
pk = {n: peak_shortcut(n) for n in range(3, LIM + 1, 2)}
pkf = {n: peak_full(n) for n in range(3, LIM + 1, 2)}
check(len(pk) == 49999, "49,999 odd seeds")
check(all(p % 6 == 2 for p in pk.values()) and all(p % 6 == 4 for p in pkf.values()), "residues of peaks")
check(pk[5] == 8 and pk[7] == 26 and pk[23] == 80 and pkf[5] == 16 and pkf[7] == 52 and pkf[23] == 160, "peaks 5,7,23")
# boundary seeds
check(peak_shortcut(1) == 1 and peak_full(1) == 1, "seed 1 peak 1 (odd; excluded from Lemma 8.1 by n>1)")
check(peak_shortcut(3) == 8 and peak_shortcut(9) == 26, "n=3 -> 8, n=9 -> 26")
# proof of Lemma 8.1 needs n>1 odd: for the seed itself to be the peak it must not be odd>1
check(all(pk[n] > n for n in pk), "odd seed >1 is never its own shortcut peak")
best = max(pk, key=lambda n: Fraction(pk[n], n))
print("  max peak/seed: %.3f at seed %d (peak %d)" % (pk[best] / best, best, pk[best]))
check(best == 77671 and pk[best] == 785412368, "max ratio witness")
check(2 * 785412368 == pkf[77671], "full-map peak is twice the shortcut peak here")
# ladder
for k in (5, 10, 20):
    n = 2**k - 1
    for j in range(k):
        n = C(n)
    check(n == 3**k - 1, "ladder k=%d" % k)
# general ladder for a not necessarily odd: value after j<k steps is odd regardless of a
for a in (1, 2, 3, 5, 7, 12):
    for k in (1, 2, 3, 6):
        n = a * 2**k - 1
        for j in range(k):
            n = C(n)
        check(n == a * 3**k - 1, "ladder a=%d k=%d" % (a, k))
print("  ladder C^k(a 2^k-1) = a 3^k-1 confirmed for a in {1,2,3,5,7,12}, k in {1,2,3,6} (a odd is NOT needed)")


def window(p):
    for q, nm in ((p + 1, "-1"), (p - 1, "+1")):
        s = math.isqrt(q)
        if s * s == q:
            return nm, s
    return None


hits = sum(1 for p in pk.values() if window(p) and window(p)[1] % 4 == 1)
hits_any = sum(1 for p in pk.values() if window(p))
byr = Counter(window(p)[1] % 4 for p in pk.values() if window(p))
hits_full = sum(1 for p in pkf.values() if window(p) and window(p)[1] % 4 == 1)
print("  seed-weighted: 4k+1 window %d; any square +-1 %d; by base mod 4 %s; full-map 4k+1 window %d" % (hits, hits_any, dict(byr), hits_full))
check(hits == 248 and hits_any == 535 and byr[1] == 248 and byr[3] == 287 and byr[0] == 0 and byr[2] == 0 and hits_full == 0, "seed-weighted counts")
cnt = Counter(pk.values())
check(len(cnt) == 20341, "distinct peaks")
top = cnt.most_common(4)
print("  distinct peaks %d; top multiplicities %s" % (len(cnt), top))
check(top == [(125252, 448), (4616, 408), (638468, 347), (172772, 161)], "top multiplicities")


def v3(n):
    k = 0
    while n % 3 == 0:
        n //= 3
        k += 1
    return k


print("  v_3(p+1) of the top-multiplicity peaks:", [(p, v3(p + 1)) for p, _ in top])
cells = Counter()
cseeds = Counter()
for p, c in cnt.items():
    w = window(p)
    if w:
        cells[(w[0], w[1] % 4)] += 1
        cseeds[(w[0], w[1] % 4)] += c
        if w[0] == "-1":
            check(w[1] % 3 == 0, "3|s in -1 window")
        else:
            check(w[1] % 3 != 0, "3 !| s in +1 window")
print("  cells (window, base mod 4): distinct %s ; seeds %s" % (dict(cells), dict(cseeds)))
check(cells[("-1", 1)] == 46 and cells[("-1", 3)] == 47 and cells[("+1", 1)] == 49 and cells[("+1", 3)] == 50, "distinct cells")
check(cseeds[("-1", 1)] == 152 and cseeds[("-1", 3)] == 195 and cseeds[("+1", 1)] == 96 and cseeds[("+1", 3)] == 92, "seed cells")


# null model, independent formulation: for each dyadic range, count window values q = 2 mod 6 directly
def window_values(lo, hi, nm, r, extra=lambda q, s: True):
    out = 0
    s = 1
    while s * s - 1 < hi:
        q = s * s + (1 if nm == "+1" else -1)
        if lo <= q < hi and q % 6 == 2 and s % 4 == r and extra(q, s):
            out += 1
        s += 1
    return out


def count_class(lo, hi, mod, res):
    # integers in [lo,hi) congruent to res mod mod
    return len(range(lo + ((res - lo) % mod), hi, mod))


dy = Counter(p.bit_length() - 1 for p in cnt)
null = {}
for nm in ("-1", "+1"):
    for r in (1, 3):
        e = Fraction(0)
        for j, c in dy.items():
            lo, hi = 2**j, 2 ** (j + 1)
            e += Fraction(c * window_values(lo, hi, nm, r), count_class(lo, hi, 6, 2))
        null[(nm, r)] = float(e)
print("  null expected (uniform on 2 mod 6 per dyadic range): %s" % {k: round(v, 1) for k, v in null.items()})
check(abs(null[("-1", 1)] - 20.4) < 0.06 and abs(null[("-1", 3)] - 20.1) < 0.06 and abs(null[("+1", 1)] - 40.8) < 0.06 and abs(null[("+1", 3)] - 40.9) < 0.06, "null values")
plus = cells[("+1", 1)] + cells[("+1", 3)]
minus = cells[("-1", 1)] + cells[("-1", 3)]
pn = null[("+1", 1)] + null[("+1", 3)]
mn = null[("-1", 1)] + null[("-1", 3)]
z = (plus - pn) / math.sqrt(pn)
print("  +1 window %d vs %.1f (z=%.2f); -1 window %d vs %.1f (z=%.2f)" % (plus, pn, z, minus, mn, (minus - mn) / math.sqrt(mn)))
check(abs(z - 1.91) < 0.02 and abs(pn - 81.8) < 0.1 and abs(mn - 40.6) < 0.1, "z and null sums")

# ---- the tautology: every distinct peak p is 3a-1 with ladder seed 2a-1 = (2p-1)/3 whose peak is p
taut = 0
for p in cnt:
    a = (p + 1) // 3
    seed = 2 * a - 1
    check(seed == (2 * p - 1) // 3 and seed % 2 == 1 and C(seed) == p, "predecessor")
    if peak_shortcut(seed) == p:
        taut += 1
print("  TAUTOLOGY CHECK: distinct peaks p for which the k=1 ladder seed (2p-1)/3 has peak p: %d of %d" % (taut, len(cnt)))
check(taut == len(cnt), "the lane's 'explained' criterion with k=1 holds for EVERY peak, hence explains nothing")
# more generally: for every k <= v_3(p+1) the ladder seed has peak exactly p
gen_ok = 0
for p in cnt:
    for k in range(1, v3(p + 1) + 1):
        seed = ((p + 1) // 3**k) * 2**k - 1
        if peak_shortcut(seed) != p:
            raise AssertionError("ladder seed peak mismatch")
    gen_ok += 1
print("  PROVED+checked: for every peak p and every 1<=k<=v_3(p+1), the ladder seed ((p+1)/3^k)2^k-1 has peak p")
print("    (the ladder rises monotonically to p and the orbit of p never exceeds p).  So the lane's")
print("    'explained 93, unexplained none' is a tautology; it does not explain the -1 excess.")

# ---- the real mechanism: v_3(p+1) and the seed box
print()
print("  v_3(p+1) distribution of distinct peaks (null for the class 2 mod 6: P(v>=2)=1/3, P(v>=3)=1/9):")
for label, sel in (("p < 10^5 (box-free: every such p is reached from its own predecessor)", lambda p: p < LIM),
                   ("p >= 10^5 (must be reached from a seed < 10^5)", lambda p: p >= LIM)):
    ps = [p for p in cnt if sel(p)]
    vs = Counter(min(v3(p + 1), 6) for p in ps)
    n_ = len(ps)
    frac2 = sum(c for v, c in vs.items() if v >= 2) / n_
    frac3 = sum(c for v, c in vs.items() if v >= 3) / n_
    print("    %-70s n=%5d  P(v>=2)=%.3f  P(v>=3)=%.3f  hist %s" % (label, n_, frac2, frac3, sorted(vs.items())))
    if sel(LIM - 1):
        check(abs(frac2 - 1 / 3) < 0.03, "box-free peaks: v_3(p+1) is at the null (fraction 1/3 within 0.03)")
        box_free_frac2 = frac2
    else:
        check(frac2 > 0.5, "peaks above the box are strongly biased to 9 | p+1")
        box_frac2 = frac2
# -1 window peaks: how many lie above the box?
minus_ps = sorted(p for p in cnt if window(p) and window(p)[0] == "-1")
above = sum(1 for p in minus_ps if p >= LIM)
print("  -1 window distinct peaks: %d, of which %d are >= 10^5 (above the seed box) and %d below" % (len(minus_ps), above, len(minus_ps) - above))
plus_ps = sorted(p for p in cnt if window(p) and window(p)[0] == "+1")
above_p = sum(1 for p in plus_ps if p >= LIM)
print("  +1 window distinct peaks: %d, of which %d are >= 10^5" % (len(plus_ps), above_p))
# conditional null: uniform on the class 8 mod 18 (9 | p+1), per dyadic range, using the observed number of
# distinct peaks with 9 | p+1 in that range
dy9 = Counter(p.bit_length() - 1 for p in cnt if (p + 1) % 9 == 0)
cond = 0.0
for j, c in dy9.items():
    lo, hi = 2**j, 2 ** (j + 1)
    cond += c * (window_values(lo, hi, "-1", 1) + window_values(lo, hi, "-1", 3)) / count_class(lo, hi, 18, 8)
n9 = sum(dy9.values())
print("  distinct peaks with 9 | p+1: %d (null share 1/3 of %d = %.0f)" % (n9, len(cnt), len(cnt) / 3))
print("  -1 window conditional on 9 | p+1 (uniform on 8 mod 18 per dyadic range): observed %d, expected %.1f, z=%.2f"
      % (minus, cond, (minus - cond) / math.sqrt(cond)))
COND_NULL = cond
COND_Z = (minus - cond) / math.sqrt(cond)
# split the conditional comparison below/above the box
for label, sel in (("p < 10^5", lambda p: p < LIM), ("p >= 10^5", lambda p: p >= LIM)):
    dyx = Counter(p.bit_length() - 1 for p in cnt if (p + 1) % 9 == 0 and sel(p))
    e = sum(c * (window_values(2**j, 2**(j + 1), "-1", 1) + window_values(2**j, 2**(j + 1), "-1", 3)) / count_class(2**j, 2**(j + 1), 18, 8) for j, c in dyx.items())
    obs = sum(1 for p in minus_ps if sel(p))
    print("    %-10s -1 window: observed %d, conditional expected %.1f" % (label, obs, e))
# does the -1 window also prefer LARGE v_3?  compare v_3(p+1) among -1 peaks with the conditional null (v even>=2)
vminus = Counter(v3(p + 1) for p in minus_ps)
print("  v_3(p+1) among -1 window peaks (always even): %s" % sorted(vminus.items()))
# seeds sharing a peak: multiplicity vs v_3(p+1)
mult_by_v = defaultdict(list)
for p, c in cnt.items():
    mult_by_v[min(v3(p + 1), 7)].append(c)
print("  mean multiplicity (seeds per distinct peak) by v_3(p+1):")
for v in sorted(mult_by_v):
    L = mult_by_v[v]
    print("    v=%d%s: n=%5d  mean seeds %.2f  max %d" % (v, "+" if v == 7 else " ", len(L), sum(L) / len(L), max(L)))
check(sum(mult_by_v[1]) / len(mult_by_v[1]) < sum(mult_by_v[3]) / len(mult_by_v[3]), "multiplicity grows with v_3(p+1)")

# ---- the 2-adic side: s odd (p even) forces 8 | s^2-1, while s^2+1 = 2 mod 8
def v2(n):
    k = 0
    while n % 2 == 0:
        n //= 2
        k += 1
    return k


print()
print("  v_2(p) among -1 window peaks: %s ; among +1 window peaks: %s"
      % (sorted(Counter(min(v2(p), 6) for p in minus_ps).items()), sorted(Counter(v2(p) for p in plus_ps).items())))
check(all(p % 8 == 0 for p in minus_ps) and all(p % 8 == 2 for p in plus_ps), "8 | s^2-1 and s^2+1 = 2 mod 8 for odd s")
# PROVED: no shortcut peak is 6 mod 8 (p -> p/2 = 3 mod 4 -> (3p+2)/4 odd -> (9p+10)/8 > p)
check(all(p % 8 != 6 for p in cnt), "no peak is 6 mod 8")
m24 = sorted(Counter(p % 24 for p in cnt).items())
print("  PROVED: no shortcut peak is 6 mod 8, so peaks are 2, 8, 20 mod 24.  Census by p mod 24: %s" % m24)
check([r for r, _ in m24] == [2, 8, 20], "peaks mod 24 in {2,8,20}")
# below the box the peak set is exact
below_set = set(q for q in range(8, LIM, 6) if peak_shortcut(q) == q)
check(below_set == set(p for p in cnt if p < LIM), "distinct peaks below the box = q = 2 mod 6 with orbit never exceeding q")
print("  distinct peaks below the box = %d = #{q = 2 mod 6, q < 10^5 : orbit(q) never exceeds q} (exact)" % len(below_set))
cand24 = Counter(q % 24 for q in range(8, LIM, 6))
pk24 = Counter(p % 24 for p in below_set)
print("  below-box peak rate by q mod 24: %s" % {r: "%d/%d=%.3f" % (pk24[r], cand24[r], pk24[r] / cand24[r]) for r in sorted(cand24)})
check(pk24[14] == 0 and abs(pk24[8] / cand24[8] - 0.889) < 0.002 and abs(pk24[2] / cand24[2] - 0.688) < 0.002 and abs(pk24[20] / cand24[20] - 0.685) < 0.002, "rates mod 24")
# multiplicity lemma: unique seed iff 9 !| p+1 ; ladder = odd-preimage chain of length v_3(p+1)
check(all(c == 1 for p, c in cnt.items() if (p + 1) % 9 != 0), "9 !| p+1 => unique seed")
m_in = [p for p in cnt if (p + 1) % 9 == 0 and (2 * p - 1) // 3 <= LIM]
check(all(cnt[p] >= 2 for p in m_in), "9 | p+1 and m <= 10^5 => at least two seeds")
for p in cnt:
    chain, x = [], p
    while (2 * x - 1) % 3 == 0:
        x = (2 * x - 1) // 3
        chain.append(x)
    check(len(chain) == v3(p + 1) and all(chain[k - 1] == ((p + 1) // 3**k) * 2**k - 1 for k in range(1, len(chain) + 1)), "odd-preimage chain = ladder")
lad_in = sum(sum(1 for k in range(1, v3(p + 1) + 1) if ((p + 1) // 3**k) * 2**k - 1 <= LIM) for p in cnt)
print("  PROVED+checked: unique seed iff 9 !| p+1 (%d such peaks); the odd-preimage chain of p is exactly the ladder" % sum(1 for p in cnt if (p + 1) % 9 != 0))
print("    seeds ((p+1)/3^k)2^k-1, of length v_3(p+1); %d of the %d census seeds are ladder seeds, %d are not." % (lad_in, sum(cnt.values()), sum(cnt.values()) - lad_in))
check(lad_in == 28483 and sum(cnt.values()) - lad_in == 21516, "ladder seed split")


# stratified null (a): (dyadic range, min(v_2,4), min(v_3(q+1),4)) -- lumps 2 mod 8 with 6 mod 8: a documented artefact
def strat_a(q):
    return (q.bit_length() - 1, min(v2(q), 4), min(v3(q + 1), 4))


# stratified null (b): (dyadic range, q mod 24, min(v_3(q+1),4)) -- the correct one
def strat_b(q):
    return (q.bit_length() - 1, q % 24, min(v3(q + 1), 4))


def build_rates(strat, M, key_of_r):
    cand_s = Counter()
    cls = {r: key_of_r(r) for r in range(M) if r % 6 == 2}
    for j in range(3, 31):
        lo, hi = 2**j, 2 ** (j + 1)
        for r, key in cls.items():
            t_lo = max(0, math.ceil((lo - r) / M))
            t_hi = math.ceil((hi - r) / M)
            if t_hi > t_lo:
                cand_s[(j,) + key] += t_hi - t_lo
    check(sum(cand_s.values()) == sum(count_class(2**j, 2**(j + 1), 6, 2) for j in range(3, 31)), "stratum sizes")
    peak_s = Counter(strat(p) for p in cnt)
    return {k: peak_s[k] / cand_s[k] for k in cand_s if cand_s[k]}


rate_a = build_rates(strat_a, 2**5 * 3**5, lambda r: (min(v2(r), 4), min(v3(r + 1), 4)))
rate_b = build_rates(strat_b, 24 * 81, lambda r: (r % 24, min(v3(r + 1), 4)))


def strat_expect(rate, strat, nm, r_mod4=None, sel=lambda q: True):
    e, s = 0.0, 1
    while s * s - 1 < 2**31:
        q = s * s + (1 if nm == "+1" else -1)
        if q >= 8 and q % 6 == 2 and (r_mod4 is None or s % 4 == r_mod4) and sel(q):
            e += rate.get(strat(q), 0.0)
        s += 1
    return e


print("  stratified null (a) by (range, v_2, v_3) -- WRONG because v_2=1 lumps 2 mod 8 (rate ~0.69) with 6 mod 8 (rate 0):")
za = {}
for nm, obs in (("-1", minus_ps), ("+1", plus_ps)):
    e = strat_expect(rate_a, strat_a, nm)
    za[nm] = (len(obs) - e) / math.sqrt(e)
    print("    %s window: observed %d, expected %.1f, z=%.2f" % (nm, len(obs), e, za[nm]))
check(abs(za["-1"]) < 3 and za["+1"] > 3, "null (a): -1 at chance, +1 spuriously in excess")
print("  STRATIFIED NULL (b) by (dyadic range, q mod 24, v_3(q+1) capped at 4), exact stratum sizes:")
ZB = {}
for nm, obs in (("-1", minus_ps), ("+1", plus_ps)):
    e_all = strat_expect(rate_b, strat_b, nm)
    e_lo = strat_expect(rate_b, strat_b, nm, sel=lambda q: q < LIM)
    e_hi = strat_expect(rate_b, strat_b, nm, sel=lambda q: q >= LIM)
    o_lo = sum(1 for p in obs if p < LIM)
    o_hi = len(obs) - o_lo
    ZB[nm] = (len(obs), e_all, (len(obs) - e_all) / math.sqrt(e_all))
    print("    %s window: observed %d (below box %d, above %d); expected %.1f (%.1f, %.1f); z=%.2f"
          % (nm, len(obs), o_lo, o_hi, e_all, e_lo, e_hi, ZB[nm][2]))
    for r in (1, 3):
        o = sum(1 for p in obs if math.isqrt(p + (1 if nm == "-1" else -1)) % 4 == r)
        e = strat_expect(rate_b, strat_b, nm, r_mod4=r)
        print("      base %d mod 4: observed %d, expected %.1f" % (r, o, e))
check(abs(ZB["-1"][2]) < 1 and abs(ZB["+1"][2]) < 1, "both windows at chance under null (b) (|z|<1)")
check(abs(ZB["-1"][1] - 88.7) < 0.1 and abs(ZB["+1"][1] - 100.9) < 0.1, "null (b) expectations 88.7 / 100.9")
# below the box, coarser and finer stratifications agree (q mod 24 already suffices)
Wlo = {nm: [q for q in range(8, LIM, 6) if window(q) and window(q)[0] == nm] for nm in ("-1", "+1")}
for M in (24, 1152, 27648):
    cand = Counter(q % M for q in range(8, LIM, 6))
    pkc = Counter(p % M for p in below_set)
    rt = {r: pkc[r] / cand[r] for r in cand}
    out = []
    for nm in ("-1", "+1"):
        e = sum(rt[q % M] for q in Wlo[nm])
        o = sum(1 for q in Wlo[nm] if q in below_set)
        out.append("%s: %d vs %.1f" % (nm, o, e))
        check(abs(o - e) < 2 * math.sqrt(e), "below-box window at chance mod %d" % M)
    print("  below the box, null by q mod %5d: %s" % (M, "; ".join(out)))

# shifted windows (independent code)
row = []
for s in range(-12, 13, 2):
    h = 0
    for p in pk.values():
        for q in (p - s - 1, p - s, p - s + 1):
            if q >= 1 and math.isqrt(q) ** 2 == q and math.isqrt(q) % 4 == 1:
                h += 1
                break
    row.append((s, h))
print("  shifted windows:", row)
check(row == [(-12, 23), (-10, 0), (-8, 24), (-6, 115), (-4, 91), (-2, 152), (0, 248), (2, 96), (4, 0), (6, 150), (8, 150), (10, 21), (12, 21)], "shifted windows")
# ladder multiplicities and 196/169/160
check(cnt[6560] == 21 and cnt[59048] == 48 and cnt[164024] == 17 and cnt[4782968] == 17, "ladder multiplicities")
check(6560 == 3**8 - 1 and 59048 == 3**10 - 1 and 164024 == 25 * 3**8 - 1 == 405**2 - 1 and 4782968 == 3**14 - 1, "ladder values")
pr = list(sympy.primerange(2, 38))
check(len(pr) == 12 and sum(pr) == 197 and sum(pr[:11]) == 160, "first 12 primes sum 197 (first 11 sum 160)")
print("  first 12 primes sum to %d (first 11 to %d)" % (sum(pr), sum(pr[:11])))
o196 = [196]
while o196[-1] != 7:
    o196.append(C(o196[-1]))
check(o196 == [196, 98, 49, 74, 37, 56, 28, 14, 7], "196 orbit")
check(sorted(n for n in pk if pk[n] == 170) == [75, 113] and not [n for n in pk if pk[n] in (168, 160)], "170/168/160")
check(sorted(n for n in pkf if pkf[n] == 160) == [15, 23, 35, 53], "full-map peak 160")
print("  196 orbit, peak-170 seeds [75,113], full-map-peak-160 seeds [15,23,35,53]: confirmed")

# ----------------------------------------------------------------------------- A8
hdr("A8. Descent certificate table, and the stopping-time terminology")


def syr(n, L):
    ks = []
    for _ in range(L):
        m = 3 * n + 1
        k = 0
        while m % 2 == 0:
            m //= 2
            k += 1
        ks.append(k)
        n = m
    return ks, n


exp = {(27, 5): (6, -5120), (7, 2): (2, -40), (23, 3): (7, 2304), (1, 1): (2, 0), (31, 5): (6, -5760), (2**20 - 1, 20): (24, -3638566269747200)}
for (n, L), (K_e, m_e) in exp.items():
    ks, fin = syr(n, L)
    K = sum(ks)
    Bv = sum(3 ** (L - 1 - j) * 2 ** sum(ks[:j]) for j in range(L))
    check(2**K * fin == 3**L * n + Bv, "identity n=%d L=%d" % (n, L))
    check(K == K_e and 2**K * n - (3**L * n + Bv) == m_e, "margin n=%d L=%d" % (n, L))
ks, fin = syr(27, 41)
check(sum(ks) == 70 and fin == 1, "27: L=41 odd steps, K=70, reaches 1")
ks37, f37 = syr(27, 37)
check(sum(ks37) == 59 and f37 == 23, "27: L=37, K=59, value 23")
# first L with positive margin
ksall, _ = syr(27, 41)
firstL = next(L for L in range(1, 42) if 2 ** sum(ksall[:L]) * 27 > 3**L * 27 + sum(3 ** (L - 1 - j) * 2 ** sum(ksall[:j]) for j in range(L)))
check(firstL == 37, "first positive margin at L=37")
# full-map stopping time (first index with value < 27) and total stopping time (first index with value 1)
x, i, st, tot = 27, 0, None, None
while x != 1:
    x = x // 2 if x % 2 == 0 else 3 * x + 1
    i += 1
    if st is None and x < 27:
        st = i
tot = i
print("  27 under the full map: stopping time (first value < 27) = %d, total stopping time (reach 1) = %d" % (st, tot))
check(st == 96 and tot == 111, "stopping time 96, total stopping time 111")
print("  The lane note calls 96 the 'classical total stopping time': WRONG LABEL.  96 = sigma(27), the stopping time;")
print("  the total stopping time sigma_inf(27) = 111 = 41 + 70.")
print("  margins table (K, margin): confirmed for all seven rows.")

# ----------------------------------------------------------------------------- A9
hdr("A9. mod-30 wheel")
check([r for r in range(30) if math.gcd(r, 30) == 1] == [1, 7, 11, 13, 17, 19, 23, 29], "wheel")
print("  confirmed")

# ----------------------------------------------------------------------------- summary
hdr("AUDIT SUMMARY")
print("  CONFIRMED: A1 (identity/repair/flows), A2 (simplex laws), A3 (g-operator; 'symmetric' load-bearing),")
print("    A4 (census; note 35,39 also absent at n=6), A5 (Bang/Zsigmondy/Cipolla; wild_typing lane already has it),")
print("    A6 (3-cycle), A7 counts (248/535/287/20341/cells/null/z/shifted windows/ladder values), A8 table, A9.")
print("  REFUTED AS STATED: section 8's 'the excess is entirely Lemma 8.2 ... explained 93, unexplained none' --")
print("    the criterion holds for all %d distinct peaks (k=1 is the immediate predecessor), so it is vacuous." % len(cnt))
print("  TRUE MECHANISM (FINITE-EXACT): distinct peaks below the seed box have 9 | p+1 at the null rate;")
print("    peaks above the box are biased to 9 | p+1 (ladder ancestors are small); conditional on 9 | p+1 alone the")
print("    -1 square window is still observed %d vs expected %.1f (z=%.2f), because s odd forces 8 | s^2-1 and no" % (minus, COND_NULL, COND_Z))
print("    peak is 6 mod 8 (peak rates 0.688 / 0.889 / 0 / 0.685 on 2 / 8 / 14 / 20 mod 24).  Stratifying by")
print("    (dyadic range, q mod 24, v_3(q+1)) puts BOTH windows at chance: -1 window %d vs %.1f (z=%.2f), +1 window %d vs %.1f (z=%.2f)."
      % (ZB["-1"][0], ZB["-1"][1], ZB["-1"][2], ZB["+1"][0], ZB["+1"][1], ZB["+1"][2]))
print("  NEW EXACT CONTENT (elementary): peaks are 2,8,20 mod 24; unique seed iff 9 !| p+1; ladder = odd-preimage chain of length v_3(p+1).")
print("  WRONG LABEL: '37+59=96 is the classical total stopping time' -> 96 is the stopping time; 111 is total.")
print()
print("FAILS: %d" % len(FAILS))
print("ALL AUDIT CHECKS PASSED" if not FAILS else "SOME AUDIT CHECKS FAILED")
