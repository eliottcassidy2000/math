"""Independent referee re-computation (audit of opus-2026-10-28-S67) -- THM-4838 numerics.

von Staudt-Clausen sets, orbifold Euler characteristics, triangle groups (2,3,p), genus of X(p) and X_0(p),
Monster / Baby Monster orders (from their DECIMAL orders, factored here), J-homomorphism denominators,
the M24 eta list, the vSC sinc-product deficits (by the BB 2001 sign-sum formula, not by the session's
probabilistic code), sigma^{-1}(144) and the prime triples with prod(p+1) = 144.

Run:  python audit_vsc_monster.py
"""
import math
from fractions import Fraction as Fr
from itertools import combinations, product as iproduct

from sympy import bernoulli, isprime, primerange, factorint, divisors, divisor_sigma, zeta, Rational, pi, N as Nnum

FAIL = []


def say(*a):
    print(*a, flush=True)


def check(cond, label):
    say(("  [PASS] " if cond else "  [FAIL] ") + label)
    if not cond:
        FAIL.append(label)


def B(n):
    b = bernoulli(n)
    return Fr(int(b.p), int(b.q))


def vsc(n2k):
    """{p prime : (p-1) | n2k}"""
    return [d + 1 for d in divisors(n2k) if isprime(d + 1)]


# ======================================================================== 1
say("=" * 96)
say("1. vSC(2k) = {p : (p-1) | 2k}; chi_k = 2 - |vSC| + sum 1/p; N_k = B_2k + sum 1/p (2k <= 200)")
say("=" * 96)
exact, Nk = [], {}
mod_ok = True
growth_ok = True
for k in range(1, 101):
    v = vsc(2 * k)
    s = sum(Fr(1, p) for p in v)
    chi = 2 - len(v) + s
    Nk[2 * k] = s + B(2 * k)
    if Nk[2 * k].denominator != 1 or (chi + B(2 * k)).denominator != 1:
        mod_ok = False
    if chi == -B(2 * k):
        exact.append(2 * k)
    if 2 * k >= 20 and not (abs(B(2 * k)) > 2 * math.sqrt(2 * k) + 2):
        growth_ok = False
    if len(v) > len(divisors(2 * k)) or len(divisors(2 * k)) > 2 * math.sqrt(2 * k):
        growth_ok = False
check(mod_ok, "von Staudt-Clausen and chi_k = -B_2k (mod 1) for every 2k <= 200")
check(exact == [4, 6, 8, 10], "exact equality chi_k = -B_2k only for 2k in {4, 6, 8, 10} (2k <= 200)")
check(all(Nk[2 * k] == 1 for k in range(1, 7)), "N_k = 1 for 2k <= 12")
say("  N_k for 2k = 14, 16, 18, 20: %s" % [int(Nk[j]) for j in (14, 16, 18, 20)])
check(growth_ok, "|B_2k| > 2 sqrt(2k) + 2 >= d(2k) + 2 >= |vSC(2k)| + 2 for 20 <= 2k <= 200")
# analytic step: |B_{2k+2}/B_{2k}| = (2k+1)(2k+2)/(4 pi^2) * zeta(2k+2)/zeta(2k) >= 11 for 2k >= 20
rat = min(float(abs(Rational(B(2 * k + 2).numerator, B(2 * k + 2).denominator) /
                    Rational(B(2 * k).numerator, B(2 * k).denominator))) for k in range(10, 100))
say("  min_{20 <= 2k <= 198} |B_(2k+2)/B_2k| = %.3f  (>= (21*22/(4 pi^2)) * zeta(22)/zeta(20) = %.3f)"
    % (rat, float(Nnum(Rational(21 * 22) / (4 * pi ** 2) * zeta(22) / zeta(20)))))
say("  -> the growth step needed by THM-4838 (1) and (2) is a one-line ratio argument; it is NOT in the theorem file.")

# ======================================================================== 2
say("")
say("=" * 96)
say("2. Triangle groups (2,3,p): chi = 1/p - 1/6 vs -B_(p-1)")
say("=" * 96)
ex, md, safe_like = [], [], []
for p in primerange(3, 2000):
    chi = Fr(1, p) - Fr(1, 6)
    if chi == -B(p - 1):
        ex.append(p)
    if (chi + B(p - 1)).denominator == 1:
        md.append(p)
    if sorted(set(vsc(p - 1))) == sorted({2, 3, p}):
        safe_like.append(p)
check(ex == [5, 7, 11], "chi(2,3,p) = -B_(p-1) exactly iff p in {5, 7, 11} (p < 2000)")
say("  congruence mod 1 holds for p in %s ..." % md[:16])
say("  vSC(p-1) = {2,3,p} as SETS for p in %s ..." % safe_like[:16])
check([p for p in md if p >= 5] == [p for p in safe_like if p >= 5], "for p >= 5: congruence iff vSC(p-1) = {2,3,p}")
check(3 in safe_like and 3 not in md,
      "p = 3 is a COUNTEREXAMPLE to THM-4838 (2) as scoped 'p >= 3': vSC(2) = {2,3} = {2,3,p} but chi(2,3,3) + B_2 = 1/3")
check(md[:8] == [5, 7, 11, 23, 47, 59, 83, 107], "list 5, 7, 11, 23, 47, 59, 83, 107 (note line 151)")

# ======================================================================== 3
say("")
say("=" * 96)
say("3. X(p) -> X(1) branched (2,3,p): 2g-2 = |PSL2(F_p)| (1/6 - 1/p); genus of X_0(p)")
say("=" * 96)
for p, g in [(5, 0), (7, 3), (11, 26)]:
    G = p * (p * p - 1) // 2
    two_g_minus_2 = G * (Fr(1, 6) - Fr(1, p))
    check(two_g_minus_2 == 2 * g - 2 == G * B(p - 1),
          "p = %d: |PSL2| = %d, 2g - 2 = %s = |PSL2| * B_%d (g = %d)" % (p, G, two_g_minus_2, p - 1, g))


def genus_X0_prime(p):
    mu = p + 1
    if p == 2:
        n2, n3 = 1, 0
    elif p == 3:
        n2, n3 = 0, 1
    else:
        n2 = 2 if p % 4 == 1 else 0
        n3 = 2 if p % 3 == 1 else 0
    g = 1 + Fr(mu, 12) - Fr(n2, 4) - Fr(n3, 3) - 1     # nu_infinity = 2 for prime level
    assert g.denominator == 1
    return int(g)


gz = [p for p in primerange(2, 5000) if genus_X0_prime(p) == 0]
check(gz == [2, 3, 5, 7, 13] == [p for p in primerange(2, 5000) if 24 % (p - 1) == 0]
      == [p for p in primerange(2, 5000) if 12 % (p - 1) == 0] == sorted(vsc(12)),
      "genus X_0(p) = 0 iff p in {2,3,5,7,13} = {p:(p-1)|24} = {p:(p-1)|12} = vSC(12) (p < 5000)")
check(genus_X0_prime(11) == 1 and genus_X0_prime(37) == 2, "sanity: g(X_0(11)) = 1, g(X_0(37)) = 2")

# ======================================================================== 4
say("")
say("=" * 96)
say("4. Monster and Baby Monster (decimal orders from Wikipedia, factored here)")
say("=" * 96)
M = 808017424794512875886459904961710757005754368000000000
BM = 4154781481226426191177580544000000
fM, fB = factorint(M), factorint(BM)
say("  |M|  = %s" % " ".join("%d^%d" % (p, e) for p, e in sorted(fM.items())))
say("  |B|  = %s" % " ".join("%d^%d" % (p, e) for p, e in sorted(fB.items())))
check(fM == {2: 46, 3: 20, 5: 9, 7: 6, 11: 2, 13: 3, 17: 1, 19: 1, 23: 1, 29: 1, 31: 1, 41: 1, 47: 1, 59: 1, 71: 1},
      "|M| factorisation")
check(fB == {2: 41, 3: 13, 5: 6, 7: 2, 11: 1, 13: 1, 17: 1, 19: 1, 23: 1, 31: 1, 47: 1}, "|B| factorisation")
check(sorted(p for p, e in fM.items() if e >= 3) == [2, 3, 5, 7, 13], "p^3 | |M| iff p in {2,3,5,7,13}")
q = M // BM
assert q * BM == M
fq = factorint(q)
say("  |M|/|B| = %s" % " ".join("%d^%d" % (p, e) for p, e in sorted(fq.items())))
core = {p: e for p, e in fq.items() if p <= 13}
check(core == {2: 5, 3: 7, 5: 3, 7: 4, 11: 1, 13: 2} and 2 ** 5 * 3 ** 7 * 5 ** 3 * 7 ** 4 * 11 * 13 ** 2 == 39046339332000,
      "small-prime (p <= 13) core of |M|/|B| = 2^5 3^7 5^3 7^4 11 13^2 = 39,046,339,332,000 (Session43 N)")
say("  primes > 13 in |M|/|B|: %s (the 'core' is a choice: it drops 29, 41, 59, 71)" % sorted(p for p in fq if p > 13))
check(sorted(p for p, e in core.items() if e >= 2) == [2, 3, 5, 7, 13], "exponent >= 2 primes of N = {2,3,5,7,13}")
check(all(fM[p] >= 1 + 24 // (p - 1) for p in (2, 3, 5, 7, 13)) and fM[13] == 3,
      "room for p^(1+24/(p-1)) in |M| (13^(1+2) is a full Sylow 13-subgroup)")

# ======================================================================== 5
say("")
say("=" * 96)
say("5. J-homomorphism / Eisenstein denominators; {p:(p+1)|24}; sigma = 144")
say("=" * 96)
dens = [(B(2 * k) / (4 * k)).denominator for k in range(1, 7)]
say("  denominators of B_2k/(4k), k = 1..6: %s (|im J| in pi_(4k-1)^s)" % dens)
check(dens == [24, 240, 504, 480, 264, 65520], "264 (pi_19) and 65520 (pi_23) as stated")
check(-Fr(20) / B(10) == -264, "E_10 = 1 - 264 sum sigma_9(n) q^n")
p24 = [p for p in primerange(2, 1000) if 24 % (p + 1) == 0]
p12 = [p for p in primerange(2, 1000) if 12 % (p + 1) == 0]
m12 = [p for p in primerange(2, 1000) if 12 % (p - 1) == 0]
check(p24 == [2, 3, 5, 7, 11, 23] and p12 == [2, 3, 5, 11] and m12 == [2, 3, 5, 7, 13],
      "{p:(p+1)|24} = {2,3,5,7,11,23}; {p:(p+1)|12} = {2,3,5,11}; {p:(p-1)|12} = {2,3,5,7,13}")
d12 = sorted(d + 1 for d in divisors(12))
check(d12 == [2, 3, 4, 5, 7, 13], "{d+1 : d | 12} = {2,3,4,5,7,13}")
pre = [n for n in range(1, 200) if divisor_sigma(n) == 144]
check(pre == [66, 70, 94, 115, 119], "sigma^{-1}(144) = {66, 70, 94, 115, 119} (sigma(n) > n, so n < 144 suffices)")
cands = [p for p in primerange(2, 145) if 144 % (p + 1) == 0]
sets = [c for r in range(1, len(cands) + 1) for c in combinations(cands, r)
        if math.prod(p + 1 for p in c) == 144]
say("  prime sets with prod(p+1) = 144: %s" % sets)
check(sorted(c for c in sets if len(c) == 3) == [(2, 3, 11), (2, 5, 7)], "only 3-prime sets: {2,3,11}, {2,5,7}")
check([c for c in sets if set(c) <= {2, 3, 5, 7, 13}] == [(2, 5, 7)], "unique subset of vSC(12) with prod(p+1) = 144")

# ======================================================================== 6
say("")
say("=" * 96)
say("6. vSC sinc products by the BB 2001 Theorem 2 sign-sum formula (independent of the session's code)")
say("=" * 96)


def bb_value_over_pi(a0, avals):
    """int_0^oo sinc(a0 x) prod sinc(a_k x) dx / pi = C_n/(2 a0),
    C_n = (2^n n! prod a)^-1 sum_gamma eps_gamma b_gamma^n sgn(b_gamma)  (BB 2001 Thm 2; Wikipedia)."""
    n = len(avals)
    S = Fr(0)
    for g in iproduct((1, -1), repeat=n):
        b = a0 + sum(gi * ai for gi, ai in zip(g, avals))
        eg = 1
        for gi in g:
            eg *= gi
        S += eg * b ** n * ((b > 0) - (b < 0))
    pa = Fr(1)
    for ai in avals:
        pa *= ai
    return S / (2 ** n * math.factorial(n) * pa) / (2 * a0)


for k in range(1, 7):
    v = sorted(vsc(2 * k))
    val = bb_value_over_pi(Fr(1), [Fr(1, p) for p in v])
    deficit = Fr(1, 2) - val
    say("  2k = %-3d vSC = %-18s excess = %-10s deficit/pi = %s" % (2 * k, v, sum(Fr(1, p) for p in v) - 1, deficit))
    if 2 * k in (2, 6, 10):
        check(deficit == 0, "2k = %d: exactly pi/2" % (2 * k))
    if 2 * k in (4, 8):
        check(deficit == Fr(1, 43200), "2k = %d: deficit pi/43200" % (2 * k))
    if 2 * k == 12:
        check(deficit == Fr(691 ** 5 - 271 ** 5, 3840 * 2730 ** 4), "2k = 12: deficit pi (691^5 - 271^5)/(3840 * 2730^4)")
# BBB 2008 Example 2(b) (reciprocals of ALL primes, a0 = 1/2): pi(1 - 1/240) for N = 2
v = bb_value_over_pi(Fr(1, 2), [Fr(1, 3), Fr(1, 5)])
check(v == Fr(1) - Fr(1, 240),
      "BBB 2008 Ex. 2(b): int_0^oo sinc(x/2) sinc(x/3) sinc(x/5) dx = pi (1 - 1/240) [prior art: prime-reciprocal sinc products]")
say("  (value/pi = %s)" % v)

say("")
say("FAILURES: %s" % (FAIL if FAIL else "none"))
