#!/usr/bin/env python3
"""
Vieta, Schanuel, and the repo's rational/irrational threads        (mac-mini-S138)
=================================================================================
Owner: "consider Schanuel's conjecture and the algebraic proof by contradiction using
Vieta's formula, the quadratic polynomial proof for e and pi, and previous work in the
repo regarding the multiplication/addition duality. come to understand threads we have
done regarding rationals and irrationals."

THE VIETA ARGUMENT (three lines, classical).  If e+pi AND e*pi were both algebraic,
then e and pi would be the two roots of  x^2 - (e+pi)x + (e*pi) = 0,  a quadratic with
ALGEBRAIC coefficients -- so e and pi would themselves be algebraic, contradicting
Hermite and Lindemann.  Hence AT LEAST ONE of e+pi, e*pi is transcendental.  Which one
is unknown.  SCHANUEL'S CONJECTURE would give BOTH, and much else.

WHY THIS IS THE ADDITION/MULTIPLICATION DUALITY.  Vieta says the SUM and the PRODUCT are
the two coefficients; the argument is precisely that you cannot have both the additive
and the multiplicative combination be tame while the elements are wild.  And Schanuel is
the same duality at the transcendence level: Q-LINEAR (additive) independence of the z_i
forces the exponentials e^{z_i} to contribute genuinely new algebraic content.  exp is
the bridge -- exactly the role log plays between the repo's two arithmetics (the additive
tiling hypercube and the multiplicative H-norm, THM-1460 / death-star-S60).

WHAT IS ACTUALLY CHECKED HERE (the analogy is NOT the deliverable):
  A  Vieta IS already in canon: char-poly coefficients ARE the elementary symmetric
     functions of the spectrum.  For skew-Seidel the coefficients are INTEGERS while the
     eigenvalues are irrational -- the ALGEBRAIC side of the Vieta dichotomy, in contrast
     with e/pi which sits on the transcendental side.
  C  A rational/transcendental dichotomy inside MY OWN THM-1520: the saddle limit is
     exp(a_{D-1}/(D a_D)), which is EXACTLY 1 when a_{D-1} = 0 and TRANSCENDENTAL
     otherwise (Lindemann).  Verified numerically.
  D  The repo's own additive/multiplicative bridge, made exact: dim_Q span{log H(T)}
     equals the number of DISTINCT PRIMES dividing the H-values -- by unique
     factorization, no transcendence needed.  The elementary shadow of Baker/Schanuel.
"""
from fractions import Fraction as F
from math import factorial, exp, log, gcd
from itertools import permutations
import numpy as np

# ================================================================ PART A
print("=" * 78)
print("PART A -- Vieta is already canon: char-poly coefficients = elementary symmetric")
print("=" * 78)
print("  THM-1440's n=7 cospectral pair has  p(x) = x(x^2+7)(x^4+14x^2+17),")
print("  integer coefficients [1,0,21,0,115,0,119,0], with eigenvalues of iS equal to")
print("  0, +-sqrt(7), +-sqrt(7-4sqrt2), +-sqrt(7+4sqrt2)  -- IRRATIONAL, splitting field Q(sqrt2).")
print()
roots = [0.0, 7**0.5, -(7**0.5),
         (7 - 4*2**0.5)**0.5, -((7 - 4*2**0.5)**0.5),
         (7 + 4*2**0.5)**0.5, -((7 + 4*2**0.5)**0.5)]
coeffs = np.poly(roots)          # elementary symmetric functions, up to sign
print(f"{'k':>3} {'e_k(roots) via np.poly':>26} {'nearest integer':>16} {'|error|':>12}")
for k, c in enumerate(coeffs):
    print(f"{k:>3} {c:>26.10f} {round(c):>16} {abs(c-round(c)):>12.2e}")
print("  => every elementary symmetric function of these irrational eigenvalues is an")
print("     INTEGER.  That is the ALGEBRAIC branch of Vieta: all symmetric functions tame,")
print("     the elements themselves wild but still algebraic.")
print()
print("  THE CONTRAST WITH e AND pi.  There the Vieta argument runs the other way:")
print("    if e+pi and e*pi were BOTH algebraic, e and pi would be roots of")
print("        x^2 - (e+pi) x + (e*pi)")
print("    with algebraic coefficients, hence algebraic -- contradicting Hermite/Lindemann.")
print("    So AT LEAST ONE of e+pi, e*pi is transcendental.  WHICH ONE IS UNKNOWN.")
print("    Schanuel's conjecture would give BOTH (and the algebraic independence of e,pi).")
print(f"    numerically: e+pi = {exp(1)+np.pi:.10f},  e*pi = {exp(1)*np.pi:.10f}")
print("    Neither is known to be irrational individually.  The DISJUNCTION is a theorem;")
print("    each DISJUNCT is open.  That asymmetry is the whole flavour of the subject.")

# ================================================================ PART C
print()
print("=" * 78)
print("PART C -- a rational/transcendental dichotomy INSIDE THM-1520")
print("=" * 78)
print("  THM-1520's saddle lemma:  L(p^m)/(a_D^m (Dm)!)  ->  exp( a_{D-1} / (D a_D) ).")
print("  If p has rational (or algebraic) coefficients, the exponent is rational (algebraic),")
print("  so by LINDEMANN the limit is TRANSCENDENTAL unless the exponent is 0.  Hence:")
print()
print("      the saddle limit is RATIONAL  <=>  it equals exactly 1  <=>  a_{D-1} = 0.")
print()
def L(coeffs):
    return sum(a*factorial(k) for k, a in enumerate(coeffs))
def polymul(a, b):
    out = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): out[i+j] += x*y
    return out
def polypow(a, m):
    r = [1]
    for _ in range(m): r = polymul(r, a)
    return r
print(f"{'p':>20} {'a_{D-1}':>9} {'exponent':>10} {'predicted limit':>18} {'m=14 ratio':>13} {'kind':>15}")
for co, nm in (([0,1], "v"), ([-1,1], "v - 1"), ([3,1], "v + 3"),
               ([2,-3,1], "v^2-3v+2"), ([0,0,1], "v^2"), ([0,1,0,2], "2v^3 + v"),
               ([5,0,2,1], "v^3+2v^2+5")):
    D = len(co)-1
    r = F(co[D-1], D*co[D])
    pred = exp(float(r))
    val = L(polypow(co, 14)); ratio = float(F(val, factorial(D*14)) / F(co[D])**14)
    kind = "RATIONAL (=1)" if r == 0 else "transcendental"
    print(f"{nm:>20} {co[D-1]:>9} {str(r):>10} {pred:>18.9f} {ratio:>13.6f} {kind:>15}")
print()
print("  So the one analytic constant this session produced is e^rational -- and it is")
print("  rational EXACTLY on the codimension-one locus a_{D-1} = 0.  The dichotomy is")
print("  clean and it is Lindemann doing the work, not anything new.")

# ================================================================ PART D
print()
print("=" * 78)
print("PART D -- the repo's OWN additive/multiplicative bridge, made exact")
print("=" * 78)
print("  H is multiplicative under ordinal sum (H(T1 (+) T2) = H(T1)H(T2)), so log H is")
print("  ADDITIVE.  The transcendence-flavoured question 'are these logs independent?' has")
print("  an ELEMENTARY answer here: {log p : p prime} is Q-linearly independent by UNIQUE")
print("  FACTORIZATION (clear denominators; prod p^{n_i} = 1 forces all n_i = 0).  Hence")
print()
print("      dim_Q span{ log H(T) }  =  #{ distinct primes dividing some H(T) }.")
print()
def scaffold(n):
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    return pairs, {p: k for k, p in enumerate(pairs)}, len(pairs)
def ham_paths(w, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for S in range(1 << n):
        for v in range(n):
            c = dp[S][v]
            if not c or not (S >> v & 1): continue
            for u in range(n):
                if S >> u & 1 or not w[v][u]: continue
                dp[S | 1 << u][u] += c
    return sum(dp[(1 << n)-1])
def factor(x):
    f = {}; d = 2
    while d*d <= x:
        while x % d == 0: f[d] = f.get(d, 0)+1; x //= d
        d += 1
    if x > 1: f[x] = f.get(x, 0)+1
    return f

print(f"{'n':>3} {'#H values':>10} {'H set':>34} {'distinct primes':>16} {'dim_Q':>6}")
for n in range(3, 8):
    pairs, idx, E = scaffold(n)
    Hs = set()
    for code in range(1 << E):
        w = [[0]*n for _ in range(n)]
        for e, (i, j) in enumerate(pairs):
            if code >> e & 1: w[j][i] = 1
            else:             w[i][j] = 1
        Hs.add(ham_paths(w, n))
        if E > 15 and len(Hs) > 60: break
    primes = set()
    for h in Hs:
        if h > 1: primes |= set(factor(h))
    shown = sorted(Hs)[:6]
    print(f"{n:>3} {len(Hs):>10} {str(shown)+('...' if len(Hs)>6 else ''):>34} "
          f"{str(sorted(primes)):>16} {len(primes):>6}")
print()
print("  Every H is ODD (Redei), so 2 NEVER appears -- the prime 2 is structurally absent")
print("  from the multiplicative side, which is exactly why the ADDITIVE side (the tiling")
print("  hypercube, F_2^m) is where all the 2-adic structure lives.  The two arithmetics")
print("  partition the primes: 2 is purely additive, the odd primes purely multiplicative.")

print()
print("=" * 78)
print("SUMMARY -- what is a theorem here and what is only an analogy")
print("=" * 78)
print("  THEOREM (classical):  at least one of e+pi, e*pi is transcendental (Vieta).")
print("  THEOREM (canon):      char-poly coefficients ARE the elementary symmetric")
print("                        functions -- Vieta is already how the repo reads spectra.")
print("  THEOREM (this file):  THM-1520's saddle limit is e^rational, hence transcendental")
print("                        except on a_{D-1} = 0, where it is exactly 1.")
print("  THEOREM (this file):  dim_Q span{log H} = #distinct primes in the H-values, by")
print("                        unique factorization.  And 2 is absent by Redei.")
print("  ANALOGY ONLY:         'Schanuel is the repo's two arithmetics at the transcendence")
print("                        level.'  This is a framing.  Nothing here bears on Schanuel,")
print("                        and no repo object is known to be transcendental in a way")
print("                        that is not already Lindemann.")
