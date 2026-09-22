#!/usr/bin/env python3
"""
collatz_mod6_20260917_zsigmondy_triad.py  (session collatz-mod6-20260917, LANE = zsigmondy_triad)

Wildcard lane: the user's proposed parallel between
  (i)   x -> 2x+1 from 0            (terms 2^n-1; Bang/Zsigmondy exception n=6, 63)
  (ii)  x -> x^2-7/4 from 0         (no primitive prime divisor at n=3, f^3(0)=-7/256)
  (iii) x -> x^2-29/16, 3-cycle -7/4 -> 5/4 -> -1/4      (THM-4139 / THM-4146)
  (iv)  the integer systems x^2+0, x^2-1, x^2-2 (the user's "x^2+{0,1,2}")

Everything load-bearing is exact (Python int / Fraction / sympy over Q).
All checks use explicit `raise` (active under python -O).
Optional external cross-check: PARI/GP `thue` (skipped gracefully if gp is absent).

Inheritance (read, not re-derived):
  05-knowledge/results/arithmetic_braids_20260917_collatz.md   (rows 6j+1->9j+2 etc., ord_9(2)=6)
  05-knowledge/results/arithmetic_braids_20260917_summand.md
  05-knowledge/results/arithmetic_braids_20260917_divisors.md
  05-knowledge/results/arithmetic_braids_20260917_geometry.md  (eta-parabola c=-(eta^2+7)/4 (4); parabolic cycle at -7/4 (9);
        n=3 numerator test Q(a,b)=+-1 (11) with census b<=2000 and the Krieger citation; Gaussian squaring x=2A/C -> x^2-2 (16))
  01-canon/theorems/THM-3341-*.md, THM-3333-*.md  (Gaussian squaring of primitive triples; read for the x^2-2 conjugacy only)
  01-canon/theorems/THM-4139-*.md  (PrePer(x^2-29/16,Q), unique AP 3-cycle, Psi_6 of squaring, 63 census)
  01-canon/theorems/THM-4146-*.md  (sigma-parametrization c=-(sigma^2+sigma+2), multiplier, 3:4:5 forcing)
"""
import sys, re, time, shutil, subprocess, itertools
from fractions import Fraction as Fr
from math import gcd, isqrt
import sympy as sp
from sympy import Rational as R

T0 = time.time()
FAILS = []

def check(cond, msg):
    if not cond:
        FAILS.append(msg)
        raise RuntimeError("CHECK FAILED: " + msg)

def hr(title):
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)

def factor_str(n):
    n = abs(n)
    if n == 0:
        return "0"
    if n == 1:
        return "1"
    return " * ".join(f"{p}^{e}" if e > 1 else f"{p}" for p, e in sorted(sp.factorint(n).items()))

GP = shutil.which("gp")

def gp_run(script, timeout=300):
    """Run PARI/GP on `script`; gp -q exits at EOF, so no \\q is needed.  Returns output lines or None."""
    if not GP:
        return None
    try:
        r = subprocess.run([GP, "-q"], input=script + "\n", capture_output=True, text=True, timeout=timeout)
        return [l for l in r.stdout.strip().splitlines() if not l.startswith("  ***")]
    except Exception as e:
        print("  PARI/GP call failed:", e)
        return None

def parse_pairs(line):
    return {(int(a_), int(b_)) for a_, b_ in re.findall(r"\[(-?\d+),\s*(-?\d+)\]", line)}

def normalize_pairs(sols):
    out = set()
    for a_, b_ in sols:
        if b_ < 0 or (b_ == 0 and a_ < 0):
            a_, b_ = -a_, -b_
        out.add((a_, b_))
    return out

def primitive_part(terms, n):
    """terms[1..n] integers. Return R = part of |terms[n]| coprime to all earlier terms,
    obtained by repeated gcd stripping.  R != 1  <=>  terms[n] has a primitive prime divisor
    (a prime dividing terms[n] and no terms[m], m<n).  No factoring needed."""
    Rn = abs(terms[n])
    if Rn == 0:
        return 0
    prod = 1
    for m in range(1, n):
        prod *= abs(terms[m]) if terms[m] != 0 else 1
    while True:
        g = gcd(Rn, prod)
        if g == 1:
            return Rn
        Rn //= g

# ---------------------------------------------------------------------------------------
hr("S0. Definitions, universe, inheritance")
print("f_c(x) = x^2 + c over Q.  G_n(c) = f_c^n(0) (Gleason polynomials).  N_n = numerator of f_c^n(0).")
print("'primitive prime divisor at n' = prime p | N_n with p not dividing N_m for any m<n.")
print("Bang/Zsigmondy universe: 2^n-1, 1<=n<=40 (exact factorization via sympy).")
print("Critical-orbit census universe: c=a/b, b in {1,2,4,8,16,32,64}, |a|<=200, gcd(a,b)=1, 1<=n<=8;")
print("  the integer sub-universe c in [-60,60] is contained in b=1.")
print("Inherited and NOT re-derived: three-row typing, ord_9(2)=6 row exponents, R(n)=4n+1 fibre,")
print("  PrePer(x^2-29/16,Q) (THM-4139 sec.1 / THM-4146 sec.1), sigma-parametrization (THM-4146 (10)),")
print("  Psi_6 of squaring = Phi_9 Phi_21 Phi_63 (THM-4139 (37)),")
print("  arithmetic_braids_20260917_geometry.md: eta-parabola c=-(eta^2+7)/4 and disc (eta^2+eta+7)^2 (4), the parabolic")
print("  cycle P(x)=x^3+x^2/2-9x/4-1/8 with Phi_3=P^2 and multiplier 1 at -7/4 (9), the n=3 test f^3(0)=aQ(a,b)/b^4 with")
print("  Q=+-1 <=> no new prime (11) and its census b<=2000, Krieger's Theorem 6.1 boundary, x=2A/C -> x^2-2 (16).")
print("  Everything below is re-verified only where a check is printed; inherited statements are cited, not claimed.")

# ---------------------------------------------------------------------------------------
hr("S1. (i) 2^n-1, n<=40: primitive prime divisors, Bang/Zsigmondy exceptions")
mers = [0] + [2**n - 1 for n in range(1, 41)]
exceptions = []
print(f"{'n':>3} {'2^n-1':>15}  {'factorization':<34} {'Phi_n(2)':>10}  primitive primes")
nonprim_phi = {}
for n in range(1, 41):
    prim = []
    if mers[n] > 1:
        for p in sp.factorint(mers[n]):
            if all(mers[m] % p != 0 for m in range(1, n)):
                prim.append(p)
    phi_n_2 = int(sp.cyclotomic_poly(n, 2))
    # non-primitive primes of Phi_n(2)
    npp = [p for p in sp.factorint(phi_n_2) if p not in prim]
    nonprim_phi[n] = (phi_n_2, npp)
    if not prim:
        exceptions.append(n)
    print(f"{n:>3} {mers[n]:>15}  {factor_str(mers[n]):<34} {phi_n_2:>10}  {prim if prim else '--- NONE ---'}"
          + (f"   [non-primitive prime(s) in Phi_n(2): {npp}]" if npp else ""))
    # cross-check with gcd-stripping primitive part
    Rn = primitive_part(mers, n)
    check((Rn != 1) == bool(prim), f"primitive-part/gcd cross-check at n={n}")
print("Exceptions (no primitive prime divisor):", exceptions)
check(exceptions == [1, 6], "Bang exceptions for base 2 up to 40 must be exactly n=1,6")
# n>2: any non-primitive prime of Phi_n(2) is the largest prime factor of n and divides Phi_n(2) exactly once
for n in range(3, 41):
    phi, npp = nonprim_phi[n]
    for p in npp:
        check(p == max(sp.factorint(n)), f"non-primitive prime {p} of Phi_{n}(2) is not the largest prime of n")
        check(phi % (p * p) != 0, f"non-primitive prime {p} divides Phi_{n}(2) more than once")
print("FINITE-EXACT (n<=40): every non-primitive prime of Phi_n(2), n>2, is the largest prime of n and divides")
print("  Phi_n(2) exactly once.  (This is the classical lemma behind Bang/Zsigmondy.)")
print("63 = 2^6-1 = Phi_1(2)Phi_2(2)Phi_3(2)Phi_6(2) =", [int(sp.cyclotomic_poly(d, 2)) for d in (1, 2, 3, 6)],
      " -> Phi_6(2) = 3 = largest prime of 6, already in 2^2-1: NO primitive prime.")
print("Catalan reading (PROVED, elementary): 2^6-1 = (2^3-1)(2^3+1) = 7 * 9 and 2^3+1 = 3^2; a primitive prime of")
print("  2^6-1 would have to divide (2^3+1)/gcd(2^3+1,2^2-1) = 9/3 = 3 -> impossible.  The exception IS 3^2-2^3=1.")

# ord_9(2)=6 and the row exponents
o = sp.n_order(2, 9)
check(o == 6, "ord_9(2) must be 6")
rows = {2: [], 5: [], 8: []}
for e in range(0, 60):
    r = pow(2, e, 9)
    if r in rows:
        rows[r].append(e % 6)
check(set(rows[2]) == {1} and set(rows[5]) == {5} and set(rows[8]) == {3}, "row exponent classes mod 6")
print("ord_9(2) = 6 <=> 9 | 63.  Powers 2^e with e<60: 2^e=2 mod 9 iff e=1 mod 6; =5 iff e=5 mod 6; =8 iff e=3 mod 6")
print("  (FINITE-EXACT e<60; PROVED for all e since ord_9(2)=6 and 2^1=2, 2^5=32=5, 2^3=8 mod 9).")
# No prime p has ord_p(2)=6 (PROVED): such p would be a primitive prime divisor of 2^6-1.
for p in sp.primerange(3, 10**5):
    if sp.n_order(2, p) == 6:
        check(False, f"prime {p} with ord_p(2)=6 found")
print("PROVED: no prime p has ord_p(2)=6 (a prime with ord_p(2)=n is exactly a primitive prime divisor of 2^n-1);")
print("  checked hostile-style for all p<10^5.  Hence 9 is the least modulus of multiplicative order 6 for 2,")
print("  which is why the inherited three-row exponent law lives modulo 9 with period 6, not modulo a prime.")
# Zsigmondy for 2^n+1 (relevant later: 2^3+1=9)
exc_plus = []
plus = [0] + [2**n + 1 for n in range(1, 41)]
for n in range(1, 41):
    if primitive_part(plus, n) == 1:
        exc_plus.append(n)
print("2^n+1, n<=40, terms with no primitive prime divisor:", exc_plus)
check(exc_plus == [3], "2^n+1 exception must be exactly n=3 (2^3+1=9)")

# ---------------------------------------------------------------------------------------
hr("S2. (ii) zero orbit of x^2-7/4 to n=8, and the bounded census")

def zero_orbit_numerators(c, nmax):
    x = Fr(0)
    nums = [0]
    orb = [x]
    for n in range(1, nmax + 1):
        x = x * x + c
        nums.append(x.numerator)
        orb.append(x)
    return nums, orb

c74 = Fr(-7, 4)
nums, orb = zero_orbit_numerators(c74, 8)
print("n : f^n(0)  [numerator factorization ; primitive part]")
for n in range(1, 9):
    Rn = primitive_part(nums, n)
    if n <= 6:
        fs = factor_str(nums[n])
    else:
        # partial factorization (trial division) for display only
        m = abs(nums[n]); small = []
        for p in sp.primerange(2, 10**5):
            while m % p == 0:
                small.append(p); m //= p
        fs = " * ".join(map(str, small)) + (f" * [cofactor {len(str(m))} digits, {'prime' if sp.isprime(m) else 'composite/unfactored'}]" if m > 1 else "")
    print(f"{n} : {orb[n] if n <= 4 else str(orb[n].numerator)[:12] + '.../' + '4^' + str(2 ** (n - 1))}"
          f"   num = {'-' if nums[n] < 0 else ''}{fs} ; primitive part {Rn if len(str(Rn)) < 40 else str(Rn)[:20] + '...(' + str(len(str(Rn))) + ' digits)'}")
check(orb[3] == Fr(-7, 256), "f^3(0) at c=-7/4 must be -7/256")
check(primitive_part(nums, 3) == 1, "n=3 must have no primitive prime divisor at c=-7/4")
for n in [1, 2, 4, 5, 6, 7, 8]:
    check(primitive_part(nums, n) != 1, f"n={n} should have a primitive prime divisor at c=-7/4")
print("f^3(0)/f(0) = ", orb[3] / orb[1], " = 2^-6 = 1/(63+1)   [the '63' the user sees in the -7/4 orbit]")
check(orb[3] / orb[1] == Fr(1, 64), "f^3(0)/f(0) = 1/64")

# census
def preperiodic_zero(c, nmax=64):
    """For integer c only (rational c with b>1 can never have 0 preperiodic: v_p(f^n(0)) = 2^(n-1) v_p(c) -> -inf)."""
    seen = {}
    x = Fr(0); k = 0
    while x not in seen and k <= nmax:
        seen[x] = k
        x = x * x + c; k += 1
        if abs(x) > abs(c) + 2:
            return None
    return (seen[x], k - seen[x]) if x in seen else None

box_b = [1, 2, 4, 8, 16, 32, 64]
NMAX = 8
failures = {n: [] for n in range(1, NMAX + 1)}
preper = []
count = 0
for b in box_b:
    for a in range(-200, 201):
        if gcd(a, b) != 1 and not (b == 1):
            continue
        if b == 1 and a == 0:
            pass
        c = Fr(a, b)
        count += 1
        if b == 1:
            pp = preperiodic_zero(c)
            if pp is not None:
                preper.append((c, pp))
                continue
        nums, orb = zero_orbit_numerators(c, NMAX)
        for n in range(1, NMAX + 1):
            if primitive_part(nums, n) == 1:
                failures[n].append(c)
print(f"Census universe: {count} parameters c=a/b (b in {box_b}, |a|<=200, gcd=1), orbit to n={NMAX}.")
print("Parameters with 0 preperiodic (reported separately):", [(str(c), f"tail {pp[0]}, period {pp[1]}") for c, pp in preper])
check([c for c, _ in preper] == [Fr(-2), Fr(-1), Fr(0)], "0 preperiodic exactly for c in {-2,-1,0}")
for n in range(1, NMAX + 1):
    L = failures[n]
    print(f"n={n}: {len(L):>3} non-preperiodic parameters with NO primitive prime divisor: "
          + (", ".join(str(c) for c in L) if len(L) <= 40 else ", ".join(str(c) for c in L[:40]) + " ..."))
# exact predictions
pred1 = sorted({Fr(s, b) for b in box_b for s in (1, -1)} - {Fr(-1)})   # c=-1 is preperiodic, reported separately
pred2 = sorted({Fr(-b + s, b) for b in box_b for s in (1, -1)} - {Fr(0), Fr(-2)})
check(sorted(failures[1]) == pred1, "n=1 failures must be exactly c = +-1/b")
check(sorted(failures[2]) == pred2, "n=2 failures must be exactly c = -1 +- 1/b (non-preperiodic ones)")
check(failures[3] == [Fr(-7, 4)], "n=3: the only non-preperiodic failure in the box must be -7/4")
for n in range(4, NMAX + 1):
    check(failures[n] == [], f"n={n}: unexpected failure {failures[n]}")
print("PROVED (all c=a/b, b>=1, gcd(a,b)=1):")
print("  n=1: N_1=a, so no primitive prime iff |a|=1, i.e. c=+-1/b (includes the cusp c=1/4).")
print("  n=2: N_2=a(a+b); primes of a+b are coprime to a, so no primitive prime iff |a+b|=1, i.e. c=-1+-1/b")
print("       (includes the period-doubling parameter c=-3/4=-1+1/4).  INFINITELY many n=2 failures.")
print("  n=3: N_3=a*F(a,b), F(a,b)=a^3+2a^2b+ab^2+b^3, gcd(F,a)=gcd(b^3,a)=1, gcd(F,a+b)=gcd(b^3,a+b)=1,")
print("       so every prime of F is primitive and failure <=> F(a,b)=+-1  (a cubic THUE equation; see S3).")
print("       [n=3 is inherited: geometry note (11) proves exactly this with F=Q and censuses b<=2000; new here: S3/S3b.]")
print("FINITE-EXACT (this box, n<=8): -7/4 is the ONLY non-preperiodic n=3 failure, and there are NO failures for 4<=n<=8.")

# ---------------------------------------------------------------------------------------
hr("S3. Mechanism at c=-7/4: dynatomic discriminant, parabolic 3-cycle, sigma-parabola, Thue/plastic units")
x, c, s, sig, t, X = sp.symbols('x c s sigma t X')
f = lambda z: z**2 + c
f1 = f(x); f2 = sp.expand(f(f1)); f3 = sp.expand(f(f2)); f4 = sp.expand(f(f3))
Phi3 = sp.expand(sp.cancel((f3 - x) / (f1 - x)))
print("Phi_3(x,c) =", Phi3)
disc3 = sp.factor(sp.discriminant(Phi3, x))
print("disc_x Phi_3 =", disc3)
check(sp.expand(disc3 + (4*c + 7)**3 * (16*c**2 + 4*c + 7)**2) == 0, "disc_x Phi_3 = -(4c+7)^3 (16c^2+4c+7)^2")
res3 = sp.factor(sp.resultant(Phi3, sp.diff(f3, x) - 1, x))
print("Res_x(Phi_3, (f^3)'(x)-1) =", res3)
check(sp.expand(res3 - (4*c + 7)**3 * (16*c**2 + 4*c + 7)**3) == 0, "multiplier-1 locus = (4c+7)(16c^2+4c+7)")
print("PROVED: Phi_3(x,c) has a repeated root in x  <=>  (4c+7)(16c^2+4c+7)=0; the real branch is exactly 4c+7=0.")
print("  Delta_3(c) := (4c+7)(16c^2+4c+7) = 64c^3+128c^2+56c+49 is the period-3 multiplier-1 (parabolic) locus.")
# at c=-7/4 Phi_3 is a perfect square of the parabolic cycle cubic
Pparab = 8*x**3 + 4*x**2 - 18*x - 1
Phi3_74 = sp.expand(Phi3.subs(c, R(-7, 4)))
check(sp.expand(Phi3_74 * 64 - Pparab**2) == 0, "Phi_3(x,-7/4) = (8x^3+4x^2-18x-1)^2 / 64")
rts = sp.Poly(Pparab, x).all_roots()
cf = sp.Poly(Pparab, x).all_coeffs()          # Vieta: product of the three roots = -cf[3]/cf[0]
prod_roots = -R(cf[3], 1) / cf[0]
check(prod_roots == R(1, 8), "product of parabolic cycle points = 1/8")
print("At c=-7/4: Phi_3 = (8x^3+4x^2-18x-1)^2/64 (a perfect square: two 3-cycles collide).  Multiplier of the cycle")
print("  = prod f'(x_i) = 8 * prod x_i = 8 * (1/8) = 1  -> PARABOLIC (saddle-node) 3-cycle.  PROVED (inherited: geometry (9)).")
print("  Parabolic cycle numerically:", [sp.N(r, 8) for r in rts], " (note x_1 ~ -1.75 ~ c and x_3 ~ -0.055 ~ 0)")
# sigma parabola (THM-4146 (10)) re-verified and re-parametrized by s = sigma + 1/2
Psig = X**3 - sig*X**2 - (sig**2 + 2*sig + 3)*X + (sig**3 + 2*sig**2 + 3*sig + 1)
csig = -(sig**2 + sig + 2)
check(sp.expand(Phi3.subs({x: X, c: csig}) - sp.expand(Psig * Psig.subs(sig, -1 - sig))) == 0,
      "Phi_3(X, c(sigma)) = P_sigma(X) P_{-1-sigma}(X)")
lam_sig = sp.expand(-8 * (sig**3 + 2*sig**2 + 3*sig + 1))
lam_s = sp.expand(lam_sig.subs(sig, s - R(1, 2)))
c_s = sp.expand(csig.subs(sig, s - R(1, 2)))
print("sigma-parametrization (THM-4146 (10), re-verified): c = -(sigma^2+sigma+2), multiplier lambda = 8 alpha beta gamma =", lam_sig)
print("  Put s = sigma + 1/2.  Then  c = ", c_s, "   and   lambda(s) =", lam_s)
check(c_s == -s**2 - R(7, 4), "c = -7/4 - s^2")
check(lam_s.subs(s, 0) == 1 and lam_s.subs(s, R(-1, 4)) == R(35, 8) and lam_s.subs(s, R(1, 4)) == R(-23, 8),
      "lambda(0)=1, lambda(-1/4)=35/8, lambda(1/4)=-23/8")
print("PROVED (inherited: geometry note (4) states this as c=-(eta^2+7)/4 with eta=2s; re-verified here):")
print("  the unmarked 3-cycle curve of x^2+c is the parabola c = -7/4 - s^2, a 2:1 cover of the c-line branched")
print("  exactly at the real parabolic parameter c=-7/4 (s=0, lambda=1).  A rational unmarked 3-cycle (rational cycle-sum)")
print("  exists iff -7/4 - c is a rational square.  The cycle points are rational iff Morton's t (S4) is rational.")
print("  -29/16 = -7/4 - (1/4)^2: its two 3-cycles are s=-1/4 (sigma=-3/4, the rational AP cycle, lambda=35/8 = THM-4139 (32))")
print("  and s=+1/4 (sigma=-1/4, lambda=-23/8, irrational).   -2 = -7/4 - (1/2)^2: s=+-1/2, sigma in {0,-1} (see S4, S6).")
disc_s = sp.factor(sp.discriminant(Psig.subs(sig, s - R(1, 2)), X))
print("disc_X P_s =", disc_s, "  => splitting field of any rational-s 3-cycle is Q or a cyclic cubic field (PROVED).")
check(sp.expand(disc_s - (4*s**2 + 2*s + 7)**2) == 0, "disc P_s = (4s^2+2s+7)^2")

# Gleason polynomials and the identity behind the 1/64
G = [sp.Integer(0)]
for n in range(1, 6):
    G.append(sp.expand(G[-1]**2 + c))
print("Gleason polynomials G_n(c)=f_c^n(0), factored:")
for n in range(1, 6):
    print(f"  G_{n} =", sp.factor(G[n]))
H3 = c**3 + 2*c**2 + c + 1
check(sp.expand(G[3] - c * H3) == 0, "G_3 = c H_3")
id3 = sp.factor(sp.expand(64 * G[3] - c))
print("64 G_3(c) - c =", id3)
check(sp.expand(id3 - c * (4*c + 7) * (16*c**2 + 4*c + 9)) == 0, "64G_3 - c = c(4c+7)(16c^2+4c+9)")
print("=> at c=-7/4: H_3(-7/4) = 1/64 exactly, f^3(0) = c/64, N_3 = N_1 = -7.  The parabolic factor (4c+7) divides")
print("   64G_3-c, but the complex parabolic factor (16c^2+4c+7) does NOT (cofactor is 16c^2+4c+9).")
print("Pattern 'N_n = N_1 at the real period-n parabolic parameter', i.e. Delta_n^{real} | 2^(2^n-2) G_n(c) - c :")
for n in range(1, 5):
    Pn = sp.expand(2**(2**n - 2) * G[n] - c)
    print(f"  n={n}: 2^({2**n - 2}) G_{n} - c =", sp.factor(Pn) if Pn != 0 else 0)
Phi4 = sp.expand(sp.cancel((f4 - x) / (f2 - x)))
res4 = sp.factor(sp.resultant(Phi4, sp.diff(f4, x) - 1, x))
print("Res_x(Phi_4,(f^4)'-1) =", res4)
D4prim = 64*c**3 + 144*c**2 + 108*c + 135
check(sp.expand(res4 - (4*c + 5)**4 * (16*c**2 - 8*c + 5)**4 * D4prim**4) == 0, "Delta_4 factorization")
g4 = sp.gcd(sp.expand(2**14 * G[4] - c), D4prim * (4*c + 5) * (16*c**2 - 8*c + 5))
check(g4 == 1, "gcd(2^14 G_4 - c, Delta_4) must be 1")
realroot4 = [r for r in sp.Poly(D4prim, c).all_roots() if r.is_real][0]
print(f"  n=4: real period-4 saddle-node c_4 = root of {D4prim} ~ {sp.N(realroot4, 12)};  gcd(2^14 G_4 - c, Delta_4) = {g4}")
print("  REFUTED beyond n=3: the pattern holds for n=1,2,3 (c=1/4: N_1=1; c=-3/4: f^2(0)=c/4; c=-7/4: f^3(0)=c/64) and")
print("  FAILS at n=4 (exact gcd = 1).  It is a low-degree coincidence, not a theorem about parabolic parameters.")
print("  Also 2^n-2 = 2n has the unique solution n=3, which is the only reason the exponent 6 in f^3(0)/f(0)=2^-6")
print("  equals Bang's index 6.  (2^n-2 = 0,2,6,14,... vs 2n = 2,4,6,8,...)")

# n=3 characterization: Thue equation and plastic-number units
print()
print("n=3 non-Zsigmondy parameters <=> F(a,b) = a^3+2a^2b+ab^2+b^3 = +-1 (Thue).  F(a,b)=N(a - b theta)=N(a + b rho^2)")
print("  where theta^3+2theta^2+theta+1=0, theta=-rho^2, rho = plastic number (rho^3=rho+1), field disc -23.")
rho_poly = t**3 - t - 1
theta_min = sp.expand(sp.minimal_polynomial(-sp.CRootOf(rho_poly, 0)**2, t))
check(theta_min == t**3 + 2*t**2 + t + 1, "-rho^2 has minimal polynomial t^3+2t^2+t+1 (the Gleason H_3)")
print("  PROVED: the airplane centre (real root of H_3) is c = -rho^2 =", sp.N(-sp.CRootOf(rho_poly, 0)**2, 12))
print("  and -7/4 = -1.75 is within 0.0049 of it.")
# rho^k = A_k + B_k rho + C_k rho^2, find k with B_k = 0 (units with vanishing middle coefficient)
def rho_pow(k):
    # represent in basis 1, rho, rho^2; rho^3 = rho + 1 ; rho^-1 = rho^2 - 1
    v = (1, 0, 0)
    mul = (0, 1, 0) if k >= 0 else (-1, 0, 1)
    for _ in range(abs(k)):
        A, B, C = v
        a0, a1, a2 = mul
        # (A + B r + C r^2)(a0 + a1 r + a2 r^2), reduce r^3 = r+1, r^4 = r^2 + r
        c0 = A*a0; c1 = A*a1 + B*a0; c2 = A*a2 + B*a1 + C*a0; c3 = B*a2 + C*a1; c4 = C*a2
        v = (c0 + c3, c1 + c3 + c4, c2 + c4)
    return v
unit_hits = []
for k in range(-300, 301):
    A, B, C = rho_pow(k)
    if B == 0:
        unit_hits.append((k, A, C))
print("  Units +-rho^k = A + C rho^2 with vanishing rho-coefficient, |k|<=300:", unit_hits)
check([(k, A, C) for k, A, C in unit_hits] == [(-14, -7, 4), (-5, 2, -1), (-1, -1, 1), (0, 1, 0), (2, 0, 1)],
      "vanishing-middle-coefficient units")
print("  -> c = A/C in {-7/4 (k=-14), -2 (k=-5), -1 (k=-1), 0 (k=2)} and b=0 (k=0: c=infinity, not a parameter).")
print("     So -7/4 = ratio in rho^-14 = 4 rho^2 - 7.  (Z[rho] is the full ring of integers since disc -23 is squarefree.)")
# direct Thue search
rho_num = sp.N(sp.CRootOf(rho_poly, 0), 30)
r2 = float(rho_num**2)
thue_sols = []
for b in range(1, 10**6 + 1):
    a0 = round(-r2 * b)
    for a in (a0 - 2, a0 - 1, a0, a0 + 1, a0 + 2):
        if a*a*a + 2*a*a*b + a*b*b + b*b*b in (1, -1) and gcd(a, b) == 1:
            thue_sols.append((a, b))
print("  Direct search F(a,b)=+-1, 1<=b<=10^6, a within 2 of -rho^2 b (all others have |F|>1 since |a/b+rho^2| >= 1/(2b) there):", thue_sols)
check(thue_sols == [(-2, 1), (-1, 1), (0, 1), (-7, 4)], "Thue solutions with b<=10^6")
out = gp_run("tnf=thueinit(x^3+2*x^2+x+1,1);print(thue(tnf,1));print(thue(tnf,-1));")
PARI_OK = False
if out is not None and len(out) >= 2:
    print("  PARI/GP thueinit(flag=1, unconditional) thue(F,+1):", out[0], "  thue(F,-1):", out[1])
    plus, minus = parse_pairs(out[0]), parse_pairs(out[1])
    check(minus == {(-a_, -b_) for a_, b_ in plus}, "thue(F,-1) = -thue(F,+1) (F has odd degree)")
    check(normalize_pairs(plus | minus) == {(-7, 4), (-1, 1), (0, 1), (1, 0), (-2, 1)}, "PARI Thue solution set")
    print("  => COMPLETE (CITED: PARI thue, Bilu-Hanrot algorithm, flag 1 = no GRH): all solutions of F(a,b)=+-1 are")
    print("     +-(1,0), +-(0,1), +-(-1,1), +-(-2,1), +-(-7,4).  Hence over ALL of Q the n=3 non-Zsigmondy parameters are")
    print("     exactly c in {0,-1,-2} (0 preperiodic) and c = -7/4.")
    PARI_OK = True
else:
    print("  PARI/GP not available or failed; completeness of the Thue list is then OPEN here (FINITE-EXACT to b<=10^6).")

# ---------------------------------------------------------------------------------------
hr("S3b. Gleason product G_n = prod_{d|n} H_d, pairwise resultants, and the n<=6 primitive-divisor reduction")
Gl = {n: G[n] for n in range(1, 6)}
Gl[6] = sp.expand(G[5]**2 + c)
Hd = {}
for n in range(1, 7):
    q = Gl[n]
    for d in sp.divisors(n):
        if d < n:
            q = sp.quo(q, Hd[d], c)
    Hd[n] = sp.expand(q)
    prod_ = sp.Integer(1)
    for d in sp.divisors(n):
        prod_ *= Hd[d]
    check(sp.expand(prod_ - Gl[n]) == 0, f"G_{n} = prod_(d|{n}) H_d exactly (division is exact)")
    check(sp.Poly(Hd[n], c).LC() == 1, f"H_{n} monic")
    check(sp.gcd(Hd[n], sp.diff(Hd[n], c)) == 1, f"H_{n} squarefree")
    print(f"  H_{n}: degree {int(sp.degree(Hd[n], c)):>2}, irreducible over Q: {sp.Poly(Hd[n], c).is_irreducible}, "
          f"H_{n}(0) = {Hd[n].subs(c, 0)}, H_{n}(-1) = {Hd[n].subs(c, -1)}, H_{n}(-2) = {Hd[n].subs(c, -2)}")
check(sp.expand(Hd[3] - H3) == 0, "H_3 = c^3+2c^2+c+1")
check(sp.expand(Hd[4] - (c**6 + 3*c**5 + 3*c**4 + 3*c**3 + 2*c**2 + 1)) == 0, "H_4 as displayed in S3")
def sylvester_resultant(f_, g_, var):
    """Res(f,g) as the Sylvester determinant, i.e. the standard convention Res(f,g) = lc(f)^deg g * prod_{f(a)=0} g(a)
    (PARI polresultant agrees).  SymPy's resultant() (subresultant PRS) returns the opposite sign for some pairs below
    (audit of 2026-09-21: e.g. Res(H_1,H_3) = H_3(0) = +1, where sympy.resultant gives -1); only |Res| = 1 is load-bearing."""
    pf, pg = sp.Poly(f_, var), sp.Poly(g_, var)
    m_, n_ = pf.degree(), pg.degree()
    Mx = sp.zeros(m_ + n_, m_ + n_)
    for i in range(n_):
        for j, v in enumerate(pf.all_coeffs()):
            Mx[i, i + j] = v
    for i in range(m_):
        for j, v in enumerate(pg.all_coeffs()):
            Mx[n_ + i, i + j] = v
    return Mx.det()
print("  Res_c(H_d, H_n) for 1 <= d < n <= 6 (exact Sylvester determinants, standard sign convention):")
for n in range(2, 7):
    row = []
    for d in range(1, n):
        r_ = sylvester_resultant(Hd[d], Hd[n], c)
        check(r_ in (1, -1), f"Res(H_{d},H_{n}) must be +-1, got {r_}")
        check(abs(sp.resultant(Hd[d], Hd[n], c)) == 1, f"|sympy PRS resultant| = 1 for (d,n)=({d},{n})")
        row.append(f"Res(H_{d},H_{n})={int(r_):>2}")
    print("    " + "; ".join(row))
check(sylvester_resultant(Hd[1], Hd[3], c) == Hd[3].subs(c, 0) == 1, "Res(H_1,H_3) = H_3(0) = +1 (sign sanity of the convention)")
print("PROVED (n<=6, from the exact resultants above).  Write H_n(a,b) = b^deg H_n(a/b) (monic in a).  For c=a/b in lowest terms:")
print("  N_n = prod_(d|n) H_d(a,b) exactly (total degree 2^(n-1); each H_d(a,b) = a^deg mod b is coprime to b).")
print("  If a prime p divided H_d(a,b) and H_n(a,b) with d<n, then p does not divide b, and a/b mod p would be a common root")
print("  of H_d and H_n mod p, so p | Res(H_d,H_n) = +-1: impossible.  Hence the factors H_d(a,b), d|n, are pairwise coprime,")
print("  every prime of H_n(a,b) is primitive at n, and every other prime of N_n divides some N_d, d|n, d<n.  Therefore")
print("  N_n has a primitive prime divisor  <=>  |H_n(a,b)| >= 2,   i.e.  non-Zsigmondy at n  <=>  H_n(a,b) = +-1  (n<=6),")
print("  a Thue equation of degree 1,1,3,6,15,27 for n=1..6.  (This is the quadratic-critical-orbit analogue of")
print("  'non-primitive primes of Phi_n(2) divide n'; for n<=6 it is even cleaner: there are none at all.)")
# small-height searches for H_4, H_5, H_6 = +-1
def hom_coeffs(Pc):
    return [int(v) for v in sp.Poly(Pc, c).all_coeffs()]
def hom_eval(coeffs, a_, bpow):
    r_ = 0
    for i, co in enumerate(coeffs):
        r_ = r_ * a_ + co * bpow[i]
    return r_
BFULL, BROOT = 100, 5000
for n in (4, 5, 6):
    co = hom_coeffs(Hd[n]); deg = len(co) - 1
    rr_ = [float(r_) for r_ in sp.Poly(Hd[n], c).all_roots() if r_.is_real]
    hits = set()
    for b in range(1, BROOT + 1):
        bpow = [b**i for i in range(deg + 1)]
        if b <= BFULL:
            cand = range(-3 * b, 3 * b + 1)
        else:
            cand = {a_ for r_ in rr_ for a_ in range(round(r_ * b) - 2, round(r_ * b) + 3)}
        for a_ in cand:
            if gcd(a_, b) == 1 and abs(hom_eval(co, a_, bpow)) == 1:
                hits.add((a_, b))
    hits = sorted(hits)
    print(f"  H_{n}(a,b) = +-1, all |a|<=3b for b<={BFULL} and a within 2 of b*(real roots {[round(v, 4) for v in rr_]}) for b<={BROOT}: {hits}")
    check(all(b == 1 and a_ in (0, -1, -2) for a_, b in hits), f"H_{n}=+-1 has no non-preperiodic solution in the searched box")
print("  FINITE-EXACT: in the searched boxes the only solutions are c in {0,-1,-2} (0 preperiodic): no non-preperiodic")
print("  exception at n=4,5,6 of that height (consistent with the S2 census box).")
THUE_CLOSED = []
for n in (4, 5):
    if not PARI_OK:
        print(f"  PARI/GP unavailable: n={n} completeness OPEN here (FINITE-EXACT boxes only).")
        continue
    sn = str(Hd[n]).replace("**", "^").replace("c", "x")
    t_ = time.time()
    out = gp_run(f"tnf=thueinit({sn},1);print(thue(tnf,1));print(thue(tnf,-1));", timeout=300)
    if out is not None and len(out) >= 2:
        print(f"  PARI/GP thueinit(flag=1, unconditional) for H_{n} (degree {int(sp.degree(Hd[n], c))}): thue(+1): {out[0]}  thue(-1): {out[1]}")
        plus, minus = parse_pairs(out[0]), parse_pairs(out[1])
        if n % 2:
            check(minus == {(-a_, -b_) for a_, b_ in plus}, f"thue(H_{n},-1) = -thue(H_{n},+1) (odd degree)")
        else:
            check(minus == set(), f"H_{n}=-1 has no solution")
        check(normalize_pairs(plus | minus) == {(1, 0), (0, 1), (-1, 1), (-2, 1)}, f"PARI Thue H_{n} = +-1 solution set")
        THUE_CLOSED.append(n)
        print(f"  => H_{n}(a,b)=+-1 only at +-(1,0),(0,1),(-1,1),(-2,1)  [b=0 is c=infinity; the rest are c=0,-1,-2, all preperiodic].")
    else:
        print(f"  PARI/GP H_{n} Thue call failed: n={n} completeness OPEN here (FINITE-EXACT boxes only).")
if THUE_CLOSED == [4, 5]:
    print("  PROVED + CITED (resultants above + PARI thue, unconditional): for EVERY rational c with infinite critical orbit,")
    print("  N_4 and N_5 have primitive prime divisors.  Together with n=3 (S3): the only non-preperiodic rational parameter")
    print("  whose critical orbit misses a primitive prime at some 3<=n<=5 is c=-7/4, at n=3 only.  The draft's n=4 OPEN item")
    print("  is closed.  n=6 (degree 27) is left at the finite box above: a side run of thueinit on H_6 did not finish within")
    print("  the time budget of this lane, so H_6(a,b)=+-1 is OPEN beyond that box.")

# ---------------------------------------------------------------------------------------
hr("S4. Families: fixed points, 2-cycles, Morton 3-cycles; cycle fields; the Chebyshev conductors 7 and 9")
r_, s_ = sp.symbols('r s_')
print("Fixed points x^2-x+c=0 rational iff 1-4c=r^2, c=(1-r^2)/4: r=1->0, r=3->-2, r=5/2->-21/16, r=9/2->-77/16 (PROVED).")
for rr, cc in [(1, 0), (3, -2), (R(5, 2), R(-21, 16)), (R(9, 2), R(-77, 16))]:
    check((1 - R(rr)**2) / 4 == cc, f"fixed-point family r={rr}")
print("2-cycles x^2+x+c+1=0 rational iff -3-4c=s^2, c=-(3+s^2)/4: s=1->-1, s=2->-7/4, s=0->-3/4 (degenerate, multiplier -1),")
print("  s=3->-3, s=5/2->-37/16, s=9/2->-93/16 (PROVED).")
for ss, cc in [(1, -1), (2, R(-7, 4)), (0, R(-3, 4)), (3, -3), (R(5, 2), R(-37, 16)), (R(9, 2), R(-93, 16))]:
    check(-(3 + R(ss)**2) / 4 == cc, f"2-cycle family s={ss}")
Dt = 2*t*(t + 1)
p0 = (t**3 + 2*t**2 + t + 1) / Dt; p1 = (t**3 - t - 1) / Dt; p2 = -(t**3 + 2*t**2 + 3*t + 1) / Dt
ct = -(t**6 + 2*t**5 + 4*t**4 + 8*t**3 + 9*t**2 + 4*t + 1) / (4*t**2*(t + 1)**2)
for a_, b_ in [(p0, p1), (p1, p2), (p2, p0)]:
    check(sp.simplify(a_**2 + ct - b_) == 0, "Morton cycle p0->p1->p2->p0")
check(ct.subs(t, 1) == R(-29, 16), "c(1) = -29/16")
check([p0.subs(t, 1), p1.subs(t, 1), p2.subs(t, 1)] == [R(5, 4), R(-1, 4), R(-7, 4)], "t=1 points")
sig_t = sp.factor(sp.cancel(p0 + p1 + p2)); s_t = sp.factor(sp.cancel(sig_t + R(1, 2)))
print("Morton-type family (THM-4139 (13), re-verified symbolically): c(t) = -(t^6+2t^5+4t^4+8t^3+9t^2+4t+1)/(4t^2(t+1)^2),")
print("  t=1 -> c=-29/16 with points (5/4,-1/4,-7/4).   sigma(t) =", sig_t, ";   s(t) = sigma+1/2 =", s_t)
check(sp.simplify(ct + R(7, 4) + s_t**2) == 0, "c(t) = -7/4 - s(t)^2")
check(sp.simplify(sig_t.subs(t, -1/(t + 1)) - sig_t) == 0, "sigma invariant under rho(t)=-1/(t+1)")
check(sp.simplify(sp.cancel((p1 - p2) / (p0 - p1)) - t) == 0, "t = (p1-p2)/(p0-p1)")
print("  PROVED: c(t) = -7/4 - s(t)^2 (Morton's t lies over the sigma-parabola), sigma is invariant under the relabeling")
print("  rho(t) = -1/(t+1), and t = (p1-p2)/(p0-p1) is the ratio of consecutive cycle differences (t=1 <=> AP <=> -29/16).")
num7 = sp.factor(sp.numer(sp.cancel(s_t + R(1, 2)))); num9 = sp.factor(sp.numer(sp.cancel(s_t - R(1, 2))))
print("  c(t) = -2  <=>  s = +-1/2:   s=-1/2 <=> ", num7, "= 0 (roots 1/(2cos(2 pi k/7)));   s=+1/2 <=> ", num9, "= 0 (roots -2cos(2 pi k/9)).")
check(num7 == t**3 + 2*t**2 - t - 1 and num9 == t**3 - 3*t - 1, "conductor-7/9 numerators")
Phi3_m2 = sp.factor(Phi3.subs(c, -2))
print("  Phi_3(x,-2) =", Phi3_m2, "  = (min poly of 2cos(2pi/9)) * (min poly of 2cos(2pi/7)).")
check(sp.expand(Phi3.subs(c, -2) - (x**3 - 3*x + 1) * (x**3 + x**2 - 2*x - 1)) == 0, "Phi_3(x,-2) factorization")
# verify the reciprocal / negation claims exactly
for cub, rel in [(t**3 + 2*t**2 - t - 1, "reciprocal"), (t**3 - 3*t - 1, "negation")]:
    q = sp.Poly(cub, t)
    m7 = sp.Poly(t**3 + t**2 - 2*t - 1, t); m9 = sp.Poly(t**3 - 3*t + 1, t)
    if rel == "reciprocal":
        check(sp.expand(sp.Poly(sp.expand(t**3 * m7.as_expr().subs(t, 1/t)), t).as_expr() + cub) == 0 or
              sp.expand(sp.Poly(sp.expand(t**3 * m7.as_expr().subs(t, 1/t)), t).as_expr() - cub) == 0, "reciprocal of conductor-7 cubic")
    else:
        check(sp.expand(m9.as_expr().subs(t, -t) + cub) == 0 or sp.expand(m9.as_expr().subs(t, -t) - cub) == 0, "negation of conductor-9 cubic")

# cycle fields at s in {0, -1/4, +1/4, -1/2, +1/2}: exact embeddings
def exact_embed(Ppoly, mpoly, var):
    """Find rational a0,a1,a2 with P(a0+a1 y+a2 y^2)=0 mod m(y), exactly verified. Returns coefficients or None."""
    import mpmath as mp
    mp.mp.dps = 60
    Pr = [mp.mpf(str(v)) for v in sp.Poly(Ppoly, var).all_coeffs()]
    Mr = [mp.mpf(str(v)) for v in sp.Poly(mpoly, var).all_coeffs()]
    Xr = sorted([mp.re(z) for z in mp.polyroots(Pr, maxsteps=200, extraprec=200)])
    Yr = sorted([mp.re(z) for z in mp.polyroots(Mr, maxsteps=200, extraprec=200)])
    for perm in itertools.permutations(range(3)):
        A = mp.matrix([[1, Yr[i], Yr[i]**2] for i in range(3)])
        bvec = mp.matrix([Xr[perm[i]] for i in range(3)])
        sol = mp.lu_solve(A, bvec)
        coef = [Fr(str(mp.nstr(sol[i], 40))).limit_denominator(10**6) for i in range(3)]
        expr = sum(sp.Rational(cf.numerator, cf.denominator) * var**i for i, cf in enumerate(coef))
        remv = sp.rem(sp.expand(Ppoly.subs(var, expr)), mpoly, var)
        if remv == 0:
            return coef
    return None

fields = [
    (0, "parabolic cycle at c=-7/4", 8*X**3 + 4*X**2 - 18*X - 1, X**3 + X**2 - 2*X - 1, "Q(2cos 2pi/7), conductor 7"),
    (R(-1, 4), "rational AP cycle at c=-29/16", sp.expand(64 * Psig.subs(sig, R(-3, 4))), None, "Q (splits)"),
    (R(1, 4), "second cycle at c=-29/16", sp.expand(64 * Psig.subs(sig, R(-1, 4))), X**3 - X**2 - 10*X + 8, "cyclic cubic of conductor 31"),
    (R(-1, 2), "cycle at c=-2, sigma=-1", sp.expand(Psig.subs(sig, -1)), X**3 + X**2 - 2*X - 1, "Q(2cos 2pi/7), conductor 7"),
    (R(1, 2), "cycle at c=-2, sigma=0", sp.expand(Psig.subs(sig, 0)), X**3 - 3*X + 1, "Q(2cos 2pi/9), conductor 9"),
]
print("Cycle fields on the parabola c=-7/4-s^2 (disc P_s = (4s^2+2s+7)^2):")
for sv, name, P, m, fname in fields:
    P = sp.Poly(P, X)
    dq = 4*sv**2 + 2*sv + 7
    fl = sp.factor_list(P.as_expr())
    if m is None:
        check(len(fl[1]) == 3, "AP cycle cubic splits into three linear factors over Q")
        print(f"  s={sv}: {name}: P = {P.as_expr()} splits over Q; 4s^2+2s+7 = {dq}.  Field: {fname}.")
    else:
        coef = exact_embed(P.as_expr(), m, X)
        check(coef is not None, f"exact embedding of {name} into Q[X]/({m})")
        check(sp.Poly(P.as_expr(), X).is_irreducible, f"{name} cubic irreducible over Q")
        print(f"  s={sv}: {name}: P = {P.as_expr()} irreducible; 4s^2+2s+7 = {dq}.  Root = {coef[0]} + ({coef[1]}) y + ({coef[2]}) y^2")
        print(f"        with y a root of {m}  [verified exactly mod m].  Field: {fname}.")
out = gp_run("print(nfdisc(x^3+x^2-2*x-1));print(nfdisc(x^3-x^2-10*x+8));print(nfdisc(x^3-3*x+1));"
             "print(nfdisc(8*x^3+4*x^2-18*x-1));print(nfdisc(64*x^3+16*x^2-164*x+23));")
if out is not None and len(out) >= 5:
    nfd = [int(v) for v in out[:5]]
    check(nfd == [49, 961, 81, 49, 961], f"PARI nfdisc cross-check, got {nfd}")
    print("  PARI nfdisc (CITED cross-check): Q(2cos 2pi/7) 49; y^3-y^2-10y+8 961=31^2; Q(2cos 2pi/9) 81; the parabolic")
    print("  cubic 8X^3+4X^2-18X-1 has field discriminant 49 and the s=1/4 cubic 64X^3+16X^2-164X+23 has 961:", nfd)
else:
    print("  (PARI nfdisc cross-check skipped: gp unavailable or failed.)")
print("  Inherited: geometry note sec.2 already proves the -7/4 parabolic cycle is the affine image L(x)=-x-1/2 of the reversed")
print("  conductor-7 cycle of x^2-2, so the field equality below is a re-verification by exact embedding, not a discovery.")
print("  FINITE-EXACT/PROVED: the parabolic 3-cycle of x^2-7/4 and the '2^3-1=7' 3-cycle of x^2-2 generate the SAME field")
print("  Q(2cos 2pi/7); the second 3-cycle at -29/16 generates the conductor-31 cyclic cubic; the AP cycle splits.")
print("  Which cycle at c=-2 lies where: sigma=-1 <=> sum 2cos(2pi k/7) = -1 (conductor 7 = 2^3-1);")
print("  sigma=0 <=> sum 2cos(2pi k/9) = 0 (conductor 9 = 2^3+1 = 3^2, the Catalan side).")

# ---------------------------------------------------------------------------------------
hr("S5. Complete rational preperiodic graphs (exact algorithm), isomorphism types, and the 'microcosm' incidence")

def beta_escape(xq, cq):
    """True iff |x| > beta(c) = (1+sqrt(1-4c))/2, i.e. x^2 - |x| + c > 0 with |x|>=1/2  (exact)."""
    ax = abs(xq)
    return ax >= Fr(1, 2) and ax * ax - ax + cq > 0

def preper_set(cq, maxsteps=None):
    """PrePer(f_c, Q) as a dict x -> f(x).  PROVED complete:
      * odd p with v_p(c)<0: v_p(f(x)) = 2 v_p(x) whenever v_p(x)<0 and 2v_p(x) != v_p(c) -> need v_p(x) = v_p(c)/2, c's
        denominator exponent must be even; if v_p(x) >= 0 and v_p(c)<0, v_p(f(x)) = v_p(c) then escape.  Same at p=2.
      * so x = u/D, D = prod p^{m_p}, den(c) = D^2; real: c<=1/4 and |x| <= beta(c) else strict escape.
      * orbit must keep denominator D and stay in [-beta, beta]; finitely many candidates so it revisits."""
    if cq > Fr(1, 4):
        return {}
    den = cq.denominator
    D = 1
    for p, e in sp.factorint(den).items():
        if e % 2:
            return {}
        D *= p ** (e // 2)
    # bound |u| <= D*beta
    U = (D * (1 + isqrt(int((1 - 4 * cq) * 4 * D * D)) // (2 * D) + 2)) // 2 + D + 2
    U = int(U) + 2
    graph = {}
    for u in range(-U, U + 1):
        xq = Fr(u, D)
        if beta_escape(xq, cq):
            continue
        orbit = []
        y = xq; ok = False
        seen = set()
        while True:
            if y in graph or y in seen:
                ok = True; break
            if y.denominator != D and not (D == 1 and y.denominator == 1):
                ok = False; break
            if beta_escape(y, cq):
                ok = False; break
            seen.add(y); orbit.append(y)
            y = y * y + cq
            if len(orbit) > 4 * U + 10:
                ok = False; break
        if ok:
            for z in orbit:
                graph[z] = z * z + cq
    return graph

def graph_type(graph):
    """Canonical description: for each cycle, cycle length and the rooted preimage-tree shape at each cycle point."""
    if not graph:
        return "empty"
    preim = {}
    for a_, b_ in graph.items():
        preim.setdefault(b_, []).append(a_)
    incycle = set()
    for v in graph:
        y = v
        for _ in range(len(graph) + 1):
            y = graph[y]
        # y is in a cycle
        z = y
        while True:
            incycle.add(z); z = graph[z]
            if z == y:
                break
    def tree(v):
        kids = [tree(w) for w in preim.get(v, []) if w not in incycle]
        return tuple(sorted(kids))
    comps = []
    done = set()
    for v in sorted(incycle):
        if v in done:
            continue
        cyc = []; z = v
        while z not in done:
            done.add(z); cyc.append(z); z = graph[z]
        shapes = tuple(sorted(tree(w) for w in cyc))
        comps.append((len(cyc), shapes))
    return str(sorted(comps))

def show_graph(cq):
    g = preper_set(cq)
    edges = ", ".join(f"{a_}->{b_}" for a_, b_ in sorted(g.items()))
    return g, edges

specials = [Fr(1, 4), Fr(3, 16), Fr(0), Fr(-3, 4), Fr(-1), Fr(-13, 16), Fr(-21, 16), Fr(-7, 4), Fr(-2), Fr(-29, 16), Fr(-37, 16), Fr(-77, 16)]
graphs = {}
print("c            |PrePer|  type [(cycle length, preimage-tree shapes per cycle point)]   edges")
for cq in specials:
    g, edges = show_graph(cq)
    graphs[cq] = g
    print(f"{str(cq):<12} {len(g):>7}  {graph_type(g):<60} {edges}")
# hard checks against inherited / hand results
check(set(graphs[Fr(-29, 16)]) == {Fr(m, 4) for m in (-7, -5, -3, -1, 1, 3, 5, 7)}, "PrePer(-29/16) = THM-4139 (2)")
check(graphs[Fr(-29, 16)][Fr(-7, 4)] == Fr(5, 4) and graphs[Fr(-29, 16)][Fr(5, 4)] == Fr(-1, 4) and graphs[Fr(-29, 16)][Fr(-1, 4)] == Fr(-7, 4), "3-cycle")
check(set(graphs[Fr(0)]) == {Fr(-1), Fr(0), Fr(1)}, "PrePer(0)")
check(set(graphs[Fr(-1)]) == {Fr(-1), Fr(0), Fr(1)}, "PrePer(-1)")
check(set(graphs[Fr(-2)]) == {Fr(k) for k in range(-2, 3)}, "PrePer(-2)")
check(set(graphs[Fr(-7, 4)]) == {Fr(1, 2), Fr(-1, 2), Fr(3, 2), Fr(-3, 2)}, "PrePer(-7/4) = {+-1/2,+-3/2}")
check(set(graphs[Fr(-3, 4)]) == {Fr(1, 2), Fr(-1, 2), Fr(3, 2), Fr(-3, 2)}, "PrePer(-3/4) = {+-1/2,+-3/2}")
check(graph_type(graphs[Fr(-7, 4)]) != graph_type(graphs[Fr(-3, 4)]), "-7/4 and -3/4: same point set, different graphs")
check(set(graphs[Fr(-21, 16)]) == {Fr(m, 4) for m in (-7, -5, -3, -1, 1, 3, 5, 7)}, "PrePer(-21/16) = quarter-odd integers |m|<=7")
print("User's naming corrected: 'x^2+0' = c=0 (0,1 fixed; -1->1); 'x^2+1' as described is x^2-1 (0<->-1, 1->0);")
print("  'x^2+2' as described is x^2-2 (2,-1 fixed; 1->-1; 0->-2->2).  All three graphs above match the user's description.")
print("PROVED: PrePer(x^2-7/4,Q) = {+-1/2,+-3/2} (2-cycle 1/2<->-3/2, -1/2->-3/2, 3/2->1/2); hence 0,-1,-2 are NOT")
print("  preperiodic for -7/4, and -7/4 is not preperiodic for itself (v_2 = -2 != -1).")
print("PROVED: -3/4 and -7/4 have the SAME rational preperiodic set but non-isomorphic graphs (two fixed points vs one 2-cycle).")
print("Isomorphism classes among the listed c (same type string => isomorphic functional graphs):")
types = {}
for cq, g in graphs.items():
    types.setdefault(graph_type(g), []).append(str(cq))
for k_, v_ in types.items():
    print("  ", v_, ":", k_)
print("  (Poonen 1998 lists 12 rational preperiodic graphs for x^2+c conjecturally exhaustive; labels not imported: UNCITED-RECOLLECTION.)")

# incidence matrix: is c_i a preperiodic point of f_{c_j}?
print()
print("Incidence: row = point x (a special parameter), column = parameter c; '*' iff x in PrePer(f_c,Q).")
hdr = "x \\ c       " + " ".join(f"{str(cq):>7}" for cq in specials)
print(hdr)
for xq in specials:
    row = " ".join(f"{'*':>7}" if xq in graphs[cq] else f"{'.':>7}" for cq in specials)
    print(f"{str(xq):<12} {row}")
check(Fr(-7, 4) in graphs[Fr(-29, 16)] and Fr(-7, 4) in graphs[Fr(-21, 16)] and Fr(-7, 4) in graphs[Fr(-37, 16)] and Fr(-7, 4) in graphs[Fr(-77, 16)],
      "-7/4 is preperiodic for -21/16 (to fixed 7/4), periodic of period 3 for -29/16, 2 for -37/16, 1 for -77/16")

def params_with_preperiodic_point(xq):
    """All c in Q with xq in PrePer(f_c,Q).  PROVED complete: den(c) = den(x)^2 exactly (valuation argument as in
    preper_set), c <= 1/4, c <= |x|-x^2 if |x|>=1/2 (else |x|>beta), and x^2+c >= -beta(c) which gives
    c >= -(1+x^2) - sqrt(1+x^2)."""
    D = xq.denominator
    den = D * D
    cmax = Fr(1, 4)
    if abs(xq) >= Fr(1, 2):
        cmax = min(cmax, abs(xq) - xq * xq)
    cmin_f = -(1 + float(xq) ** 2) - (1 + float(xq) ** 2) ** 0.5
    amin = int(cmin_f * den) - 2
    amax = int(cmax * den) + 1
    out = []
    for a in range(amin, amax + 1):
        if gcd(a, den) != 1 and den > 1:
            continue
        cq = Fr(a, den)
        if cq > cmax:
            continue
        # iterate x
        seen = set(); y = xq; ok = False
        for _ in range(4 * (2 * int(abs(cmin_f)) + 2) * D + 50):
            if y in seen:
                ok = True; break
            if beta_escape(y, cq) or (y.denominator != D and not (D == 1 and y.denominator == 1)):
                break
            seen.add(y); y = y * y + cq
        if ok:
            # exact period of the eventual cycle
            orbit = []; y = xq
            while y not in orbit:
                orbit.append(y); y = y * y + cq
            tail = orbit.index(y); per = len(orbit) - tail
            out.append((cq, tail, per))
    return out

print()
print("All c in Q for which a given x is preperiodic (complete lists, PROVED by the valuation+real-escape bounds):")
for xq in [Fr(1, 4), Fr(-3, 4), Fr(-7, 4), Fr(-29, 16), Fr(0), Fr(-1), Fr(-2), Fr(-1, 4), Fr(5, 4)]:
    L = params_with_preperiodic_point(xq)
    print(f"  x={str(xq):<7}: " + ", ".join(f"c={cq} (tail {tl}, period {pr})" for cq, tl, pr in L))
    if xq == Fr(-7, 4):
        check([(cq, tl, pr) for cq, tl, pr in L] == [(Fr(-93, 16), 1, 2), (Fr(-77, 16), 0, 1), (Fr(-37, 16), 0, 2), (Fr(-29, 16), 0, 3), (Fr(-21, 16), 1, 1)],
              "parameters with -7/4 preperiodic: -93/16 (tail 1 -> 2-cycle {-11/4, 7/4}), -77/16, -37/16, -29/16, -21/16")
    if xq == Fr(1, 4):
        check([(cq, tl, pr) for cq, tl, pr in L] == [(Fr(-29, 16), 1, 3), (Fr(-21, 16), 0, 2), (Fr(-13, 16), 1, 2), (Fr(-5, 16), 1, 1), (Fr(3, 16), 0, 1)],
              "parameters with 1/4 preperiodic")
    if xq == Fr(-3, 4):
        check([(cq, tl, pr) for cq, tl, pr in L] == [(Fr(-45, 16), 2, 1), (Fr(-37, 16), 1, 2), (Fr(-29, 16), 2, 3), (Fr(-21, 16), 0, 1), (Fr(-13, 16), 0, 2), (Fr(-5, 16), 2, 1), (Fr(3, 16), 1, 1)],
              "parameters with -3/4 preperiodic")
    if xq == Fr(-29, 16):
        print("     (-29/16 is a fixed point iff c = x - x^2 = -1305/256, in a 2-cycle iff c = -1-x-x^2 = -633/256.)")
L29 = params_with_preperiodic_point(Fr(-29, 16))
check([(cq, tl, pr) for cq, tl, pr in L29] == [(Fr(-1561, 256), 1, 2), (Fr(-1305, 256), 0, 1), (Fr(-633, 256), 0, 2), (Fr(-377, 256), 1, 1)],
      "parameters with -29/16 preperiodic")
check(all(pr <= 2 for _, _, pr in L29), "-29/16 never has exact period >= 3 over Q")
# exact period n parameters via rational roots of Phi_n(x0, c), n=1..4
print("Rational c with x0 of EXACT period n, from rational roots of Phi_n(x0,c) in c (n<=4):")
Phi1 = f1 - x; Phi2 = sp.expand(sp.cancel((f2 - x) / (f1 - x)))
for x0 in [R(1, 4), R(-3, 4), R(-7, 4), R(-29, 16), R(-2)]:
    per = {}
    for n, Ph in [(1, Phi1), (2, Phi2), (3, Phi3), (4, Phi4)]:
        poly = sp.Poly(sp.expand(Ph.subs(x, x0)), c)
        roots = [rt for rt in sp.roots(poly, filter='Q').keys()]
        # exclude degenerate (lower period) roots
        good = []
        for rt in roots:
            cq = Fr(int(rt.p), int(rt.q)); y = Fr(int(x0.p), int(x0.q)); orb_ = [y]
            for _ in range(n):
                y = y * y + cq; orb_.append(y)
            if orb_[-1] == orb_[0] and len(set(orb_[:-1])) == n:
                good.append(cq)
        per[n] = good
    print(f"  x0={x0}: " + "; ".join(f"period {n}: {per[n] if per[n] else 'none'}" for n in (1, 2, 3, 4)))
    check(per[4] == [], "no rational period-4 points (Morton 1998, CITED) - consistent")
print("  -7/4 has exact period 1,2,3 for c = -77/16, -37/16, -29/16 respectively (each unique), and no rational period 4")
print("  (Morton 1998: x^2+c has no rational points of exact period 4 - CITED; period 5: Flynn-Poonen-Schaefer 1997 - CITED;")
print("  period 6: Stoll 2008 under BSD - CITED).  So the chain 'parameter -> point of the next system' is exactly:")
print("  1/4 (period-1 parabolic) is a fixed point of c=3/16 and a 2-cycle point of c=-21/16;  -3/4 (period-2 parabolic)")
print("  is a fixed point of c=-21/16 and a 2-cycle point of c=-13/16;  -7/4 (period-3 parabolic) has periods 1,2,3 for")
print("  c=-77/16,-37/16,-29/16.  -29/16 has periods 1,2 for c=-1305/256,-633/256 and no period 3 over Q.")

# ---------------------------------------------------------------------------------------
hr("S6. (iv) Chebyshev conjugacy, 63 = 7*9, and the typed analogy Bang(n=6) vs critical-orbit(n=3)")
# x = y + 1/y conjugates y->y^2 to x->x^2-2
y = sp.symbols('y')
check(sp.simplify((y + 1/y)**2 - 2 - (y**2 + 1/y**2)) == 0, "Chebyshev conjugacy")
# orbit of x=5/2 (y=2) under x^2-2
xq = Fr(5, 2); ferm = []
for k in range(1, 6):
    ferm.append(xq.numerator); xq = xq * xq - 2
print("x=y+1/y conjugates y->y^2 to x->x^2-2.  y=2 <-> x=5/2; orbit numerators under x^2-2:", ferm,
      " = Fermat numbers F_1..F_5 (pairwise coprime: every term has a primitive prime; no Bang-type exception).")
check(ferm == [2**(2**k) + 1 for k in range(1, 6)], "Fermat numerators")
print("Where 2^n-1 really enters x^2-2: period-n points of x^2-2 are 2cos(2 pi k/N) with N | 2^n-1 or N | 2^n+1")
print("  (y^(2^n) = y^(+-1)).  For n=3: N | 7 or N | 9, i.e. Phi_3(x,-2) = (x^3+x^2-2x-1)(x^3-3x+1) [S4], and")
print("  63 = 2^6-1 = (2^3-1)(2^3+1) = 7*9 is the product of the two conductors of the two 3-cycles of x^2-2.")
print("  Bang's n=6 exception <=> 2^3+1 = 9 = 3^2 has no new prime <=> no prime-conductor 3-cycle of x^2-2 on the 2^n+1 side")
print("  <=> no prime p with ord_p(2)=6 (S1).  THM-4139 (37): the exact 6-cycles of squaring have conductors 9,21,63 only.")
# 2^n+1 side: which n<=40 have no prime conductor from 2^n+1: only n=3 (S1). Cross-check period-n counts of x^2-2 over roots of unity
def period_points_x2m2(n):
    # number of points of exact period n of x^2-2 on the real segment [-2,2] = (1/2)*#{y on unit circle, y != +-1 ..}
    # use dynatomic degree count: exact period-n points of a degree-2 poly: sum_{d|n} mu(n/d) 2^d
    return sum(sp.mobius(n // d) * 2**d for d in sp.divisors(n))
print("  Exact-period counts of x^2-2 (all points): n=1..6 ->", [period_points_x2m2(n) for n in range(1, 7)],
      " (n=3: 6 = 3 [conductor 7] + 3 [conductor 9]).")
print()
print("TYPED ANALOGY.")
print("  Source: Bang/Zsigmondy for 2^n-1, index n=6, primitive part Phi_6(2)=3 (non-primitive since 3|6, 3|2^2-1).")
print("  Target: critical orbit of x^2+c, index n=3, primitive part F(a,b)=a^3+2a^2b+ab^2+b^3 at (a,b)=(-7,4): value 1.")
print("  Map: NONE on objects (multiplicative group orbit vs polynomial critical orbit).  Shared predicate: 'the n-th term")
print("  has no primitive prime divisor' <=> 'no prime p sees the base point with exact period n mod p' (both PROVED:")
print("  p primitive for 2^n-1 <=> ord_p(2)=n; p primitive for N_n <=> 0 has exact period n under x^2+c mod p).")
print("  Preserved structure: a 'cyclotomic-type' factorization of the n-th term over d|n: 2^n-1 = prod Phi_d(2) versus")
print("  G_n(c) = prod_{d|n} H_d(c) (Gleason: G_1=c, G_2=c(c+1), G_3=c H_3, G_4=c(c+1)H_4, verified n<=6 in S3b), with the")
print("  primitive part = the d=n factor.  Lost: in Bang, non-primitive primes of Phi_n(a) divide n (uniform theorem, all n);")
print("  for H_n(a,b) the factors are pairwise coprime for n<=6 (S3b) but no uniform statement is proved here, and the")
print("  quadratic side has infinitely many n=2 exceptions (c=-1+-1/b) where Zsigmondy has finitely many.  Sidecar at n=3: the Thue")
print("  equation F(a,b)=+-1 = units of Z[rho] with vanishing rho-coefficient (S3).  Cheapest decisive test: Phi_6(2)=3")
print("  (a prime dividing 6) versus F(-7,4)=1 (a unit) - computed above, exact.")
print("  Verdict on '63 <-> -7/4': REFUTED as a structural identification.  Exact facts only: (a) f^3(0)/f(0)=2^-6=1/(63+1)")
print("  with exponent 2^3-2=6 (equal to Bang's 6 only because 2^n-2=2n iff n=3); (b) 63=7*9 is attached to c=-2, the")
print("  Chebyshev parameter, as the product of the conductors of its two 3-cycles (geometry note (6), THM-4139 (36)-(38)); c=-2 = -7/4 - (1/2)^2 and")
print("  -29/16 = -7/4 - (1/4)^2 are the s=1/2 and s=1/4 points of the 3-cycle parabola centred at -7/4; (c) the parabolic")
print("  3-cycle at -7/4 and the conductor-7 3-cycle at -2 generate the same field Q(2cos 2pi/7) (geometry note sec.2, inherited).")
print("  Gaussian-squaring face (inherited, geometry (13)-(16), THM-3341 sec.4, THM-3333): x=2A/C of a primitive triple obeys")
print("  x -> x^2-2 exactly, with denominators C^(2^n) and pairwise coprime odd legs, i.e. every step of such an x^2-2 orbit has")
print("  a primitive numerator prime; the Fermat orbit of 5/2 is the |y|=2 (non-unit-circle) instance of the same conjugacy.")
print("  SCOPE: no map found from THM-3341's Pell-hypotenuse selector to the critical-orbit exception; not claimed.")

# ---------------------------------------------------------------------------------------
hr("S7. Status ledger")
print("PROVED     : S1 Catalan reading of n=6; no prime of order 6; S2 n=1,2,3 characterizations (n=3 <=> Thue, inherited geometry (11));")
print("             S3b G_n = prod H_d, Res(H_d,H_n)=+-1 (d<n<=6) => non-Zsigmondy at n<=6 iff H_n(a,b)=+-1;")
print("             S3 disc_x Phi_3 = -(4c+7)^3(16c^2+4c+7)^2, parabolic multiplier 1 at -7/4, sigma-parabola c=-7/4-s^2,")
print("             lambda(s), 64G_3-c identity, airplane centre = -rho^2, unit correspondence; S4 Morton identities,")
print("             t=(p1-p2)/(p0-p1), c(t)=-7/4-s(t)^2, Phi_3(x,-2) factorization; S5 all PrePer sets and incidence lists.")
print("FINITE-EXACT: Bang table n<=40; 2^n+1 table; census box (n<=8); Thue search b<=10^6; unit search |k|<=300;")
print(f"             H_4,H_5,H_6 = +-1 boxes (b<={BFULL} full, b<={BROOT} near real roots);")
print("             exact-embedding certificates for the cycle fields; n=4 pattern gcd.")
print("CITED      : Bang 1886 / Zsigmondy 1892 (a^n-b^n exceptions: n=1 with a-b=1; n=2 with a+b a power of 2; (2,1,6));")
print("             Mihailescu 2004 (Catalan); Morton 1998 (no rational period 4); Flynn-Poonen-Schaefer 1997 (period 5);")
print("             Stoll 2008 (period 6, BSD); PARI thue (Bilu-Hanrot) for the complete F, H_4 and H_5 Thue lists; PARI nfdisc;")
print("             Krieger, Primitive prime divisors in the critical orbit of z^d+c (via the geometry note, read there 2026-09-17);")
print("             THM-4139/THM-4146; geometry note (4),(9),(11),(16).")
print("UNCITED-RECOLLECTION: Poonen 1998 graph list labels; Doerksen-Haensch 2012 primitive-divisor theorem; Gleason's")
print("             2-adic simple-root theorem as the only known input toward a uniform Res(H_d,H_n)=+-1.")
print("REFUTED    : 'N_n=N_1 at every real parabolic parameter' (fails n=4, exact gcd); '63 <-> -7/4 structural map';")
print(f"             'the -7/4 n=3 failure is the only non-preperiodic failure in the box' (n=1,2 have {len(failures[1])} and {len(failures[2])} more).")
print("OPEN       : the Thue equation H_6(a,b)=+-1 (degree 27) beyond the finite box (H_4, H_5 are closed by PARI when gp is")
print("             present); any n>=7 statement (no resultant computation, no reduction); a uniform 'Res(H_d,H_n)=+-1'")
print("             theorem for all d<n (only n<=6 computed).")
print()
print(f"ALL CHECKS PASSED ({time.time() - T0:.1f}s).  Failures: {FAILS}")
