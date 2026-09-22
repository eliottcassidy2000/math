#!/usr/bin/env python3
"""
collatz_mod6_20260917_zsigmondy_triad_audit.py   (session collatz-mod6-20260917, lane zsigmondy_triad, ADVERSARIAL AUDIT)

Independent recomputation of every load-bearing number in
  05-knowledge/results/collatz_mod6_20260917_zsigmondy_triad.md
with code that does NOT import the lane script.  Polynomial algebra is done in PARI/GP (polresultant, poldisc, thue,
nfdisc, nfisisom, factor) where the lane used SymPy, and in SymPy where the lane used hand identities; orbits and
preperiodic enumeration use Fractions with a different (pigeonhole) completeness argument.  Hostile boxes go beyond
the lane's universe (all denominators b<=64, not only powers of two).  Every check is an explicit `raise` (python -O safe).

Output: 05-knowledge/results/collatz_mod6_20260917_zsigmondy_triad_audit.out
"""
import re, time, shutil, subprocess, itertools
from fractions import Fraction as Fr
from math import gcd, isqrt
import sympy as sp

T0 = time.time()
NCHK = 0

def chk(cond, msg):
    global NCHK
    NCHK += 1
    if not cond:
        raise RuntimeError("AUDIT CHECK FAILED: " + msg)
    print("  ok  " + msg)

def hr(t):
    print("\n" + "=" * 78 + "\n" + t + "\n" + "=" * 78)

GP = shutil.which("gp")
def gp(script, timeout=300):
    """Run a PARI/GP script (fresh process, no gprc); return the list of printed lines."""
    if not GP:
        raise RuntimeError("PARI/GP (gp) is required for this audit and was not found")
    full = "default(parisizemax, 600000000);\n" + script + "\n"
    r = subprocess.run([GP, "-q", "-f"], input=full, capture_output=True, text=True, timeout=timeout)
    if r.returncode != 0 and not r.stdout.strip():
        raise RuntimeError("gp failed: " + r.stderr[:500])
    return [l for l in r.stdout.splitlines() if l.strip() and not l.startswith("  ***")]

def gpint(s):
    return int(s.strip())

def pairs(line):
    return {(int(a), int(b)) for a, b in re.findall(r"\[(-?\d+),\s*(-?\d+)\]", line)}

def norm(S):
    out = set()
    for a, b in S:
        if b < 0 or (b == 0 and a < 0):
            a, b = -a, -b
        out.add((a, b))
    return out

def prim_primes(terms, n):
    """primes of terms[n] dividing no terms[m], 1<=m<n; by FULL factorization (SymPy factorint).  terms[n] != 0."""
    ps = sorted(sp.factorint(abs(terms[n])).keys()) if abs(terms[n]) > 1 else []
    return [p for p in ps if all(terms[m] % p != 0 for m in range(1, n))]

def prim_part(terms, n):
    """gcd-stripped part of |terms[n]| coprime to all earlier terms (own implementation)."""
    r = abs(terms[n])
    g = 1
    for m in range(1, n):
        if terms[m] != 0:
            g = g * abs(terms[m])
    while True:
        d = gcd(r, g)
        if d == 1:
            return r
        r //= d

print("Audit of lane zsigmondy_triad.  gp =", GP, "| PARI version:", gp("print(version());")[0])

# ---------------------------------------------------------------------------------------
hr("A. Bang/Zsigmondy for 2^n-1 and 2^n+1 (n<=40): recomputed by full factorization")
M = [0] + [2**n - 1 for n in range(1, 41)]
P = [0] + [2**n + 1 for n in range(1, 41)]
exc_m = [n for n in range(1, 41) if not prim_primes(M, n)]
exc_p = [n for n in range(1, 41) if not prim_primes(P, n)]
chk(exc_m == [1, 6], f"2^n-1 without primitive prime, n<=40: {exc_m} == [1,6]")
chk(exc_p == [3], f"2^n+1 without primitive prime, n<=40: {exc_p} == [3]")
# non-primitive primes of Phi_n(2) and the classical lemma, n<=40
nonprim = {}
for n in range(3, 41):
    phi = gpint(gp(f"print(polcyclo({n},2));")[0])
    ps = prim_primes(M, n)
    fac = gp(f"f=factor({phi}); print(f[,1]~); print(f[,2]~);")
    primes_phi = [int(v) for v in re.findall(r"-?\d+", fac[0])]
    exps = [int(v) for v in re.findall(r"-?\d+", fac[1])]
    for p, e in zip(primes_phi, exps):
        if p not in ps:
            nonprim.setdefault(n, []).append((p, e))
            chk(p == max(sp.primefactors(n)) and e == 1, f"Phi_{n}(2)={phi}: non-primitive prime {p} is max prime of n, exponent {e}==1")
chk({n: [p for p, _ in v] for n, v in nonprim.items()} == {6: [3], 18: [3], 20: [5], 21: [7]},
    f"non-primitive primes of Phi_n(2), 3<=n<=40: {nonprim}")
chk([gpint(gp(f"print(polcyclo({d},2));")[0]) for d in (1, 2, 3, 6)] == [1, 3, 7, 3], "Phi_d(2), d|6 = 1,3,7,3")
chk(all(prim_part(M, n) != 1 for n in range(2, 41) if n != 6) and prim_part(M, 6) == 1 and prim_part(M, 1) == 1,
    "gcd-stripping primitive part agrees with factorization for 2^n-1, n<=40")
# orders
chk(gpint(gp("print(znorder(Mod(2,9)));")[0]) == 6, "ord_9(2) = 6")
mods6 = [m for m in range(3, 10001, 2) if pow(2, 6, m) == 1 and all(pow(2, k, m) != 1 for k in (1, 2, 3))]
chk(mods6 == [9, 21, 63], f"odd moduli m<=10^4 with ord_m(2)=6: {mods6} (all composite; least is 9)")
chk(all(not sp.isprime(m) for m in mods6), "no prime has ord_p(2)=6 (the moduli of order 6 are exactly the divisors of 63 not dividing 7 or 3: 9,21,63)")
chk(all((pow(2, e, 9) == 2) == (e % 6 == 1) and (pow(2, e, 9) == 5) == (e % 6 == 5) and (pow(2, e, 9) == 8) == (e % 6 == 3) for e in range(0, 200)),
    "2^e = 2,5,8 (mod 9) iff e = 1,5,3 (mod 6), e<200")
chk(2**6 - 1 == 7 * 9 and 3**2 - 2**3 == 1 and 9 // gcd(9, 3) == 3, "63 = 7*9, 3^2-2^3 = 1, (2^3+1)/gcd(2^3+1,2^2-1) = 3")
chk([n for n in range(1, 61) if 2**n - 2 == 2 * n] == [3], "2^n-2 = 2n iff n=3 (n<=60)")

# ---------------------------------------------------------------------------------------
hr("B. Critical orbit of x^2-7/4 to n=8")
def orbit(c, nmax):
    x = Fr(0); xs = [x]
    for _ in range(nmax):
        x = x * x + c; xs.append(x)
    return xs
O = orbit(Fr(-7, 4), 8)
N = [o.numerator for o in O]
chk(O[1:5] == [Fr(-7, 4), Fr(21, 16), Fr(-7, 256), Fr(-114639, 65536)], "f^n(0), n<=4: -7/4, 21/16, -7/256, -114639/65536")
chk(all(O[n].denominator == 4 ** (2 ** (n - 1)) for n in range(1, 9)), "denominators 4^(2^(n-1)) (no cancellation), n<=8")
facs = {n: gp(f"print(factor({abs(N[n])}));")[0] for n in range(1, 7)}
for n in range(1, 7):
    print(f"   N_{n} = {'-' if N[n] < 0 else ''}{facs[n]}")
chk(abs(N[4]) == 3 * 7 * 53 * 103 and abs(N[5]) == 7 * 419 * 563 * 3407 and abs(N[6]) == 3 * 7 * 47237 * 636069519847,
    "N_4=-3*7*53*103, N_5=7*419*563*3407, N_6=-3*7*47237*636069519847 (PARI factor)")
chk(all(sp.isprime(p) for p in (53, 103, 419, 563, 3407, 47237, 636069519847)), "listed factors are prime")
pp = [prim_part(N, n) for n in range(1, 9)]
chk(pp[:6] == [7, 3, 1, 5459, 803701079, 30046015909012739], f"primitive parts n<=6: {pp[:6]}")
chk([bool(prim_primes(N, n)) for n in range(1, 7)] == [True, True, False, True, True, True], "primitive prime exists at n=1,2,4,5,6 and not at n=3 (full factorization)")
q7 = abs(N[7]) // 7
chk(abs(N[7]) % 7 == 0 and abs(N[7]) % 49 != 0 and len(str(q7)) == 38 and gp(f"print(isprime({q7}));")[0].strip() == "1",
    f"N_7 = -7 * (38-digit prime {str(q7)[:10]}...), PARI isprime")
chk(N[7] < 0, "N_7 < 0")
chk(len(str(pp[7])) == 73 and pp[6] == q7, f"primitive part at n=7 is that prime; at n=8 it has {len(str(pp[7]))} digits (== 73)")
chk(pp[7] != 1 and gp(f"print(gcd({pp[7]},{7*3*53*103*419*563*3407*47237*636069519847*q7}));")[0].strip() == "1", "n=8 primitive part is coprime to all earlier numerators")
chk(O[3] / O[1] == Fr(1, 64), "f^3(0)/f(0) = 1/64 = 1/(63+1)")

# ---------------------------------------------------------------------------------------
hr("C. Census: lane box (b in {1,2,...,64} powers of two, |a|<=200) and a hostile box with ALL b<=64")
def zero_preperiodic(c):
    """exact: 0 is preperiodic for x^2+c (c rational) iff c in {0,-1,-2} (geometry note sec.3, PROVED there); recomputed:
    for b>=2 the denominators grow strictly; for integers iterate with escape |x|>|c|+1 => x^2+c > |x|."""
    if c.denominator > 1:
        return None
    seen = {}; x = Fr(0); k = 0
    while x not in seen:
        seen[x] = k; x = x * x + c; k += 1
        if abs(x) > abs(c) + 1:
            return None
    return (seen[x], k - seen[x])

def census(bs, amax, nmax):
    fails = {n: [] for n in range(1, nmax + 1)}; pre = []; cnt = 0
    for b in bs:
        for a in range(-amax, amax + 1):
            if gcd(a, b) != 1:
                continue
            c = Fr(a, b); cnt += 1
            pp_ = zero_preperiodic(c)
            if pp_ is not None:
                pre.append((c, pp_)); continue
            xs = orbit(c, nmax); nums = [o.numerator for o in xs]
            for n in range(1, nmax + 1):
                if prim_part(nums, n) == 1:
                    fails[n].append(c)
    return cnt, pre, fails

cnt, pre, fails = census([1, 2, 4, 8, 16, 32, 64], 200, 8)
chk(cnt == 1601, f"lane box size {cnt} == 1601")
chk(pre == [(Fr(-2), (2, 1)), (Fr(-1), (0, 2)), (Fr(0), (0, 1))], f"0 preperiodic exactly for c=-2 (tail 2, period 1), -1 (tail 0, period 2), 0 (tail 0, period 1): {pre}")
chk([len(fails[n]) for n in range(1, 9)] == [13, 12, 1, 0, 0, 0, 0, 0], f"failure counts n=1..8: {[len(fails[n]) for n in range(1, 9)]}")
chk(sorted(fails[1]) == sorted({Fr(1, b) for b in (1, 2, 4, 8, 16, 32, 64)} | {Fr(-1, b) for b in (2, 4, 8, 16, 32, 64)}), "n=1 failures = +-1/b minus the preperiodic c=-1 (13 of them)")
chk(sorted(fails[2]) == sorted({Fr(-1) + s * Fr(1, b) for b in (2, 4, 8, 16, 32, 64) for s in (1, -1)}), "n=2 failures = -1+-1/b, b>=2 (12 of them)")
chk(fails[3] == [Fr(-7, 4)], "n=3 failure in the lane box = -7/4 only")
# cross-check n<=4 by full factorization on the lane box (independent of gcd-stripping)
bad = 0
for b in (1, 2, 4, 8, 16, 32, 64):
    for a in range(-200, 201):
        if gcd(a, b) != 1 or Fr(a, b) in (Fr(0), Fr(-1), Fr(-2)):
            continue
        xs = orbit(Fr(a, b), 4); nums = [o.numerator for o in xs]
        for n in range(1, 5):
            if (prim_part(nums, n) == 1) != (not prim_primes(nums, n)):
                bad += 1
chk(bad == 0, "gcd-stripping and full-factorization primitive tests agree on the lane box for n<=4")
# hostile box: all denominators 1..64
cnt2, pre2, fails2 = census(range(1, 65), 120, 6)
print(f"   hostile box: {cnt2} parameters (all b<=64, |a|<=120), n<=6 failure counts {[len(fails2[n]) for n in range(1, 7)]}")
chk(pre2 == pre, "hostile box: same three preperiodic parameters")
chk(fails2[3] == [Fr(-7, 4)] and all(fails2[n] == [] for n in (4, 5, 6)), "hostile box (all b<=64): n=3 fails only at -7/4; no failures at n=4,5,6")
chk(all(abs(c.numerator) == 1 for c in fails2[1]) and all(abs(c.numerator + c.denominator) == 1 for c in fails2[2]),
    "hostile box: n=1 failures have |a|=1, n=2 failures have |a+b|=1 (odd denominators included)")
# the n<=3 hand proof, checked on random (a,b)
import random
random.seed(20260917)
for _ in range(2000):
    a = random.randint(-10**6, 10**6); b = random.randint(1, 10**6)
    if gcd(a, b) != 1:
        continue
    F = a**3 + 2*a*a*b + a*b*b + b**3
    if F % a != (b**3) % a if a not in (0,) else False:
        raise RuntimeError("F mod a")
    if (F - b**3) % (a + b) != 0 if a + b != 0 else False:
        raise RuntimeError("F mod a+b")
    c = Fr(a, b); xs = orbit(c, 3)
    if xs[3] != Fr(a * F, b**4) or xs[2] != Fr(a * (a + b), b * b):
        raise RuntimeError("closed forms")
chk(True, "N_2 = a(a+b), N_3 = a*F(a,b), F = b^3 mod a and mod (a+b): 2000 random (a,b) up to 10^6")
chk((-7)**3 + 2*49*4 + (-7)*16 + 64 == 1, "F(-7,4) = 1")

# ---------------------------------------------------------------------------------------
hr("D. Gleason product G_n = prod_{d|n} H_d (PARI): degrees, irreducibility, values, resultants n<=6, and EXTENSION to n=7,8")
GL = """
G=vector(8); G[1]=x; for(n=2,8, G[n]=G[n-1]^2+x);
H=vector(8); for(n=1,8, q=G[n]; fordiv(n,d, if(d<n, q=q/H[d])); H[n]=q);
"""
out = gp(GL + """
for(n=1,8, print(n," ",type(H[n])," ",poldegree(H[n])," ",pollead(H[n])," ",polisirreducible(H[n])," ",subst(H[n],x,0)," ",subst(H[n],x,-1)," ",subst(H[n],x,-2)," ",issquarefree(H[n])));
print(H[3]); print(H[4]);
for(n=2,8, s=""; for(d=1,n-1, s=concat(s, Str(polresultant(H[d],H[n]))); s=concat(s," ")); print(n,": ",s));
""", timeout=600)
rows = [l.split() for l in out[:8]]
chk(all(r[1] == "t_POL" for r in rows), "all quotients G_n / prod_{d|n,d<n} H_d are exact polynomials (n<=8)")
chk([int(r[2]) for r in rows] == [1, 1, 3, 6, 15, 27, 63, 120], f"deg H_n, n=1..8: {[int(r[2]) for r in rows]}")
chk(all(r[3] == "1" for r in rows), "every H_n is monic")
chk(all(r[4] == "1" for r in rows), "every H_n, n<=8, is irreducible over Q (PARI polisirreducible)")
chk(all(r[8] == "1" for r in rows), "every H_n, n<=8, is squarefree")
chk([r[5] for r in rows[1:]] == ["1"] * 7, "H_n(0) = 1 for 2<=n<=8")
chk([r[6] for r in rows[2:]] == ["1"] * 6, "H_n(-1) = 1 for 3<=n<=8")
chk([int(r[7]) for r in rows[:6]] == [-2, -1, -1, 1, -1, -1], f"H_n(-2), n=1..6: {[int(r[7]) for r in rows[:6]]}")
chk(out[8].replace(" ", "") == "x^3+2*x^2+x+1" and out[9].replace(" ", "") == "x^6+3*x^5+3*x^4+3*x^3+2*x^2+1", "H_3, H_4 as displayed in the note")
res = {}
for l in out[10:]:
    n, vals = l.split(":")
    res[int(n)] = [int(v) for v in vals.split()]
out2 = gp(GL + "for(n=1,8, print(Vec(H[n])));")
Hco = [[int(v) for v in re.findall(r"-?\d+", l)] for l in out2[:8]]
want = {2: [1], 3: [1, 1], 4: [1, 1, -1], 5: [1, 1, -1, 1], 6: [1, 1, -1, -1, -1]}
chk({n: res[n] for n in range(2, 7)} == want, f"resultant table Res(H_d,H_n), d<n<=6 (standard Sylvester convention, PARI): {want}")
# the standard convention, independently: Sylvester determinant in SymPy; and the SymPy-PRS sign quirk that the lane's draft table inherited
cc = sp.Symbol('c')
Hs = {n: sp.Poly(sum(v * cc**(len(Hco[n-1]) - 1 - i) for i, v in enumerate(Hco[n-1])), cc) for n in range(1, 7)}
def sylvester(f, g):
    m, n_ = f.degree(), g.degree(); M = sp.zeros(m + n_, m + n_)
    for i in range(n_):
        for j, v in enumerate(f.all_coeffs()): M[i, i + j] = v
    for i in range(m):
        for j, v in enumerate(g.all_coeffs()): M[n_ + i, i + j] = v
    return int(M.det())
syl = {n: [sylvester(Hs[d], Hs[n]) for d in range(1, n)] for n in range(2, 7)}
chk(syl == want, "SymPy Sylvester determinants agree with PARI (standard convention)")
prs = {n: [int(sp.resultant(Hs[d].as_expr(), Hs[n].as_expr(), cc)) for d in range(1, n)] for n in range(2, 7)}
print(f"   SymPy {sp.__version__} resultant() (PRS) signs: {prs}")
chk(prs == {2: [1], 3: [-1, -1], 4: [1, 1, -1], 5: [-1, -1, 1, 1], 6: [-1, -1, 1, -1, 1]},
    "the lane's original table = SymPy-PRS signs, which differ from the standard convention at 9 of 15 entries (e.g. Res(H_1,H_3) = H_3(0) = +1, not -1); only |Res|=1 is load-bearing")
print(f"   EXTENSION: Res(H_d,H_7), d<7: {res[7]};  Res(H_d,H_8), d<8: {res[8]}")
chk(all(v in (1, -1) for v in res[7] + res[8]), "Res(H_d,H_n) = +-1 also for n=7 (deg 63) and n=8 (deg 120), all d<n: the n<=6 reduction extends to n<=8 (FINITE-EXACT)")
# numerical check of the reduction theorem on random (a,b): N_n = prod H_d(a,b), pairwise coprime, n<=8
def Hhom(n, a, b):
    co = Hco[n - 1]; d = len(co) - 1; r = 0
    for i, v in enumerate(co):
        r = r * a + v * b ** i
    return r
random.seed(7)
ntested = 0
for _ in range(300):
    a = random.randint(-3000, 3000); b = random.randint(1, 3000)
    if gcd(a, b) != 1 or Fr(a, b) in (Fr(0), Fr(-1), Fr(-2)):
        continue
    ntested += 1
    xs = orbit(Fr(a, b), 8)
    for n in range(1, 9):
        hs = [Hhom(d, a, b) for d in sp.divisors(n)]
        prod_ = 1
        for h in hs: prod_ *= h
        if prod_ != xs[n].numerator:
            raise RuntimeError(f"N_n != prod H_d(a,b) at {a}/{b}, n={n}")
        for i in range(len(hs)):
            for j in range(i):
                if gcd(hs[i], hs[j]) != 1:
                    raise RuntimeError("factors not coprime")
        if (prim_part([o.numerator for o in xs], n) == 1) != (abs(Hhom(n, a, b)) == 1):
            raise RuntimeError("reduction theorem mismatch")
chk(ntested >= 150, f"N_n = prod_{{d|n}} H_d(a,b) with pairwise coprime factors and 'no primitive prime at n <=> H_n(a,b)=+-1': {ntested} random reduced c=a/b (|a|,b<=3000), all n<=8")

# ---------------------------------------------------------------------------------------
hr("E. Thue equations H_3=F, H_4, H_5 (PARI thue, flag 1 = unconditional) and the plastic units")
t_ = time.time()
out = gp(GL + """
for(n=3,5, tnf=thueinit(H[n],1); print(n," ",thue(tnf,1)," | ",thue(tnf,-1)));
""", timeout=600)
print(f"   (thueinit/thue for H_3,H_4,H_5 took {time.time() - t_:.1f}s)")
sol = {}
for l in out:
    n, rest = l.split(" ", 1)
    plus, minus = rest.split("|")
    sol[int(n)] = (pairs(plus), pairs(minus))
chk(norm(sol[3][0] | sol[3][1]) == {(-7, 4), (-1, 1), (0, 1), (1, 0), (-2, 1)} and sol[3][1] == {(-a, -b) for a, b in sol[3][0]},
    f"F(a,b)=+-1 complete solution set (mod sign) = {{(-7,4),(-2,1),(-1,1),(0,1),(1,0)}}: {sorted(norm(sol[3][0] | sol[3][1]))}")
chk(norm(sol[4][0]) == {(1, 0), (0, 1), (-1, 1), (-2, 1)} and len(sol[4][0]) == 8 and sol[4][1] == set(),
    f"H_4(a,b)=+1: 8 sign-variants of (1,0),(0,1),(-1,1),(-2,1); H_4=-1: none.  {sorted(sol[4][0])}")
chk(norm(sol[5][0] | sol[5][1]) == {(1, 0), (0, 1), (-1, 1), (-2, 1)} and sol[5][1] == {(-a, -b) for a, b in sol[5][0]},
    f"H_5(a,b)=+-1 solutions (mod sign) = (1,0),(0,1),(-1,1),(-2,1): {sorted(sol[5][0])}")
# plastic number
out = gp("""
r=Mod(x, x^3-x-1); print(poldisc(x^3-x-1)); print(minpoly(-r^2));
hits=[]; for(k=-300,300, u=lift(r^k); if(polcoeff(u,1)==0, hits=concat(hits,[[k,polcoeff(u,0),polcoeff(u,2)]]))); print(hits);
print(lift(r^-14)); print(lift(r^-5)); print(lift(r^-1)); print(lift(r^2));
print(nfinit(x^3-x-1).disc); K=bnfinit(x^3-x-1); u=K.fu[1]; print(u); print(#K.fu," ",(u==r)||(u==-r)||(u==r^-1)||(u==-r^-1)," ",K.no);
""")
chk(out[0].strip() == "-23", "disc(t^3-t-1) = -23")
chk(out[1].replace(" ", "") == "x^3+2*x^2+x+1", "minimal polynomial of -rho^2 is H_3 (PARI minpoly)")
hits = [tuple(int(v) for v in m) for m in re.findall(r"\[(-?\d+),\s*(-?\d+),\s*(-?\d+)\]", out[2])]
chk(hits == [(-14, -7, 4), (-5, 2, -1), (-1, -1, 1), (0, 1, 0), (2, 0, 1)], f"units rho^k = A + C rho^2 (zero rho-coefficient), |k|<=300: {hits}")
chk(out[3].replace(" ", "") == "4*x^2-7" and out[4].replace(" ", "") == "-x^2+2" and out[5].replace(" ", "") == "x^2-1" and out[6].replace(" ", "") == "x^2",
    "rho^-14 = 4rho^2-7, rho^-5 = 2-rho^2, rho^-1 = rho^2-1, rho^2 = rho^2")
chk(out[7].strip() == "-23", "field discriminant of Q(rho) is -23 (Z[rho] is the maximal order)")
print("   fundamental unit(s) of Q(rho) per PARI bnfinit:", out[8].strip())
chk(out[9].split() == ["1", "1", "1"], "PARI bnfinit: unit rank 1, fundamental unit = +-rho^{+-1} (unit group = +-rho^Z), class number 1  (CITED: PARI bnfinit, unconditional only after bnfcertify; used here as a cross-check, the Thue closure does not depend on it)")
rho = float(sp.N(sp.CRootOf(sp.Symbol('t')**3 - sp.Symbol('t') - 1, 0), 30))
chk(abs(-rho**2 - (-1.75487766625)) < 1e-10 and abs(-7/4 - (-rho**2)) < 0.0049, f"airplane centre -rho^2 = {-rho**2:.11f}, |-7/4 - centre| = {abs(-1.75 + rho**2):.6f} < 0.0049")
r2 = rho**2; hits2 = []
for b in range(1, 10**5 + 1):
    a0 = round(-r2 * b)
    for a in range(a0 - 3, a0 + 4):
        if gcd(a, b) == 1 and a**3 + 2*a*a*b + a*b*b + b**3 in (1, -1):
            hits2.append((a, b))
chk(sorted(hits2) == sorted([(-7, 4), (-2, 1), (-1, 1), (0, 1)]), f"direct search F=+-1, b<=10^5, a within 3 of -rho^2 b: {hits2}")
# boxes for H_6 (and H_4,H_5) = +-1
out = gp(GL + "for(n=4,6, print(n,\" \",polrootsreal(H[n])));")
realroots = {}
for l in out:
    n, rest = l.split(" ", 1)
    realroots[int(n)] = sorted(float(v) for v in re.findall(r"-?\d+\.\d+", rest))
print("   real roots:", {n: [round(v, 4) for v in realroots[n]] for n in (4, 5, 6)})
chk([round(v, 4) for v in realroots[4]] == [-1.9408, -1.3107] and [round(v, 4) for v in realroots[5]] == [-1.9854, -1.8608, -1.6254]
    and [round(v, 4) for v in realroots[6]] == [-1.9964, -1.9668, -1.9073, -1.7729, -1.476], "real roots of H_4, H_5, H_6 as in the note")
for n in (4, 5, 6):
    hits3 = set()
    for b in range(1, 5001):
        cand = range(-3 * b, 3 * b + 1) if b <= 100 else {a for r in realroots[n] for a in range(round(r * b) - 2, round(r * b) + 3)}
        for a in cand:
            if gcd(a, b) == 1 and abs(Hhom(n, a, b)) == 1:
                hits3.add((a, b))
    chk(hits3 == {(0, 1), (-1, 1), (-2, 1)}, f"H_{n}(a,b)=+-1 in the box (|a|<=3b, b<=100; near real roots b<=5000): {sorted(hits3)}")

# ---------------------------------------------------------------------------------------
hr("F. Dynatomic algebra (PARI): disc_x Phi_3, multiplier-1 locus, the parabolic square, 64G_3-c, the n=4 refutation")
DYN = """
f1=x^2+c; f2=subst(f1,x,f1); f3=subst(f1,x,f2); f4=subst(f1,x,f3);
Phi3=(f3-x)/(f1-x); Phi4=(f4-x)/(f2-x);
"""
out = gp(DYN + """
print(type(Phi3)," ",type(Phi4));
D3=poldisc(Phi3,x); print(D3 + (4*c+7)^3*(16*c^2+4*c+7)^2);
R3=polresultant(Phi3, deriv(f3,x)-1, x); print(R3 - ((4*c+7)*(16*c^2+4*c+7))^3);
P=8*x^3+4*x^2-18*x-1; print(64*subst(Phi3,c,-7/4) - P^2);
print(polrootsreal(P));
G3=subst(f3,x,0); G2=subst(f2,x,0); G4=subst(f4,x,0);
print(64*G3 - c - c*(4*c+7)*(16*c^2+4*c+9)); print(4*G2 - c - c*(4*c+3));
R4=polresultant(Phi4, deriv(f4,x)-1, x); D4=64*c^3+144*c^2+108*c+135;
print(R4 - ((4*c+5)*(16*c^2-8*c+5)*D4)^4);
print(gcd(2^14*G4 - c, (4*c+5)*(16*c^2-8*c+5)*D4));
print(polrootsreal(D4)); print(polisirreducible(D4));
""")
chk(out[0].split() == ["t_POL", "t_POL"], "Phi_3, Phi_4 are exact polynomial quotients")
chk(out[1].strip() == "0", "disc_x Phi_3 = -(4c+7)^3 (16c^2+4c+7)^2  (PARI poldisc, sign included)")
chk(out[2].strip() == "0", "Res_x(Phi_3, (f^3)'-1) = ((4c+7)(16c^2+4c+7))^3")
chk(out[3].strip() == "0", "64 Phi_3(x,-7/4) = (8x^3+4x^2-18x-1)^2")
roots = sorted(float(v) for v in re.findall(r"-?\d+\.\d+", out[4]))
chk(len(roots) == 3 and [round(v, 7) for v in roots] == [-1.7469796, -0.0549581, 1.3019377], f"parabolic cycle points {roots}")
chk(abs(8 * roots[0] * roots[1] * roots[2] - 1) < 1e-9 and Fr(1, 8) == Fr(1, 8), "multiplier 8*prod(x_i) = 1 (Vieta: product of roots = 1/8)")
chk(out[5].strip() == "0" and out[6].strip() == "0", "64G_3 - c = c(4c+7)(16c^2+4c+9);  4G_2 - c = c(4c+3)")
chk(out[7].strip() == "0", "Res_x(Phi_4,(f^4)'-1) = ((4c+5)(16c^2-8c+5)(64c^3+144c^2+108c+135))^4")
chk(out[8].strip() == "1", "gcd(2^14 G_4 - c, Delta_4) = 1: 'N_n = N_1 at the real parabolic parameter' FAILS at n=4 (REFUTED, minimal witness n=4)")
c4 = float(re.findall(r"(-?\d+\.\d+)", out[9])[0])
chk(abs(c4 - (-1.94055078898)) < 1e-10 and out[10].strip() == "1", f"real period-4 saddle node c_4 = {c4:.11f} (irrational: cubic irreducible)")
# 1/4 and -3/4 trivial cases of the pattern
chk(orbit(Fr(1, 4), 1)[1] == Fr(1, 4) and orbit(Fr(-3, 4), 2)[2] == Fr(-3, 16), "c=1/4: N_1=1; c=-3/4: f^2(0) = c/4")

# ---------------------------------------------------------------------------------------
hr("G. The 3-cycle parabola c = -7/4 - s^2, lambda(s), cycle fields (PARI nfdisc / nfisisom)")
s, sig, X, c = sp.symbols('s sigma X c')
Psig = X**3 - sig*X**2 - (sig**2 + 2*sig + 3)*X + (sig**3 + 2*sig**2 + 3*sig + 1)      # THM-4146 (11)
csig = -(sig**2 + sig + 2)                                                            # THM-4146 (10)
Ps = sp.expand(Psig.subs(sig, s - sp.Rational(1, 2)))
chk(sp.expand(csig.subs(sig, s - sp.Rational(1, 2)) + sp.Rational(7, 4) + s**2) == 0, "c(sigma) with sigma = s-1/2 is c = -7/4 - s^2")
lam = sp.expand(-8 * (sig**3 + 2*sig**2 + 3*sig + 1)).subs(sig, s - sp.Rational(1, 2))
chk(sp.expand(lam - (1 - 14*s - 4*s**2 - 8*s**3)) == 0, "lambda(s) = 1 - 14s - 4s^2 - 8s^3")
chk([lam.subs(s, v) for v in (0, sp.Rational(-1, 4), sp.Rational(1, 4))] == [1, sp.Rational(35, 8), sp.Rational(-23, 8)], "lambda(0)=1, lambda(-1/4)=35/8, lambda(1/4)=-23/8")
chk(sp.expand(sp.discriminant(Ps, X) - (4*s**2 + 2*s + 7)**2) == 0, "disc_X P_s = (4s^2+2s+7)^2  (= THM-4146 (12) (4sigma^2+6sigma+9)^2 = geometry (4) (eta^2+eta+7)^2, eta=2s)")
chk(sp.expand((4*sig**2 + 6*sig + 9).subs(sig, s - sp.Rational(1, 2)) - (4*s**2 + 2*s + 7)) == 0, "4sigma^2+6sigma+9 -> 4s^2+2s+7")
# Phi_3(X,c(s)) = P_s * P_{-s}
x = sp.Symbol('x')
f3 = x; 
for _ in range(3): f3 = sp.expand(f3**2 + c)
Phi3 = sp.cancel((f3 - x) / (x**2 + c - x))
chk(sp.expand(Phi3.subs({x: X, c: -sp.Rational(7, 4) - s**2}) - sp.expand(Ps * Ps.subs(s, -s))) == 0, "Phi_3(X, -7/4-s^2) = P_s(X) P_{-s}(X): the parabola double-covers the c-line, branched at s=0")
cubics = {}
for sv in (0, sp.Rational(-1, 4), sp.Rational(1, 4), sp.Rational(-1, 2), sp.Rational(1, 2)):
    Pv = sp.Poly(Ps.subs(s, sv), X)
    den = sp.ilcm(*[sp.Rational(v).q for v in Pv.all_coeffs()])
    cubics[sv] = sp.Poly(sp.expand(Pv.as_expr() * den), X)
chk(cubics[0].as_expr() == 8*X**3 + 4*X**2 - 18*X - 1, "s=0 cubic 8X^3+4X^2-18X-1")
chk(cubics[sp.Rational(-1, 4)].as_expr() == 64*X**3 + 48*X**2 - 132*X - 35, "s=-1/4 cubic 64X^3+48X^2-132X-35")
chk(cubics[sp.Rational(1, 4)].as_expr() == 64*X**3 + 16*X**2 - 164*X + 23, "s=1/4 cubic 64X^3+16X^2-164X+23")
chk(cubics[sp.Rational(-1, 2)].as_expr() == X**3 + X**2 - 2*X - 1 and cubics[sp.Rational(1, 2)].as_expr() == X**3 - 3*X + 1, "s=-1/2: X^3+X^2-2X-1 (2cos 2pi/7); s=1/2: X^3-3X+1 (2cos 2pi/9)")
chk(sorted(sp.Poly(cubics[sp.Rational(-1, 4)], X).ground_roots().keys()) == [sp.Rational(-7, 4), sp.Rational(-1, 4), sp.Rational(5, 4)], "s=-1/4 cubic splits: roots -7/4, -1/4, 5/4 (the AP cycle)")
chk([4*sv**2 + 2*sv + 7 for sv in (0, sp.Rational(-1, 4), sp.Rational(1, 4), sp.Rational(-1, 2), sp.Rational(1, 2))] == [7, sp.Rational(27, 4), sp.Rational(31, 4), 7, 9], "4s^2+2s+7 column: 7, 27/4, 31/4, 7, 9")
out = gp("""
A=8*x^3+4*x^2-18*x-1; B=64*x^3+16*x^2-164*x+23; C7=x^3+x^2-2*x-1; C9=x^3-3*x+1; C31=x^3-x^2-10*x+8;
print(nfdisc(C7)," ",nfdisc(C31)," ",nfdisc(C9)," ",nfdisc(A)," ",nfdisc(B));
print(polisirreducible(A)," ",polisirreducible(B)," ",polisirreducible(C31));
print(nfisisom(polredabs(A),C7)!=0," ",nfisisom(polredabs(B),C31)!=0," ",nfisisom(polredabs(A),C9)==0," ",nfisisom(C7,C9)==0," ",nfisisom(polredabs(B),C7)==0);
print(polgalois(C7)," ",polgalois(C9)," ",polgalois(C31)," ",polgalois(polredabs(B)));
print(polredabs(A)," ",polredabs(B));
""")
chk(out[0].split() == ["49", "961", "81", "49", "961"], f"nfdisc [C7, C31, C9, parabolic cubic, s=1/4 cubic] = {out[0].split()} == [49, 961, 81, 49, 961]")
chk(out[1].split() == ["1", "1", "1"], "parabolic cubic, s=1/4 cubic, y^3-y^2-10y+8 irreducible")
chk(out[2].split() == ["1", "1", "1", "1", "1"], "PARI nfisisom: parabolic cubic field = Q(2cos 2pi/7); s=1/4 cubic field = Q[y]/(y^3-y^2-10y+8); the 7-, 9-, 31-fields are pairwise distinct")
chk(all("3, 1" in g or "[3, 1" in g for g in re.findall(r"\[[^\]]*\]", out[3])), f"Galois groups all cyclic C3: {out[3]}")
print("   polredabs of the two non-reference cubics:", out[4])

# ---------------------------------------------------------------------------------------
hr("H. Families and Morton's chart (SymPy, own derivation)")
r, t = sp.symbols('r t')
chk([sp.Rational(1 - rv**2, 4) if isinstance(rv, int) else (1 - rv**2) / 4 for rv in (1, 3, sp.Rational(5, 2), sp.Rational(9, 2))] == [0, -2, sp.Rational(-21, 16), sp.Rational(-77, 16)],
    "fixed points rational iff c=(1-r^2)/4: r=1,3,5/2,9/2 -> 0,-2,-21/16,-77/16")
chk([-(3 + sv**2) / sp.Integer(4) for sv in (0, 1, 2, 3, sp.Rational(5, 2), sp.Rational(9, 2))] == [sp.Rational(-3, 4), -1, sp.Rational(-7, 4), -3, sp.Rational(-37, 16), sp.Rational(-93, 16)],
    "2-cycle x^2+x+c+1=0 rational iff c=-(3+s^2)/4: s=0,1,2,3,5/2,9/2 -> -3/4,-1,-7/4,-3,-37/16,-93/16")
chk(sp.discriminant(x**2 + x + c + 1, x) == -4*c - 3 and sp.discriminant(x**2 - x + c, x) == 1 - 4*c, "discriminants -3-4c and 1-4c")
chk(sp.expand((x**2 + x + c + 1).subs(c, sp.Rational(-3, 4)) - (x + sp.Rational(1, 2))**2) == 0, "c=-3/4: the 2-cycle degenerates to the double fixed point -1/2 (multiplier -1)")
Dt = 2*t*(t + 1)
p0 = (t**3 + 2*t**2 + t + 1) / Dt; p1 = (t**3 - t - 1) / Dt; p2 = -(t**3 + 2*t**2 + 3*t + 1) / Dt
ct = -(t**6 + 2*t**5 + 4*t**4 + 8*t**3 + 9*t**2 + 4*t + 1) / (4*t**2*(t + 1)**2)
chk(all(sp.simplify(u**2 + ct - v) == 0 for u, v in ((p0, p1), (p1, p2), (p2, p0))), "THM-4139 (13): p0->p1->p2->p0 under x^2+c(t)")
sig_t = sp.cancel(p0 + p1 + p2); s_t = sp.cancel(sig_t + sp.Rational(1, 2))
chk(sp.simplify(s_t - (t**3 + t**2 - 2*t - 1) / (2*t*(t + 1))) == 0 and sp.simplify(sig_t - (t**3 - 3*t - 1) / (2*t*(t + 1))) == 0, "s(t) = (t^3+t^2-2t-1)/(2t(t+1)), sigma(t) = (t^3-3t-1)/(2t(t+1))")
chk(sp.simplify(ct + sp.Rational(7, 4) + s_t**2) == 0, "c(t) = -7/4 - s(t)^2")
chk(sp.simplify(sig_t.subs(t, -1 / (t + 1)) - sig_t) == 0 and sp.simplify(sp.cancel((p1 - p2) / (p0 - p1)) - t) == 0, "sigma invariant under t -> -1/(t+1); t = (p1-p2)/(p0-p1)")
chk(ct.subs(t, 1) == sp.Rational(-29, 16) and [p0.subs(t, 1), p1.subs(t, 1), p2.subs(t, 1)] == [sp.Rational(5, 4), sp.Rational(-1, 4), sp.Rational(-7, 4)], "t=1: c=-29/16, points 5/4,-1/4,-7/4")
n7 = sp.factor(sp.numer(sp.together(s_t + sp.Rational(1, 2)))); n9 = sp.factor(sp.numer(sp.together(s_t - sp.Rational(1, 2))))
chk(sp.expand(n7 - (t**3 + 2*t**2 - t - 1)) == 0 and sp.expand(n9 - (t**3 - 3*t - 1)) == 0, "c(t)=-2 iff t^3+2t^2-t-1=0 (s=-1/2) or t^3-3t-1=0 (s=+1/2)")
chk(sp.expand(t**3 * (x**3 + x**2 - 2*x - 1).subs(x, 1/t) + (t**3 + 2*t**2 - t - 1)) == 0 and sp.expand((x**3 - 3*x + 1).subs(x, -t) + (t**3 - 3*t - 1)) == 0,
    "t^3+2t^2-t-1 is the (negated) reciprocal of the 2cos(2pi/7) cubic; t^3-3t-1 is the (negated) negation of the 2cos(2pi/9) cubic")
chk(sp.expand(Phi3.subs(c, -2) - (x**3 - 3*x + 1) * (x**3 + x**2 - 2*x - 1)) == 0, "Phi_3(x,-2) = (x^3-3x+1)(x^3+x^2-2x-1)")
chk(sp.minimal_polynomial(2*sp.cos(2*sp.pi/7), x) == x**3 + x**2 - 2*x - 1 and sp.minimal_polynomial(2*sp.cos(2*sp.pi/9), x) == x**3 - 3*x + 1, "minimal polynomials of 2cos(2pi/7), 2cos(2pi/9)")
chk(sp.expand(sp.discriminant(x**3 + x**2 - 2*x - 1, x)) == 49 and sp.discriminant(x**3 - 3*x + 1, x) == 81, "disc 49 and 81")

# ---------------------------------------------------------------------------------------
hr("I. Rational preperiodic sets: own enumerator (pigeonhole on the finite candidate set), graphs, incidence")
def beta_bound(cq):
    """smallest integer U with U/D >= beta(c) = (1+sqrt(1-4c))/2 (c<=1/4), D = sqrt(den c); exact via isqrt."""
    D = isqrt(cq.denominator)
    num = (1 - 4 * cq)                  # rational >= 0
    # sqrt(num) <= (isqrt(num.n * num.d) + 1) / num.d
    sq = Fr(isqrt(num.numerator * num.denominator) + 1, num.denominator)
    beta_up = (1 + sq) / 2
    return D, int(beta_up * D) + 1

def preper(cq):
    """PrePer(f_c, Q) as dict x -> f(x).  Complete: preperiodic x has den(x)=D with den(c)=D^2 (valuation argument, every
    prime), and every orbit point y satisfies |y| <= beta (else y^2+c > |y| and the orbit increases strictly), so the
    orbit lies in the finite set S = {u/D : |u| <= U}; after |S|+1 steps a repeat must occur (pigeonhole)."""
    if cq > Fr(1, 4):
        return {}
    for p, e in sp.factorint(cq.denominator).items():
        if e % 2:
            return {}
    D, U = beta_bound(cq)
    S = [Fr(u, D) for u in range(-U, U + 1)]
    Sset = set(S)
    g = {}
    for x0 in S:
        y = x0; seen = []
        for _ in range(len(S) + 2):
            if y in g or y in seen:
                for z in seen: g[z] = z * z + cq
                break
            if y not in Sset:
                break
            seen.append(y); y = y * y + cq
    return g

def canon(g):
    if not g: return "()"
    pre = {}
    for a, b in g.items(): pre.setdefault(b, []).append(a)
    cyc = set()
    for v in g:
        y = v
        for _ in range(len(g)): y = g[y]
        z = y
        while True:
            cyc.add(z); z = g[z]
            if z == y: break
    def tree(v):
        return "(" + "".join(sorted(tree(w) for w in pre.get(v, []) if w not in cyc)) + ")"
    comps = []; done = set()
    for v in sorted(cyc):
        if v in done: continue
        L = []; z = v
        while z not in done:
            done.add(z); L.append(z); z = g[z]
        rot = min("".join(tree(L[(i + j) % len(L)]) + "|" for j in range(len(L))) for i in range(len(L)))
        comps.append(f"C{len(L)}:{rot}")
    return " ".join(sorted(comps))

specials = [Fr(1, 4), Fr(3, 16), Fr(0), Fr(-3, 4), Fr(-1), Fr(-13, 16), Fr(-21, 16), Fr(-7, 4), Fr(-2), Fr(-29, 16), Fr(-37, 16), Fr(-77, 16), Fr(-93, 16), Fr(-45, 16), Fr(-5, 16)]
G = {cq: preper(cq) for cq in specials}
for cq in specials:
    print(f"   c={str(cq):<7} |PrePer|={len(G[cq]):>2}  {canon(G[cq])}   edges: " + ", ".join(f"{a}->{b}" for a, b in sorted(G[cq].items())))
chk(set(G[Fr(0)]) == {Fr(-1), Fr(0), Fr(1)} and G[Fr(0)][Fr(-1)] == 1 and G[Fr(0)][Fr(1)] == 1 and G[Fr(0)][Fr(0)] == 0, "PrePer(x^2): {-1,0,1}, 0 fixed, 1 fixed, -1->1")
chk(set(G[Fr(-1)]) == {Fr(-1), Fr(0), Fr(1)} and G[Fr(-1)][Fr(0)] == -1 and G[Fr(-1)][Fr(-1)] == 0 and G[Fr(-1)][Fr(1)] == 0, "PrePer(x^2-1): 0<->-1, 1->0")
chk(set(G[Fr(-2)]) == {Fr(k) for k in range(-2, 3)} and G[Fr(-2)][Fr(0)] == -2 and G[Fr(-2)][Fr(-2)] == 2 and G[Fr(-2)][Fr(2)] == 2 and G[Fr(-2)][Fr(1)] == -1 and G[Fr(-2)][Fr(-1)] == -1, "PrePer(x^2-2): {-2..2}, 0->-2->2 fixed, 1->-1 fixed")
chk(set(G[Fr(-7, 4)]) == {Fr(1, 2), Fr(-1, 2), Fr(3, 2), Fr(-3, 2)} and G[Fr(-7, 4)][Fr(1, 2)] == Fr(-3, 2) and G[Fr(-7, 4)][Fr(-3, 2)] == Fr(1, 2) and G[Fr(-7, 4)][Fr(-1, 2)] == Fr(-3, 2) and G[Fr(-7, 4)][Fr(3, 2)] == Fr(1, 2),
    "PrePer(x^2-7/4) = {+-1/2,+-3/2}: 2-cycle 1/2<->-3/2, -1/2->-3/2, 3/2->1/2")
chk(set(G[Fr(-3, 4)]) == set(G[Fr(-7, 4)]) and G[Fr(-3, 4)][Fr(-1, 2)] == Fr(-1, 2) and G[Fr(-3, 4)][Fr(3, 2)] == Fr(3, 2) and canon(G[Fr(-3, 4)]) != canon(G[Fr(-7, 4)]),
    "PrePer(x^2-3/4) is the same set, fixed points -1/2 and 3/2, non-isomorphic graph")
chk(set(G[Fr(-21, 16)]) == {Fr(m, 4) for m in range(-7, 8, 2)} and G[Fr(-21, 16)][Fr(7, 4)] == Fr(7, 4) and G[Fr(-21, 16)][Fr(-3, 4)] == Fr(-3, 4) and G[Fr(-21, 16)][Fr(1, 4)] == Fr(-5, 4) and G[Fr(-21, 16)][Fr(-5, 4)] == Fr(1, 4) and G[Fr(-21, 16)][Fr(-7, 4)] == Fr(7, 4),
    "PrePer(-21/16) = odd m/4, |m|<=7: fixed 7/4, -3/4; 2-cycle {1/4,-5/4}; -7/4 -> 7/4")
chk(set(G[Fr(-29, 16)]) == {Fr(m, 4) for m in range(-7, 8, 2)} and G[Fr(-29, 16)][Fr(-7, 4)] == Fr(5, 4) and G[Fr(-29, 16)][Fr(5, 4)] == Fr(-1, 4) and G[Fr(-29, 16)][Fr(-1, 4)] == Fr(-7, 4), "PrePer(-29/16) = THM-4139 (2), 3-cycle (3)")
chk(set(G[Fr(-77, 16)]) == {Fr(7, 4), Fr(-7, 4), Fr(11, 4), Fr(-11, 4)} and G[Fr(-77, 16)][Fr(-7, 4)] == Fr(-7, 4) and G[Fr(-77, 16)][Fr(11, 4)] == Fr(11, 4), "PrePer(-77/16) = {+-7/4,+-11/4}, fixed -7/4 and 11/4")
chk(set(G[Fr(1, 4)]) == {Fr(1, 2), Fr(-1, 2)}, "PrePer(1/4) = {+-1/2} (cusp: 1/2 fixed)")
chk({str(cq): len(G[cq]) for cq in (Fr(3, 16), Fr(-13, 16), Fr(-37, 16), Fr(-93, 16), Fr(-45, 16), Fr(-5, 16))} == {'3/16': 4, '-13/16': 6, '-37/16': 6, '-93/16': 4, '-45/16': 6, '-5/16': 6},
    "|PrePer| for 3/16, -13/16, -37/16, -93/16, -45/16, -5/16 = 4, 6, 6, 4, 6, 6 (the last two: {+-9/4,+-5/4,+-3/4} and {+-5/4,+-3/4,+-1/4}, beta = 9/4 and 5/4)")
chk(set(G[Fr(-45, 16)]) == {Fr(m, 4) for m in (9, -9, 5, -5, 3, -3)} and G[Fr(-45, 16)][Fr(-3, 4)] == Fr(-9, 4) and G[Fr(-45, 16)][Fr(-9, 4)] == Fr(9, 4) == G[Fr(-45, 16)][Fr(9, 4)], "-45/16: -3/4 -> -9/4 -> 9/4 fixed (tail 2)")
chk(G[Fr(-93, 16)][Fr(-7, 4)] == Fr(-11, 4) and G[Fr(-93, 16)][Fr(-11, 4)] == Fr(7, 4) and G[Fr(-93, 16)][Fr(7, 4)] == Fr(-11, 4), "-93/16: -7/4 -> -11/4 <-> 7/4 (tail 1 into the 2-cycle {-11/4, 7/4})")
classes = {}
for cq in specials[:12]:
    classes.setdefault(canon(G[cq]), []).append(str(cq))
iso = sorted(v for v in classes.values() if len(v) > 1)
chk(iso == [['-13/16', '-37/16'], ['3/16', '-3/4', '-77/16']], f"isomorphism classes among the twelve: {iso}; all other graphs distinct")
chk(all(v not in G[Fr(-7, 4)] for v in (Fr(0), Fr(-1), Fr(-2), Fr(-7, 4))), "0,-1,-2,-7/4 are not preperiodic for x^2-7/4")
# incidence among the twelve
inc = {str(xq): [str(cq) for cq in specials[:12] if xq in G[cq]] for xq in specials[:12]}
chk(inc['-7/4'] == ['-21/16', '-29/16', '-37/16', '-77/16'] and inc['1/4'] == ['3/16', '-13/16', '-21/16', '-29/16'] and inc['-3/4'] == ['3/16', '-13/16', '-21/16', '-29/16', '-37/16'],
    f"incidence among the twelve: -7/4 in {inc['-7/4']}, 1/4 in {inc['1/4']}, -3/4 in {inc['-3/4']}")

def params_for_point(xq):
    """All rational c with xq preperiodic.  Own bounds: den(c) = den(x)^2 = D^2; c <= 1/4; if |x|>=1/2 then |x|<=beta(c)
    gives c <= |x|-x^2; and x^2+c = f(x) >= -beta(c) >= -(1+sqrt(|c|)) (beta <= 1 + sqrt|c| for c<0) gives, for |c|>=4,
    |c|/2 <= |c|-sqrt|c| <= 1+x^2, so c >= -max(4, 2(1+x^2)).  Then test each candidate by the pigeonhole enumerator."""
    D = xq.denominator; den = D * D
    cmax = Fr(1, 4)
    if abs(xq) >= Fr(1, 2):
        cmax = min(cmax, abs(xq) - xq * xq)
    cmin = -max(Fr(4), 2 * (1 + xq * xq))
    out = []
    for a in range(int(cmin * den) - 1, int(cmax * den) + 2):
        if den > 1 and gcd(a, den) != 1:
            continue
        cq = Fr(a, den)
        if cq > cmax or cq < cmin:
            continue
        g = preper(cq)
        if xq in g:
            orb = [xq]; y = xq
            while True:
                y = g[y]
                if y in orb: break
                orb.append(y)
            tail = orb.index(y); out.append((cq, tail, len(orb) - tail))
    return out

L = {xq: params_for_point(xq) for xq in (Fr(1, 4), Fr(-3, 4), Fr(-7, 4), Fr(-29, 16))}
for xq, v in L.items():
    print(f"   x={str(xq):<6}: " + ", ".join(f"c={cq} ({tl},{pr})" for cq, tl, pr in v))
chk(L[Fr(-7, 4)] == [(Fr(-93, 16), 1, 2), (Fr(-77, 16), 0, 1), (Fr(-37, 16), 0, 2), (Fr(-29, 16), 0, 3), (Fr(-21, 16), 1, 1)], "parameters with -7/4 preperiodic: -93/16 (1,2), -77/16 (0,1), -37/16 (0,2), -29/16 (0,3), -21/16 (1,1)")
chk(L[Fr(1, 4)] == [(Fr(-29, 16), 1, 3), (Fr(-21, 16), 0, 2), (Fr(-13, 16), 1, 2), (Fr(-5, 16), 1, 1), (Fr(3, 16), 0, 1)], "parameters with 1/4 preperiodic: -29/16 (1,3), -21/16 (0,2), -13/16 (1,2), -5/16 (1,1), 3/16 (0,1)")
chk(L[Fr(-3, 4)] == [(Fr(-45, 16), 2, 1), (Fr(-37, 16), 1, 2), (Fr(-29, 16), 2, 3), (Fr(-21, 16), 0, 1), (Fr(-13, 16), 0, 2), (Fr(-5, 16), 2, 1), (Fr(3, 16), 1, 1)], "parameters with -3/4 preperiodic: seven, as in the note")
chk(L[Fr(-29, 16)] == [(Fr(-1561, 256), 1, 2), (Fr(-1305, 256), 0, 1), (Fr(-633, 256), 0, 2), (Fr(-377, 256), 1, 1)], "parameters with -29/16 preperiodic: -1561/256 (1,2), -1305/256 (0,1), -633/256 (0,2), -377/256 (1,1); never period 3")
chk(Fr(-29, 16) - Fr(29, 16)**2 == Fr(-1305, 256) and -1 + Fr(29, 16) - Fr(29, 16)**2 == Fr(-633, 256), "-1305/256 = x-x^2 and -633/256 = -1-x-x^2 at x=-29/16")
# exact periods via rational roots of Phi_n(x0, c) (PARI factor, linear factors)
f2 = sp.expand((x**2 + c)**2 + c); f4 = sp.expand(f3**2 + c)
Phi = {1: x**2 + c - x, 2: sp.cancel((f2 - x) / (x**2 + c - x)), 3: Phi3, 4: sp.cancel((f4 - x) / (f2 - x))}
def rational_c_roots(x0, n):
    poly = sp.Poly(sp.expand(Phi[n].subs(x, x0)), c)
    den = sp.ilcm(*[sp.Rational(v).q for v in poly.all_coeffs()])
    pstr = str(sp.expand(poly.as_expr() * den)).replace("**", "^")
    out = gp(f"f=factor({pstr}); v=[]; for(i=1,#f~, if(poldegree(f[i,1])==1, v=concat(v,[-polcoeff(f[i,1],0)/polcoeff(f[i,1],1)]))); print(v);")[0]
    return sorted(Fr(v) for v in re.findall(r"-?\d+(?:/\d+)?", out))
per = {}
for x0 in (sp.Rational(-7, 4), sp.Rational(1, 4), sp.Rational(-3, 4), sp.Rational(-29, 16)):
    xq = Fr(int(x0.p), int(x0.q)); per[xq] = {}
    for n in (1, 2, 3, 4):
        good = []
        for cq in rational_c_roots(x0, n):
            y = xq; orb = [y]
            for _ in range(n):
                y = y * y + cq; orb.append(y)
            if orb[-1] == orb[0] and len(set(orb[:-1])) == n:
                good.append(cq)
        per[xq][n] = good
    print(f"   x0={xq}: exact periods " + "; ".join(f"{n}: {per[xq][n] or 'none'}" for n in (1, 2, 3, 4)))
chk(per[Fr(-7, 4)] == {1: [Fr(-77, 16)], 2: [Fr(-37, 16)], 3: [Fr(-29, 16)], 4: []}, "-7/4 has exact period 1,2,3 for c=-77/16,-37/16,-29/16 (unique) and no rational period 4")
chk(per[Fr(1, 4)] == {1: [Fr(3, 16)], 2: [Fr(-21, 16)], 3: [], 4: []} and per[Fr(-3, 4)] == {1: [Fr(-21, 16)], 2: [Fr(-13, 16)], 3: [], 4: []}, "1/4: periods 1,2 at 3/16,-21/16; -3/4: periods 1,2 at -21/16,-13/16; no rational period 3 or 4")
chk(per[Fr(-29, 16)] == {1: [Fr(-1305, 256)], 2: [Fr(-633, 256)], 3: [], 4: []}, "-29/16: periods 1,2 at -1305/256, -633/256; no rational period 3 or 4")

# ---------------------------------------------------------------------------------------
hr("J. Chebyshev side: Fermat numerators, exact-period counts, the 3-cycles of x^2-2 and their conductors")
xq = Fr(5, 2); fer = []
for k in range(1, 6):
    fer.append(xq.numerator); xq = xq * xq - 2
chk(fer == [2**(2**k) + 1 for k in range(1, 6)] == [5, 17, 257, 65537, 4294967297], f"orbit numerators of 5/2 under x^2-2 are F_1..F_5: {fer}")
chk(all(gcd(fer[i], fer[j]) == 1 for i in range(5) for j in range(i)), "Fermat numbers pairwise coprime")
y = sp.Symbol('y')
chk(sp.simplify((y + 1/y)**2 - 2 - (y**2 + 1/y**2)) == 0, "x=y+1/y conjugates y->y^2 to x->x^2-2")
counts = [sum(sp.mobius(n // d) * 2**d for d in sp.divisors(n)) for n in range(1, 7)]
chk(counts == [2, 2, 6, 12, 30, 54], f"exact-period point counts of a degree-2 polynomial, n=1..6: {counts}")
# explicit 3-cycles of x^2-2 on 2cos(2 pi k/N): k -> 2k mod N up to sign
def cycles_mod(Nn):
    seen = set(); cyc = []
    for k in range(0, Nn // 2 + 1):
        kk = min(k % Nn, (-k) % Nn)
        if kk in seen: continue
        L = []; z = kk
        while z not in L:
            L.append(z); seen.add(z); z = min((2 * z) % Nn, (-2 * z) % Nn)
        cyc.append((L, len(L) - L.index(z)))
    return cyc
c7 = cycles_mod(7); c9 = cycles_mod(9)
chk(sorted(p for _, p in c7) == [1, 3] and sorted(p for _, p in c9) == [1, 1, 3], f"doubling on k mod 7 (mod sign): periods {[p for _, p in c7]}; mod 9: {[p for _, p in c9]} (k=3: 2cos(2pi/3)=-1 is a FIXED point)")
chk(63 == 7 * 9 == (2**3 - 1) * (2**3 + 1), "63 = 7*9 = (2^3-1)(2^3+1)")
# the 2cos values really form 3-cycles of x^2-2
import mpmath as mp
mp.mp.dps = 30
for Nn in (7, 9):
    v = 2 * mp.cos(2 * mp.pi / Nn); w = v
    for _ in range(3): w = w * w - 2
    if abs(w - v) > mp.mpf(10)**-25 or abs((v * v - 2) - v) < mp.mpf(10)**-5:
        raise RuntimeError("3-cycle check")
chk(True, "2cos(2pi/7) and 2cos(2pi/9) have exact period 3 under x^2-2 (numerically, 30 digits)")
chk(sp.simplify(2 * sp.cos(2 * sp.pi / 3)) == -1 and (-1)**2 - 2 == -1, "2cos(2pi/3) = -1 is fixed by x^2-2: the prime conductor 3 on the 2^3+1 side gives no 3-cycle")

print(f"\nALL {NCHK} AUDIT CHECKS PASSED ({time.time() - T0:.1f}s)")
