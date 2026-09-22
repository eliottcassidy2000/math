#!/usr/bin/env python3
"""
Independent audit (proof-audit lens) of lane zsigmondy_triad.
Does NOT import the explorer's script.  Recomputes every load-bearing number with
its own code, uses PARI/GP (not sympy) for discriminants/resultants/nfdisc where
possible, and probes hostile boundaries the explorer did not.
"""
import subprocess, shutil, time
from fractions import Fraction as Fr
from math import gcd, isqrt
import sympy as sp

T0 = time.time()
BAD = []
def chk(cond, msg):
    print(("  ok   " if cond else "  FAIL ") + msg)
    if not cond:
        BAD.append(msg)

def gp(cmd):
    g = shutil.which("gp")
    if not g:
        return None
    r = subprocess.run([g, "-q", "-f"], input=cmd + "\n", capture_output=True, text=True, timeout=300)
    return r.stdout.strip()

def prim_primes(terms, n):
    """primes p | terms[n] with p not dividing terms[m], 1<=m<n (terms[n]!=0)."""
    return sorted(p for p in sp.factorint(abs(terms[n])) if all(terms[m] % p != 0 for m in range(1, n)))

def prim_part(terms, n):
    """part of |terms[n]| coprime to all nonzero earlier terms (own loop)."""
    r = abs(terms[n])
    for m in range(1, n):
        t = abs(terms[m])
        if t == 0:
            continue
        g = gcd(r, t)
        while g > 1:
            r //= g
            g = gcd(r, t)
    return r

print("=== A. Bang/Zsigmondy, base 2, n<=40 ===")
M = [0] + [2**n - 1 for n in range(1, 41)]
exc = [n for n in range(1, 41) if M[n] > 1 and not prim_primes(M, n)] + ([1] if M[1] == 1 else [])
chk(sorted(exc) == [1, 6], f"exceptions = {sorted(exc)}")
# lemma: non-primitive prime of Phi_n(2), n>2, is max prime of n and divides Phi_n(2) once
for n in range(3, 41):
    phi = int(sp.cyclotomic_poly(n, 2))
    pp = prim_primes(M, n)
    for p in sp.factorint(phi):
        if p not in pp:
            chk(p == max(sp.factorint(n)) and phi % (p*p) != 0, f"n={n}: non-primitive {p} of Phi_n(2)")
P = [0] + [2**n + 1 for n in range(1, 41)]
chk([n for n in range(1, 41) if prim_part(P, n) == 1] == [3], "2^n+1 exception index = 3 only")
mods = [m for m in range(2, 64) if gcd(m, 2) == 1 and sp.n_order(2, m) == 6]
chk(mods == [9, 21, 63], f"moduli with ord_m(2)=6: {mods} (all composite, least 9)")
chk(all(not sp.isprime(m) for m in mods), "no prime has ord_p(2)=6")
chk(int(sp.cyclotomic_poly(6, 2)) == 3, "Phi_6(2)=3")

print("=== B. critical orbit at c=-7/4 ===")
def orbit_nums(c, N):
    x = Fr(0); L = [0]
    for _ in range(N):
        x = x*x + c; L.append(x.numerator)
    return L
nums = orbit_nums(Fr(-7, 4), 8)
chk(nums[1:5] == [-7, 21, -7, -114639], f"numerators {nums[1:5]}")
pp = [prim_part(nums, n) for n in range(1, 9)]
chk(pp[0] == 7 and pp[1] == 3 and pp[2] == 1 and pp[3] == 5459 and all(v > 1 for v in pp[3:]), f"primitive parts {pp[:4]}...")
chk(Fr(nums[3], 4**4) / Fr(nums[1], 4) == Fr(1, 64), "f^3(0)/f(0) = 1/64")
chk(sp.factorint(5459) == {53: 1, 103: 1}, "5459 = 53*103")

print("=== C. census (explorer box, and an extended box) ===")
def census(bs, A, N):
    fails = {n: [] for n in range(1, N+1)}; cnt = 0; pre = []
    for b in bs:
        for a in range(-A, A+1):
            if gcd(a, b) != 1:
                continue
            c = Fr(a, b); cnt += 1
            L = orbit_nums(c, N)
            if 0 in L[1:]:
                pre.append(c); continue
            for n in range(1, N+1):
                if prim_part(L, n) == 1:
                    fails[n].append(c)
    return cnt, pre, fails
cnt, pre, fails = census([1, 2, 4, 8, 16, 32, 64], 200, 8)
chk(cnt == 1601, f"box size {cnt}")
chk(sorted(pre) == [Fr(-2), Fr(-1), Fr(0)], f"0 preperiodic for {sorted(pre)}")
chk([len(fails[n]) for n in range(1, 9)] == [13, 12, 1, 0, 0, 0, 0, 0], f"failure counts {[len(fails[n]) for n in range(1,9)]}")
chk(fails[3] == [Fr(-7, 4)], "n=3 failure = -7/4 only")
chk(all(abs(c.numerator) == 1 for c in fails[1]), "n=1 failures are +-1/b")
chk(all(abs(c.numerator + c.denominator) == 1 for c in fails[2]), "n=2 failures are -1+-1/b")
cnt2, pre2, fails2 = census([1, 2, 4, 8, 16, 32, 64, 128, 256], 400, 6)
print(f"  extended box: {cnt2} params, n<=6 failure counts {[len(fails2[n]) for n in range(1,7)]}")
chk(fails2[3] == [Fr(-7, 4)] and all(fails2[n] == [] for n in (4, 5, 6)), "extended box: n=3 only -7/4, none at n=4..6")
# hostile: odd denominators (b=3,5,9,25) -- explorer's box has only powers of 2
cnt3, pre3, fails3 = census([3, 5, 7, 9, 25, 27], 150, 6)
print(f"  odd-denominator box: {cnt3} params, n<=6 failure counts {[len(fails3[n]) for n in range(1,7)]}; n=3 fails {fails3[3]}")
chk(fails3[3] == [] and all(fails3[n] == [] for n in (4, 5, 6)), "odd denominators: no n>=3 failures")

print("=== D. n<=3 characterization: exact identities ===")
a, b, c, x, s, t, X, y = sp.symbols('a b c x s t X y')
G1 = c; G2 = sp.expand(G1**2 + c); G3 = sp.expand(G2**2 + c); G4 = sp.expand(G3**2 + c)
H3 = c**3 + 2*c**2 + c + 1
chk(sp.expand(G3 - c*H3) == 0, "G_3 = c H_3")
F = a**3 + 2*a**2*b + a*b**2 + b**3
chk(sp.expand(b**3*H3.subs(c, a/b)) == sp.expand(F), "F = b^3 H_3(a/b)")
chk(sp.expand(F.subs(a, 0) - b**3) == 0 and sp.expand(F.subs(a, -b) - b**3) == 0, "F = b^3 mod a and mod (a+b)")
H4 = sp.factor(G4 / (c*(c+1)))
chk(sp.expand(H4*c*(c+1) - G4) == 0 and sp.Poly(H4, c).degree() == 6, f"G_4 = c(c+1)H_4, H_4 = {sp.expand(H4)}")
# can primes of H_4(a,b) fail to be primitive?  Res_c(H_3,H_4):
res34 = sp.resultant(H3, sp.expand(H4), c)
print(f"  Res_c(H_3,H_4) = {res34}")
H4h = sp.expand(b**6*sp.expand(H4).subs(c, a/b))
chk(sp.expand(H4h.subs(a, 0) - b**6) == 0 and sp.expand(H4h.subs(a, -b) - b**6) == 0, "H_4(a,b) coprime to a and a+b")
# search for (a,b) with gcd(H4(a,b), F(a,b)) > 1  (would make a prime of H_4 non-primitive at n=4)
wit = []
for bb in range(1, 40):
    for aa in range(-120, 121):
        if gcd(aa, bb) != 1:
            continue
        g = gcd(int(H4h.subs({a: aa, b: bb})), int(F.subs({a: aa, b: bb})))
        if g > 1:
            wit.append((aa, bb, g))
print(f"  (a,b) with gcd(H_4(a,b),F(a,b))>1, b<40,|a|<=120: {wit[:8]}{' ...' if len(wit)>8 else ''} (count {len(wit)})")
chk((res34 in (1, -1)) == (len(wit) == 0), "resultant vs witness consistency")

print("=== E. dynatomic discriminant / resultant via PARI (independent of sympy) ===")
f3 = sp.expand(((x**2 + c)**2 + c)**2 + c); f1 = x**2 + c
Phi3 = sp.expand(sp.cancel((f3 - x)/(f1 - x)))
out = gp(f"Phi3={str(Phi3).replace('**','^')}; print(factor(poldisc(Phi3,x))); print(factor(polresultant(Phi3, deriv({str(f3).replace('**','^')},x)-1, x)));")
print("  PARI:", out.replace("\n", " | ") if out else "gp absent")
if out:
    L = out.splitlines()
    chk("4*c + 7, 3" in L[0] and "16*c^2 + 4*c + 7, 2" in L[0] and L[0].startswith("[-1"), "disc_x Phi_3 = -(4c+7)^3(16c^2+4c+7)^2 (PARI)")
    chk("4*c + 7, 3" in L[1] and "16*c^2 + 4*c + 7, 3" in L[1], "Res(Phi_3,(f^3)'-1) = ((4c+7)(16c^2+4c+7))^3 (PARI)")
Pp = 8*x**3 + 4*x**2 - 18*x - 1
chk(sp.expand(64*Phi3.subs(c, sp.Rational(-7, 4)) - Pp**2) == 0, "Phi_3(x,-7/4) = P^2/64")
rts = sp.Poly(Pp, x).all_roots()
mult = sp.nsimplify(sp.prod([2*r for r in rts]))
chk(sp.simplify(mult - 1) == 0, f"multiplier of parabolic cycle = {mult}")
chk(sp.expand(64*G3 - c - c*(4*c + 7)*(16*c**2 + 4*c + 9)) == 0, "64G_3 - c = c(4c+7)(16c^2+4c+9)")
chk(sp.expand(4*G2 - c - c*(4*c + 3)) == 0, "4G_2 - c = c(4c+3)")
# 2-cycle multiplier = 4(c+1): parabolic at -3/4
chk(sp.expand(sp.resultant(x**2 + x + c + 1, 4*x*(x**2 + c) - 1, x)) == sp.expand((4*c + 3)**2), "2-cycle multiplier-1 locus = 4c+3")
f4 = sp.expand(f3**2 + c); f2 = sp.expand(f1**2 + c)
Phi4 = sp.expand(sp.cancel((f4 - x)/(f2 - x)))
out4 = gp(f"Phi4={str(Phi4).replace('**','^')}; R=polresultant(Phi4, deriv({str(f4).replace('**','^')},x)-1, x); print(factor(R)); print(gcd({str(sp.expand(2**14*G4 - c)).replace('**','^')}, R)); print(polroots(64*c^3+144*c^2+108*c+135)[1]);")
print("  PARI:", out4.replace("\n", " | ") if out4 else "gp absent")
if out4:
    L4 = out4.splitlines()
    chk("64*c^3 + 144*c^2 + 108*c + 135, 4" in L4[0] and "4*c + 5, 4" in L4[0] and "16*c^2 - 8*c + 5, 4" in L4[0], "Delta_4 factorization (PARI)")
    chk(L4[1].strip() == "1", "gcd(2^14 G_4 - c, Delta_4) = 1 (PARI)")
    chk(L4[2].startswith("-1.9405507"), f"real period-4 saddle node {L4[2][:14]}")
chk(sp.Poly(64*c**3 + 144*c**2 + 108*c + 135, c).is_irreducible, "the period-4 saddle-node cubic is irreducible (c_4 irrational)")

print("=== F. sigma-parabola from THM-4146 (10)/(12) ===")
sig = sp.symbols('sigma')
chk(sp.expand(-(sig**2 + sig + 2) - (-(sig + sp.Rational(1, 2))**2 - sp.Rational(7, 4))) == 0, "c = -7/4 - (sigma+1/2)^2 is completing the square of THM-4146 (10)")
chk(sp.expand((4*sig**2 + 6*sig + 9).subs(sig, s - sp.Rational(1, 2)) - (4*s**2 + 2*s + 7)) == 0, "THM-4146 (12) disc factor becomes 4s^2+2s+7")
lam = sp.expand(-8*(sig**3 + 2*sig**2 + 3*sig + 1)).subs(sig, s - sp.Rational(1, 2))
chk(sp.expand(lam - (1 - 14*s - 4*s**2 - 8*s**3)) == 0, "lambda(s) = 1-14s-4s^2-8s^3")
chk(lam.subs(s, 0) == 1 and lam.subs(s, -sp.Rational(1, 4)) == sp.Rational(35, 8) and lam.subs(s, sp.Rational(1, 4)) == -sp.Rational(23, 8), "lambda values")
Psig = X**3 - sig*X**2 - (sig**2 + 2*sig + 3)*X + (sig**3 + 2*sig**2 + 3*sig + 1)
chk(sp.expand(Phi3.subs({x: X, c: -(sig**2 + sig + 2)}) - Psig*Psig.subs(sig, -1 - sig)) == 0, "Phi_3 = P_sigma P_{-1-sigma}")
chk(sp.expand(Psig.subs(sig, -sp.Rational(1, 2))*8 - Pp.subs(x, X)) == 0, "P_{sigma=-1/2} = parabolic cubic /8")
# cycle-field certificates (polynomial remainders), and nfdisc via PARI
certs = [("s=0", 8*X**3 + 4*X**2 - 18*X - 1, sp.Rational(3, 2) - y**2, y**3 + y**2 - 2*y - 1, 49),
         ("s=1/4", 64*X**3 + 16*X**2 - 164*X + 23, -sp.Rational(1, 4) + y/2, y**3 - y**2 - 10*y + 8, 961),
         ("s=-1/2", X**3 + X**2 - 2*X - 1, y, y**3 + y**2 - 2*y - 1, 49),
         ("s=1/2", X**3 - 3*X + 1, y, y**3 - 3*y + 1, 81)]
for name, Pc, emb, m, dexp in certs:
    chk(sp.rem(sp.expand(Pc.subs(X, emb)), m, y) == 0 and sp.Poly(Pc, X).is_irreducible, f"{name}: embedding certificate + irreducible")
    o = gp(f"print(nfdisc({str(Pc).replace('**','^').replace('X','x')}));print(polgalois({str(Pc).replace('**','^').replace('X','x')})[1]);")
    if o:
        chk(o.splitlines()[0] == str(dexp) and o.splitlines()[1] == "3", f"{name}: nfdisc={o.splitlines()[0]}, Galois order {o.splitlines()[1]}")
chk(sp.expand(64*Psig.subs({sig: -sp.Rational(3, 4)})) == sp.expand(64*(X + sp.Rational(7, 4))*(X - sp.Rational(5, 4))*(X + sp.Rational(1, 4))), "s=-1/4 cubic splits with roots -7/4,5/4,-1/4")
chk(sp.expand(Psig.subs(sig, sp.Rational(1, 4) - sp.Rational(1, 2))*64 - (64*X**3 + 16*X**2 - 164*X + 23)) == 0, "s=1/4 cubic as stated")
# complex parabolic parameters: P_s degenerate exactly when 4s^2+2s+7=0
for sv in sp.solve(4*s**2 + 2*s + 7, s):
    cv = sp.simplify(-sp.Rational(7, 4) - sv**2)
    chk(sp.simplify(16*cv**2 + 4*cv + 7) == 0, f"degenerate s={sv} sits over complex parabolic c={cv}")

print("=== G. Thue equation F(a,b)=+-1: units, direct search, PARI ===")
rho = sp.CRootOf(t**3 - t - 1, 0)
chk(sp.expand(sp.minimal_polynomial(-rho**2, t)) == t**3 + 2*t**2 + t + 1, "minpoly(-rho^2) = H_3")
chk(sp.discriminant(t**3 - t - 1, t) == -23, "disc(t^3-t-1) = -23 (squarefree => Z[rho] maximal, unit rank 1)")
# rho^-14 = 4 rho^2 - 7 exactly: check (4 rho^2 - 7) * rho^14 == 1 in Q[t]/(t^3-t-1)
prod14 = sp.rem(sp.expand((4*t**2 - 7)*t**14), t**3 - t - 1, t)
chk(prod14 == 1, f"(4rho^2-7)*rho^14 = {prod14}")
chk(sp.rem(sp.expand((t**2 - 2)*t**5), t**3 - t - 1, t) == 1, "(rho^2-2) rho^5 = 1, i.e. -2+rho^2 = rho^-5 (sign: -2+rho^2 = +rho^-5)")
chk(sp.rem(sp.expand((t**2 - 1)*t), t**3 - t - 1, t) == 1, "rho^2-1 = rho^-1")
# units +-rho^k with zero rho-coefficient, |k|<=400, own arithmetic via sympy rem
hits = []
for k in range(-400, 401):
    e = sp.rem(sp.expand(t**k) if k >= 0 else sp.expand((t**2 - 1)**(-k)), t**3 - t - 1, t)
    pl = sp.Poly(e, t)
    if pl.coeff_monomial(t) == 0:
        hits.append((k, pl.coeff_monomial(1), pl.coeff_monomial(t**2)))
chk(hits == [(-14, -7, 4), (-5, 2, -1), (-1, -1, 1), (0, 1, 0), (2, 0, 1)], f"unit hits {hits}")
# direct search with a correct completeness argument for the window:
# |F(a,b)| = b^3 |a/b-th1| |a/b-th2|^2 >= b^3 |a/b-th1| Im(th2)^2, Im(th2)^2 = 1/rho^2 - (Re th2)^2
th = sp.Poly(t**3 + 2*t**2 + t + 1, t).all_roots()
im2 = min(float(abs(sp.im(r)))**2 for r in th if not r.is_real)
print(f"  Im(theta_2)^2 = {im2:.4f}; window |a + rho^2 b| <= 1/(im2 b^2) < 1 for b>=2")
r2 = float(rho**2)
sols = set()
for bb in range(1, 10**5 + 1):
    a0 = round(-r2*bb)
    for aa in range(a0 - 3, a0 + 4):
        if gcd(aa, bb) == 1 and aa**3 + 2*aa**2*bb + aa*bb**2 + bb**3 in (1, -1):
            sols.add((aa, bb))
chk(sols == {(-2, 1), (-1, 1), (0, 1), (-7, 4)}, f"direct Thue search b<=1e5: {sorted(sols)}")
o = gp("tnf=thueinit(x^3+2*x^2+x+1,1);print(thue(tnf,1));print(thue(tnf,-1));")
print("  PARI (flag=1):", o.replace("\n", " | ") if o else "gp absent")
if o:
    S = set()
    for line in o.splitlines():
        for pr in line.strip("[]").split("], ["):
            u, v = [int(z) for z in pr.strip("[]").split(",")]
            if v < 0 or (v == 0 and u < 0):
                u, v = -u, -v
            S.add((u, v))
    chk(S == {(1, 0), (0, 1), (-1, 1), (-2, 1), (-7, 4)}, f"PARI complete set (mod sign) {sorted(S)}")

print("=== H. Morton chart identities (THM-4139 (13)) ===")
Dt = 2*t*(t + 1)
p0 = (t**3 + 2*t**2 + t + 1)/Dt; p1 = (t**3 - t - 1)/Dt; p2 = -(t**3 + 2*t**2 + 3*t + 1)/Dt
ct = -(t**6 + 2*t**5 + 4*t**4 + 8*t**3 + 9*t**2 + 4*t + 1)/(4*t**2*(t + 1)**2)
chk(all(sp.simplify(u**2 + ct - v) == 0 for u, v in [(p0, p1), (p1, p2), (p2, p0)]), "cycle p0->p1->p2->p0")
st = sp.cancel(p0 + p1 + p2 + sp.Rational(1, 2))
chk(sp.simplify(ct + sp.Rational(7, 4) + st**2) == 0, "c(t) = -7/4 - s(t)^2")
chk(sp.factor(sp.numer(st)) == t**3 + t**2 - 2*t - 1 and sp.denom(st) == 2*t**2 + 2*t or sp.simplify(st - (t**3 + t**2 - 2*t - 1)/(2*t*(t + 1))) == 0, "s(t) as stated")
chk(sp.simplify((p1 - p2)/(p0 - p1) - t) == 0, "t = (p1-p2)/(p0-p1)")
sol29 = sp.solve(sp.numer(sp.cancel(ct + sp.Rational(29, 16))), t)
print(f"  c(t) = -29/16 at t in {sol29}")
chk(set(sp.nsimplify(v) for v in sol29 if v.is_rational) == {1, -2, sp.Rational(-1, 2)}, "rational t with c(t)=-29/16 are {1,-1/2,-2}; the irrational ones are the conductor-31 cycle")
chk(sp.factor(sp.numer(sp.cancel(st + sp.Rational(1, 2)))) == t**3 + 2*t**2 - t - 1 and sp.factor(sp.numer(sp.cancel(st - sp.Rational(1, 2)))) == t**3 - 3*t - 1, "c(t)=-2 numerators")
chk(sp.expand(Phi3.subs(c, -2) - (x**3 - 3*x + 1)*(x**3 + x**2 - 2*x - 1)) == 0, "Phi_3(x,-2) factorization")
chk(sp.expand(sp.minimal_polynomial(2*sp.cos(2*sp.pi/7), x)) == x**3 + x**2 - 2*x - 1 and sp.expand(sp.minimal_polynomial(2*sp.cos(2*sp.pi/9), x)) == x**3 - 3*x + 1, "minpolys of 2cos(2pi/7), 2cos(2pi/9)")

print("=== I. rational preperiodic graphs: own brute-force ===")
def preper(cq):
    """Brute force: den(x)^2 = den(c) forced (valuation), |x| <= beta(c) (real escape); iterate."""
    if cq > Fr(1, 4):
        return {}
    dn = cq.denominator; D = isqrt(dn)
    if D*D != dn:
        return {}
    disc = 1 - 4*cq
    # beta = (1+sqrt(disc))/2 ; |u/D| <= beta  <=>  (2|u| - D)^2 <= D^2 disc when 2|u| >= D
    cand = []
    U = D*3 + int(D*float(disc)**0.5) + 3
    for u in range(-U, U + 1):
        if gcd(u, D) != 1 and D > 1:
            continue
        au = abs(u)
        if 2*au >= D and Fr((2*au - D)**2, D*D) > disc:
            continue
        cand.append(Fr(u, D))
    cs = set(cand); g = {}
    for x0 in cand:
        z = x0; seen = []
        while z in cs and z not in seen:
            seen.append(z); z = z*z + cq
        if z in seen:
            for w in seen:
                g[w] = w*w + cq
    return g
expected = {
    Fr(1, 4): {Fr(-1, 2): Fr(1, 2), Fr(1, 2): Fr(1, 2)},
    Fr(0): {Fr(-1): Fr(1), Fr(0): Fr(0), Fr(1): Fr(1)},
    Fr(-1): {Fr(-1): Fr(0), Fr(0): Fr(-1), Fr(1): Fr(0)},
    Fr(-2): {Fr(k): Fr(k)**2 - 2 for k in range(-2, 3)},
    Fr(-3, 4): {Fr(-3, 2): Fr(3, 2), Fr(-1, 2): Fr(-1, 2), Fr(1, 2): Fr(-1, 2), Fr(3, 2): Fr(3, 2)},
    Fr(-7, 4): {Fr(-3, 2): Fr(1, 2), Fr(-1, 2): Fr(-3, 2), Fr(1, 2): Fr(-3, 2), Fr(3, 2): Fr(1, 2)},
    Fr(-21, 16): {Fr(m, 4): Fr(m, 4)**2 - Fr(21, 16) for m in (-7, -5, -3, -1, 1, 3, 5, 7)},
    Fr(-29, 16): {Fr(m, 4): Fr(m, 4)**2 - Fr(29, 16) for m in (-7, -5, -3, -1, 1, 3, 5, 7)},
}
for cq, e in expected.items():
    g = preper(cq)
    chk(g == e, f"PrePer({cq}) = {sorted(g)}")
g = preper(Fr(-29, 16))
chk(g[Fr(-7, 4)] == Fr(5, 4) and g[Fr(5, 4)] == Fr(-1, 4) and g[Fr(-1, 4)] == Fr(-7, 4), "THM-4139 (3) 3-cycle")
sizes = {str(cq): len(preper(cq)) for cq in [Fr(3, 16), Fr(-13, 16), Fr(-37, 16), Fr(-77, 16), Fr(-93, 16)]}
print(f"  |PrePer| for 3/16,-13/16,-37/16,-77/16,-93/16: {sizes}")
chk(sizes == {'3/16': 4, '-13/16': 6, '-37/16': 6, '-77/16': 4, '-93/16': 6}, "sizes as in explorer table (+ -93/16)")
# incidence: which listed parameters are preperiodic points of which
specials = [Fr(1, 4), Fr(3, 16), Fr(0), Fr(-3, 4), Fr(-1), Fr(-13, 16), Fr(-21, 16), Fr(-7, 4), Fr(-2), Fr(-29, 16), Fr(-37, 16), Fr(-77, 16)]
PP = {cq: preper(cq) for cq in specials}
inc = {str(xq): [str(cq) for cq in specials if xq in PP[cq]] for xq in specials}
chk(inc['-7/4'] == ['-21/16', '-29/16', '-37/16', '-77/16'], f"-7/4 in PrePer of {inc['-7/4']}")
chk(inc['1/4'] == ['3/16', '-13/16', '-21/16', '-29/16'], f"1/4 in PrePer of {inc['1/4']}")
chk(inc['-3/4'] == ['3/16', '-13/16', '-21/16', '-29/16', '-37/16'], f"-3/4 in PrePer of {inc['-3/4']}")
chk(inc['-29/16'] == [] and inc['-13/16'] == [] and inc['-21/16'] == [] and inc['-37/16'] == [] and inc['-77/16'] == [] and inc['3/16'] == [], "no other incidences")

print("=== J. all c with a given x preperiodic (own bounds) ===")
def params_for(xq):
    D = xq.denominator; den = D*D
    x2 = xq*xq
    cmax = Fr(1, 4) if abs(xq) < Fr(1, 2) else min(Fr(1, 4), abs(xq) - x2)
    # c >= -(1+x^2) - sqrt(1+x^2): use a safe integer bound
    lo = -(1 + x2) - Fr(isqrt(int((1 + x2)*10**12)) + 1, 10**6)
    out = []
    a_lo = (lo*den).numerator // (lo*den).denominator - 1
    a_hi = (cmax*den).numerator // (cmax*den).denominator + 1
    for aa in range(a_lo, a_hi + 1):
        if gcd(aa, den) != 1:
            continue
        cq = Fr(aa, den)
        if cq > cmax:
            continue
        g = preper(cq)
        if xq in g:
            orb = [xq]; z = xq
            while True:
                z = z*z + cq
                if z in orb:
                    break
                orb.append(z)
            tail = orb.index(z); out.append((cq, tail, len(orb) - tail))
    return out
L74 = params_for(Fr(-7, 4))
chk(L74 == [(Fr(-93, 16), 1, 2), (Fr(-77, 16), 0, 1), (Fr(-37, 16), 0, 2), (Fr(-29, 16), 0, 3), (Fr(-21, 16), 1, 1)], f"x=-7/4: {L74}")
L14 = params_for(Fr(1, 4))
chk([v[0] for v in L14] == [Fr(-29, 16), Fr(-21, 16), Fr(-13, 16), Fr(-5, 16), Fr(3, 16)], f"x=1/4: {L14}")
L34 = params_for(Fr(-3, 4))
chk([v[0] for v in L34] == [Fr(-45, 16), Fr(-37, 16), Fr(-29, 16), Fr(-21, 16), Fr(-13, 16), Fr(-5, 16), Fr(3, 16)], f"x=-3/4: {L34}")
L29 = params_for(Fr(-29, 16))
chk([v[0] for v in L29] == [Fr(-1561, 256), Fr(-1305, 256), Fr(-633, 256), Fr(-377, 256)] and max(v[2] for v in L29) == 2, f"x=-29/16: {L29}")
chk(Fr(-29, 16) - Fr(29, 16)**2 == Fr(-1305, 256) and -1 - Fr(-29, 16) - Fr(29, 16)**2 == Fr(-633, 256), "fixed/2-cycle formulas at -29/16")
# exact-period parameters for -7/4 by direct formulas: fixed c=x-x^2, 2-cycle c=-1-x-x^2
chk(Fr(-7, 4) - Fr(49, 16) == Fr(-77, 16) and -1 + Fr(7, 4) - Fr(49, 16) == Fr(-37, 16), "-7/4 fixed at -77/16, 2-cycle at -37/16")
# is -7/4 of exact period 3 for any c other than -29/16?  Phi_3(-7/4, c) rational roots
r3 = sp.roots(sp.Poly(Phi3.subs(x, sp.Rational(-7, 4)), c), filter='Q')
chk(set(r3) == {sp.Rational(-29, 16)}, f"rational c with Phi_3(-7/4,c)=0: {list(r3)}")
r4 = sp.roots(sp.Poly(Phi4.subs(x, sp.Rational(-7, 4)), c), filter='Q')
chk(len(r4) == 0, "no rational c with Phi_4(-7/4,c)=0")

print("=== K. Chebyshev side ===")
xq = Fr(5, 2); fn = []
for k in range(1, 6):
    fn.append(xq.numerator); xq = xq*xq - 2
chk(fn == [2**(2**k) + 1 for k in range(1, 6)], f"Fermat numerators {fn}")
chk(sp.expand(sp.cyclotomic_poly(7, y)*sp.cyclotomic_poly(9, y)*(y - 1)*(y + 1) * 1) != 0, "placeholder")
# period-3 points of y->y^2: y^8=y => y^7=1 (with y^2=y^-1 side: y^9=1); exact period 3 excludes y^1=y
chk(sp.expand((y**8 - y)/(y**2 - y)) == sp.expand(sp.cyclotomic_poly(7, y)) and sp.expand((y**8 - 1/y)*y/(y**2 - 1)) == sp.expand(sp.cyclotomic_poly(9, y)*sp.cyclotomic_poly(3, y)), "period-3 conductors of squaring: 7 and 9 (3 gives the fixed point -1)")
chk(2**3 - 1 == 7 and 2**3 + 1 == 9 and 7*9 == 63, "63 = 7*9")
chk([n for n in range(1, 8) if 2**n - 2 == 2*n] == [3], "2^n-2 = 2n iff n=3")

print()
print(f"ALL DONE in {time.time()-T0:.1f}s; failures: {BAD}")
