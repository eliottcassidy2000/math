#!/usr/bin/env python3
"""
collatz_mod6_20260917_zsigmondy_triad_audit_recompute.py
Adversarial verifier (lane zsigmondy_triad, lens recompute).  INDEPENDENT recomputation: does not import
the explorer's script.  Exact integer / Fraction / sympy arithmetic; sympy used only as an algebra engine.
Every load-bearing check is an explicit `raise`-based assertion (active under -O).
"""
import sys, time
from fractions import Fraction as Fr
from math import gcd, isqrt
import sympy as sp
from sympy import Rational as R, symbols, expand, factor, discriminant, resultant, Poly, gcd as spgcd

T0 = time.time()
BAD = []
def chk(cond, msg):
    print(("  OK   " if cond else "  FAIL ") + msg)
    if not cond:
        BAD.append(msg)

def primes_of(n):
    return sorted(sp.factorint(abs(n)).keys()) if abs(n) > 1 else []

def prim_primes(terms, n):
    """primitive primes of terms[n]: primes of terms[n] dividing no terms[m], 1<=m<n (terms[m]=0 treated as divisible by all)."""
    if terms[n] == 0:
        return None
    out = []
    for p in primes_of(terms[n]):
        if all(terms[m] % p != 0 for m in range(1, n)):
            out.append(p)
    return out

def has_prim_gcdfree(terms, n):
    """factorization-free: does terms[n] have a prime coprime to all earlier terms?"""
    x = abs(terms[n])
    if x == 0:
        return None
    for m in range(1, n):
        t = abs(terms[m])
        if t == 0:
            return False
        while True:
            g = gcd(x, t)
            if g == 1:
                break
            x //= g
    return x != 1

print("=== A. Bang/Zsigmondy for 2^n-1 and 2^n+1, n<=40 ===")
mers = [0] + [2**n - 1 for n in range(1, 41)]
exc = [n for n in range(1, 41) if not prim_primes(mers, n)]
chk(exc == [1, 6], f"2^n-1 exceptions n<=40 = {exc}")
# cyclotomic values Phi_n(2) by Moebius product (independent of sympy.cyclotomic_poly)
def phi_n_2(n):
    num = 1; den = 1
    for d in sp.divisors(n):
        mu = sp.mobius(n // d)
        if mu == 1: num *= 2**d - 1
        elif mu == -1: den *= 2**d - 1
    assert num % den == 0
    return num // den
nonprim = {}
for n in range(3, 41):
    ph = phi_n_2(n)
    pp = prim_primes(mers, n)
    npp = [p for p in primes_of(ph) if p not in pp]
    if npp:
        nonprim[n] = npp
        for p in npp:
            chk(p == max(primes_of(n)) and ph % (p*p) != 0, f"n={n}: non-primitive prime {p} of Phi_n(2) is largest prime of n, exponent 1")
chk(nonprim == {6: [3], 18: [3], 20: [5], 21: [7]}, f"non-primitive primes of Phi_n(2), 3<=n<=40: {nonprim}")
chk([phi_n_2(d) for d in (1, 2, 3, 6)] == [1, 3, 7, 3], "Phi_d(2), d|6 = 1,3,7,3")
# consistency of gcd-free test with factorization for all n
chk(all(has_prim_gcdfree(mers, n) == bool(prim_primes(mers, n)) for n in range(1, 41)), "gcd-free primitive test agrees with factorization, 2^n-1")
# ord_p(2)=6 impossible: proof-level (p|63, p not | 3, 7) + sweep
chk([p for p in sp.primerange(3, 10**5) if sp.n_order(2, p) == 6] == [], "no prime p<10^5 with ord_p(2)=6")
chk([m for m in range(1, 200, 2) if sp.n_order(2, m) == 6] == [9, 21, 63], "odd moduli m<200 with ord_m(2)=6 are exactly 9,21,63")
plus = [0] + [2**n + 1 for n in range(1, 41)]
chk([n for n in range(1, 41) if not prim_primes(plus, n)] == [3], "2^n+1 exceptions n<=40 = [3]")
chk([n for n in range(1, 10) if 2**n - 2 == 2*n] == [3], "2^n-2=2n iff n=3 (n<10; for n>=4 LHS>RHS by induction)")

print("=== B. Critical orbit of x^2-7/4; census recount; hostile extended census ===")
def orbit_nums(c, nmax):
    x = Fr(0); nums = [0]; orb = [Fr(0)]
    for n in range(nmax):
        x = x*x + c; nums.append(x.numerator); orb.append(x)
    return nums, orb
nums, orb = orbit_nums(Fr(-7, 4), 8)
chk(orb[1:5] == [Fr(-7, 4), Fr(21, 16), Fr(-7, 256), Fr(-114639, 65536)], "orbit values n=1..4")
chk(orb[3] / orb[1] == Fr(1, 64), "f^3(0)/f(0) = 1/64")
chk([abs(nums[n]) for n in range(1, 5)] == [7, 21, 7, 114639], "numerators n<=4")
chk(sp.factorint(114639) == {3: 1, 7: 1, 53: 1, 103: 1}, "114639 = 3*7*53*103")
chk(sp.factorint(abs(nums[5])) == {7: 1, 419: 1, 563: 1, 3407: 1}, f"N_5 = 7*419*563*3407 (N_5={nums[5]})")
chk(sp.factorint(abs(nums[6])) == {3: 1, 7: 1, 47237: 1, 636069519847: 1}, "N_6 = 3*7*47237*636069519847")
hp = [has_prim_gcdfree(nums, n) for n in range(1, 9)]
chk(hp == [True, True, False, True, True, True, True, True], f"primitive prime present n=1..8: {hp}")
chk(prim_primes(nums, 4) == [53, 103] and prim_primes(nums, 2) == [3], "primitive primes at n=2,4")
# denominators: 4^(2^(n-1))
chk(all(orb[n].denominator == 4**(2**(n-1)) for n in range(1, 9)), "denominators 4^(2^(n-1)), no cancellation")

# census recount, same box
box_b = [1, 2, 4, 8, 16, 32, 64]
fails = {n: [] for n in range(1, 9)}; pre = []; cnt = 0
for b in box_b:
    for a in range(-200, 201):
        if gcd(a, b) != 1:
            continue
        cnt += 1
        c = Fr(a, b)
        nm, ob = orbit_nums(c, 8)
        if 0 in nm[1:] or len(set(ob)) < 9:
            pre.append(c); continue
        for n in range(1, 9):
            if not has_prim_gcdfree(nm, n):
                fails[n].append(c)
chk(cnt == 1601, f"box size {cnt}")
chk(sorted(pre) == [Fr(-2), Fr(-1), Fr(0)], f"0 preperiodic (orbit revisits within 8 steps) exactly for {sorted(pre)}")
chk(len(fails[1]) == 13 and len(fails[2]) == 12 and fails[3] == [Fr(-7, 4)] and all(fails[n] == [] for n in range(4, 9)),
    f"census counts: {[len(fails[n]) for n in range(1, 9)]}")
chk(sorted(fails[1]) == sorted(set(Fr(s, b) for b in box_b for s in (1, -1)) - {Fr(-1)}), "n=1 failures = +-1/b minus preperiodic -1")
chk(sorted(fails[2]) == sorted(set(Fr(-b + s, b) for b in box_b for s in (1, -1)) - {Fr(0), Fr(-2)}), "n=2 failures = -1+-1/b minus preperiodic")

# hostile extension: denominators not powers of 2, larger |a|, n<=6
box2 = [3, 5, 6, 7, 9, 10, 12, 15, 25, 27, 49, 100, 128, 256]
ext = {n: [] for n in range(1, 7)}
for b in box2:
    for a in range(-400, 401):
        if gcd(a, b) != 1:
            continue
        c = Fr(a, b)
        nm, ob = orbit_nums(c, 6)
        for n in range(1, 7):
            if not has_prim_gcdfree(nm, n):
                ext[n].append(c)
chk(all(c.numerator in (1, -1) for c in ext[1]) and len(ext[1]) == 2*len(box2), "extended box n=1 failures are exactly +-1/b")
chk(all(abs(c.numerator + c.denominator) == 1 for c in ext[2]) and len(ext[2]) == 2*len(box2), "extended box n=2 failures are exactly -1+-1/b")
chk(ext[3] == [] and all(ext[n] == [] for n in range(4, 7)), f"extended box (b in {box2}, |a|<=400): no n=3..6 failures: {[len(ext[n]) for n in range(3,7)]}")
# hostile: b=4 with large |a|, n=3 (Thue) - up to |a|<=5000
f3 = [Fr(a, 4) for a in range(-5000, 5001, 2) if gcd(a, 4) == 1 and not has_prim_gcdfree(orbit_nums(Fr(a, 4), 3)[0], 3)]
chk(f3 == [Fr(-7, 4)], f"b=4, |a|<=5000: n=3 failures {f3}")

print("=== C. Symbolic dynatomic / Gleason identities (sympy, own code) ===")
x, c, s, sg, t, X = symbols('x c s sigma t X')
f = lambda z: z**2 + c
F1 = f(x); F2 = expand(f(F1)); F3 = expand(f(F2)); F4 = expand(f(F3))
q, r = sp.div(Poly(F3 - x, x), Poly(F1 - x, x))
chk(r.is_zero, "(f^3-x)/(f-x) exact")
Phi3 = q.as_expr()
chk(Poly(Phi3, x).degree() == 6, "deg_x Phi_3 = 6")
d3 = factor(discriminant(Phi3, x))
chk(expand(d3 + (4*c + 7)**3*(16*c**2 + 4*c + 7)**2) == 0, f"disc_x Phi_3 = -(4c+7)^3(16c^2+4c+7)^2   [got {d3}]")
r3 = factor(resultant(Phi3, sp.diff(F3, x) - 1, x))
chk(expand(r3 - ((4*c + 7)*(16*c**2 + 4*c + 7))**3) == 0, f"Res_x(Phi_3,(f^3)'-1) = ((4c+7)(16c^2+4c+7))^3  [got {r3}]")
# 16c^2+4c+7 = 0 is the satellite (rotation 1/3) locus: c = mu/2 - mu^2/4 with mu^3=1, mu != 1
mu = symbols('mu')
csat = mu/2 - mu**2/4
chk(sp.rem(expand((16*csat**2 + 4*csat + 7)), mu**2 + mu + 1, mu) == 0, "16c^2+4c+7=0 <=> c = mu/2-mu^2/4 with mu a primitive cube root of unity (satellite, NOT saddle-node)")
Pp = 8*x**3 + 4*x**2 - 18*x - 1
chk(expand(64*Phi3.subs(c, R(-7, 4)) - Pp**2) == 0, "Phi_3(x,-7/4) = Pp^2/64")
chk(Poly(Pp, x).is_irreducible, "Pp irreducible over Q")
cf = Poly(Pp, x).all_coeffs(); chk(-cf[3]/cf[0] == R(1, 8), "prod roots Pp = 1/8, multiplier 8*(1/8)=1")
rts = sorted(Poly(Pp, x).nroots(n=12), key=lambda z: sp.re(z))
print("   parabolic cycle points sorted:", rts)
o = [rts[0]]
for _ in range(3):
    o.append(sp.N(o[-1]**2 - R(7, 4), 12))
print("   dynamical order from x_min:", [sp.N(v, 8) for v in o])
chk(abs(o[3] - o[0]) < 1e-8 and abs(o[1] - rts[2]) < 1e-8 and abs(o[2] - rts[1]) < 1e-8, "cycle order is -1.747 -> 1.302 -> -0.055 -> (note's tuple is SORTED, not dynamical order)")
# sigma parabola
Ps = X**3 - sg*X**2 - (sg**2 + 2*sg + 3)*X + (sg**3 + 2*sg**2 + 3*sg + 1)
csg = -(sg**2 + sg + 2)
chk(expand(Phi3.subs({x: X, c: csg}) - expand(Ps*Ps.subs(sg, -1 - sg))) == 0, "Phi_3(X,c(sigma)) = P_sigma P_{-1-sigma}")
lam = expand(-8*(sg**3 + 2*sg**2 + 3*sg + 1))
# own check of multiplier: product of 2*roots = 8*(-const term of P_sigma with sign) via Vieta
chk(expand(lam - 8*(-(sg**3 + 2*sg**2 + 3*sg + 1))) == 0, "lambda = 8*prod roots by Vieta")
lam_s = expand(lam.subs(sg, s - R(1, 2))); c_s = expand(csg.subs(sg, s - R(1, 2)))
chk(c_s == -s**2 - R(7, 4) and lam_s == -8*s**3 - 4*s**2 - 14*s + 1, f"c=-7/4-s^2, lambda(s)={lam_s}")
chk(lam_s.subs(s, 0) == 1 and lam_s.subs(s, -R(1, 4)) == R(35, 8) and lam_s.subs(s, R(1, 4)) == -R(23, 8), "lambda(0)=1, lambda(-1/4)=35/8, lambda(1/4)=-23/8")
ds = factor(discriminant(Ps.subs(sg, s - R(1, 2)), X))
chk(expand(ds - (4*s**2 + 2*s + 7)**2) == 0, f"disc P_s = (4s^2+2s+7)^2  [got {ds}]")
eta = symbols('eta')
chk(expand((4*s**2 + 2*s + 7).subs(s, eta/2) - (eta**2 + eta + 7)) == 0 and expand(c_s.subs(s, eta/2) + (eta**2 + 7)/4) == 0,
    "s = eta/2: parabola and disc are geometry.md (4) [c=-(eta^2+7)/4, disc=(eta^2+eta+7)^2] verbatim (INHERITED)")
# Gleason
G = [sp.Integer(0)]
for n in range(1, 5):
    G.append(expand(G[-1]**2 + c))
H3 = c**3 + 2*c**2 + c + 1
chk(expand(G[3] - c*H3) == 0 and expand(G[2] - c*(c + 1)) == 0, "G_2=c(c+1), G_3=cH_3")
chk(sp.rem(G[4], c*(c + 1), c) == 0, "c(c+1) | G_4")
H4 = sp.quo(G[4], c*(c + 1), c)
chk(Poly(H4, c).is_irreducible and Poly(H3, c).is_irreducible, "H_3, H_4 irreducible")
chk(expand(64*G[3] - c - c*(4*c + 7)*(16*c**2 + 4*c + 9)) == 0, "64G_3-c = c(4c+7)(16c^2+4c+9)")
chk(expand(4*G[2] - c - c*(4*c + 3)) == 0, "4G_2-c = c(4c+3)")
q4, r4 = sp.div(Poly(F4 - x, x), Poly(F2 - x, x)); chk(r4.is_zero, "(f^4-x)/(f^2-x) exact")
Phi4 = q4.as_expr()
res4 = factor(resultant(Phi4, sp.diff(F4, x) - 1, x))
D4p = 64*c**3 + 144*c**2 + 108*c + 135
chk(expand(res4 - ((4*c + 5)*(16*c**2 - 8*c + 5)*D4p)**4) == 0, f"Res_x(Phi_4,(f^4)'-1) = ((4c+5)(16c^2-8c+5)(64c^3+144c^2+108c+135))^4  [got {res4}]")
Delta4 = (4*c + 5)*(16*c**2 - 8*c + 5)*D4p
g4 = spgcd(Poly(expand(2**14*G[4] - c), c), Poly(expand(Delta4), c))
chk(g4.degree() == 0, f"gcd(2^14 G_4 - c, Delta_4) = {g4.as_expr()}")
# also directly: 2^14 G_4(c) - c evaluated at the three real parabolic period-4 params is nonzero
chk(expand(2**14*G[4] - c).subs(c, -R(5, 4)) != 0, "2^14 G_4(-5/4) != -5/4 (satellite period-4 parameter also fails)")
rr4 = [z for z in Poly(D4p, c).nroots(n=15) if abs(sp.im(z)) < 1e-12]
chk(len(rr4) == 1 and abs(rr4[0] + 1.940550788982) < 1e-9, f"real saddle-node c_4 = {rr4}")
chk(expand(Phi3.subs(c, -2) - (x**3 - 3*x + 1)*(x**3 + x**2 - 2*x - 1)) == 0, "Phi_3(x,-2) = q9*q7 (INHERITED geometry.md (6))")
# Morton chart (THM-4139 (13)) identities
Dt = 2*t*(t + 1)
p0 = (t**3 + 2*t**2 + t + 1)/Dt; p1 = (t**3 - t - 1)/Dt; p2 = -(t**3 + 2*t**2 + 3*t + 1)/Dt
ct = -(t**6 + 2*t**5 + 4*t**4 + 8*t**3 + 9*t**2 + 4*t + 1)/(4*t**2*(t + 1)**2)
chk(all(sp.cancel(a_**2 + ct - b_) == 0 for a_, b_ in [(p0, p1), (p1, p2), (p2, p0)]), "Morton chart is a 3-cycle")
sigt = sp.cancel(p0 + p1 + p2); st = sp.cancel(sigt + R(1, 2))
chk(sp.cancel(ct + R(7, 4) + st**2) == 0, "c(t) = -7/4 - s(t)^2")
chk(sp.cancel(sigt.subs(t, -1/(t + 1)) - sigt) == 0, "sigma(t) invariant under t -> -1/(t+1)")
chk(sp.cancel((p1 - p2)/(p0 - p1) - t) == 0 and sp.cancel(p0 + p1 - t) == 0, "t = (p1-p2)/(p0-p1) = p0+p1 (INHERITED geometry.md (5a))")
chk(sp.factor(sp.numer(sp.cancel(st + R(1, 2)))) == t**3 + 2*t**2 - t - 1 and sp.factor(sp.numer(sp.cancel(st - R(1, 2)))) == t**3 - 3*t - 1, "s(t)=-+1/2 numerators")
chk(ct.subs(t, 1) == -R(29, 16) and [p0.subs(t, 1), p1.subs(t, 1), p2.subs(t, 1)] == [R(5, 4), -R(1, 4), -R(7, 4)], "t=1 -> -29/16, (5/4,-1/4,-7/4)")
# rational-s cycle field: square disc => Q or cyclic cubic; the s=-1/4 cubic splits, s=1/4 irreducible
P14 = expand(64*Ps.subs(sg, -R(1, 4))); Pm14 = expand(64*Ps.subs(sg, -R(3, 4)))
chk(P14 == 64*X**3 + 16*X**2 - 164*X + 23 and Poly(P14, X).is_irreducible, "s=1/4 cubic = 64X^3+16X^2-164X+23 irreducible")
chk(sorted(sp.roots(Poly(Pm14, X), filter='Q').keys()) == [-R(7, 4), -R(1, 4), R(5, 4)], "s=-1/4 cubic splits with roots -7/4,-1/4,5/4")
chk(expand(Ps.subs(sg, -1)) == X**3 + X**2 - 2*X - 1 and expand(Ps.subs(sg, 0)) == X**3 - 3*X + 1, "sigma=-1,0 cubics = q7, q9")
# embeddings (exact remainder)
y = symbols('y')
chk(sp.rem(expand(Pp.subs(x, R(3, 2) - y**2)), y**3 + y**2 - 2*y - 1, y) == 0, "parabolic root = 3/2 - y^2, y=2cos(2pi/7)")
chk(sp.rem(expand(P14.subs(X, -R(1, 4) + y/2)), y**3 - y**2 - 10*y + 8, y) == 0, "s=1/4 root = -1/4 + y/2, y^3-y^2-10y+8=0")
chk(expand((y**3 + y**2 - 10*y - 8).subs(y, -y)*(-1)) == y**3 - y**2 - 10*y + 8, "y^3-y^2-10y+8 is the conductor-31 cubic x^3+x^2-10x-8 under y->-y")
chk(sp.rem(expand(Pp.subs(x, -y - R(1, 2))), y**3 + y**2 - 2*y - 1, y) == 0, "parabolic cycle = L(q7 cycle), L=-y-1/2 (INHERITED geometry.md (8))")
# plastic number
rho = sp.CRootOf(t**3 - t - 1, 0)
chk(sp.minimal_polynomial(-rho**2, t) == t**3 + 2*t**2 + t + 1, "minpoly(-rho^2) = H_3")
chk(abs(float(-rho**2) + 1.7548776662466927) < 1e-12 and abs(float(-rho**2) + 1.75) < 0.0049, "airplane centre -rho^2 = -1.75487766...")
a_, b_ = symbols('a b')
Fab = a_**3 + 2*a_**2*b_ + a_*b_**2 + b_**3
# norm of a + b rho^2 = prod over conjugates
conj = [sp.CRootOf(t**3 - t - 1, i) for i in range(3)]
Nexpr = sp.expand(sp.simplify(sp.prod([a_ + b_*rr**2 for rr in conj])))
chk(sp.simplify(sp.nsimplify(Nexpr, rational=True) - Fab) == 0 or all(abs(complex(Nexpr.subs({a_: A, b_: B})) - Fab.subs({a_: A, b_: B})) < 1e-9 for A in range(-3, 4) for B in range(-3, 4)),
    "F(a,b) = N(a + b rho^2) (checked on the box |a|,|b|<=3 numerically, exact rational disc -23)")
chk(sp.discriminant(t**3 - t - 1, t) == -23, "disc(t^3-t-1) = -23")
chk(sp.rem(Fab.subs(a_, 0), b_, b_) == 0 and expand(Fab.subs(a_, -b_)) == b_**3, "F = b^3 mod (a+b); F(0,b)=b^3 i.e. F = b^3 mod a")
# rho powers in basis (1, rho, rho^2): own recurrence via companion matrix
M = sp.Matrix([[0, 0, 1], [1, 0, 1], [0, 1, 0]])   # multiplication by rho on coords (1,rho,rho^2): rho*rho^2 = rho^3 = 1 + rho
Minv = M.inv()
def rhopow(k):
    v = sp.Matrix([1, 0, 0])
    Mk = (M if k >= 0 else Minv)**abs(k)
    return list(Mk*v)
hits = [(k, v[0], v[2]) for k in range(-600, 601) for v in [rhopow(k)] if v[1] == 0]
chk(hits == [(-14, -7, 4), (-5, 2, -1), (-1, -1, 1), (0, 1, 0), (2, 0, 1)], f"rho^k with vanishing rho-coefficient, |k|<=600: {hits}")
chk(rhopow(-14) == [-7, 0, 4], "rho^-14 = -7 + 4 rho^2")
# direct Thue search (independent bound): |F(a,b)| = |b|^3 |prod (a/b + rho_i^2)|; for |a/b + rho^2| >= 1/b, |F| >= b^2 * |a/b+rho'^2|^2 ... just brute force a in wide window
sols = set()
r2 = float(rho**2)
for b in range(1, 200001):
    a0 = round(-r2*b)
    for a in range(a0 - 3, a0 + 4):
        if gcd(a, b) == 1 and a**3 + 2*a*a*b + a*b*b + b**3 in (1, -1):
            sols.add((a, b))
chk(sols == {(-7, 4), (-1, 1), (0, 1), (-2, 1)}, f"Thue F=+-1, b<=2e5, a near -rho^2 b: {sorted(sols)}")
# window justification: other conjugates are complex with |a/b + rho_i^2| >= Im(rho_i^2) > 0.5 (numeric), so |F| >= b^3 * |a/b+rho^2| * 0.25 >= ... for |a - (-rho^2 b)| >= 4 this exceeds 1 for b>=1
im2 = abs(float(sp.im(conj[1]**2)))
chk(im2 > 0.5, f"|Im(rho_i^2)| = {im2:.4f} > 0.5 so |F(a,b)| >= |a + rho^2 b| * im2^2 * b^2 >= 4*0.25*b^2 > 1 outside the window a within 3 of -rho^2 b")

print("=== D. Complete rational preperiodic sets: independent brute force with WIDE search windows ===")
def preper_bruteforce(cq, mult=3, steps=300):
    """Candidates: all u/D with den(c)=D^2 (valuation lemma), |u| <= mult*D*(beta+1); a point is preperiodic iff its
    orbit revisits; escape when |x| > beta (exact) -- but we also run a float magnitude fuse as a hostile cross-check."""
    if cq > Fr(1, 4):
        return {}
    den = cq.denominator; D = 1
    for p, e in sp.factorint(den).items():
        if e % 2: return {}
        D *= p**(e//2)
    beta = (1 + (1 - 4*float(cq))**0.5)/2
    U = int(mult*D*(beta + 1)) + 5
    g = {}
    for u in range(-U, U + 1):
        x0 = Fr(u, D); seen = []; z = x0; ok = False
        for _ in range(steps):
            if z in seen: ok = True; break
            seen.append(z); z = z*z + cq
            if abs(z) > 1e6: break
        if ok:
            for w in seen: g[w] = w*w + cq
    return g
specials = [Fr(1,4), Fr(3,16), Fr(0), Fr(-3,4), Fr(-1), Fr(-13,16), Fr(-21,16), Fr(-7,4), Fr(-2), Fr(-29,16), Fr(-37,16), Fr(-77,16)]
PP = {cq: preper_bruteforce(cq) for cq in specials}
expected = {
    Fr(1,4): {Fr(1,2): Fr(1,2), Fr(-1,2): Fr(1,2)},
    Fr(0): {Fr(0): Fr(0), Fr(1): Fr(1), Fr(-1): Fr(1)},
    Fr(-1): {Fr(0): Fr(-1), Fr(-1): Fr(0), Fr(1): Fr(0)},
    Fr(-2): {Fr(k): Fr(k*k - 2) for k in range(-2, 3)},
    Fr(-7,4): {Fr(1,2): Fr(-3,2), Fr(-3,2): Fr(1,2), Fr(-1,2): Fr(-3,2), Fr(3,2): Fr(1,2)},
    Fr(-3,4): {Fr(1,2): Fr(-1,2), Fr(-1,2): Fr(-1,2), Fr(3,2): Fr(3,2), Fr(-3,2): Fr(3,2)},
    Fr(-29,16): {Fr(m,4): Fr(m*m - 29, 16) for m in (-7,-5,-3,-1,1,3,5,7)},
    Fr(-21,16): {Fr(m,4): Fr(m*m - 21, 16) for m in (-7,-5,-3,-1,1,3,5,7)},
    Fr(-77,16): {Fr(m,4): Fr(m*m - 77, 16) for m in (-11,-7,7,11)},
    Fr(-37,16): {Fr(m,4): Fr(m*m - 37, 16) for m in (-7,-5,-3,3,5,7)},
    Fr(-13,16): {Fr(m,4): Fr(m*m - 13, 16) for m in (-5,-3,-1,1,3,5)},
    Fr(3,16): {Fr(m,4): Fr(m*m + 3, 16) for m in (-3,-1,1,3)},
}
for cq in specials:
    chk(PP[cq] == expected[cq], f"PrePer({cq}) = {sorted(PP[cq])}")
# graph-shape facts
def cycles_and_trees(g):
    incyc = set()
    for v in g:
        z = v
        for _ in range(len(g)+1): z = g[z]
        w = z
        while True:
            incyc.add(w); w = g[w]
            if w == z: break
    pre = {}
    for a, b in g.items(): pre.setdefault(b, []).append(a)
    def tree(v): return tuple(sorted(tree(w) for w in pre.get(v, []) if w not in incyc))
    comps = []; done = set()
    for v in sorted(incyc):
        if v in done: continue
        cyc = []; z = v
        while z not in done: done.add(z); cyc.append(z); z = g[z]
        comps.append((len(cyc), tuple(sorted(tree(w) for w in cyc))))
    return tuple(sorted(comps))
T = {cq: cycles_and_trees(PP[cq]) for cq in specials}
chk(T[Fr(-7,4)] != T[Fr(-3,4)] and set(PP[Fr(-7,4)]) == set(PP[Fr(-3,4)]), "-7/4 vs -3/4: same set, non-isomorphic graphs")
chk(T[Fr(3,16)] == T[Fr(-3,4)] == T[Fr(-77,16)] and T[Fr(-13,16)] == T[Fr(-37,16)], "type coincidences {3/16,-3/4,-77/16}, {-13/16,-37/16}")
others = [Fr(1,4), Fr(0), Fr(-1), Fr(-21,16), Fr(-7,4), Fr(-2), Fr(-29,16)]
chk(len({T[cq] for cq in others} | {T[Fr(3,16)], T[Fr(-13,16)]}) == 9, "all other listed types distinct")
chk(Fr(0) not in PP[Fr(-7,4)] and Fr(-1) not in PP[Fr(-7,4)] and Fr(-2) not in PP[Fr(-7,4)] and Fr(-7,4) not in PP[Fr(-7,4)], "0,-1,-2,-7/4 not preperiodic for -7/4")

print("=== E. Parameters for which a given point is preperiodic: WIDE brute force over c=a/D^2 ===")
def params_bruteforce(xq, arange=(-6000, 200), steps=400):
    D = xq.denominator; den = D*D; out = []
    for a in range(arange[0], arange[1] + 1):
        if den > 1 and gcd(a, den) != 1: continue
        cq = Fr(a, den); seen = []; z = xq; ok = False
        for _ in range(steps):
            if z in seen: ok = True; break
            seen.append(z); z = z*z + cq
            if abs(z) > 1e9: break
        if ok:
            tail = seen.index(z); out.append((cq, tail, len(seen) - tail))
    return out
E = {}
for xq in [Fr(1,4), Fr(-3,4), Fr(-7,4), Fr(-29,16), Fr(0), Fr(-1), Fr(-2)]:
    E[xq] = params_bruteforce(xq)
    print("   x=%s: %s" % (xq, ", ".join(f"c={cq}(t{tl},p{pr})" for cq, tl, pr in E[xq])))
chk(E[Fr(-7,4)] == [(Fr(-93,16),1,2),(Fr(-77,16),0,1),(Fr(-37,16),0,2),(Fr(-29,16),0,3),(Fr(-21,16),1,1)], "params with -7/4 preperiodic")
chk([cq for cq,_,_ in E[Fr(1,4)]] == [Fr(-29,16),Fr(-21,16),Fr(-13,16),Fr(-5,16),Fr(3,16)], "params with 1/4 preperiodic")
chk([cq for cq,_,_ in E[Fr(-3,4)]] == [Fr(-45,16),Fr(-37,16),Fr(-29,16),Fr(-21,16),Fr(-13,16),Fr(-5,16),Fr(3,16)], "params with -3/4 preperiodic")
chk([cq for cq,_,_ in E[Fr(-29,16)]] == [Fr(-1561,256),Fr(-1305,256),Fr(-633,256),Fr(-377,256)] and max(pr for _,_,pr in E[Fr(-29,16)]) == 2, "params with -29/16 preperiodic; max period 2")
chk([cq for cq,_,_ in E[Fr(0)]] == [Fr(-2),Fr(-1),Fr(0)], "params with 0 preperiodic")
# exact-period parameters by rational roots of Phi_n(x0,c)
Phi1 = F1 - x; Phi2 = sp.quo(Poly(F2 - x, x), Poly(F1 - x, x)).as_expr()
def exact_period_params(x0, n, Ph):
    poly = Poly(expand(Ph.subs(x, x0)), c)
    good = []
    for rt in sp.roots(poly, filter='Q'):
        cq = Fr(int(rt.p), int(rt.q)); z = Fr(int(x0.p), int(x0.q)); ob = [z]
        for _ in range(n): z = z*z + cq; ob.append(z)
        if ob[-1] == ob[0] and len(set(ob[:-1])) == n: good.append(cq)
    return sorted(good)
per74 = [exact_period_params(-R(7,4), n, Ph) for n, Ph in [(1,Phi1),(2,Phi2),(3,Phi3),(4,Phi4)]]
chk(per74 == [[Fr(-77,16)], [Fr(-37,16)], [Fr(-29,16)], []], f"-7/4 exact periods 1..4 params: {per74}")
per29 = [exact_period_params(-R(29,16), n, Ph) for n, Ph in [(1,Phi1),(2,Phi2),(3,Phi3),(4,Phi4)]]
chk(per29 == [[Fr(-1305,256)], [Fr(-633,256)], [], []], f"-29/16 exact periods 1..4 params: {per29}")
chk(exact_period_params(R(1,4), 1, Phi1) == [Fr(3,16)] and exact_period_params(R(1,4), 2, Phi2) == [Fr(-21,16)] and exact_period_params(R(1,4), 3, Phi3) == [], "1/4: periods 1,2 at 3/16,-21/16; no period 3")
chk(exact_period_params(-R(3,4), 1, Phi1) == [Fr(-21,16)] and exact_period_params(-R(3,4), 2, Phi2) == [Fr(-13,16)] and exact_period_params(-R(3,4), 3, Phi3) == [], "-3/4: periods 1,2 at -21/16,-13/16; no period 3")
# families
chk(all(Fr(1,4)*(1 - Fr(rr)**2) == cc for rr, cc in [(1,0),(3,-2),(Fr(5,2),Fr(-21,16)),(Fr(9,2),Fr(-77,16))]), "fixed-point family c=(1-r^2)/4")
chk(all(-(3 + Fr(ss)**2)/4 == cc for ss, cc in [(1,-1),(2,Fr(-7,4)),(0,Fr(-3,4)),(3,-3),(Fr(5,2),Fr(-37,16)),(Fr(9,2),Fr(-93,16))]), "2-cycle family c=-(3+s^2)/4")

print("=== F. Chebyshev side ===")
xq = Fr(5, 2); fer = []
for k in range(1, 6): fer.append(xq.numerator); xq = xq*xq - 2
chk(fer == [2**(2**k) + 1 for k in range(1, 6)], "x^2-2 orbit of 5/2 has Fermat numerators")
chk(all(gcd(fer[i], fer[j]) == 1 for i in range(5) for j in range(i)), "Fermat numerators pairwise coprime")
chk([sum(sp.mobius(n//d)*2**d for d in sp.divisors(n)) for n in range(1, 7)] == [2, 2, 6, 12, 30, 54], "exact-period counts of a degree-2 map")
chk(sp.n_order(2, 7) == 3 and sp.n_order(2, 9) == 6 and sp.n_order(2, 21) == 6 and sp.n_order(2, 63) == 6, "ord_7(2)=3; ord_{9,21,63}(2)=6")
print("=== G. Boundary hostiles for ZT4 wording ===")
chk(orbit_nums(Fr(-1), 2)[0][2] == 0, "c=-1: N_2 = 0 (a+b=0), not covered by '|a+b|=1'")
chk(orbit_nums(Fr(0), 1)[0][1] == 0, "c=0: N_1 = 0 (a=0), not covered by '|a|=1'")
chk(orbit_nums(Fr(-2), 3)[0][3] == 2 and has_prim_gcdfree(orbit_nums(Fr(-2), 3)[0], 3) is False, "c=-2: N_3 = 2 = N_1... F(-2,1)=-1, no primitive prime at n=3 (preperiodic)")

print()
print(f"DONE in {time.time()-T0:.1f}s.  FAILED CHECKS: {BAD if BAD else 'none'}")
if BAD:
    sys.exit(1)
