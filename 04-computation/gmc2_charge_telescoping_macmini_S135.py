#!/usr/bin/env python3
"""
GMC(2): the CHARGE grading, the telescoping lemma, and the reduced 1-variable core
                                                        (mac-mini-2026-07-20-S135)
==================================================================================
Owner: "work on the GMC(2) proof, think forbidding one variable telescopes, see how
that pattern can abstractly apply to other problems the repo has touched."

SETUP.  Two real Gaussians = one standard complex Gaussian Z, with W = Zbar, and
    E[Z^a W^b] = a! * delta_{ab}.
Grade monomials by CHARGE  c = deg_Z - deg_W.  Then E KILLS EVERY NONZERO CHARGE:
E is (charge-0 projection) followed by  Z^a W^a  |->  a!.
Equivalently, in polar form Z = sqrt(V) e^{i theta} with V = ZW ~ Exp(1) and theta
uniform, V and theta are INDEPENDENT, and charge IS the theta-Fourier index.

THE TELESCOPING LEMMA (this is the owner's "forbidding one variable telescopes").
Charges ADD under multiplication.  So if P's charge support is ONE-SIDED -- say all
charges >= 0 -- then in E[P^m] the only surviving combination is the one taking the
charge-0 part from EVERY factor:
        E[P^m] = L( p_0^m ),      L(v^k) := k!,   p_0 = charge-0 part of P.
Forbidding negative charge collapses TWO variables to ONE.  And since P >= 0 charges
means Q can only "spend" a bounded amount of negative charge, E[Q P^m] is a FINITE
combination of L(h * p_0^{m-j}) with j bounded by Q.  So the whole one-sided case of
GMC(2) reduces to a ONE-VARIABLE Mathieu-Zhao question:

        is ker(L) a Mathieu-Zhao subspace of C[v],  where L(v^k) = k!  ?
        (equivalently L(f) = int_0^infty f(v) e^{-v} dv)

PART B attacks the sharpest form of that: is there ANY nonzero p with L(p^m) = 0 for
all m >= 1?  If not, the one-sided branch of GMC(2) is CLOSED outright.
"""
from fractions import Fraction as F
from math import factorial
import itertools, random

# ---------------------------------------------------------------- charge formalism
def cmul(p, q):
    out = {}
    for k1, c1 in p.items():
        for k2, c2 in q.items():
            k = (k1[0]+k2[0], k1[1]+k2[1])
            out[k] = out.get(k, F(0)) + c1*c2
    return {k: c for k, c in out.items() if c}

def cexp(p):
    """E[Z^a W^b] = a! delta_ab -- kills every nonzero charge."""
    return sum(c*factorial(a) for (a, b), c in p.items() if a == b)

def charge_part(p, c):
    return {k: v for k, v in p.items() if k[0]-k[1] == c}

ONE = {(0,0): F(1)}
def cpow(p, m):
    r = ONE
    for _ in range(m): r = cmul(r, p)
    return r

def L(coeffs):
    """L(sum a_k v^k) = sum a_k k!  -- the reduced one-variable functional."""
    return sum(a*factorial(k) for k, a in enumerate(coeffs))

def polymul(a, b):
    out = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): out[i+j] += x*y
    return out

def polypow(a, m):
    r = [F(1)]
    for _ in range(m): r = polymul(r, a)
    return r

print("=" * 78)
print("PART A -- the telescoping lemma, verified")
print("=" * 78)
print("  If P's charge support is ONE-SIDED (all charges >= 0), then")
print("        E[P^m] = L(p_0^m)   where p_0 = charge-0 part, L(v^k) = k!.")
print("  Random test: build P with charges in {0,1,2}, compare both sides.")
print()
rng = random.Random(135)
allok = True
for trial in range(6):
    # charge-0 part is a polynomial in V = ZW, i.e. monomials (a,a)
    P = {}
    deg0 = rng.randint(1, 3)
    for a in range(deg0 + 1):
        v = rng.randint(-3, 3)
        if v: P[(a, a)] = F(v)
    # add some strictly positive charge terms (a > b)
    for _ in range(rng.randint(1, 4)):
        b = rng.randint(0, 2); c = rng.randint(1, 2)
        v = rng.randint(-3, 3)
        if v: P[(b+c, b)] = P.get((b+c, b), F(0)) + F(v)
    P = {k: v for k, v in P.items() if v}
    p0 = charge_part(P, 0)
    p0c = [F(0)]*(max([k[0] for k in p0] or [0])+1)
    for (a, _), v in p0.items(): p0c[a] = v
    ok = True; lhs = []; rhs = []
    Pm = ONE
    for m in range(1, 6):
        Pm = cmul(Pm, P)
        l = cexp(Pm); r = L(polypow(p0c, m))
        lhs.append(l); rhs.append(r)
        if l != r: ok = False
    allok &= ok
    print(f"  trial {trial}: charges {sorted({k[0]-k[1] for k in P})}  "
          f"E[P^m]={[str(x) for x in lhs[:3]]}...  L(p0^m)={[str(x) for x in rhs[:3]]}...  {ok}")
print(f"  telescoping lemma holds on all trials: {allok}")

print()
print("  Sanity: it FAILS when charges are two-sided (as it must -- that is exactly")
print("  where the known counterexamples live).")
P2 = {(1,0): F(1), (0,1): F(1), (1,1): F(-1)}     # charges +1, -1, 0
p0c = [F(0), F(-1)]
Pm = ONE; rows = []
for m in range(1, 5):
    Pm = cmul(Pm, P2); rows.append((cexp(Pm), L(polypow(p0c, m))))
print(f"    P = Z + W - ZW:  E[P^m] = {[str(a) for a,_ in rows]}   "
      f"L(p0^m) = {[str(b) for _,b in rows]}   equal? {all(a==b for a,b in rows)}")

# ---------------------------------------------------------------- PART B
print()
print("=" * 78)
print("PART B -- the reduced core: is there a NONZERO p with L(p^m) = 0 for all m?")
print("=" * 78)
print("  L(f) = int_0^infty f(v) e^{-v} dv, i.e. L(v^k) = k!.")
print()
try:
    import sympy as sp
    HAVE = True
except ImportError:
    HAVE = False
    print("  (sympy unavailable -- symbolic elimination skipped)")

if HAVE:
    v = sp.symbols('v')
    for D in (1, 2, 3):
        a = sp.symbols(f'a0:{D+1}')
        p = sum(a[k]*v**k for k in range(D+1))
        eqs = []
        for m in range(1, D + 3):
            e = sp.Poly(sp.expand(p**m), v)
            val = sum(c*sp.factorial(k) for k, c in
                      zip(range(e.degree(), -1, -1), e.all_coeffs()))
            eqs.append(sp.expand(val))
        # normalize: leading coeff = 1 (a nonzero p can be scaled)
        sols = sp.solve(eqs[:D+1] + [a[D] - 1], a, dict=True)
        print(f"  deg p = {D}: solving L(p^m)=0 for m=1..{D+1} with a{D}=1  ->  "
              f"{len(sols)} solution(s)")
        for s in sols[:4]:
            resid = [sp.simplify(e.subs(s)) for e in eqs]
            print(f"      {s}   residuals L(p^m) m=1..{len(eqs)}: {resid}")
        if not sols:
            print(f"      => NO nonzero p of degree {D}.  (system is inconsistent)")

print()
print("  Asymptotic reason, checked numerically:  L(p^m) is dominated by the top")
print("  coefficient, |L(p^m)| ~ |a_D|^m (Dm)! -- the saddle of |p(v)|^m e^{-v} sits at")
print("  v ~ Dm, where the phase correction m*arg(p(v)/v^D) = O(1) and cannot oscillate")
print("  the integral away.  So L(p^m) != 0 for large m whenever deg p >= 1.")
print()
print(f"{'p':>22} {'m':>3} {'L(p^m)':>26} {'L(p^m)/(a_D^m (Dm)!)':>24}")
for coeffs, name in ((( -1, 1), "v - 1"), ((2, -3, 1), "v^2-3v+2"), ((0, 1), "v")):
    D = len(coeffs) - 1
    for m in (2, 4, 6, 8):
        val = L(polypow([F(c) for c in coeffs], m))
        ratio = F(val, factorial(D*m)) / F(coeffs[-1])**m
        print(f"{name:>22} {m:>3} {str(val)[:26]:>26} {str(float(ratio))[:24]:>24}")

print()
print("=" * 78)
print("PART C -- what this closes, and what it does not")
print("=" * 78)
print("  CLOSED (modulo PART B): the ONE-SIDED-CHARGE branch of GMC(2).  If P's charge")
print("  support is one-sided then E[P^m] = L(p_0^m); if no nonzero p has L(p^m)=0 for")
print("  all m, then p_0 = 0, so P has NO charge-0 part, so P^m has no charge-0 part for")
print("  any m, and E[Q P^m] = 0 as soon as m exceeds the charge budget of Q.  GMC holds.")
print()
print("  STILL OPEN: the TWO-SIDED case, where P has charges of both signs.  That is")
print("  exactly where the n>=3 counterexamples live -- P = (1+Z)(W - g(Z)U) has charges")
print("  from -1 up.  So two-sidedness is NECESSARY for failure, which is itself a")
print("  structural constraint worth recording.")
