#!/usr/bin/env python3
"""
klein-2026-07-20-S363 -- THE GAMMA BRIDGE'S REAL MECHANISM IS POSITIVITY, NOT DOMINATION;
and the top-r-coefficient identity that makes the bridge's algebraic core rigorous.

STATE OF THE THREAD (read before working, to avoid duplication):
  * TNC (constant-coeff toral nullcone) = Duistermaat-van der Kallen 1998, Thm 2 + Rmk 3.
    A CITATION.  My THM-1530/1550 etc are rediscoveries (mac-mini THM-1630).
  * NC2 => GMC(2): FORMALIZED in Lean (death-star mathieuZhao_of_charge_pos, kernel-pure).
  * TNC => NC2 (the r-dependent lift = the "Gamma bridge", klein-S351): THE LIVE GAP.
    death-star MISTAKE-202 REFUTED the ">50% mass" domination (top fraction -> 0.068 at m=24).
    That is the SAME mechanism my S351 bridge invoked ("top term dominates the r-average").
    So I must check my own claim, and I do below.

FOUR PARTS, each with a control:
  1. THE TOP-r-COEFFICIENT IDENTITY, proved and exact-tested: [r^{mD}]CT_u(P^m) = CT_u(Lam^m),
     Lam = the top-r-degree symbol.  Pure algebra, no analysis -- the rigorous core.
  2. DOMINATION FAILS -- independent confirmation of MISTAKE-202, and hence RETRACTION of my
     own S351 mechanism.  Control: a case with E[P^m] != 0 known.
  3. POSITIVITY is the real mechanism on the radial {-1,+1} span: alpha = r|c|^2 >= 0 makes
     E_r[alpha^j] > 0 termwise -- no domination needed.  This matches opus/boxeph Hankel-PD.
  4. WHY THE BOTH-SIGNS {-1,0,1} CASE ESCAPES POSITIVITY: the beta (charge-0) terms make L_m
     sign-indefinite, exhibited exactly.  That is the precise remaining gap.
"""
import sympy as sp
from fractions import Fraction as Fr
import itertools

r, u, t = sp.symbols('r u t')

# ---------------------------------------------------------------- PART 1
print("=" * 84)
print("PART 1 -- THE TOP-r-COEFFICIENT IDENTITY (the rigorous algebraic core)")
print("=" * 84)
print(" CLAIM: for P(r,u) = sum_q g_q(r) u^q, D = max_q deg_r g_q, and leading symbol")
print("        Lam(u) = sum_{deg g_q = D} lc(g_q) u^q,  we have  [r^{mD}] CT_u(P^m) = CT_u(Lam^m).")
print(" PROOF: CT_u(P^m) = sum_{q_1+...+q_m=0} prod_i g_{q_i}(r); a product has r-degree mD iff")
print("        every q_i is a top-degree charge, and then its [r^{mD}] coeff is prod lc(g_{q_i}).\n")

def CT_of_power(gdict, M, m):
    """CT_u(P^m) as a sympy poly in r; gdict: charge q -> g_q(r) (sympy expr)"""
    cur = {0: sp.Integer(1)}
    for _ in range(m):
        nxt = {}
        for e, v in cur.items():
            for q, g in gdict.items():
                nxt[e + q] = sp.expand(nxt.get(e + q, 0) + v * g)
        cur = nxt
    return sp.expand(cur.get(0, sp.Integer(0)))

def top_symbol_CT(gdict, m):
    """CT_u(Lam^m) where Lam = top-r-degree part"""
    D = max(sp.Poly(g, r).degree() for g in gdict.values())
    lam = {q: sp.Poly(g, r).LC() for q, g in gdict.items() if sp.Poly(g, r).degree() == D}
    cur = {0: sp.Integer(1)}
    for _ in range(m):
        nxt = {}
        for e, v in cur.items():
            for q, c in lam.items():
                nxt[e + q] = nxt.get(e + q, 0) + v * c
        cur = nxt
    return sp.expand(cur.get(0, sp.Integer(0))), D

# control: three random-ish P's, exact
tests = [
    {-1: 1 + 2*r, 0: 3*r, 1: r - 1},                 # D=1, all charges top
    {-2: r**2, -1: 5*r, 0: 1 + r**2, 1: 2*r, 2: r**2 - r},   # D=2, top charges -2,0,2
    {-1: r**3 + 1, 0: 7*r, 1: r**3},                 # D=3, top charges -1,1
]
allok = True
for gd in tests:
    for m in (2, 3, 4):
        Lm = CT_of_power(gd, None, m)
        ctlam, D = top_symbol_CT(gd, m)
        # coefficient of r^{mD} in CT(P^m) -- ZERO if the r-degree is below mD, which is the
        # correct reading (both sides vanish together).  Do NOT require deg == mD.
        coeff_mD = sp.Poly(Lm, r).coeff_monomial(r**(m*D))
        match = sp.simplify(coeff_mD - ctlam) == 0
        allok = allok and match
    print(f"   P charges {sorted(gd)}, D={D}: [r^{{mD}}]CT(P^m) = CT(Lam^m) for m=2,3,4 : "
          f"{'all match' if allok else 'MISMATCH'}")
print(f"\n   IDENTITY HOLDS on every control: {allok}")

# ---------------------------------------------------------------- PART 2
print("\n" + "=" * 84)
print("PART 2 -- DOMINATION FAILS: independent check of MISTAKE-202, retracting my S351 mechanism")
print("=" * 84)
def Er_poly(expr):
    p = sp.Poly(sp.expand(expr), r)
    return sum(c * sp.factorial(k) for (k,), c in zip(p.monoms(), p.coeffs()))
# a genuine two-sided non-one-sided P with E[P^m] != 0: use the known {-1,0,1} shape
a, b, c = sp.Integer(1), sp.Integer(1), sp.Integer(1)       # constants -> alpha = r, beta = 1
alpha, beta = r*a*c, b
def L_leg(al, be, m):
    tot = 0
    for k in range(m//2+1):
        tot += sp.factorial(m)//(sp.factorial(k)**2*sp.factorial(m-2*k)) * al**k * be**(m-2*k)
    return sp.expand(tot)
print(" P = Zb + 1 + Z  (alpha = r, beta = 1); E[P^m] known nonzero.  Track |top term|/sum|terms|")
print(" of E_r[L_m] = sum_k ([r^k]L_m) k!  as m grows:\n")
print(f"   {'m':>4} {'E[P^m]':>14} {'|top term|/sum|terms|':>24}")
for m in (2, 6, 12, 18, 24):
    Lm = L_leg(alpha, beta, m)
    p = sp.Poly(Lm, r)
    terms = [abs(int(cf)) * sp.factorial(k) for (k,), cf in zip(p.monoms(), p.coeffs())]
    D = p.degree()
    topterm = abs(int(p.coeff_monomial(r**D))) * sp.factorial(D)
    frac = sp.Rational(topterm, sum(terms))
    Em = Er_poly(Lm)
    print(f"   {m:>4} {int(Em):>14} {float(frac):>24.4f}")
print("""
   CONFIRMED: the top term's fraction COLLAPSES (matches death-star's 0.068 at m=24).  So the
   ">50% mass / top term dominates" mechanism is FALSE -- and that was exactly the mechanism my
   klein-S351 Gamma bridge invoked.  I RETRACT THE MECHANISM.  The bridge's CONCLUSION
   (E[P^m] != 0) is still true here, but NOT because the top term dominates.  Something else
   keeps it nonzero.  Part 3 identifies what.
""")

# ---------------------------------------------------------------- PART 3
print("\n" + "=" * 84)
print("PART 3 -- POSITIVITY is the real mechanism on the radial {-1,+1} span")
print("=" * 84)
print(" On the {-1,+1} span (beta = 0, charge-0 empty), the {-1,0,1} GF collapses:")
print("   E[P^{2j}] = C(2j,j) E_r[alpha^j],  alpha = r * a(r) * c(r).")
print(" For P in C[Z,Zb] the charge +1 and -1 coefficients are conjugate, so a*c = |c|^2 >= 0,")
print(" whence alpha = r|c|^2 >= 0 pointwise and E_r[alpha^j] > 0 unless c == 0.  NO DOMINATION")
print(" NEEDED -- every moment is a positive number, so it cannot vanish.  Exhibited:\n")
for cpoly, name in [(sp.Integer(1), "c=1"), (1 + r, "c=1+r"), (r**2 - 3, "c=r^2-3 (has a sign change!)")]:
    al = r * cpoly * cpoly            # |c|^2 with c real here; a = conj(c) = c for real
    print(f"   {name}: alpha = r*c^2, E_r[alpha^j] for j=1..4 = "
          f"{[int(Er_poly(sp.expand(al**j))) for j in range(1,5)]}  (all > 0: "
          f"{all(Er_poly(sp.expand(al**j)) > 0 for j in range(1,5))})")
print("""
   Note c = r^2 - 3 CHANGES SIGN, yet alpha = r c^2 >= 0 regardless -- the SQUARE is what
   matters, not c.  So positivity holds for ALL c, no degree bound, matching opus THM-1535/1540
   (Hankel positive-definiteness) and boxeph's radial route.  This is the correct mechanism,
   and it is immune to the domination collapse of Part 2 because it never compares terms.
""")

# ---------------------------------------------------------------- PART 4
print("=" * 84)
print("PART 4 -- WHY THE BOTH-SIGNS {-1,0,1} CASE ESCAPES POSITIVITY (the precise gap)")
print("=" * 84)
print(" With beta != 0, L_m(alpha,beta) = sum_k C(...) alpha^k beta^{m-2k} is NO LONGER a")
print(" positive combination: beta^{m-2k} alternates in sign when beta < 0, and even for")
print(" beta > 0 the r-integrand L_m(r) can take both signs.  Exhibited on alpha=r, beta=1-r:\n")
albe = [(r, sp.Integer(1) - r, "alpha=r, beta=1-r")]
for al, be, name in albe:
    print(f"   {name}:")
    for m in (2, 3, 4, 5):
        Lm = L_leg(al, be, m)
        p = sp.Poly(Lm, r)
        signs = set(sp.sign(cf) for cf in p.coeffs() if cf != 0)
        Em = Er_poly(Lm)
        print(f"      m={m}: L_m(r) coefficient signs = {sorted(str(s) for s in signs)}"
              f"   E_r[L_m] = {int(Em)}")
print("""
   So on the both-signs span the integrand L_m(r) is SIGN-INDEFINITE -- positivity is
   unavailable, and (Part 2) domination is false.  THE PRECISE REMAINING GAP for GMC(2):
   show E_r[L_m] != 0 for some m when (alpha,beta) is not the trivial (r|c|^2, 0), using
   NEITHER termwise positivity NOR top-term domination.  The known routes are opus's Hankel
   argument on the charge-0 sub-block and boxeph's pinch/Watson analysis; this file's
   contribution is to have (i) made the algebraic core rigorous (Part 1), (ii) killed the
   domination mechanism cleanly (Part 2, retracting my own S351), and (iii) isolated the gap
   as sign-indefiniteness of L_m rather than anything about growth rates (Part 4).
""")
