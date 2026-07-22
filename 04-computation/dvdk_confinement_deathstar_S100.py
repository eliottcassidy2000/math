import numpy as np
from math import comb, gcd
from itertools import product
from fractions import Fraction

print("="*70)
print("CONFINING THE GMC2 DvdK DEPENDENCY -- two independent elementary reductions")
print("="*70)

# ---------------------------------------------------------------
# REDUCTION 1: EDGE FACE (2 monomials) -> elementary binomial CT, NO DvdK.
# ---------------------------------------------------------------
print("\n[1] EDGE FACE (2 straddling monomials): CT is an elementary binomial.")
def edge_ct_witness(a, b):  # charges -a<0<b
    g = gcd(a, b)
    m0 = (a + b)//g          # smallest m with an integer diagonal
    k = b//g                 # copies of the negative monomial
    # channel: k copies of u^-a, (m0-k) copies of u^b ; k*(-a)+(m0-k)*b = 0 ?
    assert k*(-a) + (m0-k)*b == 0
    return m0, k, comb(m0, k)   # CT(f^m0) = C(m0,k)*c_-^k*c_+^(m0-k) != 0
for (a,b) in [(1,1),(1,2),(2,3),(3,5),(2,7),(5,4)]:
    m0,k,C = edge_ct_witness(a,b)
    print(f"   charges (-{a},{b}): first nonzero CT at m0={m0}, k={k}, C(m0,k)={C} != 0  (elementary)")

# ---------------------------------------------------------------
# REDUCTION 2: POSITIVE-real coefficients -> CT>0 always (positive walk return), NO DvdK.
#   (E[P^m] = sum of POSITIVE weights * c^r ; if all c_i>0 then every term >0.)
# ---------------------------------------------------------------
print("\n[2] POSITIVE coefficients: CT(f^m) > 0 (a nonneg walk-return), NO cancellation.")
def CT_pow(charges, coeffs, m):  # constant term of (sum coeffs[i] u^charges[i])^m
    from collections import defaultdict
    poly = defaultdict(float); poly[0]=1.0
    base = dict(zip(charges,coeffs))
    for _ in range(m):
        new=defaultdict(float)
        for e,c in poly.items():
            for q,cc in base.items():
                new[e+q]+=c*cc
        poly=new
    return poly[0]
# resonant AP face charges (-1,0,1) with POSITIVE coeffs -> central-trinomial-like, all >0
print("   resonant face charges (-1,0,1), coeffs (1,1,1): CT(f^m) =",
      [int(round(CT_pow([-1,0,1],[1,1,1],m))) for m in range(1,8)], "(central trinomial A002426, all >0)")
print("   coeffs (2,1,3) (>0): CT(f^m) =",
      [int(round(CT_pow([-1,0,1],[2,1,3],m))) for m in range(1,7)], "(all >0)")

# ---------------------------------------------------------------
# THE HARD CORNER: RESONANT face (>=3 monomials) AND SIGNED coeffs -> low powers CAN vanish.
#   This is the ONLY place DvdK is genuinely needed.
# ---------------------------------------------------------------
print("\n[3] THE HARD CORNER = resonant (>=3 charges) AND signed coeffs: low CT can vanish.")
# charges (-2,-1,1,2) with signed coeffs; scan first nonzero CT
for coeffs in [[1,-1,-1,1],[1,1,-1,-1],[1,-2,2,-1]]:
    seq=[round(CT_pow([-2,-1,1,2],coeffs,m),6) for m in range(1,9)]
    first=next((m+1 for m,v in enumerate(seq) if abs(v)>1e-9), None)
    print(f"   charges (-2,-1,1,2) coeffs {coeffs}: CT(f^1..8)={seq}  first nonzero at m={first}")
print("   -> DvdK guarantees the search terminates; here it does (ESV: by index sum|q|). ")

# ---------------------------------------------------------------
# GENERICITY: the lowest balanced face is an EDGE unless the support is special.
#   LP {min sum a_i x_i : x>=0, sum x=1, sum q_i x_i=0} has 2 equality constraints
#   => a basic optimum has <=2 nonzero x_i (an EDGE). >=3 needs a codim>=1 coincidence.
# ---------------------------------------------------------------
print("\n[4] GENERICITY: the LP for the lowest balanced face has 2 equality constraints")
print("    (mass=1, charge=0) => basic optimum supported on <=2 monomials = an EDGE.")
print("    A >=3-monomial face requires >=3 tilted heights a_i - lambda*q_i concurrent at the")
print("    straddling lambda* -- a codimension->=1 condition on the Newton support (a RESONANCE).")
