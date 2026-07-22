#!/usr/bin/env python3
# SCOPE CORRECTION: 102/116 is a unique-channel search only through mass 40 in
# one bounded support universe. The other 14 pass degree-6 span-membership
# tests; cofactors are not retained, only {-2,-1,1,2} is formalized, and no
# uniform effective bound follows.
"""dvdk_monomial_certificate_residual_boxeph_S231.py -- boxeph-2026-07-22-S231

Eliminate DvdK for the residual (coincident-channel / symmetric) straddling supports that the
death-star-S101 unique-MINIMAL-channel criterion left "hard" (the ~16%/12%).

Two levers, both coefficient-independent and ELEMENTARY (no Galois / THM-2067):

 (1) UNIQUE CHANNEL AT ANY MASS.  My Lean dvdk1_of_uniqueChannel (S230) needs a unique balanced
     channel at SOME mass, not just the minimal one.  Re-scan the S101-hard supports for a unique
     channel at any mass <= cap; those are already DvdK-free.

 (2) MONOMIAL CERTIFICATE.  CT(f^m) = constantTermRelation(q,m) is a homogeneous degree-m polynomial
     in the coefficient variables, with coefficient multinomial(m;r) on the monomial x^r for each
     balanced composition r (sum r=m, sum r_i q_i=0).  A support is DvdK-free elementarily IFF the
     ideal  I = <CT(f^m) : m>=1>  contains a MONOMIAL mu = prod x_i^{e_i}: then on the coefficient
     torus (all c_i != 0) mu(c) != 0, so from  mu = sum_m g_m CT(f^m)  some CT(f^m)(c) != 0.  Because
     each CT(f^m) is homogeneous of degree m, "I contains a degree-D monomial" is the linear-algebra
     question: is some degree-D monomial in the Q-span of { t * CT(f^m) : m<=D, deg t = D-m }?  We
     answer it by exact rational Gaussian elimination (no numpy/sympy).

 (This monomial-in-ideal is exactly V(I) cap torus = empty = DvdK1 for the support -- TRUE by
 THM-2067 -- but here produced as an explicit finite, formalizable certificate, per support.)

The general Lean theorem GMC2DvdKMonomialCertificate.dvdk1_of_monomialCertificate discharges DvdK
from any such identity; this script finds the identities.
"""
from fractions import Fraction as Fr
from math import comb, factorial, gcd
import itertools

def multinomial(r):
    num = factorial(sum(r)); d = 1
    for k in r: d *= factorial(k)
    return num // d

def balanced_comps(ch, m):
    """All r=(r_0..r_{n-1}) >=0, sum=m, sum r_i*ch_i = 0."""
    n = len(ch); out = []
    def rec(i, rem, acc, s):
        if i == n-1:
            if s + rem*ch[i] == 0: out.append(tuple(acc+[rem]))
            return
        for k in range(rem+1):
            rec(i+1, rem-k, acc+[k], s+k*ch[i])
    rec(0, m, [], 0)
    return out

def CT(ch, m):
    """constantTermRelation as dict {exponent-tuple r : multinomial(m;r)} (homogeneous degree m)."""
    return {r: Fr(multinomial(r)) for r in balanced_comps(ch, m)}

def min_channels(ch, cap=40):
    for m in range(1, cap+1):
        cs = balanced_comps(ch, m)
        if cs: return m, cs
    return None, []

def unique_channel_mass(ch, cap=40):
    """smallest mass with a UNIQUE balanced composition, or None."""
    for m in range(1, cap+1):
        cs = balanced_comps(ch, m)
        if len(cs) == 1: return m, cs[0]
    return None, None

# ---- exact rational Gaussian elimination over monomial vectors (dicts) ----
def monoms_deg(n, D):
    if n == 1: yield (D,); return
    for k in range(D+1):
        for rest in monoms_deg(n-1, D-k):
            yield (k,)+rest

def reduce_vec(v, basis):
    """reduce dict v by echelon basis (list of (pivot_key, dict)); return residual dict."""
    v = dict(v)
    for pk, bv in basis:
        if pk in v and v[pk] != 0:
            f = v[pk] / bv[pk]
            for k, val in bv.items():
                v[k] = v.get(k, Fr(0)) - f*val
                if v[k] == 0: del v[k]
    return v

def echelon(gens, order):
    """row-reduce generator dicts; order = list of monomial keys (pivot preference)."""
    basis = []
    for g in gens:
        r = reduce_vec(g, basis)
        r = {k: val for k, val in r.items() if val != 0}
        if not r: continue
        pk = min((k for k in r), key=lambda k: order.index(k))
        basis.append((pk, r))
        # keep basis reduced enough for membership tests
    return basis

def ideal_contains_monomial(ch, Dcap=9):
    """Return (D, monomial, masses_used) for the least degree D at which the CT-ideal contains a
    monomial, else None. masses_used = the m whose generators are needed (from the support)."""
    n = len(ch)
    CTs = {}
    for m in range(1, Dcap+1):
        c = CT(ch, m)
        if c: CTs[m] = c
    for D in range(2, Dcap+1):
        keys = list(monoms_deg(n, D))
        order = keys  # any fixed order
        gens = []
        gen_tags = []  # which mass each generator came from
        for m in [mm for mm in CTs if mm <= D]:
            for t in monoms_deg(n, D-m):
                # product t * CT(ch,m): shift each exponent tuple by t
                prod = {}
                for r, coef in CTs[m].items():
                    kk = tuple(t[i]+r[i] for i in range(n))
                    prod[kk] = prod.get(kk, Fr(0)) + coef
                gens.append(prod); gen_tags.append(m)
        if not gens: continue
        basis = echelon(gens, order)
        # is any unit monomial e_mu in the span? test each degree-D monomial
        for mu in keys:
            ev = {mu: Fr(1)}
            if not reduce_vec(ev, basis):  # reduces to empty => in span
                masses = sorted(set(gen_tags))
                return D, mu, masses
    return None

# ==========================================================================
print("="*76)
print("[1] Re-scan straddling supports (size 3-4, range +-4): unique channel at ANY mass")
print("="*76)
pool = [-4,-3,-2,-1,1,2,3,4]
tot=0; free_min=0; free_any=0; residual=[]
for r in (3,4):
    for c in itertools.combinations(pool, r):
        if not (min(c) < 0 < max(c)): continue
        tot += 1
        _, mincs = min_channels(list(c))
        um, _ = unique_channel_mass(list(c))
        m0, m0cs = min_channels(list(c))
        if len(m0cs) == 1: free_min += 1
        if um is not None: free_any += 1
        else: residual.append(c)
print(f"  total straddling supports: {tot}")
print(f"  DvdK-free by UNIQUE MINIMAL channel (S101):          {free_min}/{tot} = {100*free_min/tot:.1f}%")
print(f"  DvdK-free by UNIQUE channel at ANY mass (my S230):   {free_any}/{tot} = {100*free_any/tot:.1f}%")
print(f"  RESIDUAL (no unique channel at any mass<=40):        {len(residual)}/{tot} = {100*len(residual)/tot:.1f}%")
print(f"  gain from 'any mass': {free_any-free_min} supports reclassified free")
print(f"  residual list: {residual}")

print()
print("="*76)
print("[2] MONOMIAL CERTIFICATE for each residual support (eliminates DvdK, no Galois)")
print("="*76)
cert_ok = 0
for c in residual:
    res = ideal_contains_monomial(list(c), Dcap=9)
    if res:
        D, mu, masses = res
        cert_ok += 1
        print(f"  {c}: certificate at degree {D}, monomial exps {mu}, masses used {masses}  -> DvdK-FREE")
    else:
        print(f"  {c}: NO monomial certificate up to degree 9 (needs higher degree)")
print(f"\n  {cert_ok}/{len(residual)} residual supports certified DvdK-free by a finite monomial certificate")

print()
print("="*76)
print("[3] The paradigm residual {-2,-1,1,2}: explicit hand certificate (for the Lean instance)")
print("="*76)
ch = [-2,-1,1,2]
CT2 = CT(ch,2); CT3 = CT(ch,3); CT4 = CT(ch,4)
# variables a=x_{-2}, b=x_{-1}, c=x_{1}, d=x_{2}; exps (a,b,c,d)
def show(P):
    return " + ".join(f"{v}*a^{e[0]}b^{e[1]}c^{e[2]}d^{e[3]}" for e,v in sorted(P.items()))
print("  CT(f^2) =", show(CT2))
print("  CT(f^3) =", show(CT3))
print("  CT(f^4) =", show(CT4))
# claim: 12*b^2c^2 = (3ad+9bc)*CT2 - CT4    [monomial b^2c^2 = exps (0,2,2,0)]
def scal(P,s): return {k:v*s for k,v in P.items()}
def add(*Ps):
    out={}
    for P in Ps:
        for k,v in P.items():
            out[k]=out.get(k,Fr(0))+v
    return {k:v for k,v in out.items() if v!=0}
def mulmono(P, e):  # multiply poly P by monomial with exps e
    return {tuple(k[i]+e[i] for i in range(4)):v for k,v in P.items()}
# (3ad+9bc)*CT2 : g2 = 3*a d + 9*b c  = 3*(1,0,0,1) + 9*(0,1,1,0)
g2_CT2 = add(scal(mulmono(CT2,(1,0,0,1)),Fr(3)), scal(mulmono(CT2,(0,1,1,0)),Fr(9)))
rhs = add(g2_CT2, scal(CT4,Fr(-1)))
print("  (3ad+9bc)*CT(f^2) - CT(f^4) =", show(rhs))
target = {(0,2,2,0): Fr(12)}
print("  target 12*b^2c^2 =", show(target))
print("  IDENTITY HOLDS:", rhs == target)
print("  => on the torus (b,c != 0): CT(f^2)=0 forces CT(f^4) = -12 b^2 c^2 != 0.")
print("  => DvdK-free for {-2,-1,1,2} via the TWO-MASS certificate {2,4} -- elementary, no Galois.")
