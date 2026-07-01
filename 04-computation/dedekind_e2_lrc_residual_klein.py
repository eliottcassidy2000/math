#!/usr/bin/env python3
"""
dedekind_e2_lrc_residual_klein.py  --  klein-2026-06-30-S56

The B2 / E2 / DEDEKIND package on the genuinely hard residual (the coverage crux = the global
iota-odd degree / the genus-cusp-form obstruction, OPEN-Q-108 / HYP-3587).

The bridge is the SAWTOOTH  ((x)) = x - floor(x) - 1/2  (0 at integers), which is:
  * iota-ODD:  ((-x)) = -((x))       -- the sign atom of klein-S55 (HYP-3767)
  * tied to distance:  ||x|| = 1/2 - |((x))|   -- loneliness (||st|| large) <=> ((st)) near 0
                                                <=> runner st near the HALF-INTEGER (the iota-fixed pt)
Dedekind sum  s(h,k) = sum_{i=1}^{k-1} ((i/k))((hi/k))  correlates two sawtooths; its RECIPROCITY
  s(h,k) + s(k,h) = -1/4 + (1/12)(h/k + k/h + 1/(hk))          [1/12 = B2/2, B2 = 1/6]
is a TWO-MODULUS relation -- the prototype of the multi-metric sheaf's CRT/continued-fraction GLUING.
E2 = 1 - (4/B2) sum sigma_1(n) q^n is the weight-2 quasimodular Eisenstein series; its transformation
ANOMALY is the Dedekind sum, and the residual obstruction (genus = cusp forms, HYP-3587) is what the
Eisenstein/E2 'bulk' cannot see.

This script: (1) the constant unification (zeta(2)=pi^2 B2, the witness-sum heartbeat HYP-3746);
(2) ||x||=1/2-|((x))| and loneliness = runners at the half-integer; (3) Dedekind sums of the LRC
bindings (AP, covering-min, (n+q)-witness) + reciprocity descent as the sheaf gluing.
"""
from fractions import Fraction as F
from math import gcd, pi

def frac(x): return x - (x.numerator // x.denominator) if isinstance(x,F) else x - int(x//1)
def saw(x):
    """Dedekind sawtooth ((x)) for rational x (0 at integers)."""
    if isinstance(x,F):
        if x.denominator==1: return F(0)
        f = x - (x.numerator // x.denominator)   # fractional part in [0,1)
        return f - F(1,2)
    return 0.0

def dist_nearest(x):
    """||x|| = distance to nearest integer, rational."""
    f = x - (x.numerator // x.denominator)
    return min(f, 1-f)

def dedekind(h,k):
    """s(h,k) = sum_{i=1}^{k-1} ((i/k))((hi/k)), exact."""
    tot=F(0)
    for i in range(1,k):
        tot += saw(F(i,k))*saw(F(h*i,k))
    return tot

def dedekind_recip_rhs(h,k):
    return F(-1,4) + F(1,12)*(F(h,k)+F(k,h)+F(1,h*k))

if __name__=="__main__":
    print("="*72)
    print("(1) CONSTANT UNIFICATION via B2 = 1/6")
    print("="*72)
    B2=F(1,6)
    print(f"   B2 (Bernoulli) = {B2}")
    print(f"   zeta(2) = pi^2/6 = pi^2 * B2   (my witness-sum heartbeat 1/zeta(2), HYP-3746)")
    print(f"   E2 = 1 - (4/B2) sum sigma_1 q^n = 1 - 24 q - ...   (4/B2 = {4/B2})")
    print(f"   Dedekind reciprocity constant 1/12 = B2/2 = {B2/2}")
    print(f"   => one Bernoulli B2 runs the witness-sum zeta(2), E2, and Dedekind reciprocity.")

    print()
    print("="*72)
    print("(2) SAWTOOTH BRIDGE: ||x|| = 1/2 - |((x))|; loneliness = runners at the half-integer")
    print("="*72)
    for x in [F(1,7),F(3,7),F(1,2),F(5,13),F(1,14)]:
        lhs=dist_nearest(x); rhs=F(1,2)-abs(saw(x))
        print(f"   x={x}: ||x||={lhs}={float(lhs):.4f}, 1/2-|((x))|={rhs}  match={lhs==rhs}; ((-x))=-((x))? {saw(-x)==-saw(x)}")
    print("   => a lonely observer (all ||st|| large) has every runner st near 1/2 (the iota-fixed point).")

    print()
    print("="*72)
    print("(3) DEDEKIND SUMS of LRC bindings + RECIPROCITY as the sheaf gluing (CF descent)")
    print("="*72)
    n=14
    bindings = [
        ("AP tight  (D=14,a=1)", 1, 14),
        ("covering-min constr (D=Phi6=183,a=14)", 14, 183),
        ("(n+q)-witness q=13 (D=27,a=25)", 25, 27),
        ("(n+q)-witness q=11 (D=25,a=16)", 16, 25),
        ("covering-min beater n=10 (D=37,a=?)", 10, 37),
    ]
    for name,h,k in bindings:
        if gcd(h,k)!=1:
            print(f"   {name}: gcd({h},{k})={gcd(h,k)}!=1, s not defined (non-unit rotation)"); continue
        s_hk = dedekind(h,k); s_kh = dedekind(k%h if h>1 else 0, h) if h>1 else F(0)
        recip = dedekind_recip_rhs(h,k)
        check = (s_hk + dedekind(k,h)) == recip if h>0 else None
        print(f"   {name}:  s({h},{k}) = {s_hk} = {float(s_hk):+.4f};  reciprocity holds: {s_hk+dedekind(k,h)==recip}")
    print()
    print("   RECIPROCITY = the two-modulus GLUING: s(a,D) [modulus D] <-> s(D,a) [modulus a<D], the")
    print("   Euclidean/continued-fraction DESCENT. Iterating it walks D down its CF [n-1, a(n)] --")
    print("   the SAME Farey/three-gap ladder the covering-min lives on (HYP-3732/3762). The descent")
    print("   TERMINATES at small modulus (base case) UNLESS an obstruction (E2 anomaly / genus cusp")
    print("   form, HYP-3587) blocks it -- that obstruction IS the iota-odd degree residual.")
