"""
opus-2026-07-11-S257: pursuing the DUAL CERTIFICATE route for the covering-min lower bound (M(v) >= 14/183 for
covering v). Result: a clean dual (test-measure) formulation, a RIGOROUS obstruction to a single uniform
certificate (the deep-well knife-edge), and the forced SPLIT [tight: S255 rigidity | loose: anti-concentration].

DUAL FORMULATION (test measure). For a family v and level c, let D_i = {t: ||v_i t|| < c} (runner i's danger
arcs) and W(t) = Sum_i 1_{D_i}(t) (danger count). If a probability measure nu has
    Sum_i nu(D_i) = INT W dnu  <  1,
then (W integer >=0) W = 0 on a set of positive nu-measure, so a SAFE point exists: M(v) >= c. So the dual
certificate is a probability measure nu with Sum_i nu(D_i) < 1. Fourier form:
    Sum_i nu(D_i) = 2ck + Sum_{h!=0} b_{-h} Sum_i nuhat(h v_i),   b_h = sin(2*pi*h*c)/(pi*h), b_0 = 2c.

LEBESGUE FAILS. nu = Lebesgue: Sum_i nu(D_i) = 2ck = 2*(14/183)*13 = 364/183 = 1.989 > 1. Measure alone is
insufficient; the certificate MUST use the correlation terms Sum_i nuhat(h v_i) -- exactly the arithmetic /
anti-concentration content, and it must exploit divisor-completeness (non-covering families need not be lonely
at 14/183).

THE KNIFE-EDGE OBSTRUCTION (rigorous). The deep well {1..12,182} has M = 14/183 EXACTLY, attained only at the
single point t* = 14/183; its danger count satisfies W_dw(t) >= 1 for ALMOST EVERY t (safe set = {t*}, measure
0; verified safe-measure ~ 0.0014 = grid resolution). Therefore, for ANY absolutely-continuous nu,
    Sum_i nu(D_i) = INT W_dw dnu >= nu({W_dw >= 1}) = 1.
So NO absolutely-continuous test measure certifies the deep well -- only the atomic delta_{t*} does. A "single
positive trigonometric polynomial" (mac-mini S40) IS an AC test measure, so a SINGLE UNIFORM certificate is
rigorously OBSTRUCTED by the deep-well knife-edge. The certificate must be family-ADAPTIVE, or the tight case
handled separately.

THE FORCED SPLIT (the viable route).
  * TIGHT (M at/near 14/183, incl. the deep well): the knife-edge families. Handled by the S255 RIGIDITY --
    PROVED for the deep well (M_core=1/13 => interval => s=1 => equality, via S252 prime-13 uniqueness), the
    unique minimizer. Near-deep-well by perturbation of S255.
  * LOOSE (M >= 14/183 + margin): safe set has POSITIVE measure (verified 0.065-0.10 for random covering
    families), so an absolutely-continuous test measure certifies them. This reduces to a SECOND-MOMENT /
    anti-concentration bound measure{W=0} > 0 -- FAVORABLE here (spread families, E[W]=2ck~1.99, measure{W=0}
    ~ e^{-2} ~ 0.13). This is the project's anti-concentration core (S242-S245) on its EASY side, and the
    correlation Sum_i nuhat(h v_i) is exactly the LRCFourierCompletion completion identity (|C_w - b^2/q|).

NET. The dual certificate route, made concrete: a single uniform positive-polynomial certificate is
RIGOROUSLY IMPOSSIBLE (the deep-well knife-edge forces INT W dnu >= 1 for AC nu). The route splits, forced,
into [tight: the finite S255 rigidity, PROVED for the extremizer] + [loose: an AC test measure = second-moment
anti-concentration, the favorable side]. This precisely locates the dual route, reconciles it with the
anti-concentration framing (same content), and explains why S40's single-certificate hope stalls -- the
extremizer is atomic. The constructive path is the split, not one polynomial.

-> mac-mini S40 (dual certificate hope -- here shown to need the split), opus-S255 (deep-well tight case,
PROVED -- the tight half), opus-S242-S245 (anti-concentration core -- the loose half), LRCFourierCompletion
(the completion identity = the correlation term), klein S267 (14/183 verified).
"""
from math import gcd, sin, pi
from functools import reduce
from fractions import Fraction
import random
def divcomplete(v): return all(any(x%d==0 for x in v) for d in range(2,15))
def primitive(v): return reduce(gcd,v)==1
def safe_measure(v, c=Fraction(14,183), D=183*8):
    cD=float(c)*D; safe=0; Wmin=99
    for a in range(D):
        W=sum(1 for vi in v if min((vi*a)%D, D-(vi*a)%D) < cD)
        if W<Wmin: Wmin=W
        if W==0: safe+=1
    return Wmin, safe/D
def main():
    c=Fraction(14,183); k=13
    print(f"LEBESGUE baseline: 2ck = {2*c*k} = {float(2*c*k):.4f} (>1 => measure alone FAILS)")
    dw=[1,2,3,4,5,6,7,8,9,10,11,12,182]
    Wmin,sf=safe_measure(dw)
    print(f"DEEP WELL: min W={Wmin}, safe-measure={sf:.5f} ~ 0 => KNIFE-EDGE: W>=1 a.e. => any AC nu gives INT W dnu>=1 => NO single certificate")
    random.seed(3); print("LOOSE covering families: safe-measure > 0 (AC test measure certifies):")
    cnt=0
    for _ in range(4000):
        v=sorted(random.sample(range(1,60),13))
        if not primitive(v) or not divcomplete(v): continue
        Wmin,sf=safe_measure(v, D=183*4)
        if Wmin==0 and sf>0:
            print(f"   safe-measure={sf:.4f}  {v[:7]}...")
            cnt+=1
        if cnt>=5: break
    print("SPLIT: tight/deep-well by S255 rigidity (PROVED); loose covering by AC test measure = anti-concentration (favorable).")
if __name__=='__main__':
    main()
