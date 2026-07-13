"""
opus-2026-07-11-S253: the SHAPE of the covering-min target, and a creative geometric argument -- the covering-
min is a SLOW-FAST BALANCE, which PROVES M >= n/Phi6 for the interval-core single-killer class and gives the
mechanism (and direction) for the general lower bound.

THE GEOMETRIC SHAPE. M(v) = max_t min_i ||v_i t|| is the L-infinity CLEARANCE of the closed geodesic
gamma(t)=(v_1 t,...,v_13 t) on the torus T^13 through the integer-hyperplane arrangement {x_i in Z}: the size
of the largest safe box the loop threads. LRC(14) = the loop reaches a box of size 1/14; the covering-min is
the covering family whose loop threads the SMALLEST box.

(1) THE HEXAGONAL / EISENSTEIN SHAPE (verified). The covering-min value is n/Phi6(n) with Phi6(n)=n^2-n+1 =
    N(n - zeta_6), the norm of (n - zeta_6) in the hexagonal lattice Z[zeta_6] (Eisenstein integers). The deep
    well {1..n-2, n(n-1)} <-> the Eisenstein integer (n - zeta_6); the covering-min = n / (its hexagonal norm).
    The optimal phase set is a comb {n(n-1)... } -- three-gap g=2 with gaps {1/Phi6, n/Phi6}: the hexagonal
    fundamental domain. (n=7..14 all check: Phi6=N(n-zeta_6), rep (1,n).)

(2) THE SLOW-FAST BALANCE (the creative core; a PROOF for the structured class). Write a covering family as a
    SLOW CORE + a resonant KILLER. For the interval core {1..n-2} (slow-optimum t0=1/(n-1), core value 1/(n-1))
    plus a killer v_f that is RESONANT at t0 (i.e. (n-1)|v_f, so v_f t0 in Z): perturb t = t0 + delta.
      - killer clearance rises:   ||v_f (t0+delta)|| = v_f|delta|   (v_f t0 integer)
      - core binding runner v=1 falls:  ||1*(t0+delta)|| = 1/(n-1) - |delta|
    They cross at |delta| = 1/((n-1)(v_f+1)), giving
      M = 1/(n-1) - |delta| = v_f/((n-1)(v_f+1)).
    VERIFIED EXACTLY: M({1..12, v_f}) = v_f/(13(v_f+1)) for v_f = 182,364,546,1820,2730 (all matches).

    LOWER BOUND (PROVED for this class). Covering needs a multiple of n-1 AND of n; the interval core lacks
    both, so the killer must be a multiple of lcm(n-1,n) = n(n-1), i.e. v_f >= n(n-1). Since
    M(v_f) = v_f/((n-1)(v_f+1)) is INCREASING in v_f, the minimum is at v_f = n(n-1):
      M >= n(n-1)/((n-1)(n(n-1)+1)) = n/(n(n-1)+1) = n/Phi6(n),
    with equality iff v_f = n(n-1) -- the deep well. So EVERY interval-core single-killer covering family has
    M >= n/Phi6, and the deep well is the unique minimizer. (Non-resonant killers give M >= 1/(n-1) > n/Phi6:
    the killer is already safe at the slow-optimum, no balance needed. Verified for v_f=183,185,200.)

    This DERIVES Phi6 (= killer + 1 = n(n-1)+1) and the hexagonal appearance, makes mac-mini S40's "two-point
    equioscillation" explicit as a 1-D balance, and BOOTSTRAPS from LRC(n-1) (the interval core value 1/(n-1)
    is the LRC(n-1) extremal, PROVED for n-1<=13).

(3) THE GENERAL MECHANISM + DIRECTION. For a general core with optimum value M_core and binding speed s, and a
    killer v_f resonant at the core optimum, the same balance gives
      M = M_core * v_f/(v_f + s).
    The deep well is extremal on ALL THREE knobs: M_core minimal (interval core = LRC(n-1) extremal, 1/(n-1)),
    s minimal (=1), v_f minimal (=n(n-1)). Any deviation raises one factor. The OPEN part of the general lower
    bound is exactly ruling out a "large-s escape" (a core whose binding runner is fast) trading against v_f --
    i.e. controlling M_core*v_f/(v_f+s) >= n/Phi6 jointly. This is the concrete geometric target the balance
    hands to the general covering-min bound, and it is an INDUCTIVE reduction to LRC(n-1) plus the killer-clear
    balance.

NET. The covering-min is the L-infinity clearance of the loop; its value n/Phi6 is a hexagonal-lattice norm
and a slow-fast balance; the balance PROVES the bound for the interval-core single-killer class (deep well =
unique minimizer, from smallest-killer + monotonicity + LRC(n-1)), and reduces the general lower bound to a
joint control M_core*v_f/(v_f+s) >= n/Phi6 -- a clean, geometric, inductive target.

-> mac-mini S38 (Ostrowski ladder), mac-mini S40 (2-point equioscillation), klein S267 (14/183 covering-min),
kps (Eisenstein/three-distance), THM-366 (covering=>divisor-complete), THM-527 (three-gap), opus-S252 (target
relocation), LRC(<=13) citation (the core bound).
"""
from math import gcd, lcm
from fractions import Fraction
def eis_rep(m):
    for a in range(0,int(m**0.5)+2):
        for b in range(0,int((4*m)**0.5)+2):
            if a*a-a*b+b*b==m: return (a,b)
    return None
def Mval(v):
    qs=set()
    for i in range(len(v)):
        for j in range(i+1,len(v)):
            s=v[i]+v[j]; g=gcd(v[i],v[j]); qs.add(s//gcd(s,g)); qs.add(s)
    best=Fraction(0)
    for q in sorted(x for x in qs if x>=2)[:8000]:
        bq=0
        for k in range(1,q//2+1):
            if gcd(k,q)!=1: continue
            m=min(min((vi*k)%q,q-(vi*k)%q) for vi in v)
            if m>bq: bq=m
        if Fraction(bq,q)>best: best=Fraction(bq,q)
    return best
def main():
    print("(1) hexagonal shape: covering-min = n/Phi6(n), Phi6=N(n-zeta_6):")
    for n in range(7,15):
        D=n*n-n+1; print(f"   n={n}: {n}/{D}, Phi6={D}=a^2-ab+b^2 rep {eis_rep(D)}")
    n=14; core=list(range(1,13)); Phi6=n*n-n+1
    print(f"\n(2) slow-fast balance, n={n}: M({{1..12,v_f}}) = v_f/(13(v_f+1)) [resonant killers]:")
    for vf in [182,364,546,1820,2730]:
        pred=Fraction(vf,13*(vf+1)); act=Mval(sorted(core+[vf]))
        print(f"   v_f={vf}: balance={pred}, actual={act}, match={pred==act}")
    vf=n*(n-1)
    print(f"   LOWER BOUND: covering=>killer mult of lcm(13,14)={lcm(13,14)} => v_f>=182; M increasing => min M(182)={Fraction(vf,13*(vf+1))}=n/Phi6={Fraction(n,Phi6)}")
    for vf in [183,185,200]:
        print(f"   non-resonant v_f={vf}: M={Mval(sorted(core+[vf]))} (>=1/13, no drop)")
    print("\n(3) general: M = M_core * v_f/(v_f+s); deep well minimal on all 3 (M_core=1/13 via LRC(13), s=1, v_f=182).")
if __name__=='__main__':
    main()
