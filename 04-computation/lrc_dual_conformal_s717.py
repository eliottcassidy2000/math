#!/usr/bin/env python3
"""
S717 — Dual conformal symmetry of the LRC: what can and cannot be expressed in spacetime.

TWO SPACES (the amplitudes analogy):
  SPACETIME = the cylinder S^1 x [0,1], runner worldlines x_i(t)=frac(v_i t), observer = the vertical
    strand at 0 (the LRC pure braid, oracle-S540o). Loneliness = a time t with every strand far from
    the observer strand. LOCAL/geometric.
  DUAL      = the speed/residue space; the multiplier group (Z/m)* = "perspective"=Galois (THM-403/439);
    the autocorrelation/homometry (THM-415/441). SPECTRAL/arithmetic.

Like dual conformal symmetry in N=4 SYM (a HIDDEN symmetry, invisible in spacetime, manifest in the dual
region/twistor variables), the LRC has:

(A) ORDINARY (spacetime) conformal symmetries of M(S) — visible on the worldlines:
      dilation v->lambda v (rescale time): M invariant ("conformal weight 0");
      reflection/time-reversal v->-v: M invariant; permutation S_n.
(B) DUAL conformal symmetry — HIDDEN (scrambles worldlines, no spacetime meaning), manifest in the dual:
      the MULTIPLIER a in (Z/m)*: orbit-M(a.S mod m) = orbit-M(S) EXACTLY. The INVERSION v->v^{-1} mod m
      is the SPECIAL-CONFORMAL generator: a=v_i^{-1} is the bad multiplier (THM-420) that pulls runner i
      onto the central band. The hard (transversal/Paley) core = where inversions cover ALL multipliers.

WHAT CANNOT be expressed in spacetime / the dual:
  - the OBSERVER/origin is pure spacetime: M is NOT translation-invariant, but the dual (distance
    multiset / autocorrelation) IS. So "where 0 is" is invisible to the dual.
  - M is NOT determined by the dual autocorrelation alone: homometric configs (same distance multiset)
    can have DIFFERENT M => loneliness needs genuine spacetime data beyond the dual.
  - the dual symmetry (multiplier/perspective) is invisible in spacetime: it preserves M while scrambling
    the worldlines arbitrarily.

No numpy/sympy. Exact M over Q via Fractions.
"""
from fractions import Fraction as Fr
from math import gcd
from itertools import combinations
import random

def norm(x):
    f=x-(x.numerator//x.denominator); return f if f<=Fr(1,2) else 1-f

def gap(speeds):
    """exact M(S)=max_t min_v ||v t|| over all t in Q (candidate critical times)."""
    V=[abs(v) for v in speeds]; C=set()
    for i in range(len(V)):
        vi=V[i]
        for k in range(0,2*vi+1):
            t=Fr(2*k+1,2*vi)
            if 0<t<=Fr(1,2): C.add(t)
        for j in range(i):
            for d in (vi+V[j],abs(vi-V[j])):
                if d==0: continue
                kk=1
                while Fr(kk,d)<=Fr(1,2): C.add(Fr(kk,d)); kk+=1
    best=Fr(0)
    for t in C:
        m=min(norm(v*t) for v in V)
        if m>best: best=m
    return best

def orbit_gap(speeds, m):
    """M restricted to the (Z/m)* witness orbit t=b/m, gcd(b,m)=1 (THM-403/411 witness set)."""
    best=Fr(0)
    for b in range(1,m):
        if gcd(b,m)!=1: continue
        mm=min(norm(v*Fr(b,m)) for v in speeds)
        if mm>best: best=mm
    return best

def dist_multiset(S):
    return tuple(sorted(abs(a-b) for a,b in combinations(sorted(S),2)))

def canon(S):
    """translation+reflection canonical form of a finite set (for congruence test)."""
    S=sorted(S); a=[x-S[0] for x in S]; b=sorted(a[-1]-x for x in a)
    return min(tuple(a),tuple(b))

if __name__=="__main__":
    rng=random.Random(0)
    print("="*84)
    print("S717 — Dual conformal symmetry of the LRC: what can / cannot be expressed in spacetime")
    print("="*84)

    # (A) ordinary (spacetime) conformal symmetries
    print("\n(A) ORDINARY (spacetime) conformal symmetries of M(S): dilation, reflection, permutation")
    ok_dil=ok_ref=0; T=300
    for _ in range(T):
        S=tuple(sorted(rng.sample(range(1,30),4)))
        M=gap(S); lam=rng.choice([2,3,5])
        if gap(tuple(lam*v for v in S))==M: ok_dil+=1
        if gap(tuple(-v for v in S))==M: ok_ref+=1
    print(f"  dilation M(lambda S)=M(S): {ok_dil}/{T} | reflection M(-S)=M(S): {ok_ref}/{T}  (S_n trivial)")
    print("  => these act on the WORLDLINES (rescale/reverse time, relabel strands): the visible conformal group.")

    # (B) DUAL conformal symmetry: multiplier on (Z/m)* witness orbit
    print("\n(B) DUAL conformal symmetry: multiplier a in (Z/m)* — orbit-M(a.S mod m) = orbit-M(S)")
    okmul=0; Tm=300
    for _ in range(Tm):
        n=6; m=2*n-1   # canonical shell
        S=tuple(sorted(rng.sample(range(1,3*m),n-1)))
        base=orbit_gap(S,m)
        a=rng.choice([a for a in range(2,m) if gcd(a,m)==1])
        Sa=tuple((a*v)%m for v in S)
        if orbit_gap(Sa,m)==base: okmul+=1
    print(f"  multiplier invariance orbit-M(a.S mod m)=orbit-M(S): {okmul}/{Tm} (EXACT)  [the hidden symmetry]")
    # inversion = special conformal: a=v^{-1} pulls v onto the band {1}
    m=13; print(f"  INVERSION (special conformal) v->v^-1 mod {m}: a=v^-1 sends v->1 (central band) = the bad")
    for v in (2,3,5):
        ainv=pow(v,-1,m); print(f"    v={v}: a=v^-1={ainv} mod {m}, a*v mod m = {(ainv*v)%m} (=1, on the band) => bad multiplier (THM-420)")
    print("  => the multiplier scrambles worldlines (no spacetime picture) yet preserves loneliness:")
    print("     dual conformal symmetry. Inversions are the dangerous generators; the transversal/Paley")
    print("     core is hard precisely because inversions {+-v_i^-1} cover ALL multipliers (no dodge).")

    # (C) the OBSERVER/origin is pure spacetime: M not translation-invariant, dual data IS
    print("\n(C) WHAT CANNOT be in the dual: the OBSERVER/origin. M is NOT translation-invariant,")
    print("    but the dual (distance multiset) IS.")
    diff=0; T2=300
    for _ in range(T2):
        S=tuple(sorted(rng.sample(range(1,25),4))); c=rng.choice([1,2,3,4])
        Sc=tuple(v+c for v in S)
        if dist_multiset(S)==dist_multiset(Sc) and gap(S)!=gap(Sc): diff+=1
    ex=(1,2); exc=(3,4)
    print(f"  S+c has the SAME distance multiset but DIFFERENT M in {diff}/{T2} trials")
    print(f"  e.g. S={ex} M={gap(ex)} vs S+2={exc} M={gap(exc)}; dual data identical {dist_multiset(ex)==dist_multiset(exc)}")
    print("  => the observer (origin 0, the fixed worldline) breaks translation: pure spacetime, dual-invisible.")

    # (D) M is NOT dual-determined: homometric non-congruent pairs with different M
    print("\n(D) WHAT CANNOT be in the dual: M itself. Homometric (same autocorrelation) configs differ in M.")
    groups={}
    for S in combinations(range(0,15),6):
        if S[0]!=0: continue
        if gcd(*[x for x in S if x>0]) not in (1,):  # primitive
            pass
        groups.setdefault(dist_multiset(S),[]).append(S)
    found=0
    for dm,sets in groups.items():
        reps={}
        for S in sets: reps.setdefault(canon(S),S)
        if len(reps)>=2:
            ss=list(reps.values())[:2]
            # shift +1 so all speeds positive & nonzero (same shift => fair), compare M
            A=tuple(v+1 for v in ss[0]); B=tuple(v+1 for v in ss[1])
            MA,MB=gap(A),gap(B)
            if MA!=MB:
                print(f"  homometric pair (same dist multiset, non-congruent): {ss[0]} vs {ss[1]}")
                print(f"    distance multiset {dm}")
                print(f"    M (shifted +1): {A}->{MA}  vs  {B}->{MB}   DIFFERENT => M not dual-determined")
                found+=1
            if found>=2: break
    if not found:
        print("  (no M-distinguishing homometric pair in the size-6 / range-15 search window)")

    print("\n" + "="*84)
    print("CAN be expressed in SPACETIME (worldlines on the cylinder): loneliness, the gap M, covering")
    print("  depth, the braid topology, the observer/origin. Local & geometric.")
    print("CANNOT (DUAL-only / hidden): the multiplier-perspective symmetry (preserves M, scrambles")
    print("  worldlines), homometry (configs the dual cannot tell apart), the autocorrelation flatness")
    print("  of the hard Paley/transversal core. The OBSTRUCTION lives in the dual; the OBSERVER lives")
    print("  in spacetime; M needs BOTH. Ordinary(spacetime) + dual conformal = the 'perspective' Yangian.")
    print("="*84)
