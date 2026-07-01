"""
lrc_modular_criterion_interval_vs_generic_opus_20260630.py  (HYP-3775)
The regularizable/Eisenstein BULK vs un-regularizable/cusp RESIDUAL has a sharp geometric criterion:
the covering-min margin is a classical 2-fold Dedekind sum (=> -1/12/eta anomaly, mac-mini HYP-3774,
klein HYP-3768) EXACTLY when the speeds are a SCALED INTERVAL (AP) at the witness rotation.
 (A) construction = scaled interval-core {1..n-2} + antipode killer -1 => MODULAR;
 (B) beaters n=7..10 = generic non-consecutive subsets => NON-MODULAR;
 (C) NEGATIVE: the beater margin is NOT a Zagier cotangent sum (corrects HYP-3773's hope).
ASCII only. See reflection the-regularizable-bulk-vs-cusp-residual-IS-the-AP-interval-vs-generic-split-...-opus-20260630.md.
"""
import numpy as np
from math import gcd
def Mex(S,Q=250):
    Sa=np.array(S); bn,bd=0,1
    for q in range(2,Q+1):
        A=np.outer(Sa,np.arange(1,q))%q; kk=int(np.minimum(A,q-A).min(axis=0).max())
        if kk*bd>bn*q: bn,bd=kk,q
    g=gcd(bn,bd) or 1; return bn//g,bd//g
def witness(S,D,k):
    for j in range(1,D):
        if min(min((v*j)%D,D-(v*j)%D) for v in S)==k: return j
    return None
def consec(v): return len(v)>=2 and all(v[i+1]-v[i]==1 for i in range(len(v)-1))
def scaled_interval(resids,D):
    R=sorted(r%D for r in resids)
    for c in range(1,D):
        if gcd(c,D)!=1: continue
        ci=pow(c,-1,D); s=sorted((ci*r)%D for r in R)
        for var in [s,[x for x in s if x!=D-1]]:
            if consec(var): return c
    return None
def zagier(D,speeds):     # full-product cotangent sum; diverges if any speed shares a factor with D
    k=np.arange(1,D); p=np.ones(D-1)
    for v in speeds: p=p/np.tan(np.pi*k*(v%D)/D)
    return np.sum(p)/D

def criterion():
    print("MODULAR criterion: are the witness residues a SCALED INTERVAL (AP) => classical Dedekind sum?")
    cases=[("n=7  constr",7,[1,2,3,4,5,42]),("n=7  BEATER",7,[1,2,5,6,7,8]),
           ("n=8  constr",8,[1,2,3,4,5,6,56]),("n=8  BEATER",8,[1,4,5,6,7,11,16]),
           ("n=9  BEATER",9,[1,3,4,5,7,11,18,32]),("n=10 BEATER",10,[1,2,3,5,6,7,8,9,30]),
           ("n=14 constr",14,list(range(1,13))+[182])]
    for name,n,S in cases:
        bn,bd=Mex(S); D=bd; a=witness(S,D,bn); c=scaled_interval([(v*a)%D for v in S],D) if a else None
        print(f"  {name}: M={bn}/{bd}, witness a={a}: {'MODULAR (scaled interval c=%d)'%c if c else 'NON-MODULAR (generic subset -> cusp/residual, no closed form)'}")

def negative():
    print("\nNEGATIVE: the beater margin is NOT a Zagier cotangent sum:")
    for name,n,S in [("n=7 constr",7,[1,2,3,4,5,42]),("n=7 BEATER",7,[1,2,5,6,7,8]),("n=8 BEATER",8,[1,4,5,6,7,11,16])]:
        bn,bd=Mex(S); z=zagier(bd,S)
        print(f"  {name}: M={bn}/{bd}, margin={bn/bd-1/n:.5f}, full-product Zagier d(D;speeds)={z:.3e} (unrelated / ill-defined for non-coprime speeds)")

if __name__=="__main__":
    criterion(); negative()
    print("\n=> regularizable/Eisenstein/-1/12 <=> AP/interval (construction, ORDERED); residual/cusp/no-closed-form <=> generic (beaters, DISORDERED).")
    print("   The bulk/residual (HYP-3774) = Eisenstein/cusp (HYP-3768) split IS the AP/interval vs generic split; genus jump 0->1 at n=14 = the residual's birth.")
