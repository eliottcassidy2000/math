#!/usr/bin/env python3
"""
Which sub-cases of "compact covering => M >= 1/13" are PROVABLE?  (boxeph-S86)
=============================================================================
Proved reduction (sieve-margin, kps IX): M<1/13 => V covers all of {2..13}.
So the crux = "fully covering (2..14) + M<1/13 => rho>=13 (far element)".
The deep well (M=14/183<1/13) covers 13 via a FAR multiple (182=13*14).
HYPOTHESIS: if 13 is covered by a SMALL multiple, M>=1/13.  Test the clean claims:
  H1: covering + 13 in V (smallest 13-multiple = 13) => M >= 1/13 ?
  H2: covering + smallest-multiple-of-13 <= 13*K => M >= 1/13, find the K.
  H3: M<1/13 covering => rho>=13 (the full conjecture, hunt for violations).
  H4: the bounded compact sub-case: covering + v_max <= B => M >= 1/13, find B.
Find the SHARP provable sub-case (a finite check with an explicit bound).
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random, itertools

def exact_M(V):
    best=F(0); qs=set([13,14])
    for i in range(len(V)):
        for j in range(i,len(V)):
            s=V[i]+V[j]
            for d in range(1,s+1):
                if s%d==0: qs.add(d)
    for q in qs:
        for a in range(1,q):
            if gcd(a,q)==1:
                m=min(min((v*a)%q,q-(v*a)%q) for v in V); c=F(m,q)
                if c>best: best=c
    return best
def cov(V,n=14): return all(any(v%b==0 for v in V) for b in range(2,n+1))
def prim(V): return reduce(gcd,V)==1
def rho(V):
    Vs=sorted(V); return F(Vs[-1],Vs[-2])
def smallest_mult13(V):
    ms=[v for v in V if v%13==0]
    return min(ms) if ms else None

# Build a pool of PRIMITIVE COVERING families (compact and not), including deep-well variants.
def pool():
    fams=[]
    random.seed(31)
    cnt=0; att=0
    while cnt<8000 and att<400000:
        att+=1
        V=sorted(random.sample(range(1,80),13))
        if prim(V) and cov(V):
            fams.append(V); cnt+=1
    # deep-well family: {1..12, 13*14*k}
    for k in range(1,8):
        fams.append(sorted(list(range(1,13))+[182*k]))
    # {1..12, w} with w multiple of 182 (far killers)
    for w in [182,364,546,728]:
        fams.append(sorted(list(range(1,13))+[w]))
    # dilated: d*{1..12} u {13-mult}
    for d in range(1,8):
        for kk in [13,26,39,52,182]:
            V=sorted(set([d*i for i in range(1,13)]+[kk]))
            if len(V)==13 and prim(V) and cov(V):
                fams.append(V)
    # residue-type {1..11,13, 14*m}
    for m in range(1,15):
        V=sorted(list(range(1,12))+[13,14*m])
        if len(set(V))==13 and prim(V) and cov(V):
            fams.append(V)
    # dedup
    seen=set(); out=[]
    for V in fams:
        t=tuple(sorted(V))
        if t not in seen: seen.add(t); out.append(sorted(V))
    return out

if __name__=="__main__":
    fams=pool()
    print(f"pool: {len(fams)} primitive covering families\n")
    # H3 + H1 + H2 : relationship of M<1/13 to rho and to the 13-multiple size
    below13 = []
    H1_ok=True; H1_fail=[]
    for V in fams:
        M=exact_M(V)
        if M < F(1,13):
            below13.append((M,V,float(rho(V)),smallest_mult13(V)))
        if 13 in V:  # H1: 13 itself present
            if M < F(1,13):
                H1_ok=False; H1_fail.append((M,V))
    print(f"families with M<1/13: {len(below13)}")
    print("  (M, rho, smallest-multiple-of-13):")
    below13.sort(key=lambda t:t[0])
    for M,V,r,s13 in below13[:20]:
        print(f"    M={str(M):8s}({float(M):.5f}) rho={r:6.2f} min13mult={s13}  V={V}")
    allrho13 = all(r>=13 for _,_,r,_ in below13)
    print(f"\n  H3: ALL M<1/13 families have rho>=13 ?  {allrho13}")
    all13far = all((s13 is None or s13>=169) for _,_,_,s13 in below13)
    print(f"  small-13-mult: do all M<1/13 families cover 13 via a mult >= 169 (=13*13)? {all13far}")
    print(f"\n  H1: covering + (13 in V) => M>=1/13 ?  {H1_ok}")
    if H1_fail:
        for M,V in H1_fail[:5]: print(f"      FAIL M={M} V={V}")

    # H2: smallest 13-multiple <= 13*K => M>=1/13. Find threshold K.
    print("\n  H2: does 'smallest mult of 13 is SMALL' force M>=1/13?")
    for K in [1,2,3,4,6,9,12,13,14]:
        viol=[(M,V,s13) for M,V,r,s13 in below13 if s13 is not None and s13<=13*K]
        print(f"     K={K:2d} (min13mult<=169*K/13={13*K:3d}): M<1/13 families with min13mult<=13K: {len(viol)}"
              + ("" if not viol else f"  e.g. {viol[0][2]} in {viol[0][1]}"))

    # H4: bounded compact sub-case
    print("\n  H4: covering + v_max <= B => M >= 1/13 ?  (B = provable finite-check bound)")
    for B in [26,40,60,80,120,168]:
        viol=[(M,V) for M,V,r,s13 in below13 if max(V)<=B]
        print(f"     B={B:3d}: M<1/13 families with v_max<=B: {len(viol)}"
              + ("" if not viol else f"  e.g. M={viol[0][0]} V={viol[0][1]}"))
