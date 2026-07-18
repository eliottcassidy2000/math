#!/usr/bin/env python3
"""killing_structure_kps_S128c58.py -- kind-pasteur S128 cont.58.
STRUCTURE OF THE KILLING PAIRS.  If every killing pair needs a killer residue in a tiny
set (0, or 0/+-1), then 'uncertified at every q in Q' forces each q in Q to put a killer
in that tiny set -- and an integer has boundedly many divisors, so a large enough Q gives
a contradiction.  That would CLOSE the middle band.  Measure it.  PRINT DATA ONLY."""
import sys
sys.stdout.reconfigure(line_buffering=True)
def la(r,q):
    r%=q; return min(r,q-r)
def killing(P,q):
    thr=-(-q//14)
    safeA=[a for a in range(1,q) if all(la(p*a,q)>=thr for p in P)]
    if not safeA: return None,None
    badres=[set(r for r in range(q) if la(r*a,q)<thr) for a in safeA]
    K=[]
    for r in range(q):
        for s in range(r,q):
            if all((r in B) or (s in B) for B in badres): K.append((r,s))
    return K,len(safeA)
print("### do killing pairs require a killer residue near 0? ###")
print("  drop  q   |Kill|  need-a-zero  need-la<=1  need-la<=2  max min-la over killing pairs")
rows=[]
for drop in [1,2,6,12]:
    for q in [26,27,28,39,41,42]:
        P=[x for x in range(1,13) if x!=drop]
        K,ns=killing(P,q)
        if K is None: continue
        nz=sum(1 for r,s in K if r==0 or s==0)
        n1=sum(1 for r,s in K if min(la(r,q),la(s,q))<=1)
        n2=sum(1 for r,s in K if min(la(r,q),la(s,q))<=2)
        mx=max(min(la(r,q),la(s,q)) for r,s in K)
        rows.append((drop,q,len(K),nz,n1,n2,mx))
        print("  %-5d %-3d %-7d %-13d %-11d %-11d %d"%(drop,q,len(K),nz,n1,n2,mx))
print()
print("### the killing pairs with LARGEST min-residue (the ones a divisor count cannot catch) ###")
for drop in [1,12]:
    for q in [26,28]:
        P=[x for x in range(1,13) if x!=drop]
        K,ns=killing(P,q)
        if K is None: continue
        K2=sorted(K,key=lambda rs:-min(la(rs[0],q),la(rs[1],q)))[:8]
        print("  drop=%d q=%d |SafeA|=%d : %s"%(drop,q,ns,
              [(r,s,'minla=%d'%min(la(r,q),la(s,q))) for r,s in K2]))
print()
print("### CRUCIAL: can a killing pair have BOTH residues far from 0 at MANY moduli at once? ###")
print("  scan actual killer pairs (13a,14b) and count at how many q in [15,60] they are killing")
from collections import Counter
cnt=Counter(); worst=[]
for drop in [1,12]:
    P=[x for x in range(1,13) if x!=drop]; M=max(P)
    Ks={}
    for q in range(15,61):
        K,ns=killing(P,q)
        if K is not None: Ks[q]=set(K)
    for a13 in range(13,60):
        for b14 in range(12,60):
            k1,k2=13*a13,14*b14
            if k1<=13*M or k2<=13*M or k1==k2: continue
            if len(set(P+[k1,k2]))!=13: continue
            n=0; qs=[]
            for q,K in Ks.items():
                r,s=k1%q,k2%q
                if (min(r,s),max(r,s)) in K: n+=1; qs.append(q)
            cnt[n]+=1
            if n>=len(Ks)-3: worst.append((drop,k1,k2,n,len(Ks)))
print("  distribution of '#moduli in [15,60] at which the family is killing':")
for n in sorted(cnt): print("     killing at %2d of ~%d moduli : %d families"%(n,len(Ks),cnt[n]))
print("  families killing at nearly ALL moduli:",worst[:5] if worst else "NONE")
print("DONE")
