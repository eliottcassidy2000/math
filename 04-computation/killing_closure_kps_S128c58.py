#!/usr/bin/env python3
"""killing_closure_kps_S128c58.py -- kind-pasteur S128 cont.58.
(A) STRUCTURAL CLAIM to test: every killing pair has a killer residue in the DEGENERATE
    set  Deg(q) = {r : la(r,q) <= 2}  u  {q/2}  (q even).   If true, 'killing at q' forces
    q | one of  k-2,k-1,k,k+1,k+2  (either killer) or q | 2k -- ~10 fixed integers -- and
    an integer has few divisors in [15,60], so most moduli are non-killing.
(B) Is ANY family killing at EVERY modulus in [15,60]?  (my cont.58 'NONE' had a reporting
    bug -- len(Ks) was read after the loop; redo it correctly.)
PRINT DATA ONLY."""
import sys
sys.stdout.reconfigure(line_buffering=True)
def la(r,q):
    r%=q; return min(r,q-r)
def killing_set(P,q):
    thr=-(-q//14)
    safeA=[a for a in range(1,q) if all(la(p*a,q)>=thr for p in P)]
    if not safeA: return None
    badres=[set(r for r in range(q) if la(r*a,q)<thr) for a in safeA]
    K=set()
    for r in range(q):
        for s in range(r,q):
            if all((r in B) or (s in B) for B in badres): K.add((r,s))
    return K
print("### (A) structural claim: killing => a killer residue in {la<=2} u {q/2} ###")
print("  drop  q    |Kill|  violations of the claim")
viol_tot=0
for drop in range(1,13):
    P=[x for x in range(1,13) if x!=drop]
    for q in range(15,61):
        K=killing_set(P,q)
        if K is None: continue
        bad=0
        for r,s in K:
            okr = la(r,q)<=2 or (q%2==0 and r==q//2)
            oks = la(s,q)<=2 or (q%2==0 and s==q//2)
            if not (okr or oks): bad+=1
        if bad:
            viol_tot+=bad
            if viol_tot<=40: print("  %-5d %-4d %-7d %d"%(drop,q,len(K),bad))
print("  TOTAL violations across all 12 cores x q in [15,60]: %d"%viol_tot)
print()
print("### (B) families killing at EVERY modulus in [15,60]? ###")
allmods=list(range(15,61))
for drop in range(1,13):
    P=[x for x in range(1,13) if x!=drop]; M=max(P)
    Ks={}
    for q in allmods:
        K=killing_set(P,q)
        if K is not None: Ks[q]=K
    nq=len(Ks)
    tot=0; maxn=0; argmax=None; allkill=[]
    for a13 in range(13,70):
        for b14 in range(12,70):
            k1,k2=13*a13,14*b14
            if k1<=13*M or k2<=13*M or k1==k2: continue
            if len(set(P+[k1,k2]))!=13: continue
            tot+=1
            n=0
            for q,K in Ks.items():
                r,s=k1%q,k2%q
                if (min(r,s),max(r,s)) in K: n+=1
            if n>maxn: maxn=n; argmax=(k1,k2)
            if n==nq: allkill.append((k1,k2))
    print("  drop=%-3d usable moduli=%-3d families=%-6d  max #killing-moduli=%-3d at %s  killing-at-ALL: %d"%(
        drop,nq,tot,maxn,argmax,len(allkill)))
    if allkill: print("      *** COUNTEREXAMPLES to the q<=60 criterion:",allkill[:5])
print("DONE")
