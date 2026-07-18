#!/usr/bin/env python3
"""killing_pairs_kps_S128c58.py -- kind-pasteur S128 cont.58.
CAN THE MIDDLE BAND BE CLOSED BY A COVERING SYSTEM OF MODULI?

The small-modulus criterion depends on the family ONLY through (core P, killer residues).
For core P and modulus q define
   SafeA(q,P) = { a : every core speed p has la(p*a mod q, q) >= ceil(q/14) }
   BadRes(q,a) = { r : la(r*a mod q, q) < ceil(q/14) }
   Killing(q,P) = { (r,s) : for EVERY a in SafeA(q,P), r in BadRes(q,a) or s in BadRes(q,a) }
A two-killer family with core P is uncertified at q iff its killer residues form a killing
pair.  If Killing(q,P) is EMPTY for some q, that single q certifies EVERY family with core
P -- a complete proof for that core, no search and no counting.
PRINT DATA ONLY."""
import sys
sys.stdout.reconfigure(line_buffering=True)
def la(r,q):
    r%=q; return min(r,q-r)
def analyse(P,q):
    thr=-(-q//14)
    safeA=[a for a in range(1,q) if all(la(p*a,q)>=thr for p in P)]
    if not safeA: return None
    badres=[set(r for r in range(q) if la(r*a,q)<thr) for a in safeA]
    killing=[]
    for r in range(q):
        for s in range(r,q):
            if all((r in B) or (s in B) for B in badres): killing.append((r,s))
    return len(safeA),len(killing),killing
print("### per (core, q): |SafeA| and |Killing pairs| ; q with EMPTY killing set proves that core ###")
print("  drop  q    |SafeA|  |Killing|   <- 0 killing = q alone certifies the whole core")
empties={}
for drop in range(1,13):
    P=[x for x in range(1,13) if x!=drop]
    best=None
    for q in range(15,61):
        res=analyse(P,q)
        if res is None: continue
        ns,nk,_=res
        if nk==0:
            best=(q,ns); break
    if best:
        empties[drop]=best
        print("  %-5d %-4d %-8d %-10d  *** EMPTY -- core {1..12}\{%d} fully certified by q=%d alone"%(
            drop,best[0],best[1],0,drop,best[0]))
    else:
        # report the minimum killing count achieved
        rows=[]
        for q in range(15,61):
            res=analyse(P,q)
            if res is None: continue
            rows.append((res[1],q,res[0]))
        rows.sort()
        print("  %-5d  no q<=60 with empty killing set; best 3: %s"%(
            drop,[(q,'|SafeA|=%d'%ns,'|Kill|=%d'%nk) for nk,q,ns in rows[:3]]))
print()
print("### cores fully certified by ONE modulus: %d of 12 ###"%len(empties))
print("   ",{d:q for d,(q,_) in empties.items()})
print()
print("### for the remaining cores: do TWO moduli suffice? (killing sets must not co-occur) ###")
for drop in range(1,13):
    if drop in empties: continue
    P=[x for x in range(1,13) if x!=drop]
    cands=[]
    for q in range(15,61):
        res=analyse(P,q)
        if res is None: continue
        cands.append((q,set(res[2])))
    found=None
    for i in range(len(cands)):
        for j in range(i+1,len(cands)):
            q1,K1=cands[i]; q2,K2=cands[j]
            # a killer pair (k1,k2) defeats BOTH iff (k1%q1,k2%q1) in K1 (unordered) and same for q2.
            # by CRT every combination of residues is realisable, so two moduli suffice only if
            # one of the killing sets is empty -- record the smallest product of sizes instead.
            pass
    sizes=sorted((len(K),q) for q,K in cands)[:4]
    print("  drop=%-3d smallest killing sets: %s"%(drop,[(q,n) for n,q in sizes]))
print("DONE")
