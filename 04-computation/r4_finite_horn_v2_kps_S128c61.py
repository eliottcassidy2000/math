#!/usr/bin/env python3
"""r4_finite_horn_v2_kps_S128c61.py -- kind-pasteur S128 cont.61 (background).
r=4 FINITE HORN, efficient version.  v1 died on a per-core `seen` set of tuples.
Here: sort killers by kill-fraction DESCENDING, enumerate i<j<k<l with the sound
necessary condition sum frac >= 1, and use sorted order for a real early break
(for fixed i,j,k the best possible 4th is at index k+1).  Quadruples failing the
necessary condition cannot cover the core's safe set, hence are certified automatically.
PRINT DATA ONLY."""
import sys, itertools
sys.stdout.reconfigure(line_buffering=True)
def la(r,q):
    r%=q; return min(r,q-r)
QS=[(q,a) for q in range(15,41) for a in range(1,q)]
KB=400
CORES=[sorted(c) for c in itertools.combinations(range(1,13),9)]
tested=0; cert=0; bad=[]; skipped=0
for ci,P in enumerate(CORES):
    M=max(P); lo=13*M+1
    if lo>=KB: continue
    bits=[i for i,(q,a) in enumerate(QS) if all(la(p*a,q)>=-(-q//14) for p in P)]
    n=len(bits)
    if n==0: continue
    FULL=(1<<n)-1
    ent=[]
    for k in range(lo,KB):
        km=0
        for j,i in enumerate(bits):
            q,a=QS[i]
            if la(k*a,q)<-(-q//14): km|=(1<<j)
        ent.append((bin(km).count("1")/n,k,km))
    ent.sort(key=lambda e:-e[0])
    fr=[e[0] for e in ent]; kv=[e[1] for e in ent]; km=[e[2] for e in ent]
    L=len(ent)
    for i in range(L):
        if fr[i]*4<1.0: break
        for j in range(i+1,L):
            if fr[i]+fr[j]*3<1.0: break
            uij=km[i]|km[j]
            for k in range(j+1,L):
                if fr[i]+fr[j]+fr[k]*2<1.0: break
                uijk=uij|km[k]
                for l in range(k+1,L):
                    if fr[i]+fr[j]+fr[k]+fr[l]<1.0: break
                    u=uijk|km[l]
                    if u!=FULL:
                        cert+=1; tested+=1; continue
                    K=(kv[i],kv[j],kv[k],kv[l])
                    V=P+list(K)
                    tested+=1
                    if len(set(V))!=13: skipped+=1; continue
                    if not all(any(v%q==0 for v in V) for q in range(2,15)): skipped+=1; continue
                    bad.append((tuple(P),tuple(sorted(K))))
    if ci%40==0: print("  ... core %d/%d  tested=%d  bad=%d"%(ci,len(CORES),tested,len(bad)))
print()
print("### r=4 finite horn, KB = %d ###"%KB)
print("  quadruples passing the necessary condition (sum frac >= 1): %d"%tested)
print("  of those, certified by the small-modulus criterion: %d"%cert)
print("  covering-and-distinct failures (genuine UNCERTIFIED): %d"%len(bad))
print("  non-covering / non-distinct among the coverers: %d"%skipped)
if bad:
    print("  examples:")
    for P,K in bad[:12]: print("    core=%s K=%s"%(list(P),list(K)))
print("  NOTE: quadruples failing sum frac >= 1 cannot cover, hence are certified automatically.")
print("DONE")
