#!/usr/bin/env python3
"""r4_finite_horn_kps_S128c61.py -- kind-pasteur S128 cont.61 (background).
THE r=4 FINITE HORN, made feasible by a SOUND pruning.
Brute force is ~10^9-10^10 quadruples.  But a quadruple can only be UNCERTIFIED if its four
kill-sets COVER the core's safe (q,a) set, which REQUIRES sum_i |kill(k_i)| >= |bits(P)|.
That is a necessary condition, so restricting to it keeps the check exhaustive.
Measured kill-fractions are mostly ~0.1-0.2, so quadruples with sum >= 1 are a small tail.
PRINT DATA ONLY."""
import sys, itertools
sys.stdout.reconfigure(line_buffering=True)
def la(r,q):
    r%=q; return min(r,q-r)
QS=[(q,a) for q in range(15,41) for a in range(1,q)]
KB=400
print("### r=4 finite horn, KB = %d (measure-horn samples gave 354-359) ###"%KB)
CORES=[sorted(c) for c in itertools.combinations(range(1,13),9)]
tot=0; cand=0; cert=0; bad=[]
for P in CORES:
    M=max(P); lo=13*M+1
    if lo>=KB: continue
    bits=[i for i,(q,a) in enumerate(QS) if all(la(p*a,q)>=-(-q//14) for p in P)]
    n=len(bits)
    if n==0: continue
    idx={b:j for j,b in enumerate(bits)}
    FULL=(1<<n)-1
    info={}
    for k in range(lo,KB):
        km=0
        for b in bits:
            q,a=QS[b]
            if la(k*a,q)<-(-q//14): km|=(1<<idx[b])
        info[k]=(km,bin(km).count("1")/n)
    ks=sorted(info, key=lambda k:-info[k][1])
    A=[k for k in ks if k%13==0]; B=[k for k in ks if k%14==0]
    seen=set()
    for x in A:
        for y in B:
            for i3 in range(len(ks)):
                k3=ks[i3]
                s3=info[x][1]+info[y][1]+info[k3][1]
                if s3+info[ks[0]][1]<1.0: break     # sorted desc: no k4 can reach 1
                for i4 in range(i3+1,len(ks)):
                    k4=ks[i4]
                    if s3+info[k4][1]<1.0: break
                    K=tuple(sorted({x,y,k3,k4}))
                    if len(K)!=4 or K in seen: continue
                    seen.add(K)
                    V=P+list(K)
                    if len(set(V))!=13: continue
                    if not all(any(v%q==0 for v in V) for q in range(2,15)): continue
                    cand+=1
                    u=0
                    for k in K: u|=info[k][0]
                    if u==FULL: bad.append((tuple(P),K))
                    else: cert+=1
    tot+=len(seen)
print("  quadruples passing the necessary condition (sum frac >= 1): %d"%cand)
print("  of those, certified: %d"%cert)
print("  UNCERTIFIED: %d"%len(bad))
if bad:
    print("  examples:")
    for P,K in bad[:10]: print("    core=%s K=%s"%(list(P),list(K)))
print("  (quadruples failing the necessary condition are certified automatically)")
print("DONE")
