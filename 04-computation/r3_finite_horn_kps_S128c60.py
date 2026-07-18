#!/usr/bin/env python3
"""r3_finite_horn_kps_S128c60.py -- kind-pasteur S128 cont.60 (background).
THE r=3 FINITE HORN.  Measure horn (r3_pair_removal) gives threshold 430.4 for the LARGEST
killer, so the finite horn must cover every covering family whose three killers are all
< 431.  Core = 10-subset of {1..12} (66 of them); covering forces a multiple of 13 AND a
multiple of 14 among the killers, which prunes the triple enumeration hard.
Certify with the small-modulus criterion (q <= 40) via bitmask intersection.
PRINT DATA ONLY."""
import sys, itertools
sys.stdout.reconfigure(line_buffering=True)
def la(r,q):
    r%=q; return min(r,q-r)
KB=431
QS=[(q,a) for q in range(15,41) for a in range(1,q)]
print("### r=3 finite horn ###")
print("  killer bound KB = %d (measure-horn threshold was 430.4) ; (q,a) pairs = %d"%(KB,len(QS)))
def mask(k):
    m=0
    for i,(q,a) in enumerate(QS):
        if la(k*a,q)>=-(-q//14): m|=(1<<i)
    return m
KM={k:mask(k) for k in range(1,KB)}
CORES=[sorted(c) for c in itertools.combinations(range(1,13),10)]
tot=0; cert=0; bad=[]
for P in CORES:
    M=max(P); lo=13*M+1
    if lo>=KB: continue
    cm=0
    for i,(q,a) in enumerate(QS):
        if all(la(p*a,q)>=-(-q//14) for p in P): cm|=(1<<i)
    rng=range(lo,KB)
    A=[k for k in rng if k%13==0]
    B=[k for k in rng if k%14==0]
    seen=set()
    for x in A:
        for y in B:
            for z in rng:
                K=tuple(sorted({x,y,z}))
                if len(K)!=3: continue
                if K in seen: continue
                seen.add(K)
                V=P+list(K)
                if len(set(V))!=13: continue
                if not all(any(v%q==0 for v in V) for q in range(2,15)): continue
                tot+=1
                m=cm & KM[K[0]] & KM[K[1]] & KM[K[2]]
                if m: cert+=1
                else: bad.append((tuple(P),K))
print("  covering families with all three killers < %d : %d"%(KB,tot))
print("  certified by the small-modulus criterion (q <= 40) : %d"%cert)
print("  UNCERTIFIED : %d"%len(bad))
if bad:
    print("  examples (first 10):")
    for P,K in bad[:10]: print("    core=%s K=%s"%(list(P),list(K)))
print("DONE")
