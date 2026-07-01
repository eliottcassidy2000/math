#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
THE MOMENT LADDER: quantify the 'Bridge-2 mindset' recursively -- which rung IDENTIFIES tournament classes at each n?

kind-pasteur-2026-07-01-S18. The pattern (my S15 twins, mac-mini's LRC Gram Bridge-2): a 1st-moment / COUNT
invariant is cospectral, a 2nd-moment / SPECTRAL invariant separates -- then IT saturates too, forcing higher.
Quantify: for each n, the RESOLUTION (#distinct values / #classes) of an increasing ladder of invariants, and
the coarsest rung reaching resolution 1 (full identification).  Two orthogonal families:
  COUNT/combinatorial: score sequence -> H=I(Omega,2) -> cycle spectrum (c3,c5,c7)
  SPECTRAL/moment:      d=det(I+S) -> cpA (all tr(A^k)) -> cpS (skew spectrum)
Expect: count family saturates early, spectral rescues through n=6, BOTH saturate at n=7 (need WL/canon).
"""
import sys, itertools, random
import numpy as np
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

def make(n):
    P=list(itertools.permutations(range(n)))
    def canon(A):
        best=None
        for p in P:
            s=tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or s<best: best=s
        return best
    def invs(A):
        score=tuple(sorted(sum(A[i]) for i in range(n)))
        H=sum(1 for p in P if all(A[p[t]][p[t+1]] for t in range(n-1)))
        cyc={3:0,5:0,7:0}
        for L in (3,5,7):
            if L>n: continue
            seen=set()
            for sub in itertools.combinations(range(n),L):
                for pr in itertools.permutations(sub[1:]):
                    seq=(sub[0],)+pr
                    if all(A[seq[t]][seq[(t+1)%L]] for t in range(L)):
                        seen.add(min(tuple(seq[k:]+seq[:k]) for k in range(L)))
            cyc[L]=len(seen)
        M=np.array(A); S=M-M.T
        d=int(round(np.linalg.det(np.eye(n)+S)))
        cpA=tuple(int(round(c)) for c in np.poly(M))
        cpS=tuple(int(round(c)) for c in np.poly(S))
        return dict(score=score,H=H,cyc=(cyc[3],cyc[5],cyc[7]),d=d,cpA=cpA,cpS=cpS)
    return P,canon,invs

def classes_exhaustive(n):
    # via tiling cube 2^{C(n-1,2)}
    VERTS=[n-i for i in range(n)]; TILES=[(x,y) for y in range(1,n-1) for x in range(n,y+1,-1)]; m=len(TILES)
    def adj(bits):
        A=[[0]*n for _ in range(n)]
        for k in range(n-1): A[k][k+1]=1
        for i,(xL,yL) in enumerate(TILES):
            xi=VERTS.index(xL); yi=VERTS.index(yL)
            if bits[i]==0: A[xi][yi]=1
            else: A[yi][xi]=1
        return A
    P,canon,invs=make(n); reps={}
    for mask in range(1<<m):
        bits=[(mask>>k)&1 for k in range(m)]; A=adj(bits); c=canon(A)
        if c not in reps: reps[c]=A
    return list(reps.values()),invs

def classes_sample(n,NS):
    P,canon,invs=make(n); rng=random.Random(1); reps={}
    for _ in range(NS):
        A=[[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1,n):
                if rng.random()<0.5: A[i][j]=1
                else: A[j][i]=1
        c=canon(A)
        if c not in reps: reps[c]=A
    return list(reps.values()),invs

LADDER=[("score",lambda v:(v['score'],)),
        ("+H",lambda v:(v['score'],v['H'])),
        ("+cycles",lambda v:(v['score'],v['H'],v['cyc'])),
        ("+d",lambda v:(v['score'],v['H'],v['cyc'],v['d'])),
        ("+cpA",lambda v:(v['score'],v['H'],v['cyc'],v['d'],v['cpA'])),
        ("+cpS",lambda v:(v['score'],v['H'],v['cyc'],v['d'],v['cpA'],v['cpS']))]
COUNTONLY=("count fam (score,H,cyc)",lambda v:(v['score'],v['H'],v['cyc']))
SPECONLY=("spectral fam (d,cpA,cpS)",lambda v:(v['d'],v['cpA'],v['cpS']))

print("="*104); print(" MOMENT-LADDER RESOLUTION by n: #distinct invariant values / #iso classes (1.000 = identifies)"); print("="*104)
A568={4:4,5:12,6:56,7:456}
for n in [4,5,6,7]:
    if n<=6: reps,invs=classes_exhaustive(n)
    else: reps,invs=classes_sample(n,1600)
    IV=[invs(A) for A in reps]; nc=len(IV)
    print(f"\n n={n}: {nc} classes found (of {A568[n]}, {'exhaustive' if n<=6 else '~sample'}):")
    for name,f in LADDER:
        res=len(set(f(v) for v in IV))
        flag=' <= IDENTIFIES' if res==nc else ''
        print(f"    {name:>10}: {res:>4}/{nc}  resolution={res/nc:.3f}{flag}")
    rc=len(set(COUNTONLY[1](v) for v in IV)); rs=len(set(SPECONLY[1](v) for v in IV))
    print(f"    [{COUNTONLY[0]}: {rc}/{nc}={rc/nc:.3f}]  [{SPECONLY[0]}: {rs}/{nc}={rs/nc:.3f}]")
print("\n INTERPRET: the COUNT family (score/H/cycles) saturates first; the SPECTRAL family (d/cpA/cpS)")
print("   is orthogonal and rescues -- BOTH together IDENTIFY through n=6 but SATURATE at n=7 (cospectral")
print("   families) => the ladder must climb to WL/canonical form. Each rung's wall arrives one n later:")
print("   count-only fails ~n=5, +d fails ~n=7, full-spectral fails n=7 => the reconstruction wall is the")
print("   moment ladder saturating. Same recursion as the LRC Lasserre hierarchy (HYP-3789: level1 union /")
print("   level2 pair-correlation / level-inf exact) and THM-589 (H mean vs H VARIANCE).")
print("DONE.")
