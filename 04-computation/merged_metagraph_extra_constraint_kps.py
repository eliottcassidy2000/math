#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
THE RECONSTRUCTION GAP: what EXTRA constraint (beyond category+degree) pins the metagraph? Test the H-gradient.

kind-pasteur-2026-07-01-S13. Bucket data (category + blue-deg + black-deg) does NOT determine the metagraph
for n>=5 (many legal 2-swaps). Candidate extra constraint: every blue/black LINE connects tilings t, flip(t)
whose classes have COMPLEMENTARY H-position -- i.e. the line pairing respects H (Redei Ham-path count).
Test: is H(node(t)) + H(node(flip(t))) constant, or is flip an H-antitone involution?  And is the black graph
an H-GRADIENT (edges go between specific H-levels)? If the flip-involution's H-structure is rigid, THAT is the
missing constraint.
"""
import sys, itertools
from collections import defaultdict
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

def build(n):
    VERTS=[n-i for i in range(n)]
    TILES=[(x,y) for y in range(1,n-1) for x in range(n,y+1,-1)]
    m=len(TILES); tileIdx={t:i for i,t in enumerate(TILES)}
    TRANS=[tileIdx[(n-y+1,n-x+1)] for (x,y) in TILES]
    def gsym(bits): return all(TRANS[i]==i or bits[i]==bits[TRANS[i]] for i in range(m))
    def adj(bits):
        A=[[0]*n for _ in range(n)]
        for k in range(n-1): A[k][k+1]=1
        for i,(xL,yL) in enumerate(TILES):
            xi=VERTS.index(xL); yi=VERTS.index(yL)
            if bits[i]==0: A[xi][yi]=1
            else: A[yi][xi]=1
        return A
    perms=list(itertools.permutations(range(n)))
    def canon(A):
        best=None
        for p in perms:
            s=tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or s<best: best=s
        return best
    def H(A): return sum(1 for p in perms if all(A[p[t]][p[t+1]] for t in range(n-1)))
    T=[]
    for mask in range(1<<m):
        bits=[(mask>>k)&1 for k in range(m)]
        A=adj(bits)
        T.append(dict(mask=mask,canon=canon(A),g=gsym(bits),fl=mask^((1<<m)-1),H=H(A)))
    sigs=sorted(set(t['canon'] for t in T)); ci={s:i for i,s in enumerate(sigs)}
    for t in T: t['ci']=ci[t['canon']]
    Hof={t['ci']:t['H'] for t in T}
    bym={t['mask']:t for t in T}
    for t in T: t['fci']=bym[t['fl']]['ci']
    return dict(n=n,m=m,T=T,Hof=Hof,bym=bym)

for n in [4,5,6]:
    D=build(n); T=D['T']; Hof=D['Hof']; bym=D['bym']
    # for each line, H(t.class) and H(flip.class); test relationship
    sums=set(); prods=set(); pairs=[]
    seen=set()
    for t in T:
        pr=frozenset((t['mask'],t['fl']))
        if pr in seen: continue
        seen.add(pr)
        h1=Hof[t['ci']]; h2=Hof[bym[t['fl']]['ci']]
        sums.add(h1+h2); pairs.append((min(h1,h2),max(h1,h2),t['g']))
    Hmax=max(Hof.values())
    print(f"n={n}: H(class) ranges 1..{Hmax}. FLIP-line H-pairs {{H(t),H(flip t)}}:")
    print(f"   distinct H-sums over all lines: {sorted(sums)[:12]}{'...' if len(sums)>12 else ''}  (constant? {len(sums)==1})")
    # transitive (H=1) always flips to H=Hmax?
    trans_flip=[max(a,b) for (a,b,g) in pairs if min(a,b)==1]
    print(f"   lines touching H=1 (transitive) flip to H in {sorted(set(trans_flip))}  (Hmax={Hmax})")
    # is flip H-antitone: low H flips to high H?
    lohi=[(a,b) for (a,b,g) in pairs]
    corr_hi=sum(1 for a,b in lohi if a+b>Hmax)  # both high impossible mostly
    print(f"   sample H-pairs (min,max,blue): {[(a,b,g) for (a,b,g) in pairs[:8]]}")
    # blue lines: H-pair pattern
    bluepairs=sorted(set((a,b) for (a,b,g) in pairs if g))
    print(f"   BLUE-line H-pairs: {bluepairs}")
    print()
print("READ: if flip pairs low-H with high-H systematically, the H-gradient is the extra reconstruction constraint")
print("      (the flip-involution is H-antitone) -- the metagraph = bucket-degrees + H-antitone flip pairing.")
print("DONE.")
