#!/usr/bin/env python3
"""
lrc14_base_tournament_ostrowski_klein_S269.py
=============================================
klein-2026-07-12-S269  (owner request: change-of-base -> tournament + Ostrowski)

The owner's frame (opus-S250): time t=p/q => runner i at DIGIT (v_i p mod q) on a q-cell
dial; lonely = all digits in the middle band; clean base for AP {1..13} is q=n=14 (t=1/14),
extremal covering family (deep well {1..12,182}) forced to base q=183=Phi6(14), t*=14/183.

This script exhibits the THREE FACES of that one object:
 (A) BASE (arithmetic): the FORCED-BASE SPECTRUM -- which bases q a covering family clears at
     ("blocked at 14 -> first clean at q_first -> best at argq"), value ceil(q/14)/q.
 (B) TOURNAMENT: at the optimum the runners are residues mod q; the winding tournament
     (u->v iff (s_u - s_v) mod q in {1..floor((q-1)/2)}) is a CIRCULANT. We compute the
     digit GAP structure (three-gap/Steinhaus) and the cyclic-triangle count c3 (tournament
     regularity: c3 = C(N,3) - sum_v C(outdeg,2); max for the regular circulant). Observer
     (speed 0) included as vertex 0; loneliness = observer's nearest-runner gap.
 (C) OSTROWSKI: continued fractions of the optimal times -- 14/183=[0;13,14] etc.
"""
import math
from fractions import Fraction
from itertools import combinations

def dist0(r,q): return min(r,q-r)
def exact_M(v):
    n=len(v); best=Fraction(0); argq=None; argp=None
    qs=set()
    for a in range(n):
        for b in range(a,n): qs.add(v[a]+v[b])
    for q in sorted(qs):
        for p in range(1,q):
            mn=min(dist0((vl*p)%q,q) for vl in v)
            val=Fraction(mn,q)
            if val>best: best=val; argq=q; argp=p
    return best,argq,argp
def covering(v): return all(any(x%d==0 for x in v) for d in range(2,15))
def cf(fr):  # continued fraction of a Fraction in [0,1)
    a=[]; x=fr
    for _ in range(12):
        ai=x.numerator//x.denominator if x>=1 else 0
        # for x<1 first term is 0
        ai = int(x)  # floor
        a.append(ai)
        x = x-ai
        if x==0: break
        x = 1/x
    return a
def clears(v,q):  # smallest-band clear: exists p with all residues avoiding {0,+-1,..,+-(ceil(q/14)-1)}
    m=math.ceil(q/14)-1
    for p in range(1,q):
        if all(dist0((vi*p)%q,q)>m for vi in v): return p
    return None

fams=[
    ("AP {1..13}", list(range(1,14))),
    ("deep well {1..12,182}", list(range(1,13))+[182]),
    ("kps 2-block {1,2,3,4,10..18}", [1,2,3,4,10,11,12,13,14,15,16,17,18]),
    ("compressed 2*{1..12}u{13}", sorted([2*i for i in range(1,13)]+[13])),
]

print("="*76)
print("(A) FORCED-BASE SPECTRUM: clean bases q, value ceil(q/14)/q, first vs best")
print("="*76)
for nm,v in fams:
    M,argq,argp=exact_M(v)
    cov=covering(v)
    # sweep bases 14..min(argq,200) for clearing (smallest-band)
    clean=[q for q in range(14,min(argq,200)+1) if q%14!=0 and clears(v,q) is not None]
    first=clean[0] if clean else None
    blocked14 = (min(dist0((vi*1)%14,14) for vi in v)==0) if any(x%14==0 for x in v) else False
    print(f"  {nm:30s} cov={str(cov):5s} M={str(M):>8s} argq={argq:4d} (t={argp}/{argq})")
    print(f"     base-14 blocked (a runner at digit 0)? {any(x%14==0 for x in v)};  first clean base={first};  best base(argq)={argq}")

print()
print("="*76)
print("(B) TOURNAMENT + THREE-GAP at the clean base (optimum). Observer=vertex 0.")
print("="*76)
def winding_tournament_stats(v, q, p):
    # vertices: observer (speed 0) + runners; positions = residues mod q at t=p/q
    speeds=[0]+list(v)
    pos=[(s*p)%q for s in speeds]           # digit positions on the q-dial
    N=len(speeds)
    # winding tournament: u->w iff (pos_u - pos_w) mod q in {1..floor((q-1)/2)}; antipode=tie
    half=(q-1)//2
    out=[0]*N; ties=0
    for a in range(N):
        for b in range(N):
            if a==b: continue
            d=(pos[a]-pos[b])%q
            if 1<=d<=half: out[a]+=1
            elif q%2==0 and d==q//2: ties+=1  # counted twice (a,b),(b,a)
    c3 = math.comb(N,3) - sum(math.comb(o,2) for o in out)   # cyclic triangles
    # three-gap: sorted distinct positions, gaps
    sp=sorted(set(pos))
    gaps=sorted(set((sp[(i+1)%len(sp)]-sp[i])%q for i in range(len(sp))))
    # observer isolation: nearest runner distance from position 0
    obs=pos[0]
    others=[pos[i] for i in range(1,N)]
    iso=min(dist0((o-obs)%q,q) for o in others)
    return dict(N=N,pos=pos,c3=c3,ties=ties//2,gaps=gaps,ngaps=len(gaps),iso=Fraction(iso,q),outdeg=sorted(out))
for nm,v,q,p in [("AP {1..13}",list(range(1,14)),14,1),
                 ("deep well {1..12,182}",list(range(1,13))+[182],183,None)]:
    if p is None:
        _,q2,p=exact_M(v); q=q2
    st=winding_tournament_stats(v,q,p)
    N=st['N']; c3max=math.comb(N,3)  # loose ref
    print(f"  {nm} at base q={q}, t={p}/{q}:")
    print(f"     digit positions (residues): {sorted(st['pos'])}")
    print(f"     distinct GAP lengths: {st['gaps']}  => {st['ngaps']}-gap  {'(ONE-GAP = regular!)' if st['ngaps']==1 else '(Steinhaus <=3-gap)'}")
    print(f"     winding tournament: cyclic-triangles c3={st['c3']}, apex ties={st['ties']}, outdeg profile={st['outdeg']}")
    print(f"     observer isolation (nearest runner) = {st['iso']} = {float(st['iso']):.5f} = M  ({'regular circulant' if st['ngaps']==1 else 'circulant'})")

print()
print("="*76)
print("(C) OSTROWSKI: continued fractions of the optimal times (encode n=14?)")
print("="*76)
for nm,fr in [("1/14 (AP, tight)",Fraction(1,14)),("3/41 (3rd-mediant)",Fraction(3,41)),
              ("2/27 (mediant)",Fraction(2,27)),("14/183 (deep well)",Fraction(14,183)),
              ("1/13 (compressed)",Fraction(1,13))]:
    print(f"  t*={nm:22s} = {str(fr):>7s} = {float(fr):.5f}  CF = {cf(fr)}")
print("  note: 14/183 = [0;13,14] = [0; n-1, n]  -- the extremal time's CF encodes n.")
print("  Ostrowski ladder M_k = k/(13k+1): k=1->1/14 (AP), k=14->14/183 (deep well). base = 13k+1.")
print("\ndone.")
