# lrc14_tourmap_first-return-dynamics_probe_kps-S2-wf.py
# Follow-up probe: understand WHY M3 (lap-count parity) and M5 (cross-ratio) forbid classes.
# Key question: is the forbidden signal due to LONELINESS (tau being a gap-optimal lonely
# time) or just an intrinsic property of the arc-definition regardless of tau / S?
#
# Strategy:
#  (A) M3 with FREE tau (any rational tau in (0,1/2], not just LRC-constrained) on free S.
#      If the same classes are still forbidden -> the forbidding is INTRINSIC to the map,
#      NOT a loneliness constraint (weaker claim). If loneliness forbids MORE -> real lever.
#  (B) M5 is tau-free, so its forbidden class (regular pentagon) is intrinsic to the
#      cross-ratio definition. Check whether it's a genuine antisymmetry obstruction.
#  (C) Decompose M3's realized H-set: which (S,tau) give H>1, and is the parity twist the
#      source of non-transitivity? Compare M3 to its "no-twist" base (=M1a, transitive).
#  (D) Check whether the M3 forbidden classes are forbidden EXHAUSTIVELY over a much wider
#      speed range to confirm robustness.

from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations
from functools import reduce

def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
def frac(x):
    r=x-int(x); return r+1 if r<0 else r
def g(S,t): return min(nrm(v*t) for v in S)
def cand(S):
    S=sorted(set(S)); C=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): C.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): C.add(F(k,d)); k+=1
    C.add(F(1,2)); return C
def Mgap(S):
    b=F(0); at=None
    for t in cand(S):
        v=g(S,t)
        if v>b: b=v; at=t
    return b,at
def primitive(S): return reduce(gcd,S)==1

def is_tournament(n,edges):
    for i in range(n):
        for j in range(i+1,n):
            if ((i,j) in edges)==((j,i) in edges): return False
    return True
def num_3cycles(n,edges):
    cnt=0
    arc=lambda x,y:(x,y) in edges
    for a,b,c in combinations(range(n),3):
        if arc(a,b) and arc(b,c) and arc(c,a): cnt+=1
        if arc(a,c) and arc(c,b) and arc(b,a): cnt+=1
    return cnt
def ham_paths(n,edges):
    arc=lambda x,y:(x,y) in edges; cnt=0
    for perm in permutations(range(n)):
        if all(arc(perm[k],perm[k+1]) for k in range(n-1)): cnt+=1
    return cnt
def score_seq(n,edges):
    sc=[0]*n
    for (i,j) in edges: sc[i]+=1
    return tuple(sorted(sc))
def canon(n,edges):
    best=None
    for perm in permutations(range(n)):
        mat=[[0]*n for _ in range(n)]
        for (i,j) in edges: mat[perm[i]][perm[j]]=1
        key=tuple(mat[a][b] for a in range(n) for b in range(n) if a<b)
        if best is None or key<best: best=key
    return best
def all_tournament_classes(n):
    pairs=list(combinations(range(n),2)); classes={}
    for bits in range(2**len(pairs)):
        edges=set()
        for idx,(i,j) in enumerate(pairs):
            if (bits>>idx)&1: edges.add((i,j))
            else: edges.add((j,i))
        c=canon(n,edges)
        if c not in classes: classes[c]=(num_3cycles(n,edges),ham_paths(n,edges),score_seq(n,edges))
    return classes
FREE={n:all_tournament_classes(n) for n in (3,4,5)}

def classify(verts,arcfn):
    n=len(verts); edges=set()
    for a in range(n):
        for b in range(a+1,n):
            va,vb=verts[a],verts[b]
            r=arcfn(va,vb)
            if r is None: return None
            if r: edges.add((a,b))
            else: edges.add((b,a))
    if not is_tournament(n,edges): return None
    return canon(n,edges),num_3cycles(n,edges),ham_paths(n,edges),score_seq(n,edges)

def M3_arc_factory(tau):
    def arc(vi,vj):
        d=abs(vi-vj); k=int(d*tau); odd=(k%2==1)
        return (vi>vj)^odd
    return arc
def M5_arc(vi,vj):
    a=frac(F(vi+vj,2*vi)); b=frac(F(vi+vj,2*vj))
    if a==b: return None
    return a<b

def gen_speed_sets(n,maxspeed):
    return [S for S in combinations(range(1,maxspeed+1),n) if primitive(S)]

print("PROBE: first-return-dynamics M3 & M5  (kind-pasteur-S2-wf)")
print("="*70)

# ---------------------------------------------------------------------------
# (A) M3 at FREE tau vs LRC-gap tau. Does loneliness forbid MORE than free tau?
# ---------------------------------------------------------------------------
print("\n(A) M3 lap-count-parity: realized classes under FREE tau vs LRC-gap tau")
for n in (4,5):
    maxspeed={4:14,5:12}[n]
    SS=gen_speed_sets(n,maxspeed)
    # FREE tau: scan a dense rational grid of tau in (0,1) -> denominators up to D
    D=40
    free_taus=[F(p,q) for q in range(2,D+1) for p in range(1,q) if gcd(p,q)==1]
    realized_free=set(); realized_lonely=set()
    for S in SS:
        Ss=sorted(S)
        # lonely (gap) tau
        _,tg=Mgap(Ss)
        if tg is not None:
            r=classify(Ss,M3_arc_factory(tg))
            if r: realized_lonely.add(r[0])
    # free tau: use a FIXED speed set family but ALL taus -> realized over ALL S and ALL grid taus
    for S in SS:
        Ss=sorted(S)
        for tau in free_taus:
            r=classify(Ss,M3_arc_factory(tau))
            if r: realized_free.add(r[0])
    free_all=set(FREE[n].keys())
    forb_free=free_all-realized_free
    forb_lonely=free_all-realized_lonely
    print(f"  n={n}: realized@free-tau={len(realized_free)}/{len(free_all)} "
          f"(forbidden {len(forb_free)}); realized@gap-tau={len(realized_lonely)}/{len(free_all)} "
          f"(forbidden {len(forb_lonely)})")
    extra = forb_lonely - forb_free  # forbidden by loneliness but NOT by free tau
    print(f"     classes forbidden by LONELINESS but reachable at some free tau: {len(extra)}")
    for c in sorted(extra):
        n3,H,sc=FREE[n][c]; print(f"        -> H={H} c3={n3} score={sc}")

# ---------------------------------------------------------------------------
# (B) M3 INTRINSIC obstruction: characterize which classes M3 can EVER produce
#     over ALL tau (free) and ALL speed sets. Is the forbidden set a theorem of
#     the parity map itself? Print the realized H-multiset structure.
# ---------------------------------------------------------------------------
print("\n(B) M3 intrinsic: over ALL free tau, which classes appear & their H?")
for n in (4,5):
    maxspeed={4:16,5:14}[n]
    SS=gen_speed_sets(n,maxspeed)
    D=60
    free_taus=[F(p,q) for q in range(2,D+1) for p in range(1,q) if gcd(p,q)==1]
    realized={}
    for S in SS:
        Ss=sorted(S)
        for tau in free_taus:
            r=classify(Ss,M3_arc_factory(tau))
            if r and r[0] not in realized:
                realized[r[0]]=(r[1],r[2],r[3])
    free_all=set(FREE[n].keys())
    forb=free_all-set(realized.keys())
    print(f"  n={n} (speeds 1..{maxspeed}, tau denom<= {D}): realized {len(realized)}/{len(free_all)}")
    print(f"     realized H: {sorted(set(v[1] for v in realized.values()))}")
    if forb:
        for c in sorted(forb):
            n3,H,sc=FREE[n][c]; print(f"     FORBIDDEN: H={H} c3={n3} score={sc}")

# ---------------------------------------------------------------------------
# (C) M5 forbidden class probe: is the regular pentagon (H=15,c3=5) genuinely
#     impossible for the cross-ratio map? Prove by analyzing the arc rule.
#     M5 arc(i,j): compare frac((vi+vj)/(2vi)) vs frac((vi+vj)/(2vj)).
#     (vi+vj)/(2vi) = 1/2 + vj/(2vi); (vi+vj)/(2vj)=1/2+vi/(2vj).
#     frac(1/2 + vj/(2vi)) vs frac(1/2 + vi/(2vj)). Investigate the induced order.
# ---------------------------------------------------------------------------
print("\n(C) M5 structure: simplify the arc rule")
def m5_detail(vi,vj):
    a=F(vi+vj,2*vi); b=F(vi+vj,2*vj)
    fa=frac(a); fb=frac(b)
    return fa,fb
# show that arc(i,j) depends on ratio vj/vi only
print("  arc(i,j) compares frac(1/2 + vj/(2vi)) vs frac(1/2 + vi/(2vj))")
# Is M5 transitive when restricted? It produced 11/12 at n=5. Check WHY pentagon missing.
# Pentagon = regular tournament score (2,2,2,2,2). For M5 to give regular, need each vertex
# out-degree 2. Search exhaustively whether ANY 5 speeds give M5=regular.
print("  searching speeds 1..40 for M5 regular pentagon (score 2,2,2,2,2)...")
found=False
SS=gen_speed_sets(5,30)
for S in SS:
    Ss=sorted(S)
    r=classify(Ss,M5_arc)
    if r and r[3]==(2,2,2,2,2):
        print(f"    FOUND regular pentagon at S={Ss} : c3={r[1]},H={r[2]}")
        found=True; break
if not found:
    print(f"    NONE found in {len(SS)} sets (speeds<=30). Pentagon appears M5-forbidden.")

# Is M5 forbidding the pentagon because of an antisymmetry? Check: does M5 always
# produce a tournament with a DOMINANT or DOMINATED vertex (out-degree 0 or 4)?
print("  M5 score-sequence census (n=5, speeds<=20):")
from collections import Counter
sc_count=Counter()
for S in gen_speed_sets(5,20):
    r=classify(sorted(S),M5_arc)
    if r: sc_count[r[3]]+=1
for sc,ct in sorted(sc_count.items()): print(f"     score {sc}: {ct}")

# ---------------------------------------------------------------------------
# (D) M3 ROBUSTNESS: confirm the low-cycle classes (c3=1, c3=2) are forbidden
#     under gap-tau over a wide speed range at n=5.
# ---------------------------------------------------------------------------
print("\n(D) M3 @gap-tau robustness at n=5, wider speed range")
for maxspeed in (12,16,20):
    SS=gen_speed_sets(5,maxspeed)
    realized=set()
    for S in SS:
        Ss=sorted(S); _,tg=Mgap(Ss)
        if tg is None: continue
        r=classify(Ss,M3_arc_factory(tg))
        if r: realized.add(r[0])
    free_all=set(FREE[5].keys())
    forb=free_all-realized
    forbH=sorted((FREE[5][c][1],FREE[5][c][0]) for c in forb)
    print(f"  speeds<= {maxspeed} ({len(SS)} sets): realized {len(realized)}/12, "
          f"forbidden (H,c3)={forbH}")
