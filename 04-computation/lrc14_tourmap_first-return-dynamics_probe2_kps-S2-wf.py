# lrc14_tourmap_first-return-dynamics_probe2_kps-S2-wf.py
# FAST probe: cache canonicalization; answer the 4 scientific questions about M3 & M5.
from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations
from functools import reduce, lru_cache
from collections import Counter

def frac(x):
    r=x-int(x); return r+1 if r<0 else r
def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
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
def gen_speed_sets(n,maxspeed):
    return [tuple(sorted(S)) for S in combinations(range(1,maxspeed+1),n) if primitive(S)]

# --- tournament tools, cached canon ---
_PERMS={n:list(permutations(range(n))) for n in (3,4,5)}
def edges_from_arc(verts,arcfn):
    n=len(verts); up=[]  # up[idx]=1 if a->b for a<b
    for a in range(n):
        for b in range(a+1,n):
            r=arcfn(verts[a],verts[b])
            if r is None: return None
            up.append(1 if r else 0)
    return tuple(up)

@lru_cache(maxsize=None)
def analyze_up(n, up):
    # up = upper-triangle orientation bits in order (a,b), a<b, lexicographic
    # build adjacency
    pos={}
    idx=0
    for a in range(n):
        for b in range(a+1,n):
            pos[(a,b)]=up[idx]; idx+=1
    def arc(x,y):
        if x<y: return pos[(x,y)]==1
        else: return pos[(y,x)]==0
    # 3-cycles
    n3=0
    for a,b,c in combinations(range(n),3):
        if arc(a,b) and arc(b,c) and arc(c,a): n3+=1
        if arc(a,c) and arc(c,b) and arc(b,a): n3+=1
    # ham paths
    H=0
    for perm in _PERMS[n]:
        if all(arc(perm[k],perm[k+1]) for k in range(n-1)): H+=1
    # score
    sc=[0]*n
    for a in range(n):
        for b in range(a+1,n):
            if pos[(a,b)]==1: sc[a]+=1
            else: sc[b]+=1
    score=tuple(sorted(sc))
    # canon: min over relabelings of upper-tri bits
    best=None
    for perm in _PERMS[n]:
        bits=[]
        for a in range(n):
            for b in range(a+1,n):
                pa,pb=perm[a],perm[b]
                if pa<pb: bits.append(pos[(pa,pb)])
                else: bits.append(1-pos[(pb,pa)])
        key=tuple(bits)
        if best is None or key<best: best=key
    return best,n3,H,score

def classify(verts,arcfn):
    up=edges_from_arc(verts,arcfn)
    if up is None: return None
    return analyze_up(len(verts),up)

# full free sets
def all_classes(n):
    res={}
    for bits in range(2**(n*(n-1)//2)):
        up=tuple((bits>>k)&1 for k in range(n*(n-1)//2))
        c,n3,H,sc=analyze_up(n,up)
        if c not in res: res[c]=(n3,H,sc)
    return res
FREE={n:all_classes(n) for n in (3,4,5)}

def M3_arc_factory(tau):
    def arc(vi,vj):
        d=abs(vi-vj); k=int(d*tau)
        return (vi>vj)^(k%2==1)
    return arc
def M5_arc(vi,vj):
    a=frac(F(vi+vj,2*vi)); b=frac(F(vi+vj,2*vj))
    if a==b: return None
    return a<b

print("PROBE2 (fast): first-return-dynamics M3 & M5  (kind-pasteur-S2-wf)")
print("="*70)

# (A) M3: forbidden under LRC-gap tau vs under FREE tau ?
print("\n(A) M3 lap-count-parity: realized under FREE-tau vs LRC-gap-tau")
for n in (4,5):
    maxspeed={4:14,5:12}[n]
    SS=gen_speed_sets(n,maxspeed)
    D=40
    free_taus=[F(p,q) for q in range(2,D+1) for p in range(1,q) if gcd(p,q)==1]
    realized_free=set(); realized_lonely=set()
    for S in SS:
        _,tg=Mgap(S)
        if tg is not None:
            r=classify(S,M3_arc_factory(tg))
            if r: realized_lonely.add(r[0])
        for tau in free_taus:
            r=classify(S,M3_arc_factory(tau))
            if r: realized_free.add(r[0])
    free_all=set(FREE[n].keys())
    forb_free=free_all-realized_free
    forb_lonely=free_all-realized_lonely
    print(f"  n={n}: free-tau realized {len(realized_free)}/{len(free_all)} (forb {len(forb_free)}); "
          f"gap-tau realized {len(realized_lonely)}/{len(free_all)} (forb {len(forb_lonely)})")
    extra=forb_lonely-forb_free
    print(f"     forbidden by LONELINESS but reachable at some free tau: {len(extra)}")
    for c in sorted(extra):
        n3,H,sc=FREE[n][c]; print(f"        H={H} c3={n3} score={sc}")

# (B) M3 intrinsic forbidden set over ALL free tau (moderate range)
print("\n(B) M3 intrinsic forbidden set (all free tau, speeds moderate)")
for n in (4,5):
    maxspeed={4:16,5:13}[n]
    SS=gen_speed_sets(n,maxspeed)
    D=48
    free_taus=[F(p,q) for q in range(2,D+1) for p in range(1,q) if gcd(p,q)==1]
    realized={}
    for S in SS:
        for tau in free_taus:
            r=classify(S,M3_arc_factory(tau))
            if r and r[0] not in realized: realized[r[0]]=(r[1],r[2],r[3])
    free_all=set(FREE[n].keys()); forb=free_all-set(realized)
    print(f"  n={n} speeds<= {maxspeed}, tau denom<= {D}: realized {len(realized)}/{len(free_all)}; "
          f"H seen {sorted(set(v[1] for v in realized.values()))}")
    for c in sorted(forb):
        n3,H,sc=FREE[n][c]; print(f"     FORBIDDEN(intrinsic): H={H} c3={n3} score={sc}")

# (C) M5 pentagon forbidden? exhaustive-ish search + score census
print("\n(C) M5 cross-ratio: is regular pentagon forbidden?")
found=None
for S in gen_speed_sets(5,30):
    r=classify(S,M5_arc)
    if r and r[3]==(2,2,2,2,2): found=(S,r); break
print(f"  regular pentagon (2,2,2,2,2): {'FOUND '+str(found) if found else 'NOT found in speeds<=30'}")
sc_count=Counter()
for S in gen_speed_sets(5,22):
    r=classify(S,M5_arc)
    if r: sc_count[r[3]]+=1
print("  M5 n=5 score census (speeds<=22):")
for sc,ct in sorted(sc_count.items()): print(f"     score {sc}: {ct}")

# (D) M3 @gap-tau robustness across speed ranges at n=5
print("\n(D) M3 @gap-tau robustness, n=5")
for maxspeed in (12,16,20):
    SS=gen_speed_sets(5,maxspeed)
    realized=set()
    for S in SS:
        _,tg=Mgap(S)
        if tg is None: continue
        r=classify(S,M3_arc_factory(tg))
        if r: realized.add(r[0])
    free_all=set(FREE[5].keys()); forb=free_all-realized
    forbH=sorted((FREE[5][c][1],FREE[5][c][0]) for c in forb)
    print(f"  speeds<= {maxspeed} ({len(SS)} sets): realized {len(realized)}/12; forbidden (H,c3)={forbH}")
print("\nDONE")
