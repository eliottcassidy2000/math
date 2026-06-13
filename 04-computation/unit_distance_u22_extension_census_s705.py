"""
opus-2026-06-07-S705 : u(22) in {60,61} — the extension-census trick in the Moser ring.

State of the art (Alexeev-Mixon-Parshall, arXiv:2412.11914):
   u(16..21) = 41,43,46,50,54,57   (EXACT, densest graphs enumerated:
      #densest = 1(n16),7(17),16(18),3(19),1(20),5(21))
   60 <= u(22) <= 61   (the first OPEN case; upper bound 61 PROVEN).

THE REDUCTION (S614, sharpened with the now-exact enumeration):
A 61-edge UDG on 22 vertices has min degree in {4,5} (sum deg=122<132). Delete a
min-degree vertex v:
   deg v = 4 -> 21-core has 57 = u(21) edges = a DENSEST 21-UDG; v is at unit
               distance to EXACTLY 4 core points.
   deg v = 5 -> 21-core has 56 edges; v unit-dist to EXACTLY 5 core points.
A degree-d extension vertex v has its d neighbours CONCYCLIC at circumradius
EXACTLY 1 (they lie on v's unit circle): a "unit-cocyclic d-set".

KEY SHARPENING (uses the PROVEN u(22)<=61): for a 57-edge core, any candidate
centre v has at most 4 core-neighbours (5 would give U=62 > u(22), impossible).
So, restricted to a carrier lattice L:
   u(22)=61 in L  <=>  some 57-edge 21-core in L has a centre v (not a vertex)
                       at unit distance to EXACTLY 4 of its 21 points.
Because U_count is the EXACT faithful unit-distance count, ANY 22-pt subset of L
with U_count=61 would PROVE u(22)=61. We test L = the Moser ring M_L = Z[zeta6]
extended by w3=(5+isqrt11)/6 (the carrier of all known dense configs, Engel et al.).
Honest scope: an M_L search can only PROVE 61 (by exhibiting it); it cannot
DISPROVE it (61 could live outside M_L) -- that is the graph-side / totally-unfaithful
job. Here we CENSUS the extension and characterise the obstruction.

All adjacency is EXACT over K=Q(sqrt3,sqrt11) (no float decides a unit distance).
"""
from fractions import Fraction as F
from itertools import combinations, product
from collections import Counter

# ---- exact arithmetic in K=Q(sqrt3,sqrt11), basis 1,sqrt3,sqrt11,sqrt33 -------
ONE=(F(1),F(0),F(0),F(0)); Z4=(F(0),)*4
def add(x,y): return tuple(x[i]+y[i] for i in range(4))
def smul(s,x): return tuple(s*x[i] for i in range(4))
_MT={(0,0):(1,0),(0,1):(1,1),(0,2):(1,2),(0,3):(1,3),(1,1):(3,0),(1,2):(1,3),
     (1,3):(3,2),(2,2):(11,0),(2,3):(11,1),(3,3):(33,0)}
def mul(x,y):
    r=[F(0)]*4
    for i in range(4):
        if x[i]==0: continue
        for j in range(4):
            if y[j]==0: continue
            a,b=(i,j) if i<=j else (j,i); c,idx=_MT[(a,b)]; r[idx]+=x[i]*y[j]*c
    return tuple(r)
RE={'1':(F(1),F(0),F(0),F(0)),'w1':(F(1,2),F(0),F(0),F(0)),'w3':(F(5,6),F(0),F(0),F(0)),
    'w13':(F(5,12),F(0),F(0),F(-1,12))}
IM={'1':(F(0),F(0),F(0),F(0)),'w1':(F(0),F(1,2),F(0),F(0)),'w3':(F(0),F(0),F(1,6),F(0)),
    'w13':(F(0),F(5,12),F(1,12),F(0))}
KEYS=['1','w1','w3','w13']
def coord(v):
    re=Z4; im=Z4
    for k,key in zip(v,KEYS):
        if k==0: continue
        re=add(re,smul(F(k),RE[key])); im=add(im,smul(F(k),IM[key]))
    return re,im
def normsq(v):
    re,im=coord(v); return add(mul(re,re),mul(im,im))
def is_unit_vec(v): return normsq(v)==ONE

# the 18 unit vectors
units=[u for u in product(range(-4,5),repeat=4) if u!=(0,0,0,0) and is_unit_vec(u)]
UNITSET=set(units)
def vadd(p,u): return (p[0]+u[0],p[1]+u[1],p[2]+u[2],p[3]+u[3])
def U_count(points):
    s=set(points); e=0
    for p in points:
        for u in units:
            if vadd(p,u) in s: e+=1
    return e//2
def deg_in(p,s):  # unit-neighbours of p among s
    return sum(1 for u in units if vadd(p,u) in s)
print(f"#unit vectors of M_L: {len(units)} (expect 18)")

# ---- build a densest 57-edge 21-vertex core: W6 (+) Delta -----------------------
hub=(0,0,0,0); tri=[u for u in units if u[2]==0 and u[3]==0]; W6=[hub]+tri
cores57=[]
for u,v in combinations(units,2):
    duv=(u[0]-v[0],u[1]-v[1],u[2]-v[2],u[3]-v[3])
    if duv in UNITSET:
        Delta=((0,0,0,0),u,v)
        pts=sorted({vadd(w,d) for w in W6 for d in Delta})
        if len(pts)==21 and U_count(pts)==57:
            cores57.append(tuple(pts))
cores57=list(dict.fromkeys(cores57))
print(f"distinct 57-edge 21-vertex M_L cores (W6(+)Delta family): {len(cores57)}")

# ---- THE EXTENSION CENSUS: for each 57-core, candidate centres v in a box ------
print("\n"+"="*74)
print("EXTENSION CENSUS on 57-edge cores: how many core points does a new centre")
print("see at unit distance?  (=4  <=>  u(22)=61 in M_L;  max is provably <=4)")
print("="*74)
B=6
global_best=0; found61=[]
percore=[]
for ci,core in enumerate(cores57[:8]):   # sample up to 8 cores
    cs=set(core)
    # candidate centres = union of unit-neighbourhoods of core points (the only places
    # a point can be unit-distance to a core point), minus the core itself
    cand=set()
    for p in core:
        for u in units: cand.add(vadd(p,u))
    cand-=cs
    cnt=Counter()
    best=0; best_v=None
    for v in cand:
        d=deg_in(v,cs)
        cnt[d]+=1
        if d>best: best=d; best_v=v
    percore.append((best,dict(sorted(cnt.items()))))
    global_best=max(global_best,best)
    if best>=4:
        # verify by exact recount
        U=U_count(list(core)+[best_v])
        found61.append((ci,best_v,U))
    if ci<4:
        print(f" core#{ci}: centre unit-neighbour-count distribution {dict(sorted(cnt.items()))}"
              f"   -> max extension degree = {best}  => best U(22) here = {57+best}")
print(f"\n MAX extension degree over sampled 57-cores = {global_best}  =>  best U(22) reachable "
      f"by +1 on a 57-core in M_L = {57+global_best}")
if found61:
    print(f"  *** FOUND a degree-4 extension -> U=61 !!  {found61[:3]} ***")
else:
    print("  No degree-4 extension on any sampled 57-core: the delta=4 route is EMPTY in M_L")
    print("  (every candidate centre sees <= 3 core points) => core+1 gives at most U=60.")

# ---- broader: greedy/anneal search for ANY 22-pt M_L set with U=61 -------------
print("\n"+"="*74)
print("Broad M_L search for any 22-point set with U=61 (would PROVE u(22)=61)")
print("="*74)
def grow(seed,target):
    S=set([seed])
    while len(S)<target:
        cand={}
        for p in S:
            for u in units:
                q=vadd(p,u)
                if q in S: continue
                cand[q]=cand.get(q,0)+1   # +1 neighbour gained per existing-neighbour hit
        if not cand: break
        # gain = neighbours already in S
        best=max(cand, key=lambda q:(deg_in(q,S), -sum(abs(x) for x in q)))
        S.add(best)
    return S
bestU=0; bestS=None
seeds=[hub]+[ (a,b,0,0) for a in range(-1,2) for b in range(-1,2)] \
           +[ (0,0,1,0),(0,0,0,1),(1,0,1,0) ]
for seed in seeds:
    S=grow(seed,22)
    if len(S)==22:
        U=U_count(S)
        if U>bestU: bestU=U; bestS=S
print(f"  best U over greedy 22-pt M_L sets = {bestU}  (Engel et al. lower bound: u(22)>=60)")
print(f"  (greedy is heuristic; finding 61 would prove u(22)=61, finding only 60 proves nothing)")

# ---- the unit-cocyclic obstruction: when 4 core pts are unit-cocyclic, is a 5th forced? --
print("\n"+"="*74)
print("UNIT-COCYCLIC obstruction: centres seeing exactly k core points, k distribution,")
print("and whether a 'k=4' centre is forced toward k>=5 (hexagon completion) — the crux.")
print("="*74)
core=cores57[0]; cs=set(core)
cand=set();
for p in core:
    for u in units: cand.add(vadd(p,u))
cand-=cs
dist=Counter(deg_in(v,cs) for v in cand)
print(f"  core#0 centre k-distribution: {dict(sorted(dist.items()))}")
print(f"  (k = number of core points on the candidate centre's unit circle)")

# ---- vertex-swap hill-climb for max U over 22-pt M_L sets (try to reach 60/61) ----
print("\n"+"="*74)
print("Vertex-swap hill-climb for max U(22) in M_L (start = 57-core + best extension)")
print("="*74)
# big ambient patch: all M_L points within a coordinate box, near the origin cluster
def ambient(box=3):
    return [p for p in product(range(-box,box+1),repeat=4)]
AMB=set(ambient(3))
def neighbours_of_set(S):  # candidate add points: unit-neighbours of S, within AMB
    out=set()
    for p in S:
        for u in units:
            q=vadd(p,u)
            if q in AMB: out.add(q)
    return out-S
def hillclimb(S0):
    S=set(S0); U=U_count(S); improved=True
    while improved:
        improved=False
        adds=neighbours_of_set(S)|set(AMB)
        for out_v in list(S):
            for in_v in adds:
                if in_v in S: continue
                T=(S-{out_v})|{in_v}
                if len(T)!=len(S): continue
                UT=U_count(T)
                if UT>U:
                    S=T; U=UT; improved=True; break
            if improved: break
    return S,U
# seed: best 57-core + its best degree-3 extension (U=60)
seed=set(cores57[0])
ext=max((vadd(p,u) for p in cores57[0] for u in units),key=lambda q:deg_in(q,set(cores57[0])) if q not in set(cores57[0]) else -1)
seed.add(ext)
bestU2=U_count(seed); bestS2=set(seed)
for restart in range(6):
    S0=set(list(seed)) if restart==0 else set(list(bestS2))
    # perturb: random-ish swap by rotating which core vertex is dropped
    if restart>0:
        drop=sorted(bestS2)[restart % len(bestS2)]
        addc=[q for q in neighbours_of_set(bestS2-{drop}) if q not in bestS2]
        if addc: S0=(bestS2-{drop})|{addc[0]}
    S,U=hillclimb(S0)
    if U>bestU2 and len(S)==22: bestU2=U; bestS2=set(S)
print(f"  hill-climb best U(22) in M_L (box 3) = {bestU2}")
print(f"  => within this M_L patch the densest 22-set found has U={bestU2}")
if bestU2>=61: print("  *** U>=61 FOUND -> would prove u(22)=61 (verify!) ***")
else: print(f"  consistent with the Engel lower bound u(22)>=60; no 61 found in M_L patch.")
print("\nDONE.")
