# lrc14_tourmap_loneliness_lever_kps-S2-wf.py
# Hunt for a map where LONELINESS genuinely forbids MORE than free-tau does
# (i.e. the loneliness constraint, not just the arc-rule algebra, removes iso classes).
#
# We compare, for each candidate map: realized@gap-tau  vs  realized@all-free-tau.
# A "loneliness lever" exists iff   forbidden@gap-tau  STRICTLY contains  forbidden@free-tau.
#
# New maps under the first-return theme:
#  M7  BINDING-PARTNER tournament: at the lonely optimum tau*, the gap M is attained;
#      a runner v is "binding" if ||v tau*|| == M (distance exactly M from 0). THM-524 says
#      >=2 runners bind. Orient i->j by which of the pair reaches its NEXT wall (||v t||=M
#      again, the next return to the binding circle) FIRST as t increases past tau*. Genuine
#      first-return at the lonely time. (Only defined at the lonely tau* -> pure loneliness map.)
#  M8  GAP-PHASE tournament: arc i->j iff the signed phase frac(vi tau*) is "ahead" of
#      frac(vj tau*) in the first-return-to-0 race weighted by speed: next time runner v hits
#      EXACT 0 (integer) after tau* is t_v = ceil(vi tau*)/vi ; orient by t_v. At a lonely
#      time none sits AT 0, so all t_v>tau* well-defined. Pure-loneliness map.
#  M9  control: same M8 rule but at an ARBITRARY (free) tau, to get the free baseline.

from fractions import Fraction as F
from math import gcd, ceil
from itertools import combinations, permutations
from functools import reduce, lru_cache

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

_PERMS={n:list(permutations(range(n))) for n in (3,4,5)}
@lru_cache(maxsize=None)
def analyze_up(n,up):
    pos={}; idx=0
    for a in range(n):
        for b in range(a+1,n): pos[(a,b)]=up[idx]; idx+=1
    def arc(x,y): return pos[(x,y)]==1 if x<y else pos[(y,x)]==0
    n3=0
    for a,b,c in combinations(range(n),3):
        if arc(a,b) and arc(b,c) and arc(c,a): n3+=1
        if arc(a,c) and arc(c,b) and arc(b,a): n3+=1
    H=0
    for perm in _PERMS[n]:
        if all(arc(perm[k],perm[k+1]) for k in range(n-1)): H+=1
    sc=[0]*n
    for a in range(n):
        for b in range(a+1,n):
            if pos[(a,b)]==1: sc[a]+=1
            else: sc[b]+=1
    best=None
    for perm in _PERMS[n]:
        bits=[]
        for a in range(n):
            for b in range(a+1,n):
                pa,pb=perm[a],perm[b]
                bits.append(pos[(pa,pb)] if pa<pb else 1-pos[(pb,pa)])
        key=tuple(bits)
        if best is None or key<best: best=key
    return best,n3,H,tuple(sorted(sc))
def edges_from_arc(verts,arcfn):
    n=len(verts); up=[]
    for a in range(n):
        for b in range(a+1,n):
            r=arcfn(verts[a],verts[b])
            if r is None: return None
            up.append(1 if r else 0)
    return tuple(up)
def classify(verts,arcfn):
    up=edges_from_arc(verts,arcfn)
    if up is None: return None
    return analyze_up(len(verts),up)
def all_classes(n):
    res={}
    for bits in range(2**(n*(n-1)//2)):
        up=tuple((bits>>k)&1 for k in range(n*(n-1)//2))
        c,n3,H,sc=analyze_up(n,up)
        if c not in res: res[c]=(n3,H,sc)
    return res
FREE={n:all_classes(n) for n in (3,4,5)}

# ---- maps ----
def M7_arc_factory(tau, M):
    # next time runner v's distance-to-integer returns to value M, increasing past tau.
    # ||v t|| = M  <=> v t = integer +/- M. Going forward, the next such t.
    def nexttime(v):
        best=None
        # candidate t = (m + M)/v and (m - M)/v ... we want frac(v t)=M or 1-M
        for targ in (M, 1-M):
            # v t = k + targ  => t=(k+targ)/v ; smallest with t>tau
            base = v*tau - targ
            k=int(base)
            t=(F(k)+targ)/v
            while t<=tau:
                k+=1; t=(F(k)+targ)/v
            if best is None or t<best: best=t
        return best
    def arc(vi,vj):
        ti,tj=nexttime(vi),nexttime(vj)
        if ti==tj: return None
        return ti<tj
    return arc

def M8_arc_factory(tau):
    # next time runner v hits EXACT integer (0 mod 1) after tau: t_v = ceil(v*tau)/v ... but
    # if v*tau is already integer (parked) skip. Use smallest t>tau with v t in Z.
    def nexttime(v):
        x=v*tau
        k=int(x)
        if F(k)==x: k+=1
        else: k+=1
        t=F(k)/v
        while t<=tau:
            k+=1; t=F(k)/v
        return t
    def arc(vi,vj):
        ti,tj=nexttime(vi),nexttime(vj)
        if ti==tj: return None
        return ti<tj
    return arc

print("LONELINESS-LEVER hunt (M7 binding-return, M8 hit-zero-return)  kps-S2-wf")
print("="*70)

for n in (4,5):
    maxspeed={4:14,5:12}[n]
    SS=gen_speed_sets(n,maxspeed)
    free_all=set(FREE[n].keys())
    D=40
    free_taus=[F(p,q) for q in range(2,D+1) for p in range(1,q) if gcd(p,q)==1]

    # ---- M7: only defined at lonely tau* (uses the gap M). Compare to M7@free-tau where
    #      we feed a free tau and the gap-value-at-that-tau as the "M".
    real_M7_lonely=set(); real_M7_free=set()
    for S in SS:
        Mg,tg=Mgap(S)
        if tg is not None and Mg>0:
            r=classify(S,M7_arc_factory(tg,Mg))
            if r: real_M7_lonely.add(r[0])
    for S in SS:
        for tau in free_taus:
            val=g(S,tau)
            if val==0: continue
            r=classify(S,M7_arc_factory(tau,val))
            if r: real_M7_free.add(r[0])
    f7l=free_all-real_M7_lonely; f7f=free_all-real_M7_free
    print(f"\nM7 binding-return  n={n}:")
    print(f"  realized@lonely {len(real_M7_lonely)}/{len(free_all)} (forb {len(f7l)}); "
          f"realized@free {len(real_M7_free)}/{len(free_all)} (forb {len(f7f)})")
    print(f"  H@lonely {sorted(set(FREE[n][c][1] for c in real_M7_lonely))}")
    lever7=f7l-f7f
    print(f"  LONELINESS-LEVER classes (forb@lonely minus forb@free): {len(lever7)}")
    for c in sorted(lever7):
        n3,H,sc=FREE[n][c]; print(f"     H={H} c3={n3} score={sc}")

    # ---- M8: hit-zero-return ----
    real_M8_lonely=set(); real_M8_free=set()
    for S in SS:
        _,tg=Mgap(S)
        if tg is not None:
            r=classify(S,M8_arc_factory(tg))
            if r: real_M8_lonely.add(r[0])
        for tau in free_taus:
            r=classify(S,M8_arc_factory(tau))
            if r: real_M8_free.add(r[0])
    f8l=free_all-real_M8_lonely; f8f=free_all-real_M8_free
    print(f"M8 hit-zero-return  n={n}:")
    print(f"  realized@lonely {len(real_M8_lonely)}/{len(free_all)} (forb {len(f8l)}); "
          f"realized@free {len(real_M8_free)}/{len(free_all)} (forb {len(f8f)})")
    print(f"  H@lonely {sorted(set(FREE[n][c][1] for c in real_M8_lonely))}")
    lever8=f8l-f8f
    print(f"  nontrivial@lonely: {any(FREE[n][c][1]>1 for c in real_M8_lonely)}")
    print(f"  LONELINESS-LEVER classes (forb@lonely minus forb@free): {len(lever8)}")
    for c in sorted(lever8):
        n3,H,sc=FREE[n][c]; print(f"     H={H} c3={n3} score={sc}")

print("\nDONE")
