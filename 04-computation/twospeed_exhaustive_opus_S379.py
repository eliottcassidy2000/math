from fractions import Fraction as F
from itertools import combinations
LAM=F(1,14)
def teeth01(x):
    w=LAM/x; out=[]
    for j in range(0,x+1):
        a,b=max(F(j,x)-w,F(0)), min(F(j,x)+w,F(1))
        if a<b: out.append((a,b))
    return out
def union(ivs):
    ivs=sorted(ivs); out=[]
    for a,b in ivs:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def complement(u):
    out=[]; prev=F(0)
    for a,b in u:
        if a>prev: out.append((prev,a))
        prev=b
    if prev<1: out.append((prev,F(1)))
    return out
def circ_lmax(c):
    if not c: return F(0)
    if len(c)>1 and c[0][0]==0 and c[-1][1]==1:
        wrap=(c[0][1]-c[0][0])+(c[-1][1]-c[-1][0])
        inner=max((b-a) for a,b in c[1:-1]) if len(c)>2 else F(0)
        return max(wrap,inner)
    return max(b-a for a,b in c)
def ess(V,drop):
    allv=[]
    for x in V:
        if x not in drop: allv.extend(teeth01(x))
    return complement(union(allv))
def diff_arcs(E,r):
    D=union(teeth01(r)); out=[]
    for a,b in E:
        cur=[(a,b)]
        for c,d in D:
            nxt=[]
            for x,y in cur:
                if d<=x or c>=y: nxt.append((x,y)); continue
                if x<c: nxt.append((x,min(c,y)))
                if y>d: nxt.append((max(d,x),y))
            cur=nxt
        out.extend([(x,y) for x,y in cur if y>x])
    return out
def contained(E,D):
    Du=union(D)
    return all(any(c<=a and b<=d for c,d in Du) for a,b in E)
def uncovered(V):
    allv=[]
    for x in V: allv.extend(teeth01(x))
    return 1-sum(b-a for a,b in union(allv))

base=list(range(1,14))
print("(3) STEP 3 -- bound the SECOND speed, then run the EXHAUSTIVE check")
tot_s=0; maxs=0; tight=[]; combos=0
for i,j in combinations(base,2):
    E=ess(base,{i,j}); L=circ_lmax(E)
    rb=int(F(2,5)/L)
    for r in range(1,rb+1):
        Ep=diff_arcs(E,r)
        if not Ep: continue
        Lp=circ_lmax(Ep)
        sb=int((2*LAM)/Lp)
        maxs=max(maxs,sb); tot_s+=max(0,sb-r+1)
        for s in range(r,sb+1):
            combos+=1
            V=sorted(set([x for x in base if x not in (i,j)]+[r,s]))
            if len(V)!=13: continue
            if contained(E, union(teeth01(r)+teeth01(s))):
                tight.append((i,j,r,s,tuple(V)))
print(f"    max bound on second speed: {maxs}")
print(f"    total (i,j,r,s) combinations enumerated: {combos}")
print(f"    TIGHT double substitutions found: {len(tight)}")
seen=set()
for i,j,r,s,V in tight:
    if V in seen: continue
    seen.add(V)
    u=uncovered(list(V))
    tag = "= {1..13}" if V==tuple(base) else ("= {1..11,13,24}" if V==(1,2,3,4,5,6,7,8,9,10,11,13,24) else "*** NEW ***")
    print(f"      ({i},{j}) -> ({r},{s}):  {list(V)}   uncovered={u}  [{tag}]")
print(f"    distinct tight families: {len(seen)}")
