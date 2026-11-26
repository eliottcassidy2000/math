"""Replay the selected physical board cut, with exact polynomial checks."""
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
from subprocess import check_output
from sys import executable
from random import Random

here=Path(__file__).resolve().parent
txt=check_output([executable,str(here/'nested_core_check.py'),'--board'],text=True)
it=iter(map(int,txt.split()))
n=next(it)
# The board generator emits the branches in SLOTS order except the source.
from three_stars import SLOTS
branches=[h for h in SLOTS if h!=1]
assert len(branches)==n and branches[0]==2
options,ref={},{}
for h in branches:
    old,count=next(it),next(it)
    ref[h]=old
    options[h]={next(it):tuple(56*c+next(it) for c in range(3)) for _ in range(count)}
assert next(it,None) is None
marker=options.pop(2)[7]; ref.pop(2)
for h,row in options.items():
    options[h]={y:p for y,p in row.items() if set(p).isdisjoint(marker)}
assert len(options)==52 and sum(map(len,options.values()))==1741
assert sorted(h for h,v in ref.items() if v==-1)==[13,33]
ref.update({13:16,33:19})

witness=[tuple(map(int,l.split())) for l in (here/'marked_cut_switches.tsv').read_text().splitlines()]
left,right=[(h,a,b) for h,a,b,*_ in witness],[tuple(w[3:]) for w in witness]
N=len(left)
assert N==17 and len({h for h,a,b in left+right})==34
def direction(entry):
    h,a,b=entry; assert ref[h]==a and a!=b
    z=Counter(options[h][b]);z.subtract(options[h][a])
    return {i:j for i,j in z.items() if j}
def dot(u,v):
    return sum(a*v.get(i,0) for i,a in u.items())
U,V=list(map(direction,left)),list(map(direction,right))
assert all(dot(u,v)==(i==j) for i,u in enumerate(U) for j,v in enumerate(V))

def score(u,v):
    new=ref.copy()
    for bits,side in [(u,left),(v,right)]:
        for i,(h,a,b) in enumerate(side):
            if bits>>i&1:new[h]=b
    c=Counter(j for h,v in new.items() for j in options[h][v])
    return sum(x*(x-1)//2 for x in c.values())
base=score(0,0)
assert base==3
rng=Random(57)
for _ in range(200):
    u,v=(rng.randrange(1<<N) for _ in range(2))
    assert score(u,v)-score(u,0)-score(0,v)+base==(u&v).bit_count()
print("complete selections 52; options 1741; 289 mixed coefficients I_17; generic rank",1<<N)

# Row and column exponent profiles from their exact quadratic expansions.
def side_scores(D):
    origin=Counter(j for h,v in ref.items() for j in options[h][v])
    # Change of sum binom(count,2) under one signed vector with sum zero.
    linear=[dot(origin,d)+dot(d,d)//2 for d in D]
    vals=[0]*(1<<N)
    for m in range(1,1<<N):
        bit=m&-m;i=bit.bit_length()-1;k=m^bit
        v=linear[i]; j=k
        while j:
            t=j&-j;j^=t;v+=dot(D[i],D[t.bit_length()-1])
        vals[m]=v+vals[k]
    return vals
a,b=side_scores(U),side_scores(V)
for _ in range(200):
    u,v=(rng.randrange(1<<N) for _ in range(2))
    assert score(u,v)==base+a[u]+b[v]+(u&v).bit_count()
groups=defaultdict(list)
full=(1<<N)-1
for u in range(1<<N):
    groups[u.bit_count(),a[u],b[full^u]].append(u)
big=max(groups,key=lambda k:len(groups[k]))
k0=base+big[1]+big[2]
print("coefficient minor",big,"squares",k0,"edits",N,"rank floor",len(groups[big]))
for top in (2,3,4,6,8,10,12,16,24,32):
    picks=[(len(vals),key) for key,vals in groups.items() if base+key[1]+key[2]<=top]
    print("coefficient degree at most",top,"best",max(picks,default=None))

# Exhaustion of only this binary box, not the 1,741-contact search.
fixed={h for h in ref}- {h for h,a,b in left+right}
load=Counter(j for h in fixed for j in options[h][ref[h]])
m0=sum(x*(x-1)//2 for x in load.values())
stack=[(options[h][a],options[h][b]) for h,a,b in left+right]
def dfs(k,c,cap,mask=0):
    if not k:
        if c<=cap: hits.append((mask & full,mask>>N,c))
        return 1
    best=min(range(len(k)),key=lambda i:(sum(sum(load[j] for j in o)+c<=cap for o in k[i][1]),
                                       -sum(load[j] for o in k[i][1] for j in o)))
    ans=0
    jbit, choices = k[best]
    for flip,o in enumerate(choices):
        cost=sum(load[j] for j in o)
        if c+cost>cap:continue
        for j in o:load[j]+=1
        ans+=dfs(k[:best]+k[best+1:],c+cost,cap,mask|(flip<<jbit))
        for j in o:load[j]-=1
    return ans
hits=[]
assert dfs(list(enumerate(stack)),m0,12)==7667
print("binary box cut by square cap: cap, configurations, GF(2) rank")
for cap in range(13):
    rows=defaultdict(int)
    subset=[(u,v) for u,v,c in hits if c<=cap]
    col={v:i for i,v in enumerate(sorted({v for _,v in subset}))}
    for u,v in subset:rows[u]|=1<<col[v]
    piv={}
    for row in rows.values():
        while row:
            k=row.bit_length()
            if k not in piv:
                piv[k]=row;break
            row^=piv[k]
    print(cap,len(subset),len(piv))
