from itertools import permutations
from replay import outer,hams
old_order=(4,5,7,6,0,1,2,3)
edges=((1,7),(1,8),(3,4),(6,8),(4,6))
def profile(mask):
    h=[]
    for m in (1,2):
        # outer(TT_m,7), mapping the exterior indices by the shared cycle.
        tt=tuple((1<<m)-(2<<i) for i in range(m))
        a=list(outer(tt,7))
        labels={x:i for i,x in enumerate(old_order)}
        labels[8]=8+m
        for k,(x,y) in enumerate(edges):
            if mask>>k&1:
                i,j=labels[x],labels[y];a[i]^=1<<j;a[j]^=1<<i
        h.append(hams(a)[-1])
    t=(h[1]-h[0])//2
    return t,h[0]-2*t
D=[profile(i) for i in range(32)]
chains=[]
for seq in permutations(range(5)):
    sums=[0]
    for k in seq:sums.append(sums[-1]|1<<k)
    v=[D[i] for i in sums]
    if all(b>0 for a,b in v) and all(a1>=a and (a1-a)*2+b1-b>0 for (a,b),(a1,b1) in zip(v,v[1:])):
        chains.append((seq,v))
print("positive coefficient routes",len(chains))
for x in chains[:5]:print(x)
pos=mono=both=0; best=None
for seq in permutations(range(5)):
    s=0;v=[D[s]]
    for k in seq:s|=1<<k;v.append(D[s])
    positive=all(b>0 for a,b in v)
    monotone=all(a1>=a and (a1-a)*2+b1-b>0 for (a,b),(a1,b1) in zip(v,v[1:]))
    concrete=all((a1-a)*8*5960+(b1-b)*h>0 for (a,b),(a1,b1) in zip(v,v[1:]) for h in (383,389))
    pos+=positive;mono+=monotone;both+=positive and concrete
    if monotone:best=max(best if best is not None else -1000,min(b for a,b in v))
if (pos,mono,both,best)!=(36,29,0,-3):raise ValueError((pos,mono,both,best))
print("120 orders: 36 keep horizontal sign, 29 monotone for every source, 0 do both even for the prime square; best minimum coefficient -3")
