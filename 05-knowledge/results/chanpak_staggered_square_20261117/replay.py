from itertools import permutations
from base_replay import outer,hams
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
def arm(order):
    mask=0; path=[mask]
    for k in order:
        mask|=1<<k;path.append(mask)
    return path
def val(m,h):
    a,b=D[m];return a*47680+b*h
a,b=arm((0,1,2,3,4)),arm((0,3,2,1,4))
def qualify(x,y):
    X=[val(m,383) for m in x];Y=[val(n,389) for n in y]
    return all(i<j for i,j in zip(X,Y)) and all(i<j for A in (X,Y) for i,j in zip(A,A[1:]))
if not qualify(a,b) or tuple((x^y).bit_count()+1 for x,y in zip(a,b))!=(1,1,3,3,1,1):
    raise ValueError("staggered witness")
for i,(x,y) in enumerate(zip(a,b)):
    print(i,D[x],D[y],val(x,383),val(y,389),(x^y).bit_count()+1)
paths=[arm(p) for p in permutations(range(5))]
cost=[sum(x!=y for x,y in zip(a,b)) for a in paths for b in paths if qualify(a,b)]
if (len(cost),min(cost))!=(322,2):raise ValueError("full order-pair census")
print("322 increasing aligned order-pairs; at least two different-exterior rows")
