from functools import lru_cache
from itertools import product

def count(a,b,c):
    sizes=(a,b,c)
    labels=sum(([i]*n for i,n in enumerate(sizes)),[])
    N=len(labels)
    def arc(i,j):
        x,y=labels[i],labels[j]
        return i<j if x==y else (y-x)%3==1
    @lru_cache(None)
    def dp(mask,v):
        if mask==(1<<N)-1:return 1
        return sum(dp(mask | (1<<w),w) for w in range(N) if not (mask>>w)&1 and arc(v,w))
    return sum(dp(1<<v,v) for v in range(N))
for a,b,c in product(range(1,5), repeat=3):
    H=count(a,b,c); t=sum(x%2==0 for x in (a,b,c)); expected=(2*t+t*(t-1)//2)%3
    if H%3 != expected or H%3==1:raise RuntimeError((a,b,c,H,t,expected))
print('PASS: all 64 labelled cyclic block triples (1..4)^3, exact subset DP')
