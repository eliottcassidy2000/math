from functools import lru_cache
from itertools import combinations

def req(t, *a):
    if not t: raise AssertionError(a)

T4=(14,12,8,0)
SC=(14,4,8,2)
CS=(10,12,9,0)
S4=(6,12,8,1)
R5=(6,12,24,17,3)

def converse(R):
    return tuple(sum(1<<j for j,b in enumerate(R) if b>>i & 1) 
                 for i in range(len(R)))

def attach(R,u):
    return tuple(r if i==u else r+(1<<len(R)) for i,r in enumerate(R))+(1<<u,)

def profile(R):
    n=len(R)
    full=(1<<n)-1
    @lru_cache(None)
    def path(s,u):
        r=s^(1<<u)
        return sum(path(r,v) for v in range(n) if r & R[u] & (1<<v)) if r else 1
    @lru_cache(None)
    def h(s):
        return sum(path(s,u) for u in range(n) if s>>u&1) if s else 1
    E=[0]*n; O=[0]*n
    for s in range(1,full+1):
        a=O if s.bit_count()%2 else E
        for u in range(n):
            if s>>u&1: a[u]+=path(s,u)*h(full^s)
    return h(full), E, O, h

def reindex_deleted(R,x):
    I=[j for j in range(len(R)) if j!=x]
    return tuple(sum(((R[i]>>j)&1)<<k for k,j in enumerate(I)) for i in I),I

for R in [T4,SC,CS,S4,R5]:
 h,E,O,_=profile(R)
 print('anchor',R,h,E,O)

# Complete by degree-one deletion and the unique regular five.
F=[attach(T4,0),attach(SC,0),attach(CS,0)]+[attach(S4,u) for u in range(4)]+[R5]
pal={}
S_min=100
for R in F:
 h,E,O,_=profile(R)
 h0,E0,O0,_=profile(converse(R))
 req(h==h0)
 S_min=min(S_min,*map(sum,zip(E,O)),*map(sum,zip(E0,O0)))
 pal.setdefault(h,set()).update(E+E0)
 print('rooted attachment',R,h,'E',E,'E converse',E0)
req(pal=={9:{3,4,5,7,10,11,12},11:{4,6,7,9,13},13:{5,6,8,12},15:{6,8,9}})
print('W5 strong h -> E support',[(h,sorted(v))for h,v in sorted(pal.items())])
req(S_min==9)
print('minimum weighted S over rooted strong five:',S_min)
for h,E in pal.items(): req((21-h)//2 not in E)

R8=(4,5,152,3,11,31,63,123)
W,I=reindex_deleted(R8,0)
h,E,O,hsub=profile(W)
print('R8 delete vertex 0',W,h,E,O,'at root',I.index(2))
print('R8 full h',profile(R8)[0])
req(h+2*E[I.index(2)]==63)

def make_wedge(a,b):
    n=a+b+1
    # 0..a-1 high row, then a..a+b-1 low row; root r=n-1.
    R=[sum(1<<v for v in range(i+1,a+b)) for i in range(a+b)]
    R += [1|(sum(1<<v for v in range(a,a+b-1)))]
    for i in range(a+b):
        if not R[-1]>>i&1: R[i]|=1<<(n-1)
    return tuple(R)

tests=0
for a in range(2,7):
    for b in range(2,7):
        R=make_wedge(a,b)
        h,E,O,_=profile(R)
        X,Y=1<<(a-1),1<<(b-1)
        req(h==(1+X)*(1+Y),a,b,h)
        req(E[-1]==X*Y+Y-1,a,b,E[-1])
        req(profile(attach(R,len(R)-1))[0]==(1+X)*(1+3*Y)-2,a,b)
        tests+=1
print('two-wing emptying family through orders 14:', tests)
