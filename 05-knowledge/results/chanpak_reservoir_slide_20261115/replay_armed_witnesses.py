"""Literal subset DP for both directions at the count-5 and count-59 seeds."""
from functools import cache

def require(ok):
    if not ok:
        raise AssertionError("two-ended witness check failed")

@cache
def profile(a):
    n=len(a); full=(1<<n)-1
    table=[[0]*n for _ in range(full+1)]
    H=[1]+[0]*full
    for s in range(1,full+1):
        for i in range(n):
            if not (s>>i&1): continue
            t=s^(1<<i)
            table[s][i] = sum(table[t][j] for j in range(n) if t&a[i]&(1<<j)) if t else 1
        H[s]=sum(table[s])
    D=tuple(2*sum(table[s][i]*H[full^s] for s in range(1,full+1)
                  if s.bit_count()%2==0)-H[full] for i in range(n))
    incoming_clear=tuple(all(not table[full^(1<<i)][j]
                             for j in range(n) if i!=j and not (a[i]>>j&1))
                         for i in range(n))
    return H[-1],D,incoming_clear,all(table[-1])

def converse(a):
    full=(1<<len(a))-1
    return tuple(full^x^(1<<i) for i,x in enumerate(a))

def menus(a):
    h,p,s,strong=profile(a)
    h2,q,t,strong2=profile(converse(a))
    require(h==h2 and strong==strong2)
    return h,strong,tuple((1,i,X[i]) for X,E in ((p,t),) for i in range(len(a)) if E[i] and abs(X[i])==1)+tuple((0,i,X[i]) for X,E in ((q,s),) for i in range(len(a)) if E[i] and abs(X[i])==1)

def forward(a,q):
    return tuple(r|((i!=q)<<len(a)) for i,r in enumerate(a))+(1<<q,)

def reverse(a,q):
    return converse(forward(converse(a),q))

def history(seed,seq,expected):
    a=(0,); labels=[seed]
    for j,q in seq:
        a=forward(a,labels.index(q));labels.append(j)
    restored=[0]*len(labels)
    for i,row in enumerate(a):
        restored[labels[i]]=sum(1<<k for j,k in enumerate(labels) if row>>j&1)
    require(tuple(restored)==expected)

A=(12,1,2,6)
history(2,((0,2),(3,2),(1,0)),A)
require(menus(A)==(5,True,((1,0,1),(1,1,-1),(0,0,-1),(0,1,1))))
require(menus(forward(A,0))==(11,True,((1,0,1),(0,0,1))))
require(menus(reverse(A,0))==(9,True,((1,1,-1),(0,0,-1))))
print("H5: same-site stable F+1 and R−1; both branch totals and child menus verified.")

B=(114,16,3,7,12,30,62)
history(4,((0,4),(6,4),(5,4),(3,0),(2,0),(1,4)),B)
require(menus(B)==(59,True,((1,4,-1),(0,0,1))))
for a in range(3):
    Y=B
    for _ in range(a): Y=forward(Y,4)
    grid=[]
    for b in range(3):
        H,D,Er,_=profile(Y)
        H2,G,Ef,_=profile(converse(Y))
        require(H==H2 and Ef[4] and Er[0])
        x,y=(1<<a)-1,(1<<b)-1
        require((H,D[4],G[0])==(59+58*x+60*y+56*x*y,-1-4*y,1-2*x))
        grid.append(H);Y=reverse(Y,0)
    print("H59 row",a,":",*grid)
require(menus(forward(B,4))==(117,True,((1,4,-1),(0,0,-1))))
require(menus(reverse(B,0))==(119,True,((0,0,1),)))
print("H59: exact two-ended grid, root readiness, and first child menus verified.")

def join(a,b):
    n=len(a); mask=((1<<len(b))-1)<<n
    return tuple(x|mask for x in a)+tuple(x<<n for x in b)

def check_raw_components(parts):
    T=()
    factor=1
    offsets=[]
    for P in parts:
        offsets.append(len(T))
        factor*=profile(P)[0]
        T=join(T,P)
    H,D,_,_=profile(T)
    require(H==factor and H>1)
    last=parts[-1]; d=profile(last)[1]
    require(all(x>=H for x in D[:offsets[-1]]))
    require(D[offsets[-1]:]==tuple((H//profile(last)[0])*q for q in d))
    R=converse(T)
    _,E,_,_=profile(R)
    first=parts[0]; e=profile(converse(first))[1]
    require(all(x>=H for x in E[len(first):]))
    require(E[:len(first)]==tuple((H//profile(first)[0])*q for q in e))
    require(not any(abs(x)==1 for x in D) or not any(abs(x)==1 for x in E))

from itertools import product
for k in (2,3):
    for parts in product(((0,),(2,4,1),A),repeat=k):
        if any(len(P)>1 for P in parts):
            check_raw_components(parts)
print("32 reducible ordinal sums verify the raw-unit confinement and component factor laws.")
