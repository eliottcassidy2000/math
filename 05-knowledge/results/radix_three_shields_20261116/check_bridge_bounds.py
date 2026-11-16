from check_compensation import calc
from itertools import combinations
from random import Random
def adj(n,code):
    r=[0]*n
    for k,(i,j) in enumerate(combinations(range(n),2)):
        if code>>k&1:i,j=j,i
        r[i]|=1<<j
    return tuple(r)
def rev(a):
    return tuple(((1<<len(a))-1)^x^(1<<i) for i,x in enumerate(a))
def restrict(a,S):
    return tuple(sum((a[i]>>j&1)<<k for k,j in enumerate(S)) for i in S)
def p2(a):
    n=len(a);F=(1<<n)-1
    hs=[]
    for m in range(F+1):
        u=[i for i in range(n) if m>>i&1]
        hs.append(calc(restrict(a,u)) if u else 1)
    return sum(hs[i]*hs[F^i] for i in range(F+1))
def join(a,x,b,y):
    n=len(a);m=len(b);T=((1<<m)-1)<<n
    c=[r|(T^(1<<(n+y) if i==x else 0)) for i,r in enumerate(a)]
    c.extend((r<<n)|((1<<x) if j==y else 0) for j,r in enumerate(b))
    return tuple(c)
seen=set()
for n in range(1,6):
    for w in range(1<<(n*(n-1)//2)):
        a=adj(n,w)
        for x in range(n):
            h,g,e,o=calc(a,x)
            if e==1:seen.add((h,e,o))
assert seen=={(1,1,1),(3,1,2)}
rng=Random(1768)
for j in range(120):
    n,m=rng.randrange(1,6),rng.randrange(1,6)
    a,b=adj(n,rng.randrange(1<<(n*(n-1)//2))),adj(m,rng.randrange(1<<(m*(m-1)//2)))
    x,y=rng.randrange(n),rng.randrange(m)
    h,g,E,O=calc(a,x)
    k,p,Es,Os=calc(b,y)
    k,p,Ee,Oe=calc(rev(b),y)
    L=p2(b)//2
    assert calc(join(a,x,b,y),x)==(h*k+2*(E*Oe+O*Ee),g*k,L*(E+O)-g*Os,L*(E+O)-g*Es)
print('Factor-one signatures through order five and corrected retained-root bridge laws passed.')
tri=(2,4,1)
V0=tuple(x|24 for x in tri)+(16,0)
V1=(30,28,8,16,4)
for v,y,ee in [(V0,0,(1,2)),(V1,1,(3,3))]:
    assert calc(v,y)==(3,1 if v==V0 else 3,6,6)
    assert calc(rev(v),y)[2:]==ee and p2(v)==48
low=[join(tri,0,v,y) for v,y in [(V0,0),(V1,1)]]
assert [calc(a,0) for a in low]==[(17,3,66,66),(27,3,66,66)]
assert [calc(join(a,0,v,y)) for v,y in [(V0,0),(V1,1)] for a in low]==[447,477,843,873]
print('Independent exact ×3 deposit from rooted 17/27 and second-step four settings passed.')
