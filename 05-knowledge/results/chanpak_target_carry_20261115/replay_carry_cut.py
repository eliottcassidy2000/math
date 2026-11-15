from functools import lru_cache
from itertools import product
from random import Random

B,Q,T,S,h=2855,5760,4034,3914,45

def check(condition):
    if not condition:
        raise AssertionError("carry cut replay")

def half(x):
    ans=0
    while x:
        x,b=divmod(x,2)
        if not b:
            ans+=(x+1)//2
    return ans

def full(x,y):
    return B+Q*x*y+T*x+S*y+2*h*(half(x)+half(y))

def w(d): return 1<<max(0,d-2)
def bits(x,r): return tuple((x>>i)&1 for i in range(r-1,-1,-1))
def addr(bs):
    x=0
    for b in bs: x=2*x+b
    return x

def allsame(a): return addr(t>=1 for t in a),addr(t==2 for t in a)
def pair(a,b):
    xx,yy=[],[]
    source=iter(b)
    for d in a:
        u=next(source) if d==1 else d//2
        xx.append(u)
        yy.append(d-u)
    return addr(xx),addr(yy)

def cut(I,z):
    return sum(w(j-i) for jj,j in enumerate(I) for kk,i in enumerate(I[:jj]) if z[jj]!=z[kk])

def selector(a, ex, target):
    modulus=1<<ex
    keep=(1<<(ex+1))-1
    @lru_cache(None)
    def can(i,hist,need):
        if i==len(a):
            return () if need==0 else None
        for bit in ((0,1) if a[i]==1 else (0,)):
            inc=0 if a[i]!=1 else sum(w(d) for d in range(1,min(ex+1,i)+1)
                                      if a[i-d]==1 and bit!=((hist>>(d-1))&1))
            tail=can(i+1,((hist<<1)|bit)&keep,(need-inc)%modulus)
            if tail is not None:
                return ((bit,) if a[i]==1 else ())+tail
        return None
    return can(0,0,target%modulus)

total=0
for r in range(1,9):
    for x,y in product(range(1<<r),repeat=2):
        xx,yy=bits(x,r),bits(y,r)
        a=tuple(i+j for i,j in zip(xx,yy))
        I=tuple(i for i,d in enumerate(a) if d==1)
        z=tuple(xx[i] for i in I)
        u,v=allsame(a)
        check(x+y==u+v)
        check(half(x)+half(y)==half(u)+half(v)+cut(I,z))
        deg=[sum(1 for j in I if 0<abs(i-j)<=2) for i in I]
        runs=[]
        for i in I:
            if not runs or i-runs[-1][-1]>2:
                runs.append([])
            runs[-1].append(i)
        frozen=all(len(run)==1 or (len(run)==3 and run[2]-run[0]==2) for run in runs)
        check(frozen==all(d%2==0 for d in deg))
        check((full(x,y)-full(u,v))%4==2*(sum(d*b for d,b in zip(deg,z))%2))
        total+=1
print('full carry and parity checks',total)

for r in range(5,17):
    limit=(1<<r)-1
    hmin=min(half(v)+half(limit-v) for v in range(1,limit))
    check(hmin==(1<<(r//2)) if r%2 else hmin==3*(1<<(r//2-2)))
print('central minimum exact ranges',5,16)

rng=Random(145)
for _ in range(1200):
    a=rng.randrange(1,100)
    R=rng.getrandbits(a)
    u=rng.getrandbits(a)
    pos=2*a+3
    A=1<<pos
    L=sum(1<<(pos-j-2) for j in range(a) if (R>>j)&1)
    high=u<<(pos+1)
    f0=full(high+A+L,high)
    f1=full(high+L,high+A)
    check(f1-f0==Q*A*L+(S-T)*A+90*R)
    check((f1-f0-90*R)%(1<<(a+1))==0)
print('adaptive star differential checks',1200)

for _ in range(600):
    a=tuple(rng.randrange(3) for i in range(rng.randrange(1,14)))
    I=tuple(i for i,d in enumerate(a) if d==1)
    ex=rng.randrange(1,7)
    t=rng.randrange(1<<ex)
    target=next((z for z in product((0,1),repeat=len(I)) if cut(I,z)%(1<<ex)==t),None)
    got=selector(a,ex,t)
    check(got==target)
    if got is not None:
        u,v=allsame(a)
        x,y=pair(a,got)
        check((half(x)+half(y)-half(u)-half(v)-t)%(1<<ex)==0)
print('canonical local selectors',600)
