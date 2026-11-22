from collections import deque
from itertools import product
from random import Random

def accepted(N,B,edges,seeds,finals,U,V):
    seeds={(state,U*x+V) for state,x in seeds}
    E=[(i,j,k,U*t+(1-B**k)*V) for i,j,k,t in edges]
    C=1+(max(abs(t) for _,_,_,t in E)+B-2)//(B-1)
    waiting=deque((s,N,0) for s in finals);seen=set(waiting)
    while waiting:
        s,x,c=waiting.popleft()
        if (s,x+c) in seeds:return True
        for i,j,k,t in E:
            if j!=s:continue
            M=B**k
            u,d=divmod(x,M)
            a,r=divmod(d+c-t,M)
            if r:continue
            assert abs(a)<=C
            v=(i,u,a)
            if v not in seen:
                seen.add(v);waiting.append(v)
    return False

def actual(N,B,edges,seeds,finals,U,V):
    wait=deque(seeds);seen=set(seeds)
    out=set()
    while wait:
        i,x=wait.popleft()
        if i in finals and U*x+V<=N:out.add(U*x+V)
        for j,k,s,t in edges:
            if j!=i:continue
            z=B**s*x+t
            assert z>x   # chosen tests grow strictly, permitting independent height truncation
            pair=k,z
            if U*z+V<=N and pair not in seen:
                seen.add(pair)
                wait.append(pair)
    return out

r=Random(9417)
for _ in range(160):
    B=r.randint(2,5);U=r.randint(1,3);V=r.randint(-7,7)
    E=[(r.randrange(4),r.randrange(4),r.randint(1,2),r.randint(-3,7)) for i in range(12)]
    S=[(0,5),(1,7)]
    F={2,3}
    values=actual(200,B,E,S,F,U,V)
    assert {n for n in range(201) if accepted(n,B,E,S,F,U,V)}==values

B=490032;slope=63618252914778017089536;K=47124685853114258307047
E=[(0,0,1,d) for d in (360532,362956)]
for depth in range(7):
    for digits in product((360532,362956),repeat=depth):
        x=0
        for d in digits:x=x*B+d
        for k in (K,K+2):
            y=slope*x+k
            assert accepted(y,B,E,[(0,0)],{0},slope,k)
            assert not accepted(y+1,B,E,[(0,0)],{0},slope,k)
            assert not accepted(y+slope,B,E,[(0,0)],{0},slope,k)
print("160 signed synthetic controller exhaustions plus all 254 short literal endpoint words")
