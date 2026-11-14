from collections import defaultdict
from itertools import product
from random import Random
from replay_bilateral import k41, rooted, value

def carry_graph(f,h,P,Q):
    A,B=P+h,Q+h
    start=(0,0)
    g=defaultdict(list)
    seen={start}
    states=[start]
    reverse=defaultdict(set)
    for u,v in states:
        for d,e in product(range(2),repeat=2):
            x=u+f*d+h*e-A
            y=2*v-f*e-h*d+B
            if x%2 or not -A<=x//2<=f+h-A or not -B<=y<=f+h-B: continue
            w=(x//2,y)
            g[u,v].append((w,d,e));reverse[w].add((u,v))
            if w not in seen:
                seen.add(w);states.append(w)
    keep={start}
    todo=[start]
    for w in todo:
        for v in reverse[w]:
            if v not in keep:
                keep.add(v);todo.append(v)
    return {s: [w for w in g[s] if w[0] in keep] for s in keep}

C=k41
g=carry_graph(C.f,C.h,C.P,C.Q)
layers=[
    {(0,0),(-30,30)},
    {(-21,-42),(39,-57),(-36,18),(24,3)},
    {(21,-21),(51,-51)},
    {(-18,36),(42,21),(57,-39),(-3,-24)},
]
assert set(g)==set.union(*layers)
for i,L in enumerate(layers):
    assert {w for s in L for w,d,e in g[s]}==layers[(i+1)%4]

def four_steps(start):
    pts=[(start,0,0)]
    for i in range(4):
        pts=[(w,a+(d<<i),2*b+e) for s,a,b in pts for w,d,e in g[s]]
    return sorted(pts)

S=[(0,0),(-30,30)]
T={
 (0,0):[(6,9),(7,8)],
 (0,1):[(2,9),(3,8)],
 (1,0):[(6,13),(7,12)],
 (1,1):[(2,13),(3,12)]
}
for i,s in enumerate(S):
    assert four_steps(s)==sorted((S[j],a,b) for j in range(2) for a,b in T[i,j])
    for j,t in enumerate(S):
        for a,b in T[i,j]:
            for x,y in product(range(1,9),repeat=2):
                diff=value(C,16*x+a,y)-value(C,x,16*y+b)
                assert diff==s[0]*(16*x+a)+s[1]*y-t[0]*x-t[1]*(16*y+b)+6*((x%2)+(y%2)-1)+12*(i-j)

rng=Random(1729)
for _ in range(300):
    r=rng.randrange(1,31)
    modes=[0]+[rng.randrange(2) for i in range(r-1)]+[0]
    a=b=0
    for j,(s,t) in enumerate(zip(modes,modes[1:])):
        c,d=rng.choice(T[s,t])
        a+=c<<(4*j)
        b=16*b+d
        if not j: first=c%2
        last=c%2
    x,y=[rng.randrange(1,1<<32) for i in range(2)]
    l,r=(x<<(4*r))+a,(y<<(4*r))+b
    assert rooted(C,l,y)==rooted(C,x,r)
    assert value(C,l,y)-value(C,x,r)==6*(x%2+y%2-1+last-first)
print('PASS: complete bounded carry graph, base-16 boundary law, 300 composed shuttles.')

def mates(x,y):
    state,b,w,first=0,0,1,None
    while x>=16:
        x,a=divmod(x,16)
        if a not in (2,3,6,7):return
        nxt=int(a<6)
        digit=next(b for t,b in T[state,nxt] if t==a)
        b=16*b+digit
        w*=16
        d=a%2
        if first is None:first=d
        if nxt==0 and x%2+y%2-1+d-first==0:
            yield x,w*y+b
        state=nxt

assert (1,669) in set(mates(354,2))
from replay_bilateral import count,graph
profiles=[]
for x,y in [(354,2),(1,669)]:
    R=graph(C,x,y); apex=len(R)-1; Z=(1<<len(R))-1
    H,S=count(R)
    U=sum(row[apex]*H[Z^i] for i,row in enumerate(S))
    E=sum(row[apex]*H[Z^i]*((-1)**i.bit_count()) for i,row in enumerate(S))
    assert (H[-1],U,E)==(value(C,x,y),rooted(C,x,y),-15)
    profiles.append((H[-1],U,E))
assert profiles[0]==profiles[1]
print('PASS: literal nonzero internal carry, H/R+/R- =', profiles[0])
