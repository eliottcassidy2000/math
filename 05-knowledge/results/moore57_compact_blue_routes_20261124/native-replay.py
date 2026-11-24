#!/usr/bin/env python3
"""Native metagraph witnesses and identifiers (not one-tile flip paths)."""
from compact import field_host,forward,arcs

def positions(n):
    return [(i,j) for i in range(n) for j in range(i+2,n)]

def tiling(T,p):
    assert all(forward(T,u,v) for u,v in zip(p,p[1:]))
    return sum(1<<b for b,(i,j) in enumerate(positions(len(T))) if forward(T,p[j],p[i]))

def tournament(n,mask):
    k={ij:b for b,ij in enumerate(positions(n))}
    return arcs(n,lambda i,j: j==i+1 or not mask>>(k[i,j])&1)

def reflected(n,mask):
    k={ij:b for b,ij in enumerate(positions(n))}
    return sum(1<<k[n-1-j,n-1-i] for b,(i,j) in enumerate(positions(n)) if mask>>b&1)

def words(T):
    N=len(T)
    W=(1<<N)-1
    def rec(p,seen):
        if seen==W:
            yield p;return
        a=T[p[-1]]&~seen
        while a:
            b=a&-a;a^=b
            yield from rec(p+[b.bit_length()-1],seen|b)
    for i in range(N):
        yield from rec([i],1<<i)

def receipt(T):
    h=0;least=None
    for p in words(T):
        z=tiling(T,p);h+=1
        if least is None or z<least:least=z
    return least,h

def wire(n,orders,blue=False):
    T=tournament(n,0);W=(1<<len(positions(n)))-1
    M=[]
    for o in orders:
        t=tiling(T,list(map(int,o)))
        if blue:assert reflected(n,t)==t
        M.append(t)
        T=tournament(n,t^W)
    return T,M

def score_debt(T):
    return sum(abs(x-i) for i,x in enumerate(sorted(u.bit_count() for u in T)))//2

def degree_debt(scores):
    return sum(abs(i-x) for i,x in enumerate(sorted(scores)))//2

penta=arcs(5,lambda i,j:(j-i)%5 in (1,2))
T,M=wire(5,["01234","41230","14203"],True)
for W in (T,penta):
    assert (receipt(W),score_debt(W))==((22,15),3)
print("regular C5 class, height, shortest blue",receipt(T),3,"source masks",M)

old=[254,380,696,848,896,24,52,106,229,483]
ord0=["0123456789","8975634201","9512367840","2839450617","7561908342","9231458670"]
T,M=wire(10,ord0,True)
x=tiling(T,list(map(int,"4196273085")))
assert x==847909912
assert receipt(T)==receipt(old)
assert score_debt(old)==6
print("old C5 observer class, height, shortest blue",receipt(T),6,"source masks",M)

new,p=field_host(1)
end=["0123456789","8956723401","5967812304","5012367894","9123456780","7308451962","5078361294"]
T,M=wire(10,end,True)
assert tiling(T,list(map(int,"2031458697")))==tiling(new,p)
assert receipt(T)==receipt(new) and score_debt(new)==7
print("compact Petersen class, height, shortest blue",receipt(new),7,"source masks",M)
for d in (2,3,7,57):
    m=d*d+1
    assert degree_debt([m+d-i for i in range(m)]+[m-d-1+i for i in range(m)])==(d+1)*(m-d-1)
for q in (1,5):
    C,p=field_host(q)
    m=tiling(C,p)
    assert reflected(len(C),m)==m and m.bit_count()==6*q+5*q*q
    assert score_debt(C)=={1:7,5:106}[q]

# Verify the score slice and the equal-height trap at n=5 and 6.
for n, data in [(5,(5,4,9)),(6,(27,8,17))]:
    keys=range(1,n-1)
    records=set()
    for z in range(1<<(n-2)):
        K=[0]+[v for v in keys if z>>(v-1)&1]+[n-1]
        up=set(zip(K,K[1:]))
        T=arcs(n,lambda i,j: (i,j) in up)
        records.add(receipt(T)[0])
    assert len(records)==1<<(n-4)
    W=(1<<len(positions(n)))-1
    full=tournament(n,W)
    long=arcs(n,lambda i,j:(i,j)==(0,n-1))
    assert receipt(full)==(data[0],data[2]) and receipt(long)==(data[1],data[2])
    two={receipt(tournament(n,tiling(full,p)^W))[0] for p in words(full)}
    assert data[1] in two and data[0]!=data[1]
    print("equal-score equal-height trap",n,data)
