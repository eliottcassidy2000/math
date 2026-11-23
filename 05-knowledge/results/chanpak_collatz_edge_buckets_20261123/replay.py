from itertools import permutations,combinations
from random import Random

def edge(x):
    u=3*x+1
    j=(u&-u).bit_length()-1
    return u>>j,j

def solve_guard(word,m,r):
    D=1;k=0
    for j in word:k=3*k+D;D<<=j
    v=(D-k)*pow(3,-len(word),2*D)%(2*D)
    x=r+m*((v-r)*pow(m,-1,2*D)%(2*D))
    assert 0<x<m*2*D
    return x,m*2*D

def data(n,r,arcs):
    a=1
    while 2*3**(a-1)<n+2:a+=1
    m=3**a;d=2*3**(a-1)
    logs={pow(2,j,m):j for j in range(d)}
    assert len(logs)==d and len(r)==n+2 and len(set(r))==n+2
    α=[logs[(3*u+1)%m] for u in r];β=[logs[u] for u in r]
    j={(u,v):d+α[u]-β[v] for u in range(n+2) for v in range(n+2)}
    E=set(arcs)|{(0,i) for i in range(1,n+2)}|{(i,n+1) for i in range(1,n+1)}
    J=d*(n+1)+sum(α[:-1])-sum(β[1:])
    def accept(x):
        if x%m!=r[0] or not x&1:return None
        seen=[0]
        for _ in range(n+1):
            y,k=edge(x);q=y%m
            if q not in r:return None
            v=r.index(q)
            if v in seen or (seen[-1],v) not in E or k!=j[seen[-1],v]:return None
            seen.append(v);x=y
        return tuple(seen) if seen[-1]==n+1 else None
    S={}
    for P in permutations(range(1,n+1)):
        q=(0,)+P+(n+1,)
        if not all(e in E for e in zip(q,q[1:])):continue
        word=[j[e] for e in zip(q,q[1:])]
        assert sum(word)==J
        x,M=solve_guard(word,m,r[0])
        assert accept(x)==q and x not in S
        S[x]=q
    return m<<(J+1),S,accept

def orientation(n,dist):
    return {(u,v) for u in range(1,n+1) for v in range(1,n+1) if dist(u,v)}

Q=orientation(5,lambda u,v:(v-u)%5 in (1,2))
r=[1,2,4,8,7,5,10]
M,S,_=data(5,r,Q)
assert len(S)==15
e=(1,2);M1,U,_=data(5,r,(Q-{e})|{e[::-1]})
assert M==M1
assert S.keys()&U.keys()=={x for x,P in S.items() if e not in zip(P,P[1:])}

# Check *all* representatives for a small physical example, independently of enumeration.
M,S,f=data(1,[1,2,4],set())
assert set(S)=={x for x in range(1,M,2) if f(x)}

# Bounded verification of the two policy obstructions.
for m in range(3,30,2):
    O={}
    for a in range(m):
        ls=[];x=a
        while x not in ls:ls.append(x);x=2*x%m
        O[a]=min(ls)
    E={(u,v) for u in range(m) for v in range(m) if O[(3*u+1)%m]==O[v]}
    assert not any(all(((u,v) in E)^((v,u) in E) for u,v in combinations(S,2))
                   for S in combinations(range(m),4))
for P in permutations(range(1,7)):
    c=(0,)+P
    W={frozenset(c[(j+k)%7] for k in range(3)) for j in range(7)}
    assert any(all(frozenset((x+d)%7 for d in D) not in W for D in ((1,2,4),(3,5,6)))
               for x in range(7))
print('exact guards, count-preserving period, one-flip overlap and Paley window obstruction pass')

def twice_guard(u,v,m,L):
    for j in range(1,L):
        if pow(2,j,m)==(pow(2,L,m)*v-9*u-3)%m:
            return j,L-j
    raise AssertionError((u,v))

r=[x for x in range(1,10) if x%3][:5]
R={}
m=9; L=7
for P in permutations(range(1,6)):
    E=list(zip(P,P[1:]))
    if not all(e in Q for e in E): continue
    word=[j for u,v in E for j in twice_guard(r[u-1],r[v-1],m,L)]
    assert sum(word)==L*4
    x,period=solve_guard(word,m,r[P[0]-1])
    y=x; q=[y%m]
    for i,j in enumerate(word):
        y,k=edge(y)
        assert k==j
        if i%2: q.append(y%m)
    assert tuple(q)==tuple(r[u-1] for u in P)
    assert x not in R
    R[x]=P
assert len(R)==15 and period==m<<(1+L*4)

def hamp(rows):
    return sum(all(rows[u]>>v&1 for u,v in zip(p,p[1:]))
               for p in permutations(range(len(rows))))

def flip(rows,u,v):
    rows=list(rows)
    rows[u]^=1<<v
    rows[v]^=1<<u
    return tuple(rows)

def strong(rows):
    n=len(rows)
    for i in range(n):
        x=1<<i
        for _ in range(n):
            for j in range(n):
                if x>>j&1: x|=rows[j]
        if x!=(1<<n)-1: return False
    return True

A=(158,28,24,16,32,79,159,62)
B=flip(A,1,5); C=flip(A,5,7); D=flip(B,5,7)
assert list(map(hamp,(A,B,C,D)))==[61,63,65,69]
assert all(map(strong,(A,B,C,D)))
print('same-length two-step macros: 15 paths; strong literal fork: 61,63,65,69')

def orbit(n):
    out=[n]
    while n!=1:
        n,j=edge(n)
        out.append(n)
    return out

assert orbit(21)==[21,1]
assert orbit(7)==[7,11,17,13,5,1]
assert orbit(61)==[61,23,35,53,5,1]
assert orbit(63).index(61)==34
assert orbit(65)==[65,49,37,7,11,17,13,5,1]
print('numerical links: 21 -> 1, 61 -> 23, 63 -> 61 in 34, 65 -> 7 in 3')

from collections import Counter
from math import comb

def gap_counts(A,port):
    ext=[v for v in range(len(A)) if v!=port]
    counts=Counter()
    for p in permutations(ext):
        N=len(p)
        allowed=required=0
        for i in range(N+1):
            left = i==0 or (A[p[i-1]]>>port&1)
            right = i==N or (A[port]>>p[i]&1)
            if left and right:
                allowed|=1<<i
            if 0<i<N and not(A[p[i-1]]>>p[i]&1):
                required|=1<<i
        if required&~allowed:
            continue
        r=required.bit_count()
        s=(allowed^required).bit_count()
        for q in range(max(0,1-r),s+1):
            counts[r+q]+=comb(s,q)
    return dict(counts)

for A, c1, c2 in [
    ((2,8,3,5),5,2),
    ((8,17,3,6,13),13,6),
    ((4,41,2,5,15,29),15,6),
    ((8,41,3,4,15,29),17,6),
    ((2,12,9,1),5,3),
    ((2,36,57,51,35,1),15,9),
    ((2,36,41,51,39,1),17,9),
]:
    assert strong(A) and gap_counts(A,0)=={1:c1,2:c2}
assert flip((2,36,57,51,35,1),2,4)==(2,36,41,51,39,1)
for A,e,c in [
    ((50,188,25,17,32,76,31,125),(3,5),84),
    ((242,76,9,17,38,78,28,126),(3,6),132),
]:
    B=flip(A,*e)
    assert strong(A) and strong(B)
    assert gap_counts(A,7)=={1:2*c+3,2:c}
    assert gap_counts(B,7)=={1:2*c+5,2:c}
A=[8,1,3,102,15,23,55,127,127^48|(1<<7)]
for i in (4,5):
    A[i]|=1<<8
assert gap_counts(A,7)=={1:225,2:112}
assert hamp((80,9,3,5,14,31,62))==57
assert hamp((48,33,67,7,14,28,59))==113
print('all strong binary-width rail templates: exact gap counts and upper seams')
