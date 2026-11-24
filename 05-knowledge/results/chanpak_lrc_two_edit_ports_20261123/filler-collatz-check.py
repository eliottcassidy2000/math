"""Exact finite checks and a modular constructor for 27's Collatz itinerary."""
from itertools import product
from math import gcd, prod
from pathlib import Path
from runpy import run_path

bank = run_path(str(Path(__file__).with_name("third-slot-bank.py")))["data"]
entries = [(d, q, m[0]) for d, q, u, v in bank for m in (u, v)]
assert all(d % 2 == r % 2 == 1 and q % 4 == 0 and d * r % 8 in (1, 3)
           for d, q, r in entries)
m8 = {(1, 1)}
a, b, c = [(bank[i][0], bank[i][2][0]) for i in (0, 1, 2)]
while True:
    enlarged = m8 | {(d*x % 8, r*y % 8) for d, r in m8 for x, y in (a,b,c)}
    if enlarged == m8:
        break
    m8 = enlarged
assert m8 == {(d, r) for d in (1,3,5,7) for r in (1,3,5,7) if d*r%8 in (1,3)}
assert [((x*x-1)//8%2, (y*y-1)//8%2) for x,y in (a,c)] == [(0,1),(1,0)]

types = (0, 2, 6, 11, 13)
rows = [bank[i] for i in types]
move = (0, 0, 2, -1, 0)
def ratios(M):
    return tuple(v[0] * pow(u[0], -1, M) % M for _, _, u, v in rows)
def factor(exps, M):
    return prod(pow(r, j, M) for r,j in zip(ratios(M),exps)) % M
assert factor(move, 77) == 1 and factor(move, 32) == 17
anchors = {}
for a0 in product(range(5),repeat=5):
    anchors.setdefault(factor(a0,77),a0)
assert set(anchors) == {s for s in range(77) if gcd(s,77) == 1}
d8,q8,u8,_ = bank[8]
r8 = u8[0]
assert (11*pow(d8,13,64)+20*pow(r8,13,64)) % 64 == 27

def readings(T,x,M):
    D = pow(d8,13,M)*prod(pow(d,T,M) for d,_,_,_ in rows)%M
    Q = pow(q8,13,M)*prod(pow(q,T,M) for _,q,_,_ in rows)%M
    R = pow(r8,13,M)*prod(pow(u[0],T-j,M)*pow(v[0],j,M)
                           for (_,_,u,v),j in zip(rows,x))%M
    return D,Q,R

def constructor(K,u,s):
    """Return the common quota and five second-alternative copy counts."""
    assert K>=6 and u%64==27 and gcd(s,77)==1
    M=2**K
    N=2**(K-6)
    T=16*(N+1)
    base=(0,0,0,N,0)
    phase=readings(T,base,77)[2]
    a0=anchors[s*pow(phase,-1,77)%77]
    y=tuple(i+j for i,j in zip(base,a0))
    D,Q,R=readings(T,y,M)
    H=(11*D+42*Q+20*R)%M
    assert H%64==27
    target=(1+(u-H)//4*pow(5*R,-1,M//4))%(M//4)
    g=factor(move,M//4)
    k=0
    for j in range(K-6):
        B=2**(5+j)
        if pow(g,k,B)!=target%B:
            k+=2**j
        assert pow(g,k,B)==target%B
    x=tuple(i+k*j for i,j in zip(y,move))
    assert all(0<=i<=T for i in x)
    d,q,r=readings(T,x,M)
    assert (d,q)==(D,Q) and (11*d+42*q+20*r)%M==u%M
    assert readings(T,x,77)[2]==s%77
    return T,x

for K in range(6,12):
    for u in range(27,2**K,64):
        for s in (1,2,76):
            constructor(K,u,s)
for s in anchors:
    constructor(70,27,s)

def trace(n):
    a0=[n]
    while n!=1:
        n=(3*n+1)//2 if n%2 else n//2
        a0.append(n)
    return a0

def ordinary_steps(a0):
    return sum(1+n%2 for n in a0[:-1])

assert max(ordinary_steps(trace(n)) for n in range(1,27))==23
traj=trace(27)
assert (len(traj)-1,ordinary_steps(traj),max(traj))==(70,111,4616)
assert max(2*y if x%2 else y for x,y in zip(traj,traj[1:]))==9232
odd=0
for j,n in enumerate(traj[1:],1):
    odd+=traj[j-1]%2
    if j<59:
        assert n>=27 and 3**odd>=2**j
    if j==59:
        assert (n,odd)==(23,37) and 3**odd<2**j
assert odd==41
def ordinary_trace(n):
    out=[n]
    while n!=1:
        n=3*n+1 if n%2 else n//2
        out.append(n)
    return out

o27,o189=ordinary_trace(27),ordinary_trace(189)
assert len(o27)==112 and len(o189)==107 and sum(n%2 for n in o189[:-1])==38
assert o27[8]==o189[3]==142 and o27[8:]==o189[3:]
assert set(o27[:8]).isdisjoint(o189)
assert set(o189[:3]).isdisjoint(o27)

def v3(n):
    k=0
    while n%3==0:
        k+=1
        n//=3
    return k
bad=(4,9,10,12)
assert all(r%3==0 for _,_,r in entries)
assert [(v3(bank[i][0]),v3(bank[i][1])) for i in bad]==[(2,1),(2,3),(1,1),(1,1)]
assert (bank[4][0]//9)%3==(bank[4][1]//3)%3==1
for i,(d,q,u,v) in enumerate(bank):
    if i not in bad:
        assert d%3!=0 and q%3==-d%3 and d*u[0]%8==d*v[0]%8==3

indices=(0,1,2,3,5,6,7,8,11,13)
counts=(72,75,53,47,4,80,11,95,41,51)
M=2**70
D,Q,R=(pow(v,13,M) for v in (d8,q8,r8))
for i,j in zip(indices,counts):
    assert 0<=j<=96
    d,q,u,v=bank[i]
    D=D*pow(d,96,M)%M
    Q=Q*pow(q,96,M)%M
    R=R*pow(u[0],96-j,M)*pow(v[0],j,M)%M
assert (11*D+42*Q+20*R)%M==27
assert all(bank[i][0]%3!=0 for i in indices)
assert 9*(13+96*len(indices))+10==8767

print("Verified dyadic profile generators, mod77-neutral itinerary constructor, and 27's record.")
