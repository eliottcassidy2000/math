from fractions import Fraction
from itertools import product
from math import gcd, lcm, prod
from pathlib import Path
from runpy import run_path

bank = run_path(str(Path(__file__).with_name("third-slot-bank.py")))
data, markers = bank["data"], bank["markers"]
def support(n):
    p, out = 2, set()
    while p*p <= n:
        if n % p == 0:
            out.add(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        out.add(n)
    return out
D0 = prod(d for d, q, u, v in data)
Q0 = prod(q for d, q, u, v in data)
Rmin = prod(gcd(u[0],v[0]) for d,q,u,v in data)
assert gcd(gcd(D0,Q0),Rmin) == 5011875
assert all(u[0] % p == 0 and u[0] % (p*p) and v[0] % p
           for (_,_,u,v),p in zip(data,markers))
supports = [(support(d),support(q), support(u[0]) | support(v[0])) for d,q,u,v in data]
best = []
for mask in range(1,1<<len(data)):
    ii = tuple(i+1 for i in range(len(data)) if mask>>i&1)
    if best and len(ii) < len(best[0]):
        continue
    a,b,c = (set().union(*(supports[i-1][j] for i in ii)) for j in range(3))
    if a & b & c:
        continue
    if best and len(ii) > len(best[0]):
        best = []
    best.append(ii)
assert best == [(1,2,3,4,7,8,12),(1,2,3,4,8,9,12)]

period = lcm(*(p-1 for p in markers))
assert all(D0%p and Q0%p for p in markers)
# Each unarmed base factor is a unit at all markers; armed factors
# reduce to 0 at their private marker and still 1 elsewhere.
for i,(_,_,u,v) in enumerate(data):
    for j,p in enumerate(markers):
        assert pow(v[0],period,p) == 1
        assert pow(u[0],period,p) == (i!=j)
for r,corners in [(0,(53,329,113,371)),(1,(73,445,157,511))]:
    assert corners == (53+20*r,329+116*r,113+44*r,371+140*r)
    assert all(h%p for p in markers for h in corners)

L = {pow(77,j,977):j for j in range(61)}
assert len(L)==61 and pow(77,61,977)==1
logs = [L[pow((u[0]*pow(v[0],-1,977))%977,16,977)] for d,q,u,v in data]
assert logs == [17,15,48,49,1,43,50,33,30,42,40,52,26,10]
assert len(set(min(a,61-a) for a in logs))==14
def value(k):
    return min(min(a*k%61,(-a*k)%61) for a in logs)
assert value(9)==max(map(value,range(61)))==5

ids = best[1]
r = [Fraction(max(data[i-1][2][0],data[i-1][3][0]),
              min(data[i-1][2][0],data[i-1][3][0])) for i in ids]
b = r.pop(1)
assert b == Fraction(302043,301803)
phases,Cmax = [],0
for e in product(range(3),repeat=len(r)):
    A = prod(u**j for u,j in zip(r,e))
    c = 0
    while A>=b:
        A /= b
        c += 1
    Cmax = max(Cmax,c)
    phases.append(A)
assert len(set(phases))==729 and min(phases)==1 and Cmax==83
phases = sorted(phases)+[b]
G = max(y/x for x,y in zip(phases,phases[1:]))
assert 1+Fraction(4337,10**9)<G<1+Fraction(4338,10**9)
print("primitive subbanks, 14-marker wheel, 61-clock, 729 anchors pass")

P,I,T=77,[2,6,11,13],[3,1,1,3]
ds=[data[i] for i in I]
rs=[(u[0]*pow(v[0],-1,P))%P for d,q,u,v in ds]
got={prod(pow(r,e,P) for r,e in zip(rs,k))%P
     for k in product(*(range(t+1) for t in T))}
assert rs==[25,59,16,39] and got=={j for j in range(77) if gcd(j,77)==1}
D=prod(b[0]**t for b,t in zip(ds,T)); Q=prod(b[1]**t for b,t in zip(ds,T))
C=11*D+42*Q
assert (C%77,C%3,C%5)==(35,1,3)
assert all(u[0]%3==v[0]%3==0 for d,q,u,v in ds)
print("82-vertex compact mod-77 wheel passes")
