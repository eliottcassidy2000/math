#!/usr/bin/env python3
from collections import Counter
from math import gcd, isqrt, log, expm1

def require(x, msg):
    if not x: raise RuntimeError(msg)

def prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n)+1))

m=290
V={r for r in range(m) if gcd(r,m)==1}
E={(r,(r+2)%m) for r in V if (r+2)%m in V}
runs=[]
for r in sorted(V):
    if (r-2)%m not in V:
        z=[]; x=r
        while x in V:
            z.append(x); x=(x+2)%m
        runs.append(z)
require(len(V)==112 and len(E)==81, 'wheel size')
require(Counter(map(len,runs)) == Counter({4:25,3:2,2:2,1:2}), 'paths')
require(sorted(x[0] for x in runs if len(x)==1)==[27,263], 'isolates')
require(sum(map(len,runs))==len(V), 'path partition')
for h in range(1,291):
    centers=[c for c in range(m) if gcd(c-h,m)==gcd(c+h,m)==1]
    expected=(4 if h%5==0 else 3)*(28 if h%29==0 else 27)
    require(len(centers)==expected, ('center count',h))
plus=next(n for n in range(2,6092) if n%m==1 and prime(n))
minus=next(n for n in range(2,6092) if n%m==m-1 and prime(n))
twin=next((n,n+2) for n in range(2,6090) if n%m==m-1 and prime(n) and prime(n+2))
require((plus,minus,twin)==(1451,2029,(6089,6091)), 'first primes')
require(prime(3) and prime(5) and (3,5) not in E, 'exceptional prime edge')
require(gcd(1,m)==1 and not prime(1), 'unit not prime')
require(all(gcd(n,m)==1 for n in (7,9)) and not prime(9), 'edge not pair of primes')
for N in (10000,100000,1000000):
    s=sum(-expm1(log(n)/n) for n in range(1,N+1) if gcd(n,m)==1)
    t=sum((1 if n%m==1 else -1)*expm1(log(n)/n) for n in range(1,N+1) if n%m in (1,m-1))
    print('N',N,'all_units',format(s,'.12f'),'residual',format(s+28/145*log(N)**2,'.12f'),'signed_pm',format(t,'.12f'))
print('vertices',len(V),'edges',len(E),'component lengths',sorted(Counter(map(len,runs)).items()),'isolates',[27,263])
print('first +1',plus,'first -1',minus,'first boundary twin',twin)
print('center counts h divisible by 5,29:', {(a,b):(4 if a else 3)*(28 if b else 27) for a in (False,True) for b in (False,True)})
print('PASS; finite tests only; primality via trial division; no asserts')
# Complete functional-graph census, using repeated iteration independently for every state.
D={a:((a*a-29)//4)%290 for a in range(1,290,2)}
basins={}
for a in D:
    path=[]; seen={}; x=a
    while x not in seen:
        seen[x]=len(path); path.append(x); x=D[x]
    cycle=path[seen[x]:]; y=min(cycle); canon=[]
    while not canon or y!=canon[0]:
        canon.append(y); y=D[y]
    basins.setdefault(tuple(canon),[]).append(seen[x])
require(len(basins)==8 and all(len(c)==3 for c in basins),'quadratic cycles')
require(Counter(map(len,basins.values()))==Counter({5:1,20:7}),'quadratic basins')
require((5,289,283) in basins and D[9]==13 and D[11]==23,'quadratic probes')
for cyc,tails in sorted(basins.items()):
    expected={0:3,1:2} if len(tails)==5 else {0:3,1:7,2:10}
    require(dict(sorted(Counter(tails).items()))==expected,'preperiod counts')
    print('cycle',cyc,'basin',len(tails),'preperiod',expected)
print('F-CENSUS PASS; modular dynamics only')
f=lambda a:(a*a-29)//4
for a in range(-7,8,2):
    x=a
    for _ in range(8):
        if x in (-7,5,-1): break
        x=f(x)
    require(x in (-7,5,-1),('small integer basin',a))
require([9,f(9),f(f(9)),f(f(f(9)))]==[9,13,35,299],'escaping modular cycle')
print('LIFT PASS; all odd |a|<=7 enter the cycle; 9 escapes despite period-3 residue')
