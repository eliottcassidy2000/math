from itertools import permutations, product
from math import gcd

def j2(n):
    x = n
    z = n*n
    for p in range(2,n+1):
        if x % p == 0:
            while x % p == 0:
                x //= p
            z = z//(p*p)*(p*p-1)
        if x == 1: break
    return z

def genus(n):
    x=(n-6)*j2(n)
    assert x%24==0
    return 1+x//24

for n in [5,6,7,8,9,72]:
    v=j2(n)//2
    mu=n*v
    assert 2-2*genus(n)==v-mu//2+mu//3
    print(n,(mu,v,mu//2,mu//3,genus(n),2*genus(n)+1))
assert genus(72)==9505 and 1+(genus(72)-1)//2==4753
assert [(s,s%8,s%9) for s in range(72) if s%8 in (1,7) and s%9 in (1,8)]==[(1,1,1),(17,1,8),(55,7,1),(71,7,8)]

out = {3:(4,1,2),5:(8,1,3,6),11:(24,1,3,6,14),13:(18,8,3,5,14)}
for h, a in out.items():
    n=len(a)
    assert all((a[i]>>j&1)+(a[j]>>i&1)==1 for j in range(n) for i in range(j))
    assert sum(all(a[u]>>v&1 for u,v in zip(p,p[1:])) for p in permutations(range(n)))==h
units = {5**a*11**b*13**c%72 for a,b,c in product(range(6),repeat=3)}
assert units=={x for x in range(72) if gcd(x,72)==1}
assert {k*x%72 for k in [1,3,9] for x in units}==set(range(1,72,2))
assert (5**3*11%72,3*5*11%72)==(7,21)
print('all 36 odd residues; two absolute hostile products',5**3*11,3*5*11)


for a,b in product(range(1,9),range(5)):
    M=2**a*3**b
    U={1}; queue=[1]
    for x in queue:
        for h in (5,11,13):
            y=x*h%M
            if y not in U:
                U.add(y);queue.append(y)
    assert U=={x for x in range(M) if gcd(x,M)==1}
    assert {3**v*x%M for x in U for v in range(b+1)}==set(range(1,M,2))
print('all odd residues on 40 finite 2^a 3^b controls')
