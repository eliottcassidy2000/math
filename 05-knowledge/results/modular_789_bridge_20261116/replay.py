"""Independent arithmetic controls for the even table and marked clock."""
from math import gcd, lcm

def v(n,p):
    if n==0:return 10**6
    n=abs(n);e=0
    while n%p==0:
        n//=p;e+=1
    return e
def ord_(x,m):
    if m==1:return 1
    a=x%m;k=1
    while a!=1:
        a=(a*x)%m;k+=1
    return k
def actual(a,x,m):
    seen={};j=0
    while x not in seen:
        seen[x]=j
        x=(a*a*x-a+1)%m
        j+=1
    return seen[x],j-seen[x]
for h in range(1,60):
    m=2*h
    for i in range(m):
        assert i*h%m==h*(i%2)
    assert h*h%m==h*(h%2)
    R={tuple(sorted((min(i,m-i),min(j,m-j)))) for i in range(m) for j in range(m)}
    assert len(R)==(h+1)*(h+2)//2 and len([1 for i,j in R if j==h])==h+1

for m in range(2,241):
    s,b=v(m,2),v(m,3)
    u=m//(2**s*3**b)
    fixed=[]
    for x in range(m):
        z=4*x-1
        expected=((max(0,b-v(z,3))+1)//2,lcm(2**max(0,s-1),ord_(9,u//gcd(u,z))))
        assert actual(3,x,m)==expected
        if (9*x-2)%m==x:
            fixed.append(x)
    assert len(fixed)==(1 if m%2 else 2 if m%4==2 else 0)
    if m%4==2:
        assert all((3*x-1)%m==(m//2-x)%m for x in fixed)

for a in range(3,24,2):
    for s in range(1,11):
        for x in range(2**s):
            assert actual(a,x,2**s)==(0,2**max(0,s-v(a-1,2)))

for a in (3,7,9):
    m=72
    seen=set();cycles=[]
    for x in range(m):
        walk=[]
        while x not in walk and x not in seen:
            walk.append(x);x=(a*a*x-a+1)%m
        if x in walk:
            cycles.append(tuple(walk[walk.index(x):]))
        seen.update(walk)
    expect={3:[4]*2,7:[4]*6+[12]*4,9:[1]*8}[a]
    assert sorted(map(len,cycles))==expect
    print(a,'modulo 72 cycles:',sorted(map(len,cycles)))
print('even cross, stitched lifts, odd parts and dyadic periods checked')

assert ord_(49,72)==3 and ord_(49,576)==12
for x in range(72):
    y=x
    for _ in range(3):y=(49*y-6)%72
    assert y==(x-18)%72
print("Carmichael multiplier/affine order mismatch verified")
