from math import comb,isqrt
from itertools import combinations
A=(22,132,136,147,166,79,31,97)
B=(22,132,136,147,38,79,31,113)
def validate(a):
    for i,j in combinations(range(len(a)),2):
        if ((a[i]>>j)&1)+((a[j]>>i)&1)!=1:raise ValueError(a)
    for i in range(len(a)):
        seen=1<<i
        for _ in a:
            for j,r in enumerate(a):
                if seen>>j&1:seen|=r
        if seen!=(1<<len(a))-1:raise ValueError("not strong")
def hams(a):
    n=len(a); F=(1<<n)-1
    d=[[0]*n for _ in range(F+1)]; h=[1]+[0]*F
    for i in range(n):d[1<<i][i]=1
    for S in range(1,F+1):
        h[S]=sum(d[S])
        for i,c in enumerate(d[S]):
            if not c:continue
            m=a[i]&~S
            while m:
                b=m&-m;m-=b
                d[S|b][b.bit_length()-1]+=c
    return h
def coeff(a):
    h=hams(a);f=len(h)-1
    # R_j is the number of j nonempty labels covering each subset.
    prev=[1]+[0]*f;out=[]
    for _ in a:
        cur=[0]*(f+1)
        for s in range(1,f+1):
            t=s
            while t:
                cur[s]+=h[t]*prev[s^t];t=(t-1)&s
        out.append(cur[-1]);prev=cur
    return out
def prime(n):
    return n>1 and all(n%d for d in range(2,isqrt(n)+1))
def outer(B,k):
    n=len(B); t=k+1; q=n+t
    masks=[(((1<<t)-(2<<j))|(((1<<n)-1)<<t)|(int(j!=0)<<q)) for j in range(t)]
    masks.extend((r<<t)|(1<<q) for r in B)
    return tuple(masks)+(1,)
for a in (A,B):validate(a)
if tuple(x^y for x,y in zip(A,B))!=(0,0,0,0,128,0,0,16):raise ValueError("flip")
u,v=coeff(A),coeff(B)
if u[0]!=383 or v[0]!=389 or u[1]+2*u[0]!=5960 or v[1]+2*v[0]!=5960:raise ValueError("base")
if tuple(y-x for x,y in zip(u,v))!=(6,-12,-72,-48,0,0,0,0):raise ValueError("polynomial")
for a,h in ((A,383),(B,389)):
    for k in (0,1,2):
        b=outer(a,k);validate(b)
        if hams(b)[-1]!=2**k*5960+h:raise ValueError("outer")
    for k in (2,10):
        p=2**k*5960+h
        if not prime(p):raise ValueError("prime")
        print(k,p,7*p)
print("exact polynomial and literal controls verified")
