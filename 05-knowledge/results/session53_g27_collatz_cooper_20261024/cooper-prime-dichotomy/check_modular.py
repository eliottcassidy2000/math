import json
from math import comb

def req(x,s):
    if not x:raise RuntimeError(s)
N=150

def mul(a,b,m):
    out=[0]*(N+1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b[:N+1-i]):out[i+j]=(out[i+j]+x*y)%m
    return out

def power(a,k,m):
    out=[1]+[0]*N
    while k:
        if k%2:out=mul(out,a,m)
        a=mul(a,a,m);k//=2
    return out

def inv(a,m):
    out=[pow(a[0],-1,m)]+[0]*N
    for n in range(1,N+1):out[n]=-sum(a[k]*out[n-k] for k in range(1,n+1))*out[0]%m
    return out

def parameter(m):
    h=[1]+[0]*N
    for j in range(1,N+1):
        f=[0]*(N+1)
        for k in range(N//j+1):f[k*j]=comb(k+3,3)%m
        h=mul(h,f,m)
        if 7*j<=N:
            f=[0]*(N+1)
            for k in range(min(4,N//(7*j))+1):f[k*7*j]=(-1)**k*comb(4,k)%m
            h=mul(h,f,m)
    h=[0]+h[:-1];h2=mul(h,h,m)
    den=[(13*h[i]+49*h2[i]+(i==0))%m for i in range(N+1)]
    return mul(h,inv(den,m),m)
rows=[]
for p,m,corrections,bound in [(2,4,[(3,2)],24),(3,9,[],32),(7,49,[(8,14),(9,28)],140)]:
    t=parameter(m)
    lhs=[t[i//p] if i%p==0 else 0 for i in range(N+1)]
    rhs=power(t,p,m)
    for k,a in corrections:
        v=power(t,k,m);rhs=[(x+a*y)%m for x,y in zip(rhs,v)]
    req(lhs==rhs,'Frobenius lift identity')
    rows.append(dict(p=p,modulus=m,verified_through=N,required_Sturm_bound=bound))
# Constant term f=theta^3*s uses s=(1-49h²)/(1+13h+49h²), but f=theta*(qdt/dq)/t is easier.
m=49;t=parameter(m);theta=[1]+[0]*N
for d in range(1,N+1):
    chi=0 if d%7==0 else 1 if d%7 in (1,2,4) else -1
    for k in range(d,N+1,d):theta[k]+=2*chi
unit=t[1:]+[0]
deriv=[(i+1)*t[i+1]%m for i in range(N)]+[0]
f=mul(theta,mul(deriv,inv(unit,m),m),m)
A1=[0]*(N+1);A2=[0]*(N+1)
for n in range(1,N+1):
    for d in range(1,n+1):
        if n%d==0:
            ch=lambda k: 0 if k%7==0 else 1 if k%7 in (1,2,4) else -1
            A1[n]+=ch(d)*d*d;A2[n]+=ch(n//d)*d*d
expected=[((i==0)-7*pow(8,-1,m)*A1[i]-49*pow(8,-1,m)*A2[i])%m for i in range(N)]
req(f[:N]==expected,'weight3 Eisenstein identity')
print(json.dumps({'scope':'FINITE-EXACT q-series checks; all-order conclusions additionally require stated modularity and Sturm theorem','lift_checks':rows,'weight3_identity_through':N-1},indent=2,sort_keys=True))
