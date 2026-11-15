from collections import Counter
from itertools import product
from math import tau
from cmath import exp

v0 = (9,89),(5,53)
def st(v,m):
    a,b=v
    return b%m,(9*b+8*a)%m
def orbit(v,m,L):
    x,y=(z%m for z in v)
    out=[]
    for k in range(L):
        out.append((x,y))
        x,y=st((x,y),m)
    return out,(x,y)
def phase(v,m,k,l):
    a,_=orbit(v,m,k+l*(m-1)+1)
    return Counter(a[k+l*j][0] for j in range(m))
def literal(back):
    n=len(back)
    adj=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            if j-i in back[j]: adj[j]|=1<<i
            else: adj[i]|=1<<j
    d=[[0]*n for _ in range(1<<n)]
    for i in range(n): d[1<<i][i]=1
    for s in range(1,1<<n):
        for i,u in enumerate(d[s]):
            if u:
                ns=adj[i]&~s
                while ns:
                    b=ns&-ns; ns-=b
                    d[s|b][b.bit_length()-1]+=u
    return sum(d[-1])
A,B=((),(2,)),((1,),(2,))
for units,vals in ((A+B,(9,89,873)),(B+A,(5,53,517))):
    for j,h in enumerate(vals,1): assert literal(((),)+units*j)==h
for idx,v in enumerate(v0):
    p,end=orbit(v,3,4)
    assert tuple(x for x,y in p)==((0,2,0,1),(2,2,1,1))[idx] and end==p[0]
    p,end=orbit(v,5,24)
    assert set(p)==set(product(range(5),repeat=2))-{(0,0)} and end==p[0]
    assert Counter(x for x,y in p)=={0:4,1:5,2:5,3:5,4:5}
    p,end=orbit(v,8,1)
    assert p[0]==end==((1,1),(5,5))[idx]
    p,_=orbit(v,27,4)
    ds=[(79*x+97*y)%9 for x,y in p]
    assert tuple(d%3==0 for d in ds)==((False,)*4,(False,True,False,True))[idx]
    assert all(d for d in ds)
    for r in range(4,11):
        m=2**r;b=(r+2)//3
        p,end=orbit(v,m,b+max(1,2**(r-4)))
        reads=p[b:]
        assert end==reads[0]
        assert sorted(x for x,y in reads)==list(range((9,5)[idx],m,16))
    for r in range(3,7):
        m=3**r
        for i in range(4):
            q=9 if ds[i]%3 else 27
            c=phase(v,m,i,4)
            assert set(c)==set(range(p0:=orbit(v,q,i+1)[0][-1][0],m,q))
            assert len(set(c.values()))==1
    for r in range(2,6):
        m=5**r
        p,_=orbit(v,m,24*5**(r-1))
        vals=[sum(exp(tau*1j*p[i][0]/m) for i in range(j,len(p),24)) /
            5**(r-1) for j in range(24)]
        assert sum(abs(z)>1e-10 for z in vals)==4
        assert all(abs(abs(z)-5**(-(r-2)/2))<1e-10 for z in vals if abs(z)>1e-10)
print('literal starts and dyadic/ternary/quinary twists verified')
for idx,v in enumerate(v0):
    for a,b,c,wants_zero in [(8,1,1,True),(2,4+idx,1,True),(2,1,3,False)]:
        m=2**a*3**b*5**c
        from math import lcm
        period=lcm(2**max(a-4,0),4*3**max(b-2,0),24*5**max(c-1,0))
        burn=(a+2)//3
        p,end=orbit(v,m,burn+period)
        assert end==p[burn]
        s=abs(sum(exp(tau*1j*x/m) for x,y in p[burn:])/period)
        assert s<=(1e-10 if wants_zero else 5**(-(c-2)/2)/6+1e-10)
print('mixed conductors verified')
U=[0,1]
for i in range(80):U.append(9*U[-1]+8*U[-2])
for d in range(2,20):
    for m in range(2,80//d+1): assert U[d*m]>U[d]>1 and U[d*m]%U[d]==0
for j,(x,y) in enumerate(orbit(v0[0],10**100,30)[0],1): assert (x,y)==(U[j+1],U[j+2])
for j,(x,y) in enumerate(orbit(v0[1],10**100,30)[0],1): assert (x,y)==(U[j+1]-4*U[j],U[j+2]-4*U[j+1])
print('Lucas prime-index sieve verified')
