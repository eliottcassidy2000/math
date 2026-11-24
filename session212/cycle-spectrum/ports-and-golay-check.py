"""Check exact cycle seams, path corrections, and 24-contact Golay data."""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from math import comb, factorial

def comb_rows(m,x,e):
    n=m+4
    t,a,b=m+1,m+2,m+3
    g=[0]*n
    def arc(i,j):
        g[i]|=1<<j
    for i in range(m+2):
        for j in range(i+1,m+2):
            arc(i,j)
    for j in (0,t):
        arc(a,j)
        arc(j,b)
    for j in range(1,m+1):
        arc(b,j)
        if x>>(j-1)&1:
            arc(j,a)
        else:
            arc(a,j)
    arc(b,a) if e else arc(a,b)
    return g

def reads(g):
    F=(1<<len(g))-1
    @lru_cache(None)
    def dp(v,S,closure):
        if not S:
            return closure<0 or bool(g[v]>>closure&1)
        cand=g[v]&S
        t=0
        while cand:
            b=cand&-cand
            cand^=b
            t+=dp(b.bit_length()-1,S^b,closure)
        return t
    return dp(0,F^1,0),sum(dp(i,F^(1<<i),-1) for i in range(len(g)))

for m in range(7):
    for x in range(1<<m):
        C,H=reads(comb_rows(m,x,0))
        C1,H1=reads(comb_rows(m,x,1))
        assert (C,C1)==(x,x+1)
        if m:
            u=[x>>j&1 for j in range(m)]
            rev=sum(b<<(m-1-i) for i,b in enumerate(u))
            assert H1-H==2**(m+1)+1-rev-x+(u[0]==0)+(u[-1]==1)
            assert (H1-H==4)==all(u)
        else:
            assert H1-H==4

assert reads(comb_rows(3,6,0))==(6,81) and reads(comb_rows(3,7,0))==(7,81)
assert reads(comb_rows(2,2,1))==(3,31) and reads(comb_rows(2,3,1))==(4,33)

gen=sum(1<<i for i in (0,2,4,5,6,10,11))
basis=[]
for i in range(12):
    v=gen<<i
    basis.append(v | (v.bit_count()%2)<<23)
code={0}
for b in basis:
    code|={x^b for x in code}
assert len(code)==4096
assert all((a&b).bit_count()%2==0 for a in basis for b in basis)
assert Counter(w.bit_count() for w in code)=={0:1,8:759,12:2576,16:759,24:1}
assert all((v+(1<<j)) not in code for v in code for j in range(24))
assert len({v+e for v in code for e in (0,1)})==8192

N=2**24
sums=[N]
for p in range(1,9):
    z=N**(p+1)-sum(comb(p+1,i)*sums[i] for i in range(p))
    assert z%(p+1)==0
    sums.append(z//(p+1))
    diff=Fraction(sum(v**p for v in code),len(code))-Fraction(sums[p],N)
    if p<8:
        assert diff==0
    else:
        octads=sum(1<<sum(j for j in range(24) if v>>j&1)
                   for v in code if v.bit_count()==8)
        assert diff==Fraction(factorial(8)*octads,2**8)>0

print("Checked literal cycle/path seams through 10 vertices and exact Golay moments.")
