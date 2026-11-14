from decoder import *
from itertools import product
from random import Random

def compositions(n, start=()):
    if n==0: yield start; return
    for s in range(1,n+1): yield from compositions(n-s,start+(s,))

def literal_word(g):
    a=adjacency_rows(g)
    n=len(a)
    d=[{} for _ in range(1<<n)]
    for i in range(n): d[1<<i][i]=1
    for m,c in enumerate(d):
        for u,z in c.items():
            nxt=a[u]&~m
            while nxt:
                b=nxt&-nxt
                nxt^=b
                v=b.bit_length()-1
                f=d[m|b]
                f[v]=f.get(v,0)+z
    return sum(d[-1].values())

t={0:0,1:1,2:1}
for i in range(-1,-111,-1):
    t[i]=t[i+3]-t[i+2]-t[i+1]
for q in range(100):
    e,f,g=kernel(q)
    assert e*g-f*f==-t[-q-4]
assert kernel(13)==(47*47,47*56,56*56)
s=Random(337)
for _ in range(2000):
    m=matrix(s.randrange(30))
    for _ in range(s.randrange(6)):
        m=mul(mul(m,(1,0,0,1<<s.randrange(1,10))),matrix(s.randrange(30)))
    c=coeff(m,s.randrange(2),s.randrange(2))
    assert max(c)==c[-1]
    a,b=s.randrange(1,50),s.randrange(1,50)
    n=value(c,1<<a,1<<b)
    assert (a,b) in tuple(inv2(n,c))
    A,B,C,D=c
    assert (D*(1<<a)+C)*(D*(1<<b)+B)==D*n+B*C-A*D
expected={}
for n in range(3,13):
    for blocks in compositions(n-1):
        g=tuple(s-1 for s in blocks)
        if not g[0] or not g[-1]: continue
        H=direct_word(g)
        if n<=9: assert literal_word(g)==H
        expected[H]=better(expected.get(H),g)
for h,g in expected.items():
    assert solve(h)==g,(h,g,solve(h))
for h in (7,21,63,385):
    assert minimum(h) is None
for a,b in ((1,1),(1,2),(2,2)):
    assert direct_word((a,)+(0,)*13+(b,))==(47+56*(1<<a))*(47+56*(1<<b))
assert minimum(573)==(12,(1,0,0,3,0,1))
print(f'{len(expected)} canonical values in the order-12 control; 2000 asymmetric generated frames; 100 determinant checks')

def hubs(w):
    m=len(w)
    def arc(i,j):
        if i==j: return False
        if i>j: return not arc(j,i)
        if i==m or j<m: return True
        return not (int(w[i])>>(1-(j-m))&1)
    return tuple(sum(1<<j for j in range(m+2) if arc(i,j)) for i in range(m+2))

def subset_path(rows,mask):
    from functools import lru_cache
    @lru_cache(None)
    def end(s,k):
        if s==1<<k: return 1
        t=s^(1<<k)
        return sum(end(t,j) for j,arcs in enumerate(rows) if t>>j&1 and arcs>>k&1)
    return sum(end(mask,j) for j in range(len(rows)) if mask>>j&1)

w=hubs('100110')
mask=(1<<len(w))-1
assert subset_path(w,mask)==63
assert tuple(subset_path(w,mask^(1<<j)) for j in range(len(w)))==(11,37,37,33,33,25,25,1)
for digits,count in [('110',11),('1200',37),('1120',33),('1121',25),('10301',49),('1012010',343)]:
    rows=hubs(digits)
    assert subset_path(rows,(1<<len(rows))-1)==count,(digits,count)
print('H=63 and every deletion, short lower-order representatives, and 49/343 seeds verified')
