"""Standalone exact replay: python replay_compass.py"""
from fractions import Fraction as F
from itertools import product
from compass import H, orbit, bar, AUTO

INV=((3,F(-5,2),F(1,2)),(F(-3,2),2,F(-1,2)),(F(1,3),F(-1,2),F(1,6)))
DOMAIN=(1,2,3)
for i,j in product(range(3),repeat=2):
    assert sum(INV[i][k]*[1,2**(k+1),3**(k+1)][j] for k in range(3))==(i==j)
def fit(f,d,basis_change=None):
    vals={p:f(*p) for p in product(DOMAIN,repeat=d)}
    for axis in range(d):
        mat=INV if basis_change is None else [
            [sum(basis_change[i][k]*INV[k][j] for k in range(3)) for j in range(3)]
            for i in range(3)]
        vals={p:sum(mat[p[axis]-1][q-1]*vals[p[:axis]+(q,)+p[axis+1:]] for q in DOMAIN) for p in vals}
    return vals

def word(k,l):return '11'+'0'*k+'2'*l+'311'
def fixed(k,l,i,e):
    j=(0,1,2,k+2,k+l+2,k+l+3,k+l+4)[i]
    return orbit(word(k,l))[j+e*(k+l+5)]
def split(p,q,r,i,e):
    k,l=(r,p+q) if i else (p+q,r)
    j=2+p+(k if i else 0)
    return orbit(word(k,l))[j+e*(k+l+5)]

def graphs(w):
    """Rows indexed by hubs 0,1 then the ascending spine."""
    n=len(w)+2
    out=[2,0]+[sum(1<<j for j in range(i+1,n)) for i in range(2,n)]
    for i,d in enumerate(map(int,w),2):
        for h in (0,1):
            if d>>(1-h)&1:out[h]|=1<<i
            else:out[i]|=1<<h
    return out

def count(g):
    n=len(g);full=(1<<n)-1
    D=[[0]*n for _ in range(full+1)]
    for v in range(n):D[1<<v][v]=1
    for s,row in enumerate(D):
        for i,npath in enumerate(row):
            if not npath:continue
            dst=g[i]&~s
            while dst:
                bit=dst&-dst;dst-=bit
                D[s|bit][bit.bit_length()-1]+=npath
    return sum(D[-1])

def cut_graph(g,u):
    allv=(1<<len(g))-1
    return [row^(allv^u if u>>j&1 else u) for j,row in enumerate(g)]
def witnessed_spine(g):
    n=len(g);seq=sorted(range(2,n),key=lambda i:(g[i]>>2).bit_count(),reverse=True)
    if any(not(g[i]>>j&1) for u,i in enumerate(seq) for j in seq[u+1:]):return None
    return ''.join(str(2*(g[0]>>j&1)+(g[1]>>j&1)) for j in seq)

from collections import Counter
assert tuple(len(tab[0]) for tab,_ in AUTO)==(18,18)
for n in range(1,5):
    for w in map(''.join,product('0123',repeat=n)):
        assert H(w)==count(graphs(w))
        g=graphs(w)
        got=[a for i in range(1<<n) if (a:=witnessed_spine(cut_graph(g,i<<2))) is not None]
        assert Counter(got)==Counter(orbit(w))
for w in ['1102311','110022311','231122331','112233110221101','011021032231','210113301230']:
    assert H(w)==count(graphs(w))

baseline=lambda l:1+89*2**l
basis=((1,2,3),(0,1,3),(0,0,1))
for e,i in product((0,1),range(9)):
    if i<7:
        f=lambda k,l:H(fixed(k,l,i,e))-baseline(l)
        d=2
    else:
        f=lambda p,q,r:H(split(p,q,r,i-7,e))-baseline(p+q if i==8 else r)
        d=3
    co=fit(f,d,basis)
    assert all(x>=0 and F(x).denominator==1 for x in co.values())
    if (e,i) in [(0,2),(0,3),(0,7)]:assert not any(co.values())
    elif (e,i) in [(1,4),(1,5)]:
        assert co==fit(lambda k,l:(30*3**k-90)*2**l,2,basis)
    else:assert co[(1,)*d]>0
    print('half',e,'chart',i,'min coefficient',min(co.values()),'constant',co[(1,)*d])

for k,l in product(DOMAIN,repeat=2):
    assert H(fixed(k,l,0,0))==1+36*2**k*3**l
for p,q,k in product(DOMAIN,repeat=3):
    K,L,X=2**k,2**(p+q),2**p
    C,alpha,beta=47+100*K+296*K*L,23+100*K,296*K-89
    h=C-alpha*X-beta*(L//X)
    assert H(split(p,q,k,1,0))==h
    assert alpha*X*X-(C-h)*X+beta*L==0
for k,l in product(range(1,8),repeat=2):
    hs=[H(split(p,l-p,k,1,0)) for p in range(1,l)]
    assert len(hs)==len(set(hs))
    assert all((v>u)==(2*p<=l) for p,(u,v) in enumerate(zip(hs,hs[1:]),1))
scores=list(map(H,orbit(word(2,2))))
assert scores==[1297,1779,357,357,357,2147,2387,4239,5921,3883,6761,6575,6791,4175,3611,1077,1077,2619]
print('confirmed sample scores',scores)
