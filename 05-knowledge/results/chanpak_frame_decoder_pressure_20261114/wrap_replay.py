from fractions import Fraction as F
from itertools import product
from functools import cache

def subset_path(rows,mask):
    @cache
    def end(s,k):
        if s==1<<k: return 1
        t=s^(1<<k)
        return sum(end(t,j) for j,arcs in enumerate(rows) if t>>j&1 and arcs>>k&1)
    return sum(end(mask,j) for j in range(len(rows)) if mask>>j&1)

states=('', '0','1','00','01','10','11')
b=list(map(F,[1000,1000,901,1000,901,886,827]))
t=F(686,1000)
lam=F(1904,1000)
assert t**100>F(2,3)**93
assert lam**100<2**93
for i,a in enumerate(states):
    assert sum((t*b[0] if len(a)==2 and a[0]==x=='1' else b[states.index((a+x)[-2:])])
               for x in '01')<=lam*b[i],a

def wrap(rows):
    n=len(rows)
    return [sum(1<<j for j in range(1,n+1))]+[(a<<1)|(1<<(n+1)) for a in rows]+[1]

def at(rows,m):
    n=len(rows)
    H=[subset_path(rows,s) if s else 1 for s in range(1<<n)]
    total=0
    for color in product(range(m),repeat=n):
        masks=[0]*m
        for j,i in enumerate(color): masks[i]|=1<<j
        weight=1
        for x in masks: weight*=H[x]
        total+=weight
    return total

for n in range(5):
    edges=[(i,j) for i in range(n) for j in range(i+1,n)]
    for code in range(1<<len(edges)):
        rows=[0]*n
        for j,(a,c) in enumerate(edges):
            if code>>j&1: rows[a]|=1<<c
            else: rows[c]|=1<<a
        for m in range(1,5):
            assert at(wrap(rows),m)==m*m*at(rows,m)+m*(at(rows,m+1)-at(rows,m-1))
print('exact scanner certificates and all wrap color transfers through four middle vertices')
