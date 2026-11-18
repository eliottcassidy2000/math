from collections import Counter
from math import gcd, prod
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path
import csv

OGG=(2,3,5,7,11,13,17,19,23,29,31,41,47,59,71)
P=tuple(p for p in OGG if p not in (2,7))

def read():
    with (Path(__file__).resolve().parent/'ogg_witnesses.tsv').open() as f:
        return {int(r['count']):tuple(map(int,r['outgoing_masks'].split(',')))
                for r in csv.DictReader(f,delimiter='\t')}
BANK=read()

def dp(rows):
    n=len(rows); d=[[0]*n for _ in range(1<<n)]
    assert all(not m>>i&1 and m<1<<n for i,m in enumerate(rows))
    assert all((rows[i]>>j&1)+(rows[j]>>i&1)==1 for i,j in combinations(range(n),2))
    for i in range(n):d[1<<i][i]=1
    for s in range(1,1<<n):
        for i,z in enumerate(d[s]):
            if not z:continue
            m=rows[i]&~s
            while m:
                t=m&-m;m^=t;d[s|t][t.bit_length()-1]+=z
    return d
def h(rows):return sum(dp(rows)[-1])
def join(A,B):
    n=len(A)
    return tuple(a|(((1<<len(B))-1)<<n) for a in A)+tuple(b<<n for b in B)
def flip(T,i,j):
    out=list(T);out[i]^=1<<j;out[j]^=1<<i;return tuple(out)

def route(n):
    """m, known factors, unresolved exceptional leaf if any."""
    if n<=0 or n%2==0:raise ValueError("positive odd input required")
    e=Counter();m=n
    for p in (*P,7):
        while m%p==0:e[p]+=1;m//=p
    b=e.pop(7,0);factor=[]
    if b>=2:
        if b%2:factor.append(343);b-=3
        factor.extend([49]*(b//2))
    elif b==1:
        q=next((p for p in P if p!=3 and e[p]),None)
        if q is not None:factor.append(7*q);e[q]-=1
        elif e[3]>=2:factor.append(63);e[3]-=2
        else:return m,(),7*3**e[3]*m
    factor.extend(p for p in e for _ in range(e[p]))
    assert m*prod(factor)==n
    return m,tuple(factor),None

def smooth_witness(n,base=(0,)):
    m,f,leaf=route(n)
    if leaf:return None
    assert h(base)==m
    items=([base] if m>1 else [])+[BANK[k] for k in f]
    R=()
    for T in items:R=join(R,T)
    return R or (0,)

CORE=(6,12,8,1)
def main():
    want={1,49,63,343,*P,*[7*p for p in P if p!=3]}
    assert want==BANK.keys()
    assert all(h(R)==n for n,R in BANK.items())
    print('Ogg chamber: literal seeds and largest order:',len(BANK),max(map(len,BANK.values())))
    print('relative odd density of cores:',prod(Fraction(p-1,p) for p in (*P,7)))
    # Factorization coverage is proved; here test actual compiled graphs on modest controls.
    ncases=0
    for n in range(1,100000,2):
        m,f,e=route(n)
        assert (e in (None,7*m,21*m)) and (e is not None or m*prod(f)==n)
        if m==1 and e is None:
            # Verify the non-spanning witness count by exact product of proven seed counts.
            assert prod(f)==n
            if len(f)<=3 and sum(len(BANK[v]) for v in f)<=12:
                assert h(smooth_witness(n))==n;ncases+=1
    print('compiled controls checked by literal DP:',ncases)
    B=join(CORE,(2,0));D=dp(CORE);E=tuple(D[-1])
    IN=tuple(sum(1<<j for j in range(4) if CORE[j]>>i&1) for i in range(4))
    F=dp(IN); R=tuple(sum((sum(D[15^t]) if t!=15 else 1)*F[t][x]
                            for t in range(1,16) if t>>x&1) for x in range(4))
    assert E==(1,1,2,1) and R==(7,9,6,5)
    print('one common host:',B,'end/two-cover reads:',E,R)
    for i,j,k in ((1,3,1),(0,2,3),(2,0,5)):
        p,q=flip(B,i,4),flip(B,j,5)
        assert (h(p),h(q))==(14-k,14+k)
        print('interior mirror pair:',h(p),h(q),'one branch each:',p,q)
    # Complete labeled n<=6 census for the independent small absence check.
    first={}
    for n in range(1,7):
        pairs=list(combinations(range(n),2))
        for z in range(1<<len(pairs)):
            T=[0]*n
            for k,(i,j) in enumerate(pairs):
                if z>>k&1:i,j=j,i
                T[i]|=1<<j
            x=h(T);first.setdefault(x,n)
    print('middle-out least orders:',[(n,first[n]) for n in sorted(range(9,20,2),key=lambda n:(abs(n-14),n))])
    assert 7 not in first and 21 not in first
    assert [first[n] for n in range(9,20,2)]==[5,5,5,5,6,6]
    print('after the interior: 7 and 21 absent through order 6; proof gives the necessary global cutoff')

if __name__=='__main__':main()
