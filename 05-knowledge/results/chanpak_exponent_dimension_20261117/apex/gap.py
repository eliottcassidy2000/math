#!/usr/bin/env python3
import argparse,json
from math import isqrt
def F(x):
    total=1
    while x:
        p,x0=divmod(x,2)
        if not x0:total+=p+(p%2)
        x=p
    return total
def count(u,w):return u*w+F(u)+F(w)-1+(u%2)*(w%2)
def physical(u,w):
    A=tuple(map(int,bin(u)[2:])) if u else ()
    B=tuple(1-int(b) for b in reversed(bin(w)[2:])) if w else ()
    z=A+B;n=len(z)
    return tuple(((1<<n)-(1<<(i+1)))+((1-b)<<n) for i,b in enumerate(z))+(sum(k<<i for i,k in enumerate(z)),)
def pair(a,b):
    assert a>=1 and b>=1
    return (2*a,2*b+1),(2*a+1,2*b)
def inverse_gap(delta,step=0):
    if delta%2:raise ValueError('gap must be even')
    if step<0:raise ValueError('step must be nonnegative')
    d=delta//3
    residue=delta-3*d
    if residue==2:d+=1;parity=1
    elif residue==1:parity=0
    else:parity=0
    a=max(1,1-d)
    if a%2!=parity:a+=1
    a+=2*step;b=a+d
    P=pair(a,b)
    assert count(*P[1])-count(*P[0])==delta
    return {'gap':delta,'a':a,'b':b,'addresses':P,'counts':[count(*x) for x in P],'orders':[len(physical(*x)) for x in P], 'rows':[physical(*x) for x in P]}
def target_unit(N):
    if N<1:return None
    a=(isqrt(N)-1)//2
    if a<1 or not a%2:return None
    P=pair(a,a+1)
    for s,x in enumerate(P):
        if count(*x)==N:return a,s,x,physical(*x)
    return None
def H(rows):
    n=len(rows);D=[[0]*n for _ in range(1<<n)]
    for i in range(n):D[1<<i][i]=1
    for s in range(1,1<<n):
        for i,z in enumerate(D[s]):
            m=rows[i]&~s
            while z and m:
                y=m&-m;m-=y;D[s|y][y.bit_length()-1]+=z
    return sum(D[-1])
def canon(rows):
    from itertools import groupby,permutations,product
    g=[list(q) for _,q in groupby(sorted(range(len(rows)),key=lambda i:rows[i].bit_count()), lambda i:rows[i].bit_count())]
    return min(tuple(sum((rows[i]>>v&1)<<j for j,v in enumerate(x)) for i in x)
               for y in product(*(list(permutations(h)) for h in g)) for x in (sum(y,()),))
def test():
    tests=0
    for a in range(1,13):
        for b in range(1,13):
            P=pair(a,b);R=list(map(lambda p:physical(*p),P))
            assert len(R[0])==len(R[1])
            assert sum((x^y).bit_count() for x,y in zip(*R))==4
            h,g=map(H,R)
            assert (h,g)==tuple(count(*p) for p in P)
            assert g-h==3*(b-a)+(b%2)-(a%2)
            tests+=1
    for d in range(-200,201,2):
        for step in (0,1,9):
            z=inverse_gap(d,step)
            assert z['counts'][1]-z['counts'][0]==d
    for k in range(1,5):
        a=2**k-1;P=pair(a,a+1)
        assert [H(physical(*p)) for p in P]==[4**(k+1)-1,4**(k+1)+1]
        assert len(physical(*P[0]))==2*k+4
    for p in range(1,5):
        for q in range(1,5):
            low=physical(2**(p+1)-2,2**(q+1)+1)
            upper=physical(2**(p+1)-1,2**(q+1))
            assert H(upper)==2**(p+q+2)+1
            assert H(low)==2**(p+q+2)+3*(2**p-2**q)-1
    q=(126,124,248,240,96,64,128,51)
    assert q==physical(6,9)
    assert canon(q)==canon((8,1,3,198,15,31,55,119))
    assert canon(q)!=canon((8,9,3,148,7,31,63,119))
    vals={count(*p) for a in range(1,151,2) for p in pair(a,a+1)}
    assert all(bool(target_unit(n))==(n in vals) for n in range(1,max(vals)+1))
    print('PASS',tests,'literal two-arc pairs; 603 signed gap inverses, four power pairs, 16 unequal/equal dimension balances, exact 63 birth class, and unit endpoint inverse')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--gap',type=int);p.add_argument('--step',type=int,default=0);p.add_argument('--target',type=int)
    a=p.parse_args()
    if a.gap is not None:print(json.dumps(inverse_gap(a.gap,a.step)))
    elif a.target is not None:print(json.dumps(target_unit(a.target)))
    else:test()
