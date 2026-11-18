"""Keep-root binary carry; arithmetic forests; smallest paired seed boundary."""
from pathlib import Path
from math import isqrt
from itertools import combinations
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE.parent/'moonshine_prime_blocks'),str(HERE.parent/'dihedral_height')]
from rooted_pair import rooted,append
from frame_replay import join,add
from replay import profile

SEEDS=((4,9,18,21,35,15,63,127),(4,9,66,53,39,135,59,95))
# Physical order, low bit first.
PREP=(0,0,0,1)
V=35

def bits(x):return tuple((x>>j)&1 for j in range(x.bit_length()))

def build(A,seq):
    for z in seq:
        A=append(A,0,(0,),0)[0] if z else join(A,(0,))
    return A

def step(s,z):
    h,e,t,d=s
    return (h+t-e,d,2*t-d,d) if z else (h,0,2*t,d)

def state(A):
    return (*rooted(A,0)[:3],profile(tuple(j>>1 for j in A[1:]))[0])

def penalty(xs):
    J=[i for i,j in enumerate(xs) if j]
    return sum((1<<(j-i-1))+(j==i+1) for i,j in combinations(J,2))

def reach(A):
    f=(1<<len(A))-1;s=1
    while True:
        z=s
        for i,row in enumerate(A):
            if s>>i&1:z|=row
        if s==z:return z==f
        s=z

def strong(A):
    f=(1<<len(A))-1
    return reach(A) and reach(tuple(row^f^(1<<i) for i,row in enumerate(A)))

def census(path):
    S=set()
    for l in path.open():
        T=tuple(map(int,l.split()))
        for n in range(3,8):
            for I in combinations(range(7),n):
                S.add(tuple(sum((T[i]>>j&1)<<k for k,j in enumerate(I)) for i in I))
    C=set()
    H=set()
    for A in S:
        if not strong(A):continue
        h,b=profile(A)
        if h in (5,35):C.add((len(A),h,b))
        if len(A)<=6:H.add(h)
    assert C=={(4,5,32),(7,35,676)}
    assert 7 not in H and 21 not in H
    print('strong low profiles',sorted(C),'and strong spectrum through order six',sorted(H))

def run(fixtures=False):
    from functools import reduce
    for scale,T in zip((1,7),SEEDS):
        start=state(T);assert start==tuple(scale*x for x in (29,11,31,13))
        T=build(T,PREP)
        h,e,t,d=reduce(step,PREP,start)
        assert (h,e,t,d)==tuple(scale*x for x in (277,13,483,13))
        assert state(T)==(h,e,t,d)
        xs=bits(V)
        low=h+2*t*V-d*penalty(xs)
        gap=t-e-d*(V+xs[0])
        assert (low,gap)==(33749*scale,2*scale)
        for w in (0,1):
            y=reduce(step,(w,)+xs,(h,e,t,d))
            assert y[0]==low+w*gap
            A=build(T,(w,)+xs)
            if fixtures:print(len(A),0,*A)
            else:print('first bit',w,'exact output',y[0])
    assert all(all(p%d for d in range(2,isqrt(p)+1)) for p in (33749,33751))
    if fixtures:return
    for n in range(7,40001,2):
        if n in (7,21,63):continue
        if n%3==0:P=(n//3,)
        else:P=tuple((n-b)//3 for b in (-4,-2,2,4) if (n-b)%3==0)
        assert any(0<p<n and p%2 and p not in (7,21) for p in P)
    for d,T,r in (
        (0,(2,4,25,19,35,15,63,127),1),
        (1,(16,5,17,7,42,15,63,127),0),
        (2,(2,12,9,17,39,15,63,127),0)):
        h,e,t,*_=rooted(T,r)
        assert (h,t-h)==(23,d) and T[r].bit_count()==1
        U,b=append(T,r,(2,0),1)
        assert profile(U)[0]==3*h+2*d
        print('restart',h,'digit',d,'target',3*h+2*d)
    for h,T,I,K in (
        (25,(2,20,17,7,41,15,63,127),(27,9,18),(75,77,79)),
        (67,(4,9,18,37,43,71,31,127),(39,87),(203,205)),
        (67,(8,9,19,52,35,71,31,127),(10,),(201,))):
        assert profile(T)[0]==h
        assert tuple(profile(add(T,b))[0] for b in I)==K
        print('apex restart',h,K)
    # Pass an archived order-7 representative list to replay the stronger seed boundary.
    if len(sys.argv)>1:census(Path(sys.argv[1]))

if __name__=='__main__':run('--fixtures' in sys.argv)
