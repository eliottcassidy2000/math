from itertools import permutations,product,combinations
from fractions import Fraction as F
from collections import deque

def rowsflip(A,i,j):
    B=list(A);B[i]^=1<<j;B[j]^=1<<i;return tuple(B)
def relabel(A,P):
    B=[0]*len(A)
    for i in range(len(A)):
        for j in range(len(A)):
            if A[i]>>j&1:B[P[i]]|=1<<P[j]
    return tuple(B)
def H(A):
    return sum(all(A[i]>>j&1 for i,j in zip(P,P[1:])) for P in permutations(range(len(A))))
def strong(A):
    for i in range(len(A)):
        S=1<<i
        while True:
            T=S
            for j in range(len(A)):
                if S>>j&1:T|=A[j]
            if T==S:break
            S=T
        if S!=(1<<len(A))-1:return False
    return True
A=(4,9,66,37,15,23,59)
assert strong(A) and H(A)==71
children=(
 ((1,2),(2,64,11,35,15,23,61)),
 ((2,3),(8,5,33,70,15,27,55)),
 ((5,6),(4,9,34,69,15,27,55)),
)
for (i,j),C in children:
    P=list(range(7));P[i],P[j]=P[j],P[i]
    assert relabel(rowsflip(A,i,j),P)==C
    assert strong(C) and H(C)==69
for (_,B),(_,C) in combinations(children,2):
    assert not any(relabel(B,P)==C for P in permutations(range(7)))

states=((6,4,0),(2,4,1),(0,5,1))
def at(t):
    n=len(t);R=[]
    for i,k in enumerate(t):
        s=((1<<(3*n))-1)^((1<<(3*i+3))-1)
        R.extend((a<<(3*i))|s for a in states[k])
    return tuple(R)
def put(v,i,w):
    x=list(v);x[i]=w;return tuple(x)
def theta(k):
    x=(0,)*k;y=put(x,0,2)
    for j in range(1,k):
        yield x,put(x,j,1),put(x,j,2),put(put(x,j,2),0,1),put(y,j,2),put(y,j,1),y
def vs(k):
    V=[];S={(F(0),F(0))}
    for j in range(k):
        for s in range(12*5**j+1):
            u=(F(1-s*s,1+s*s),F(2*s,1+s*s))
            good=True
            for a,b in S-{(0,0)}:
                for c in (-2,-1,1,2):
                    x,y=a+c*u[0],b+c*u[1];r=x*x+y*y
                    if r in (0,1):good=False
            if good:break
        else:raise AssertionError("greedy bound")
        V.append(u);S={(a+c*u[0],b+c*u[1]) for a,b in S for c in range(-2,3)}
    return V
for k in (2,3,4):
    path=tuple(theta(k));S=set(sum(path,()))
    assert len(S)==5*(k-1)+2
    edges={frozenset((u,v)) for P in path for u,v in zip(P,P[1:])}
    V=vs(k)
    co={t:tuple(sum(t[i]*V[i][j] for i in range(k)) for j in (0,1)) for t in S}
    R={t:at(t) for t in S}
    assert len(set(co.values()))==len(S)
    for u,v in combinations(S,2):
        e=frozenset((u,v)) in edges
        assert e==(sum((a^b).bit_count() for a,b in zip(R[u],R[v]))==2)
        assert e==(sum((a-b)**2 for a,b in zip(co[u],co[v]))==1)
    assert sum((a-b)**2 for a,b in zip(co[path[0][0]],co[path[0][-1]]))==4
    if k==2:
        assert all(H(R[t])==(3 if 1 in t else 1) for t in S)
print("literal star, soft theta paths, rational unit coordinates, and exact port distance all pass")
