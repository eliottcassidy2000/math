#!/usr/bin/env python3
"""Independent literal slot-cover crosschecks of wrapped trace."""
from itertools import combinations, product
from math import factorial

def from_mask(n, mask):
    A=[[0]*n for _ in range(n)]
    for k,(a,b) in enumerate(combinations(range(n),2)):
        A[a][b]=(mask>>k)&1
        A[b][a]=1-A[a][b]
    return A

def wrap(A):
    n=len(A)
    return [row+[1,0] for row in A]+[[0]*n+[0,1],[1]*n+[0,0]]

def literal_slots(A,r):
    n=len(A); S=1<<n
    H=[0]*S; H[0]=1
    end=[[0]*n for _ in range(S)]
    for s in range(1,S):
        for v in range(n):
            if (s>>v)&1:
                t=s^(1<<v)
                end[s][v]=sum(end[t][u]*A[u][v] for u in range(n)) if t else 1
        H[s]=sum(end[s])
    C=[1]+[0]*(S-1); vals=[C[-1]]
    for _ in range(r):
        D=[0]*S
        for s in range(S):
            t=s
            while True:
                D[s]+=C[s^t]*H[t]
                if not t: break
                t=(t-1)&s
        C=D; vals.append(C[-1])
    return vals

def update(q):
    return [0]+[j*(q[j+1]+j*q[j]-q[j-1])
                for j in range(1,len(q)-1)]

def weighted_update(e,m):
    return [0]+[ ((e[j+1] if j+1<len(e) else 0)+j*j*e[j]-j*(j-1)*e[j-1]) %m
                for j in range(1,len(e))]

def require(b):
    if not b:
        raise RuntimeError("trace checker mismatch")

cases=0
for n in range(1,5):
    for mask in range(1<<(n*(n-1)//2)):
        A=from_mask(n,mask)
        q=literal_slots(A,7)
        for t in range(3):
            require(q[:5]==literal_slots(A,4))
            A=wrap(A); q=update(q); cases+=1

for m,n,mask in product(range(2,27), (2,3,4), (0,1,2,3)):
    A=from_mask(n,mask); q=literal_slots(A,40)
    k=next(j for j in range(1,41) if factorial(j)%m==0)
    e=[0]+[factorial(j-1)*q[j]%m for j in range(1,k+1)]
    for _ in range(12):
        require(q[1]%m==e[1])
        q=update(q);e=weighted_update(e,m)

def general_wrap(A,B):
    a,b=len(A),len(B)
    return ([row+[1]*b+[0] for row in A]
           + [[0]*a+row+[1] for row in B]
           + [[1]*a+[0]*(b+1)])

reservoir_cases=0
for mask_a,mask_b in product(range(8),repeat=2):
    A=from_mask(3,mask_a); B=from_mask(3,mask_b)
    q=literal_slots(A,5); u=literal_slots(B,5)
    actual=literal_slots(general_wrap(A,B),4)
    expected=[0]
    for j in range(1,5):
        d=u[j+1]-u[j]; d0=u[j]-u[j-1]
        expected.append(j*(d*q[j+1]+(u[j]-d+d0)*q[j]-d0*q[j-1]))
    require(actual==expected)
    reservoir_cases+=1

print(f"{cases} literal-profile comparisons and 3600 finite-state trace comparisons passed.")
print(f"{reservoir_cases} variable-reservoir literal-profile comparisons passed.")

def det_col_apply(M,x):
    a,b,c,d=M; h,q=x
    return ((a*h+b*q)%5,(c*h+d*q)%5)

def read_mod5_control(n,b):
    u=literal_slots(from_mask(b,(1<<(b*(b-1)//2))-1),3)
    d1=u[2]-u[1];d2=u[3]-u[2];s=(-1)**n
    return tuple(x%5 for x in (3*u[1]-u[2],d1,-2*d1,2*(u[2]+d1+(s-1)*d2)))

nonzero=set(product(range(5),repeat=2))-{(0,0)}
states=0
for n in (2,3):
    controls=[read_mod5_control(n,b) for b in ((1,) if n==2 else (1,3))]
    for source in nonzero:
        seen={source}
        while True:
            nxt=seen|{det_col_apply(M,x) for x in seen for M in controls}
            if nxt==seen: break
            seen=nxt
        require(seen==nonzero)
        states+=len(seen)
print(f"{states} mod-5 controlled reachability checks passed.")

r4=literal_slots(from_mask(5,4),3)
r12=literal_slots(from_mask(5,12),3)
require(r4[1:3]==[5,64] and r12[1:3]==[11,100])
triangle=literal_slots(from_mask(3,2),3)
R=[([0,1,2,3],-1),(r4,1),(triangle,0),([0,55,12800,0],-1)]
ops=[]
for u,g in R:
    h=u[1];q=u[2]
    M=tuple(x%5 for x in (3*h-q,q-h,-2*(q-h),2*(2*q-h)))
    ops.append((M,g))
seen={(1,1,4)}
while True:
    nxt=seen|{(g*c%3,*det_col_apply(M,(h,q))) for c,h,q in seen for M,g in ops}
    if nxt==seen: break
    seen=nxt
require(seen==set(product(range(3),range(5),range(5))))
print("75 joint mod-3/mod-5 physical register states passed.")
