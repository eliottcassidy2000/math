from near_twins import path_count
from functools import cache
from itertools import combinations

def Q_values(rows, js):
    n=len(rows);N=(1<<n)-1
    D=[[0]*n for _ in range(N+1)]
    for i in range(n):D[1<<i][i]=1
    for S in range(1,N+1):
        for j,e in enumerate(D[S]):
            if e:
                R=rows[j]&~S
                while R:
                    b=R&-R;R-=b
                    D[S|b][b.bit_length()-1]+=e
    H=[1]+list(map(sum,D[1:]))
    @cache
    def q(S,j):
        if j==1:return H[S]
        T=S; w=0
        while True:
            w+=H[T]*q(S^T,j-1)
            if not T:return w
            T=(T-1)&S
    return [q(N,j) if j else int(N==0) for j in js]

def wrap(A):
    n=len(A)
    return [a|(1<<n) for a in A]+[1<<(n+1),(1<<n)-1]

def fixed(A,r):
    if len(A)==1:return 1
    n=len(A)
    if any(bool(A[i]>>j&1)!=bool(A[r[j]]>>r[i]&1) for i in range(n) for j in range(n)):
        raise ValueError('not anti')
    m=n//2
    def walk(v,S):
        if S==(1<<n)-1:
            return A[v]>>c&1 if n%2 else A[v]>>r[v]&1
        return sum(walk(w,S|(1<<w)|(1<<r[w])) for w in range(n)
                   if not S>>w&1 and (v<0 or A[v]>>w&1))
    c=next((i for i in range(n) if r[i]==i),None)
    return walk(-1,1<<c if c is not None else 0)

pal=[sum(1<<((i+d)%7) for d in (1,2,4)) for i in range(7)]
seeds=[([0],[0],1,1), (pal,[(-i)%7 for i in range(7)],189,9),
       ([446,312,66,148,4,92,411,310,60],[4,5,6,7,0,1,2,3,8],567,9)]
for A,r,h,f in seeds:
    if (path_count(A),fixed(A,r))!=(h,f):raise ValueError((A,h,f))
    for i in range(len(A)):
        S=1<<i
        for _ in A:
            for j,a in enumerate(A):
                if S>>j&1:S|=a
        if S!=(1<<len(A))-1:raise ValueError(('not strong',A))
    B=wrap(A)
    q=Q_values(A,range(6)); w=Q_values(B,range(1,5))
    if w!=[j*j*q[j]+j*(q[j+1]-q[j-1]) for j in range(1,5)]:
        raise ValueError((q,w))
    q=Q_values(A,range(9)); C=(4*(q[1]+q[2])+q[3])%9
    for k in range(8):
        if k>=3 and q[1]%9!=(C*2**(k-3))%9:raise ValueError((A,q,k))
        q=[(i*i*q[i]+i*(q[(i+1)%9]-q[i-1]))%9 for i in range(9)]
    print('seed',h,f,'wrap evaluation at step 3 mod 9',C)
