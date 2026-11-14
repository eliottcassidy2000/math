from functools import lru_cache
from itertools import product
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent /
                      "session328_fermat_gluing_and_phase_repair_20261114"))
from profile_tools import strong

def smooth(N):
    b=0; A=[]
    while N:
        while N%3==0:
            N//=3;b+=1
        a=N.bit_length()-1
        if pow(2,a,3)!=N%3:a-=1
        assert a>=0 and (1<<a)<=N<1<<(a+2)
        A.append((a,b))
        N=(N-(1<<a))//3;b+=1
    return A

def sum_smooth(A):return sum((1<<a)*3**b for a,b in A)
for n in range(1,10001):
    A=smooth(n)
    assert sum_smooth(A)==n
    assert all(a>c and b<d for (a,b),(c,d) in zip(A,A[1:]))
print('34/45/69 exact 2-3 expansions',[(n,[(1<<a)*3**b for a,b in smooth(n)]) for n in (34,45,69)])

def marked_three(A):
    s=(len(A)-1).bit_length();p=s+A[0][0];bmax=A[-1][1]
    inner=p+1;tail0=p+bmax+2;tail1=tail0+1;size=tail1+1
    T=[(1<<size)-(1<<(v+1)) for v in range(size)]
    T[0]^=1<<tail0;T[tail0]|=1
    T[inner]^=1<<tail1;T[tail1]|=1<<inner
    assert strong(T)
    ans=set()
    for ix,(a,b) in enumerate(A):
        for U in product(range(2),repeat=a):
            for V in product(range(3),repeat=b):
                R=[[tail0],[0,tail1],[inner]]
                for j in range(s):
                    R[ix>>j&1].append(j+1)
                for j in range(A[0][0]):
                    R[U[j] if j<a else 0].append(j+s+1)
                for j in range(bmax):
                    R[V[j] if j<b else 0].append(inner+j+1)
                w=tuple(v for r in R for v in sorted(r))
                assert len(set(w))==size and all(T[x]>>y&1 for x,y in zip(w,w[1:]))
                assert [(x,y) for x,y in zip(w,w[1:]) if x>y]==[(tail0,0),(tail1,inner)]
                assert w not in ans
                ans.add(w)
    return len(ans)

def crossing(A):
    DAG=[set()]
    def new(inp):
        j=len(DAG);DAG.append(set())
        for i in inp:DAG[i].add(j)
        return j
    def mult(v,k):return new([new([v]) for _ in range(k)])
    B=[0]
    for _ in range(A[0][0]): B.append(mult(B[-1],2))
    S=new([B[A[-1][0]]])
    for (a,b),(_,d) in zip(reversed(A[:-1]),reversed(A[1:])):
        for _ in range(d-b):S=mult(S,3)
        S=new([S,B[a]])
    for _ in range(A[0][1]):S=mult(S,3)
    assert S==len(DAG)-1
    d=len(DAG)
    M=[row.copy() for row in DAG]
    for i in range(1,d):M[i].add(i)
    M[S].add(0)
    tax=sum(map(len,DAG))-d+1
    for _ in range(tax):M.append({len(M)})
    assert sum(map(len,M))==2*len(M)-1
    return M

def permanent(M):
    m=[sum(1<<j for j in x) for x in M]
    @lru_cache(None)
    def rec(r,c):
        if not r:return 1
        RR=[i for i in range(len(M)) if r>>i&1]
        i=min(RR,key=lambda i:(m[i]&c).bit_count())
        W=m[i]&c
        ans=0
        while W:
            d=W&-W;W-=d
            ans+=rec(r^(1<<i),c^d)
        return ans
    z=(1<<len(M))-1
    return rec(z,z)

def marked_matching(M):
    d=len(M);n=3*d+1;T=[0]*n
    def setarc(a,b):T[a]|=1<<b
    for a in range(n):
        for b in range(a+1,n):
            if a<=d and b>2*d:setarc(b,a)
            elif d<a<=2*d<b and b-(2*d+1) not in M[a-(d+1)]:setarc(b,a)
            else:setarc(a,b)
    assert strong(T)
    def walk(i,U,P):
        if i==d:
            assert len(set(P))==n
            assert all(T[a]>>b&1 for a,b in zip(P,P[1:]))
            return 1
        o=d+1+i
        return sum(walk(i+1,U|{j},P+[o,2*d+1+j,i+1]) for j in M[i]-U)
    return walk(0,set(),[0]),n

for n in list(range(1,121))+[243,729,2026]:
    A=smooth(n)
    assert marked_three(A)==n
    M=crossing(A)
    assert permanent(M)==n
    if n<=21 or n in (34,45,69):
        R,K=marked_matching(M)
        assert R==n
        if n in (34,45,69): print('literal matching marked n/host order',n,K)
print('absolute smooth antichain and both marked compilers: PASS')
