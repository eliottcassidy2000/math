from collections import Counter
from functools import lru_cache
from itertools import product
from math import comb
from pathlib import Path
import sys

LOW=(2,12,17,37,75,23,47)
HIGH=(2,12,17,37,43,71,31)
J=(1,2,0,0,0,1,1,1,0)

def ends(R):
    n=len(R)
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for i,c in enumerate(dp[mask]):
            if not c:continue
            q=R[i]&~mask
            while q:
                b=q&-q
                q^=b
                dp[mask|b][b.bit_length()-1]+=c
    return dp

def state(R,u):
    F=(1<<len(R))-1
    assert all(((R[i]>>j)^(R[j]>>i))&1 for i in range(len(R)) for j in range(i+1,len(R)))
    h=[sum(a) for a in ends(R)]
    h[0]=1
    start=ends(tuple(F^r^(1<<i) for i,r in enumerate(R)))
    v=[0,0]
    for q in range(1,F+1):v[q.bit_count()%2]+=start[q][u]*h[F^q]
    return h[F],*v

def seam(R,u):
    n=len(R)
    return tuple(r|((i!=u)<<n)for i,r in enumerate(R))+(1<<u,),n

def mult(M,V):
    return tuple(sum(M[3*i+j]*V[j]for j in range(3))for i in range(3))

def skeleton_counts(R,first=0,second=1):
    M=(1<<first)|(1<<second)
    F=(1<<len(R))-1
    @lru_cache(None)
    def f(last,left):
        Q=Counter({(0,0):1}) if not left else Counter()
        options=(F if last<0 else R[last])&(left|M)
        for i in range(len(R)):
            if options>>i&1:
                for (a,b),c in f(i,left if M>>i&1 else left^(1<<i)).items():
                    Q[a+(i==first),b+(i==second)]+=c
        return Q
    return Counter({k:v for k,v in f(-1,F^M).items()if min(k)>0})

C=((67,81,24),(32,86,42),(0,22,22))

def abc(n):return (1, (1<<n)-2, 3**n-3*(1<<n)+3)
def H(n,m):
    a,b=abc(n),abc(m)
    return sum(C[i][j]*a[i]*b[j] for i,j in product(range(3), repeat=2))

def graph(R,n,m):
    p=[0]*n+[1]*m+[2,3,4,5,6]
    return tuple(sum(1<<j for j,b in enumerate(p)
                 if (a==b and i<j) or R[a]>>b&1)
            for i,a in enumerate(p))

def rooted_tensor(R, root):
    if root==1:
        # Vertex 7 represents the last clone of vertex 1.
        Q=tuple(r|((((r>>1)&1) or i==1)<<7) for i,r in enumerate(R))+(R[1],)
        outside=124
        root=7
    else:
        assert root==6
        Q=R
        outside=60
    @lru_cache(None)
    def remainder(l,unused):
        A=Counter({(0,0):1}) if not unused else Counter()
        options=(127 if l<0 else Q[l])&(unused|3)
        for i in range(7):
            if not options>>i&1:continue
            for (a,b),c in remainder(i,unused if i<2 else unused^(1<<i)).items():
                A[a+(i==0),b+(i==1)]+=c
        return A
    @lru_cache(None)
    def prefix(l,used):
        A=Counter()
        for (a,b),c in remainder(-1,outside^used).items():
            A[(1+used.bit_count())%2,0,0,a,b]+=c
        for i in range(7):
            if not (Q[l]&~used&(outside|3))>>i&1:continue
            for (p,a,b,c,d),multiplier in prefix(i,used if i<2 else used|(1<<i)).items():
                A[p,a+(i==0),b+(i==1),c,d]+=multiplier
        return A
    return prefix(root,0)

def root6():
    A,B=rooted_tensor(LOW,6),rooted_tensor(HIGH,6)
    total,signed=Counter(),Counter()
    for p,a,b,c,d in A.keys()|B.keys():
        delta=B[p,a,b,c,d]-A[p,a,b,c,d]
        total[a+c,b+d]+=delta
        for i,j in product(range(a+c+1),range(b+d+1)):
            sign=(-1)**((p+a+b+a+c+b+d-i-j)%2)
            signed[i-a,j-b]+=delta*sign*comb(a+c,i)*comb(b+d,j)
    assert {k:v for k,v in total.items() if min(k)>0 and v}=={(1,1):-2,(2,1):-2,(2,2):-1}
    assert {k:v for k,v in signed.items() if v}=={(0,0):-1}
    matrices=[tuple(map(int,L.split())) for L in
              Path(__file__).with_name("base2_order5_matrices.txt").read_text().splitlines()]
    assert len(matrices)==264
    assert all(b+c>=a for M in matrices for a,b,c in (M[0:3],M[3:6],M[6:9]))
    mixed=lambda v:min(v)<0<max(v)
    first={mult(M,(2,-1,-1)) for M in matrices}
    assert {v for v in first if mixed(v)}=={(-2,1,-1),(0,-1,1)}
    images=lambda v:{w for M in matrices if mixed(w:=mult(M,v))}
    assert not images((-2,1,-1))
    assert images((0,-1,1))=={(-2,1,-1)}

def other_bank():
    certs=(
        ((8,5,9,50,7,23,63),5,4,23),
        ((2,24,3,21,37,15,63),5,2,25),
        ((2,16,3,7,45,79,31),3,2,27),
        ((2,8,3,21,39,79,31),2,0,29),
    )
    for R,a,b,h in certs:
        assert skeleton_counts(R,a,b)=={(1,1):h,(1,2):12,(2,1):12,(2,2):6}
        for n,m in [(1,1),(2,3),(3,2)]:
            P=[i for i in range(len(R)) for _ in range(n if i==a else m if i==b else 1)]
            T=tuple(sum(1<<j for j,v in enumerate(P)
                    if (u==v and i<j) or ((R[u]>>v)&1))
                    for i,u in enumerate(P))
            assert sum(ends(T)[-1])==6*(1<<(n+m))+h-24

def mixed_bank():
    for R,a,b,h in (
        ((32,5,9,19,7,94,31),4,1,69),
        ((16,5,33,7,78,27,47),3,1,71),
        ((16,5,33,7,46,75,31),3,1,73),
    ):
        T=((h,54,12),(54,45,12),(12,12,4))
        assert skeleton_counts(R,a,b)=={(i+1,j+1):T[i][j] for i,j in product(range(3),repeat=2)}
        for n,m in [(1,1),(2,3),(3,2)]:
            P=[i for i in range(len(R)) for _ in range(n if i==a else m if i==b else 1)]
            S=tuple(sum(1<<j for j,v in enumerate(P)
                    if (u==v and i<j) or ((R[u]>>v)&1))
                    for i,u in enumerate(P))
            assert sum(ends(S)[-1])==4*3**(n+m)+9*2**(n+m)+h-72

def solve(N):
    n=1
    while H(n,1)<=N:
        for d in (0,2):
            R=HIGH if d else LOW
            lo=hi=1
            while H(n,hi)+d<N: hi*=2
            while lo<=hi:
                m=(lo+hi)//2
                actual=H(n,m)+d
                if actual==N:
                    yield n,m,d,graph(R,n,m)
                    break
                if actual<N:lo=m+1
                else:hi=m-1
        n+=1

if __name__=="__main__":
    a,b=skeleton_counts(LOW),skeleton_counts(HIGH)
    assert a=={(i+1,j+1):C[i][j] for i,j in product(range(3),repeat=2) if C[i][j]}
    assert b-a=={(1,1):2} and not a-b
    A,B=rooted_tensor(LOW,1),rooted_tensor(HIGH,1)
    assert B-A=={(1,0,0,i,j):2 for i,j in product((0,1),repeat=2)} and not A-B
    root6()
    other_bank()
    mixed_bank()
    for n,m in product(range(1,5),repeat=2):
        r=n+m-1
        G,L=graph(LOW,n,m),graph(HIGH,n,m)
        a,b=state(G,r),state(L,r)
        assert a[0]==H(n,m)
        assert tuple(y-x for x,y in zip(a,b))==(2,0,2)
        a6,b6=state(G,len(G)-1),state(L,len(L)-1)
        t=1+((1<<n)-2)*(1<<(m-1))
        assert tuple(y-x for x,y in zip(a6,b6))==(2,-t,-t)
        assert (n,m,0,G) in set(solve(a[0]))
        assert (n,m,2,L) in set(solve(b[0]))
        for V,T in [(a,G),(b,L)]:
            GJ,j=seam(T,r)
            assert state(GJ,j)==mult(J,V)
    print("PASS complete ordinary and two rooted tensors, 16 literal pairs and 32 direct seam checks.")
    if len(sys.argv)>1:
        print(*list(solve(int(sys.argv[1]))),sep="\n")
