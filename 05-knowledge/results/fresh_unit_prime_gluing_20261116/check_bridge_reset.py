"""Independent subset DP for both root policies and the cross-policy witness."""
from itertools import combinations
from random import Random

def adjacency(n, code):
    rows = [0]*n
    for e,(i,j) in enumerate(combinations(range(n),2)):
        if code>>e&1: i,j=j,i
        rows[i] |= 1<<j
    return tuple(rows)

def count(rows, root=0, end=False):
    n=len(rows); f=(1<<n)-1
    if end:
        rows=tuple(f^r^(1<<i) for i,r in enumerate(rows))
    dp=[[0]*n for _ in range(f+1)]
    hs=[1]+[0]*f
    for S in range(1,f+1):
        for i in range(n):
            if (S>>i)&1:
                R=S^(1<<i)
                dp[S][i]=sum(dp[R][j] for j in range(n) if rows[i]>>j&1) if R else 1
        hs[S]=sum(dp[S])
    eo=[0,0]
    for S in range(1,f+1):
        eo[S.bit_count()%2]+=dp[S][root]*hs[f^S]
    return hs[f],hs[f^(1<<root)],*eo, sum(hs[S]*hs[f^S] for S in range(f+1))

def drop(a,t):
    vertices=[i for i in range(len(a)) if i!=t]
    return tuple(sum(1<<j for j,v in enumerate(vertices) if a[i]>>v&1) for i in vertices)

def join(a,x,b,y):
    n=len(a);m=len(b);T=((1<<m)-1)<<n
    return tuple([r|(T^(1<<(n+y) if i==x else 0)) for i,r in enumerate(a)] +
                 [(r<<n)|((1<<x) if j==y else 0) for j,r in enumerate(b)])

def verify(A,x,B,y):
    H,g,E,O,_=count(A,x)
    h,hminus,Es,Os,p=count(B,y)
    _,_,Ee,Oe,_=count(B,y,True)
    W=join(A,x,B,y)
    first=h*H+2*(Oe*E+Ee*O)
    assert E+O>=H
    assert p//2>=Os and p//2>=Es+h
    assert count(W,x)[:4]==(first,h*g,p//2*(E+O)-Os*g,p//2*(E+O)-Es*g)
    if len(B)>1:
        q=count(drop(B,y))[-1]//2
        assert count(W,len(A)+y)[:4]==(first,hminus*H,Es*H+q*(E+O),Os*H+q*(E+O))
    else:
        assert count(W,len(A)+y)[:4]==(first,H,O,H+E)

tri=(2,4,1)
f=tuple(row|24 for row in tri)+(16,0)
r=(30,28,8,16,4)
assert count(join((0,),0,f,0),0)[:4]==(5,3,18,18)
assert count(join((0,),0,r,1),2)[:4]==(9,3,18,18)
rng=Random(99)
for j in range(250):
    n,m=rng.randrange(1,6),rng.randrange(1,6)
    a=adjacency(n,rng.randrange(1<<(n*(n-1)//2)))
    b=adjacency(m,rng.randrange(1<<(m*(m-1)//2)))
    verify(a,rng.randrange(n),b,rng.randrange(m))
print("Both port-transfer laws replayed on 250 cases; synchronized +4 witness replayed.")
