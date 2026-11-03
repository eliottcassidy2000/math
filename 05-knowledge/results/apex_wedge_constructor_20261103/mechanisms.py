"""Literal two-apex sources, exact gap evaluator, wedge and native reset."""
from collections import defaultdict
from itertools import product

def need(ok,message):
    if not ok:raise RuntimeError(message)

def bits(n,length=None):
    need(n>=0,'nonnegative chart')
    return tuple(map(int,format(n,'0'+str(length or max(1,n.bit_length()))+'b')))

def source(w,u,forward=True):
    need(len(w)==len(u) and all(x in (0,1) for x in w+u),'Boolean words')
    m=len(w);a=[0]*(m+2)
    for i in range(m):
        for j in range(i+1,m):a[i]|=1<<j
    for apex,word in ((m,w),(m+1,u)):
        for i,b in enumerate(word):
            if b:a[apex]|=1<<i
            else:a[i]|=1<<apex
    a[m if forward else m+1]|=1<<(m+1 if forward else m)
    return tuple(a)

def hp(rows):
    n=len(rows);dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for i,v in enumerate(dp[mask]):
            if not v:continue
            options=rows[i]&~mask
            while options:
                b=options&-options;options-=b;j=b.bit_length()-1
                dp[mask|b][j]+=v
    return sum(dp[-1])

STATES=tuple(product(range(3),repeat=3))
INDEX={s:i for i,s in enumerate(STATES)}
def step(state,first,second):
    # Empty, nonempty-failing-right, nonempty-passing-right.
    L=(1,first,second);R=(1-first,1-second,1)
    for g in range(3):
        if state[g]==0 and not L[g]:continue
        out=list(state);out[g]=2 if R[g] else 1
        yield tuple(out)

def ordered_gap(w,u,edge):
    count={(0,0,0):1}
    for a,b in zip(w,u):
        nxt=defaultdict(int)
        for s,v in count.items():
            for t in step(s,a,b):nxt[t]+=v
        count=nxt
    return sum(v for s,v in count.items() if all(x!=1 for x in s) and (s[1]!=0 or edge))

def gap_h(w,u,forward=True):
    return ordered_gap(w,u,forward)+ordered_gap(u,w,not forward)

def kernel(m):
    return tuple(tuple(1 if j==i+1 else (1<<(i-j-1)) if i>j else 0 for j in range(m)) for i in range(m))

def adjacent_count(w,u):
    if not w:return 1
    K=kernel(len(w))
    return 1-w[-1]+u[0]+sum((1-w[i])*K[i][j]*u[j] for i in range(len(w)) for j in range(len(w)))

def jump(w,u):return adjacent_count(w,u)-adjacent_count(u,w)

def wedge(m):
    if not m:return ((0,),)
    K=kernel(m);c=[sum(K[i][j] for i in range(m))+(j==0)+(j==m-1) for j in range(m)]
    A=[[0]*(m+1) for _ in range(m+1)]
    for j in range(m):A[0][j+1]=c[j];A[j+1][0]=-c[j]
    for i in range(m):
        for j in range(m):A[i+1][j+1]=K[j][i]-K[i][j]
    need(all(x%2==0 for row in A for x in row),'integral alternating half-form')
    return tuple(tuple(x//2 for x in row) for row in A)

def pair(B,x,y):return sum(a*B[i][j]*b for i,a in enumerate(x) for j,b in enumerate(y))

def chart(P,R):
    need(P>=1 and R>=1,'positive charts retain strong endpoints')
    p=bits(P);s=tuple(1-b for b in reversed(bits(R)))
    return p+(0,1)+s,p+(1,0)+s

def build_jump(d):
    need(d>=1,'positive requested half-jump')
    return chart(1,2*d+1)

def unit_family(r):
    need(r>=2,'unit family range')
    return chart(1<<r,(1<<r)+2)

def native(rows,path):
    n=len(rows)
    need(sorted(path)==list(range(n)) and all(rows[a]>>b&1 for a,b in zip(path,path[1:])),'actual Hamilton path guard')
    U=[((1<<n)-1)^r^(1<<i) for i,r in enumerate(rows)]
    for a,b in zip(path,path[1:]):U[a]|=1<<b;U[b]&=~(1<<a)
    return tuple(U)

def reset(w,u):
    need(w[0]==u[0]==1,'common first out-neighbor')
    T=source(w,u);n=len(T);a=n-2;b=n-1
    # One temporary vertex beats only a; every old name, including a,b, stays.
    X=tuple(r|(0 if i==a else 1<<n) for i,r in enumerate(T))+(1<<a,)
    P=(a,b)+tuple(range(n-2))+(n,)
    Q=(a,n,b)+tuple(range(n-2))
    Y=native(X,P);Z=native(Y,Q)
    endpoint=tuple(r&((1<<n)-1) for r in Z[:n])
    need(endpoint==source(w,u,False),'exact old-edge reset')
    return {'ports':[a,b],'source':list(T),'target':list(endpoint),'birth_out':1<<a,'paths':[P,Q],
            'vertex_clock':[n,n+1,n+1,n+1,n],'native_moves':2,'temporary_births':1,'temporary_deletions':1}
