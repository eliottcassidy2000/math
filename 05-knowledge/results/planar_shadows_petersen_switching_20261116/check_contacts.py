from itertools import combinations

def induced(T,order):
    return tuple(sum((T[i]>>j&1)<<k for k,j in enumerate(order)) for i in order)

def reachable(a,s,mask):
    r=1<<s
    while True:
        z=r
        for j in range(len(a)):
            if r>>j&1: z|=a[j]&mask
        if z==r: return r
        r=z

def scc(a,verts):
    verts=list(verts)
    mask=sum(1<<i for i in verts)
    grouped={}
    for i in verts:grouped.setdefault(reachable(a,i,mask),[]).append(i)
    return [grouped[z] for z in sorted(grouped,key=int.bit_count,reverse=True)]

def Htable(a):
    n=len(a);F=(1<<n)-1
    path=[[0]*n for _ in range(1<<n)]
    H=[1]+[0]*F
    for s in range(1,1<<n):
        for j in range(n):
            if s>>j&1:
                u=s^(1<<j)
                path[s][j]=sum(path[u][i] for i in range(n) if a[j]>>i&1) if u else 1
        H[s]=sum(path[s])
    return H,path

def profile(a,r):
    h,path=Htable(a);F=(1<<len(a))-1
    E=sum(path[s][r]*h[F^s] for s in range(1<<len(a)) if s.bit_count()%2==0)
    O=sum(path[s][r]*h[F^s] for s in range(1<<len(a)) if s.bit_count()%2==1)
    return h[F],E,O,h[F^(1<<r)]

def pad(a,r,m):
    for _ in range(m):
        a=tuple(v+((i!=r)<<len(a)) for i,v in enumerate(a))+(1<<r,)
    return a

def reduce(a,r):
    blocks=scc(a,(v for v in range(len(a)) if v!=r))
    contacts=[i for i,B in enumerate(blocks) if any(a[r]>>j&1 for j in B)]
    if not contacts:return None
    l,u=contacts[0],contacts[-1]
    left=[B for B in blocks[:l] if len(B)>1]
    active=blocks[l:u+1]
    right=[B for B in blocks[u+1:] if len(B)>1]
    tail=sum(len(B)==1 for B in blocks[u+1:])
    C=induced(a,[v for B in left+active+right for v in B]+[r])
    return C,len(C)-1,tail,len(list(B for B in active if len(B)==1))


from math import comb
for n in range(2,16):
    good=0
    for b in range(1<<(n-2)):
        w=[1]+[(b>>i)&1 for i in range(n-2)]+[0]
        high=sum(w[i] and w[i+1:].count(0)>=3 for i in range(n))
        good+=high<=2
    assert good==n-1+comb(n-1,3)
print("Exact planar Ferrers word counts verified through length 15.")

def transitive(n):
    return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))
def blocks(A,B):
    # root last; it beats A, all A beat B, B beats root
    p,q=len(A),len(B);r=p+q
    return tuple(x+(((1<<q)-1)<<p) for x in A)+tuple((x<<p)+(1<<r) for x in B)+((1<<p)-1,)

def weights(a,r):
    D=induced(a,[i for i in range(len(a)) if i!=r])
    hd,_=Htable(D)
    n=len(D);F=(1<<n)-1
    Q=[u for u in range(n) if a[r]>>u&1]
    L=[v for v in range(n) if not(a[r]>>v&1)]
    W={}
    for u in Q:
        ways=[[0]*n for _ in range(1<<n)]
        ways[1<<u][u]=1
        for S in range(1<<u,1<<n):
            if not S>>u&1:continue
            for j,x in enumerate(ways[S]):
                if not x:continue
                for k in range(n):
                    if not(S>>k&1) and D[j]>>k&1:ways[S|1<<k][k]+=x
        for v in L:
            W[u,v]=sum(ways[S][v]*hd[F^S] for S in range(1<<n) if not S.bit_count()%2)
    return W

T=(2,4,1)
for q in range(1,8):
    A,B=blocks(T,transitive(q)),blocks(transitive(q),T)
    h=9*(1<<q)-3
    r=len(A)-1
    x=profile(A,r);y=profile(B,r)
    assert x==(h,9*(1<<(q-1)),9*(1<<(q-1))+3,3)
    assert y==(h,6*((1<<q)-1),6*((1<<q)-1)+3,3)
    for C in [A,B]:
        W=weights(C,r)
        assert len(W)==3*q and all(W.values()) and 2*sum(W.values())==h-3
    assert profile(pad(A,r,1),r)[0]==h+2*x[1]
    assert profile(pad(B,r,1),r)[0]==h+2*y[1]
print("Literal weighted cut, normalized readiness, and distinguishing extensions verified.")
