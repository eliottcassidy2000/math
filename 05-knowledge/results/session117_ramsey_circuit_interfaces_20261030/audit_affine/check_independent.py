"""Independent literal-set affine geometry audit; no producer imports."""
import json
from itertools import product

def require(b, msg):
    if not b: raise RuntimeError(msg)
def dot(a,b): return (a&b).bit_count()%2

def subspaces(n):
    spaces={frozenset([0])}
    frontier=list(spaces)
    while frontier:
        U=frontier.pop()
        for x in range(1<<n):
            if x not in U:
                W=U|frozenset(y^x for y in U)
                if W not in spaces:
                    spaces.add(W); frontier.append(W)
    return spaces

def gaussian(n,d):
    if d<0 or d>n: return 0
    a=b=1
    for i in range(d): a*= (1<<(n-i))-1; b*=(1<<(d-i))-1
    return a//b

checks=0; receipts=[]
for n in range(5):
    N=1<<n
    flats={frozenset(x^a for x in U) for U in subspaces(n) for a in range(N)}
    for d in range(n+1):
        require(sum(len(A)==1<<d for A in flats)==(1<<(n-d))*gaussian(n,d),('flat enumeration',n,d))
    for r,v,e in product(range(N),range(N),range(2)):
        c=e^dot(r,v)
        f=[e^dot(r,x)^(x==v) for x in range(N)]
        counts=[[0,0] for _ in range(n+1)]
        for A in flats:
            d=len(A).bit_length()-1
            a=next(iter(A)); U={x^a for x in A}
            actual=len({f[x] for x in A})==1
            predicted=(d==0 or (v not in A and all(dot(r,u)==0 for u in U)) or (v in A and d==1 and any(dot(r,u) for u in U)))
            require(actual==predicted,('classification',n,r,v,e,A))
            if actual: counts[d][f[a]]+=1
            checks+=1
        require(sum(counts[0])==N,('singletons',n,r,v,e))
        if r:
            require(counts[0][c]==N//2-1 and counts[0][c^1]==N//2+1,('singleton colors',n,r,v,e))
        else:
            require(counts[0][e]==N-1 and counts[0][e^1]==1,('constant singleton colors',n,r,v,e))
        for d in range(1,n+1):
            expected=((1<<(n-d))-1)*gaussian(n-(r!=0),d)+( (1<<(n-1)) if r and d==1 else 0)
            require(sum(counts[d])==expected,('total',n,r,v,e,d))
            if not r:
                require(counts[d][e]==expected and counts[d][e^1]==0,('zero-phase colors',n,r,v,e,d))
            elif d<n:
                g=gaussian(n-1,d); h=1<<(n-d-1)
                require(counts[d][c]==g*(h-1) and counts[d][c^1]==g*h+((1<<(n-1)) if d==1 else 0),('colors',n,r,v,e,d))
            else:
                require(sum(counts[d])==(1 if n==1 else 0),('full space exception',n,r,v,e))
        W=[sum((-1)**(f[x]^dot(s,x)) for x in range(N)) for s in range(N)]
        require(W==[((-1)**e)*N*(s==r)-2*((-1)**(e^dot(r^s,v))) for s in range(N)],('Walsh',n,r,v,e))
    # Every invertible linear map, exhaustively in dimensions at most 3.
    maps=0
    if n<=3:
        for columns in product(range(N),repeat=n):
            images=[0]*N
            for x in range(N):
                for i in range(n):
                    if x>>i&1: images[x]^=columns[i]
            if len(set(images))!=N: continue
            maps+=1
            for r in range(1,N):
                rp=sum(dot(r,columns[i])<<i for i in range(n))
                require(rp!=0,('dual action',n,r,columns))
                for x in range(N):
                    require((dot(r,images[x])^(images[x]==0))==(dot(rp,x)^(x==0)),('orbit action',n,r,columns,x))
    receipts.append({'n':n,'flats':len(flats),'parameter_triples':N*N*2,'GL_maps_checked':maps})
print(json.dumps({'status':'PASS','literal_flat_classifications':checks,'dimensions':receipts},sort_keys=True,indent=2))
