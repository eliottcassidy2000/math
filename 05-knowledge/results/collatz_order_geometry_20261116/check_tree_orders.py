from itertools import product,combinations
from collections import defaultdict,deque
from math import factorial,prod

def extensions(par):
    n=len(par)
    seq=[0]
    def walk(free,used):
        if not free:
            yield tuple(seq);return
        for v in free:
            f=free-{v}|{i for i in range(n) if par[i]==v}
            seq.append(v)
            yield from walk(f,used|{v})
            seq.pop()
    return walk({i for i in range(n) if par[i]==0},{0})

def tournament(par,order):
    n=len(order)
    A=[[0]*n for _ in range(n)]
    for i,j in combinations(range(n),2):
        if par[order[j]]==order[i]:A[j][i]=1
        else:A[i][j]=1
    return A

def direct_tri(A):
    return sum(A[i][j]==A[j][k]==A[k][i] for i,j,k in combinations(range(len(A)),3))

def H(A):
    n=len(A)
    d=[[0]*n for _ in range(1<<n)]
    for i in range(n):d[1<<i][i]=1
    for s in range(1,1<<n):
        for i,u in enumerate(d[s]):
            if not u:continue
            for j,v in enumerate(A[i]):
                if v and not s>>j&1:d[s|1<<j][j]+=u
    return sum(d[-1])

def law(par,order):
    rank={v:i for i,v in enumerate(order)}
    n=len(par)
    degrees=[par.count(v) for v in range(n)]
    return sum(rank[v]-rank[par[v]] for v in range(1,n))-sum(d*(d-1)//2 for d in degrees)-degrees[0]

def minimum(par):
    n=len(par)
    children=[[v for v in range(n) if par[v]==u] for u in range(n)]
    sizes=[1]*n
    b=0
    for v in range(n-1,-1,-1):
        sizes[v]+=sum(sizes[u] for u in children[v])
        b+=sum(min(sizes[x],sizes[y])-1 for x,y in combinations(children[v],2))
    def dfs(v):
        yield v
        for w in sorted(children[v],key=lambda u:sizes[u]):yield from dfs(w)
    return n-1-len(children[0])+b,tuple(dfs(0)),factorial(n)//prod(sizes)

def check():
    trees=orders=0
    first=None
    for n in range(1,8):
        for ps in product(*(range(v) for v in range(1,n))):
            par=(-1,)+ps; trees+=1
            low,dfs,hook=minimum(par)
            ext=list(extensions(par));orders+=len(ext)
            tr={q:direct_tri(tournament(par,q)) for q in ext}
            assert len(ext)==hook and low==min(tr.values())==tr[dfs]
            for q,t in tr.items():
                assert t==law(par,q)
                for i in range(1,n-1):
                    a,b=q[i:i+2];r=q[:i]+(b,a)+q[i+2:]
                    if r not in tr: continue
                    assert tr[r]-t==par.count(b)-par.count(a)
                    if first is None and tr[r]==t:
                        ha,hb=H(tournament(par,q)),H(tournament(par,r))
                        if ha!=hb:first=(par,q,r,t,ha,hb)
    print('all rooted increasing parent maps through order 7:',trees,'linear extensions:',orders)
    print('smallest equal-triangle / different-H adjacent swap:',first)

    F=lambda n:n//2 if n%2==0 else (3*n+1)//2
    basin=[1]
    for y in basin:
        if len(basin)>=100:break
        basin.append(2*y)
        if y%3==2 and y!=2:basin.append((2*y-1)//3)
    for n in (10,20,40,80):
        S=basin[:n]; inv={v:i for i,v in enumerate(S)}
        par=(-1,)+tuple(inv[F(y)] for y in S[1:])
        b=law(par,tuple(range(n))); m,r,_=minimum(par)
        print('shortcut BFS initial',n,'BFS triangles',b,'minimum',m,'order',tuple(S[i] for i in r) if n==10 else '...')
        assert m==direct_tri(tournament(par,r)) and b==direct_tri(tournament(par,tuple(range(n))))
if __name__=='__main__':check()
