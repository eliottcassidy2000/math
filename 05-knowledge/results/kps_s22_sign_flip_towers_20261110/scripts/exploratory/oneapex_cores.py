"""oneapex_cores.py - direct enumeration: one-apex tournaments Q = (apex x; binary word w on a chain z_1..z_n, x -> z_i iff w_i = 1)
and vertex pairs/triples O with W(O) = H(Q) - 1 (unit sets). Reduced words only (w_1 = 1, w_n = 0: Q strong).
Prints word, O (x or chain positions), H(Q), M = U_Q(O), odd part c."""
import sys, itertools
from functools import lru_cache
def build(w):
    n=len(w); A=[[0]*(n+1) for _ in range(n+1)]   # chain 0..n-1, apex n
    for i in range(n):
        for j in range(i+1,n): A[i][j]=1
        if w[i]: A[n][i]=1
        else: A[i][n]=1
    return A
def tables(A):
    n=len(A); full=1<<n
    out=[sum(1<<u for u in range(n) if A[v][u]) for v in range(n)]
    inn=[sum(1<<u for u in range(n) if A[u][v]) for v in range(n)]
    en=[[0]*n for _ in range(full)]; st=[[0]*n for _ in range(full)]
    for v in range(n): en[1<<v][v]=1; st[1<<v][v]=1
    for m in range(1,full):
        for v in range(n):
            c=en[m][v]
            if c:
                nx=out[v]&~m
                while nx:
                    u=(nx&-nx).bit_length()-1; nx&=nx-1; en[m|1<<u][u]+=c
            c=st[m][v]
            if c:
                nx=inn[v]&~m
                while nx:
                    u=(nx&-nx).bit_length()-1; nx&=nx-1; st[m|1<<u][u]+=c
    return en,st
NMAX=int(sys.argv[1]) if len(sys.argv)>1 else 8
rows=[]
for n in range(2,NMAX+1):
    for w in itertools.product((0,1),repeat=n):
        if w[0]!=1 or w[-1]!=0: continue
        A=build(w); N=n+1; full=(1<<N)-1
        en,st=tables(A); hs=[sum(e) for e in en]; hs[0]=1; H=hs[full]
        U=[sum(st[X][y]*hs[full^X] for X in range(1,full+1) if X>>y&1) for y in range(N)]
        W=[[0]*N for _ in range(N)]
        for X in range(1,full):
            Y=full^X
            for a in range(N):
                if X>>a&1 and en[X][a]:
                    for c in range(N):
                        if Y>>c&1 and st[Y][c]:
                            q=en[X][a]*st[Y][c]; W[a][c]+=q; W[c][a]+=q
        for k in (2,3):
            for O in itertools.combinations(range(N),k):
                if sum(W[a][b] for i,a in enumerate(O) for b in O[i+1:])==H-1:
                    M=sum(U[y] for y in O); c=M
                    while c%2==0: c//=2
                    name=tuple("x" if y==n else y+1 for y in O)
                    rows.append((n,"".join(map(str,w)),name,H,M,c))
for r in rows: print(*r)
print("total",len(rows),"distinct c:",sorted(set(r[5] for r in rows)))
