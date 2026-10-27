"""exceptional_family.py -- mac-mini-2026-10-27-potentials
The family E(M; Z) of THM-4811 (C2): doubled copies of D = C3[K1, M, Z] (K1/C3 stacks M, Z) with all arcs D -> D* except
one reversed arc.  For each member: lambda, number of minimum cuts, H, descent words and black descent words; then
two-move black certificates (a black move keeping lambda = 1 to a class with a black descent)."""
import sys, itertools
sys.setrecursionlimit(100000)
def stack(blocks):
    # blocks: list of 'K' (K1) or 'C' (C3), earlier blocks dominate later ones
    verts=[]; A={}
    idx=0; groups=[]
    for b in blocks:
        g=list(range(idx, idx+(1 if b=='K' else 3))); idx+=len(g); groups.append(g)
    n=idx; M=[[0]*n for _ in range(n)]
    for gi,g in enumerate(groups):
        if len(g)==3:
            M[g[0]][g[1]]=M[g[1]][g[2]]=M[g[2]][g[0]]=1
        for gj in range(gi+1,len(groups)):
            for u in g:
                for v in groups[gj]: M[u][v]=1
    return M
def C3sub(M):
    # D = C3[K1(x0), M, K1(z)]: x0 -> M -> z -> x0 ; x0 = index 0, M = 1..m, z = m+1
    m=len(M); n=m+2; D=[[0]*n for _ in range(n)]
    for i in range(m):
        for j in range(m): D[1+i][1+j]=M[i][j]
        D[0][1+i]=1; D[1+i][m+1]=1
    D[m+1][0]=1
    return D
def doubled(D):
    m=len(D); n=2*m; T=[[0]*n for _ in range(n)]
    for i in range(m):
        for j in range(m):
            T[i][j]=D[i][j]
            T[m+i][m+j]=D[j][i]    # converse copy: vertex v* = m+v
    for i in range(m):
        for j in range(m): T[i][m+j]=1
    # reverse the single arc x0 -> y0 (x0 = 0, y0 = m+0): y0 -> x0
    T[0][m]=0; T[m][0]=1
    return T
def analyse(T, name, maxpaths=10**7):
    n=len(T); full=(1<<n)-1
    s=[sum(r) for r in T]
    best=99; cuts=[]
    for X in range(1,full):
        f=0
        for x in range(n):
            if X>>x&1:
                for y in range(n):
                    if not X>>y&1 and T[x][y]: f+=1
            if f>best: break
        if f<best: best=f; cuts=[]
        if f==best: cuts.append(X)
    # strong?
    out=[sum(1<<j for j in range(n) if T[i][j]) for i in range(n)]
    sep=[[False]*n for _ in range(n)]
    for X in cuts:
        for a in range(n):
            if X>>a&1:
                for b in range(n):
                    if not X>>b&1: sep[a][b]=True
    cand=0; blackc=0; H=0
    def isblack(Q):
        r=[0]*n
        for i in range(n): r[Q[i]]=Q[n-1-i]
        for u in range(n):
            o=out[u]
            while o:
                v=(o&-o).bit_length()-1; o&=o-1
                if not (out[r[v]]>>r[u])&1: return True
        return False
    stack_=[]
    def dfs(p,used):
        nonlocal cand,blackc,H
        if len(p)==n:
            H+=1
            if sep[p[0]][p[-1]]:
                cand+=1
                if isblack(p): blackc+=1
            return
        o=out[p[-1]]&~used
        while o:
            w=(o&-o).bit_length()-1; o&=o-1
            p.append(w); dfs(p,used|1<<w); p.pop()
    for v in range(n):
        # only starts inside some min cut matter for candidates, but count H fully for small n
        dfs([v],1<<v)
    print(f"{name}: n={n} lambda={best} #mincuts={len(cuts)} H={H} descent candidates={cand} black candidates={blackc} scores={sorted(s)}")
for blocks in (['C'], ['C','K'], ['K','C'], ['C','C'], ['C','K','K'], ['C','K','C']):
    M=stack(blocks); D=C3sub(M); T=doubled(D)
    if len(T)<=16:
        analyse(T, 'E('+'=>'.join('C3' if b=='C' else 'K1' for b in blocks)+')')

print("---- wider family D = C3[K1, M, Z] ----")
def C3sub2(M,Z):
    m=len(M); z=len(Z); n=1+m+z; D=[[0]*n for _ in range(n)]
    for i in range(m):
        for j in range(m): D[1+i][1+j]=M[i][j]
    for i in range(z):
        for j in range(z): D[1+m+i][1+m+j]=Z[i][j]
    for i in range(m):
        D[0][1+i]=1
        for j in range(z): D[1+i][1+m+j]=1
    for j in range(z): D[1+m+j][0]=1
    return D
for Mb,Zb in ((['C'],['C']), (['C'],['K','K']), (['C','K'],['C']), (['C'],['C','K'])):
    T=doubled(C3sub2(stack(Mb),stack(Zb)))
    if len(T)<=16: analyse(T, 'E('+'=>'.join('C3' if b=='C' else 'K1' for b in Mb)+' ; '+'=>'.join('C3' if b=='C' else 'K1' for b in Zb)+')')

print("---- black distance 2 certificates for E(M;Z) ----")
import random
def tools(T):
    n=len(T); full=(1<<n)-1
    out=[sum(1<<j for j in range(n) if T[i][j]) for i in range(n)]
    return n,full,out
def lam_cuts_bits(n,full,out):
    best=99; cuts=[]
    for X in range(1,full):
        f=0
        x=X
        while x:
            v=(x&-x).bit_length()-1; x&=x-1
            f+=bin(out[v]&~X&full).count('1')
            if f>best: break
        if f<best: best=f; cuts=[]
        if f==best: cuts.append(X)
    return best,cuts
def paths_iter(n,out,starts):
    for v in starts:
        st=[(v,1<<v,[v])]
        # iterative DFS
        stack_=[[v]]
        used=1<<v
        def rec(p,used):
            if len(p)==n:
                yield tuple(p); return
            o=out[p[-1]]&~used
            while o:
                w=(o&-o).bit_length()-1; o&=o-1
                p.append(w)
                yield from rec(p,used|1<<w)
                p.pop()
        yield from rec([v],1<<v)
def isblack_bits(n,out,Q):
    r=[0]*n
    for i in range(n): r[Q[i]]=Q[n-1-i]
    for u in range(n):
        o=out[u]
        while o:
            v=(o&-o).bit_length()-1; o&=o-1
            if not (out[r[v]]>>r[u])&1: return True
    return False
def kmove_bits(n,out,Q):
    full=(1<<n)-1
    new=[(~out[u])&full&~(1<<u) for u in range(n)]
    for u,v in zip(Q,Q[1:]):
        new[u]|=1<<v; new[v]&=~(1<<u)
    return new
def has_black_descent_bits(n,out):
    full=(1<<n)-1
    l,cuts=lam_cuts_bits(n,full,out)
    if l==0: return None
    sepstarts=set()
    for X in cuts:
        for a in range(n):
            if X>>a&1: sepstarts.add(a)
    for Q in paths_iter(n,out,sorted(sepstarts)):
        if any((X>>Q[0]&1) and not (X>>Q[-1]&1) for X in cuts) and isblack_bits(n,out,Q): return True
    return False
rng=random.Random(1)
for Mb,Zb in ((['C'],['K']), (['C'],['K','K']), (['C','K'],['K']), (['C'],['C']), (['C'],['K','K','K']), (['C','K','K'],['K'])):
    T=doubled(C3sub2(stack(Mb),stack(Zb)))
    n,full,out=tools(T)
    allQ=list(paths_iter(n,out,range(n)))
    rng.shuffle(allQ)
    found=None; tried=0
    for Q in allQ:
        if not isblack_bits(n,out,Q): continue
        new=kmove_bits(n,out,Q)
        l,_=lam_cuts_bits(n,full,new)
        if l!=1: continue
        tried+=1
        if has_black_descent_bits(n,new):
            found=Q; break
        if tried>200: break
    name='E('+'=>'.join('C3' if b=='C' else 'K1' for b in Mb)+' ; '+'=>'.join('C3' if b=='C' else 'K1' for b in Zb)+')'
    print(f"{name}: n={n}: black first move to a lambda=1 class with a black descent: {'FOUND '+str(found) if found else 'not found in '+str(tried)+' tries'}")
