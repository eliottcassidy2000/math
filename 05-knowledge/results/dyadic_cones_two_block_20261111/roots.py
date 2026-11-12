"""Unit-offset roots (THM-5126): for a tournament B and a vertex set O, c0 = H(B) - sum_{pairs of O} W_ef,
M = sum_{y in O} U_B(y); unit offset iff c0 = 1; then 1 + M 2^t is a path count for every t >= 1.
Explore the family B = transitive tournament on 0..m-1 (i -> j for i < j) with a set of back arcs."""
import sys, itertools
def make(m, back):
    out=[0]*m
    for i in range(m):
        for j in range(i+1,m):
            if (j,i) in back: out[j]|=1<<i
            else: out[i]|=1<<j
    return out
def tables(out):
    """E[S][v] paths on S ending at v; St[S][v] paths on S starting at v; Hs[S]"""
    n=len(out); N=1<<n
    E=[[0]*n for _ in range(N)]
    for v in range(n): E[1<<v][v]=1
    for S in range(1,N):
        row=E[S]
        for v in range(n):
            c=row[v]
            if c:
                t=out[v]&~S
                while t:
                    b=t&-t; E[S|b][b.bit_length()-1]+=c; t^=b
    inn=[0]*n
    for v in range(n):
        for w in range(n):
            if (out[w]>>v)&1: inn[v]|=1<<w
    St=[[0]*n for _ in range(N)]
    for v in range(n): St[1<<v][v]=1
    for S in range(1,N):
        row=St[S]
        for v in range(n):
            c=row[v]
            if c:
                t=inn[v]&~S
                while t:
                    b=t&-t; St[S|b][b.bit_length()-1]+=c; t^=b
    Hs=[sum(r) for r in E]; Hs[0]=1
    return E,St,Hs
def root_data(out):
    n=len(out); N=1<<n; full=N-1
    E,St,Hs=tables(out)
    H=Hs[full]
    U=[0]*n; W=[[0]*n for _ in range(n)]
    for S in range(0,N):          # S = vertex set of P1 (possibly empty), complement = P2 (nonempty)
        C=full^S
        if C==0: continue
        hs=Hs[S]
        for y in range(n):
            if (C>>y)&1 and St[C][y]: U[y]+=hs*St[C][y]
        if S:
            for e in range(n):
                if (S>>e)&1 and E[S][e]:
                    te=E[S][e]
                    for f in range(n):
                        if (C>>f)&1 and St[C][f]: W[e][f]+=te*St[C][f]     # Q(e,f)
    Wsym=[[W[e][f]+W[f][e] for f in range(n)] for e in range(n)]
    return H,U,Wsym
def unit_roots(out, maxO=4):
    n=len(out); H,U,W=root_data(out); res=[]
    for k in range(1,maxO+1):
        for O in itertools.combinations(range(n),k):
            w=sum(W[e][f] for e,f in itertools.combinations(O,2))
            if H-w==1:
                M=sum(U[y] for y in O); res.append((O,M))
    return H,res
def odd(x):
    j=0
    while x%2==0: x//=2; j+=1
    return x,j
if __name__=="__main__":
    m=int(sys.argv[1]); nb=int(sys.argv[2])
    pairs=[(j,i) for i in range(m) for j in range(i+1,m)]
    seen={}
    for back in itertools.combinations(pairs,nb):
        out=make(m,set(back)); H,res=unit_roots(out,3)
        for O,M in res:
            c,j=odd(M)
            seen.setdefault(c,[]).append((j,H,back,O))
    for c in sorted(seen):
        best=min(seen[c]); print("c=%d  least j=%d  (H(B)=%d, back arcs %s, O=%s)  [%d roots]" % (c,best[0],best[1],best[2],best[3],len(seen[c])))
