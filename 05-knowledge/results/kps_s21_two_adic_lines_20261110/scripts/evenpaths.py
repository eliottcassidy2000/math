"""evenpaths.py - kind-pasteur-2026-11-10-S21. Check of the path-deletion forms (h(p) = H(B - V(p)), H(empty) = 1):
   W_ef = Q(e,f) + Q(f,e) = 2 * sum over directed paths p between e and f (either direction) with an EVEN number of vertices of h(p),
   U_B(y) = sum over directed paths p starting at y (any number of vertices >= 1) of h(p),
on every labelled tournament of order <= 5 and 40 random tournaments of orders 6-7. Exits 1 on failure."""
import itertools, random, sys
def H(A, keep):
    k=len(keep)
    if k==0: return 1
    idx={v:i for i,v in enumerate(keep)}
    out=[sum(1<<idx[u] for u in keep if A[v][u]) for v in keep]
    dp=[[0]*k for _ in range(1<<k)]
    for i in range(k): dp[1<<i][i]=1
    for m in range(1,1<<k):
        for i in range(k):
            c=dp[m][i]
            if c:
                nx=out[i]&~m
                while nx:
                    u=(nx&-nx).bit_length()-1; nx&=nx-1; dp[m|1<<u][u]+=c
    return sum(dp[(1<<k)-1])
def Q(A,e,f):
    n=len(A); tot=0; others=[x for x in range(n) if x not in (e,f)]
    for r in range(len(others)+1):
        for X in itertools.combinations(others,r):
            S1=list(X)+[e]; S2=[x for x in others if x not in X]+[f]
            # paths of A[S1] ending at e, paths of A[S2] starting at f
            c1=sum(1 for p in itertools.permutations(X) if all(A[a][b] for a,b in zip(p+(e,),(p+(e,))[1:])))
            rest=[x for x in others if x not in X]
            c2=sum(1 for p in itertools.permutations(rest) if all(A[a][b] for a,b in zip((f,)+p,((f,)+p)[1:])))
            tot+=c1*c2
    return tot
def paths_from(A,y):
    n=len(A); res=[]
    def rec(path):
        res.append(tuple(path))
        for u in range(n):
            if u not in path and A[path[-1]][u]: rec(path+[u])
    rec([y]); return res
def check(A):
    n=len(A); allv=list(range(n)); ok=True
    P={y:paths_from(A,y) for y in allv}
    h={}
    def hh(p):
        key=frozenset(p)
        if key not in h: h[key]=H(A,[x for x in allv if x not in key])
        return h[key]
    for y in allv:
        U=0
        for r in range(n):
            for X in itertools.combinations([x for x in allv if x!=y],r):
                rest=[x for x in allv if x!=y and x not in X]
                cnt=sum(1 for p in itertools.permutations(rest) if all(A[a][b] for a,b in zip((y,)+p,((y,)+p)[1:])))
                if cnt: U+=cnt*H(A,list(X))
        if U!=sum(hh(p) for p in P[y]): ok=False
    for e in allv:
        for f in allv:
            if e<f:
                W=Q(A,e,f)+Q(A,f,e)
                ev=sum(hh(p) for p in P[e] if p[-1]==f and len(p)%2==0)+sum(hh(p) for p in P[f] if p[-1]==e and len(p)%2==0)
                if W!=2*ev: ok=False
    return ok
ok=True; cnt=0
for n in range(2,6):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in range(1<<len(pairs)):
        A=[[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if bits>>k&1: A[i][j]=1
            else: A[j][i]=1
        ok&=check(A); cnt+=1
rng=random.Random(7)
for n in (6,7):
    for _ in range(20):
        A=[[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1,n):
                if rng.random()<0.5: A[i][j]=1
                else: A[j][i]=1
        ok&=check(A); cnt+=1
print(("PASS" if ok else "FAIL"),"even-path form of W_ef and path form of U_B(y) on %d tournaments"%cnt); sys.exit(0 if ok else 1)
