"""mean_pi_identity.py - kind-pasteur-2026-11-10-S21. Checks, for n = 3, 4 (and 5 with argument 5):
   sum over all labelled tournaments T on n vertices of pi_T(x)  =  2^C(n,2) * sum over permutations sigma of odd order
   of 2^(psi(sigma) - moved(sigma)) x^cyc(sigma),      pi_T(x) = sum_j P_j(T) C(x, j),
at x = 1, 2, 3 (psi = number of nontrivial cycles, moved = number of non-fixed points, cyc = number of cycles).
At x = 1 this is Szele's mean n!/2^(n-1). Exits 1 on failure."""
import itertools, math, sys
from fractions import Fraction
def tables(A):
    n=len(A); full=1<<n
    out=[sum(1<<u for u in range(n) if A[v][u]) for v in range(n)]
    dp=[[0]*n for _ in range(full)]
    for v in range(n): dp[1<<v][v]=1
    for m in range(1,full):
        for v in range(n):
            c=dp[m][v]
            if c:
                nx=out[v]&~m
                while nx:
                    u=(nx&-nx).bit_length()-1; nx&=nx-1; dp[m|1<<u][u]+=c
    return [sum(r) for r in dp]
def covers(A):
    n=len(A); full=(1<<n)-1; h=tables(A); h[0]=0
    layer=[1]+[0]*full; P=[0]*(n+1)
    for j in range(1,n+1):
        nxt=[0]*(full+1)
        for m in range(1,full+1):
            s=m
            while s:
                if h[s] and layer[m^s]: nxt[m]+=h[s]*layer[m^s]
                s=(s-1)&m
        layer=nxt; P[j]=layer[full]
    return P
ok=True
for n in ([3,4,5] if len(sys.argv)>1 else [3,4]):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    tot={1:0,2:0,3:0}
    for bits in range(1<<len(pairs)):
        A=[[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if bits>>k&1: A[i][j]=1
            else: A[j][i]=1
        P=covers(A)
        for x in tot: tot[x]+=sum(P[j]*math.comb(x,j) for j in range(1,n+1))
    for x in tot:
        rhs=Fraction(0)
        for perm in itertools.permutations(range(n)):
            seen=[0]*n; lens=[]
            for i in range(n):
                if not seen[i]:
                    l=0; j=i
                    while not seen[j]: seen[j]=1; j=perm[j]; l+=1
                    lens.append(l)
            if any(l%2==0 for l in lens): continue
            rhs+=Fraction(2**sum(1 for l in lens if l>1), 2**sum(l for l in lens if l>1))*x**len(lens)
        good=(tot[x]==rhs*2**len(pairs)); ok&=good
        print(("PASS" if good else "FAIL"),"n=%d x=%d: sum pi_T(x) = %d"%(n,x,tot[x]))
print("ALL PASS" if ok else "FAILURES"); sys.exit(0 if ok else 1)
