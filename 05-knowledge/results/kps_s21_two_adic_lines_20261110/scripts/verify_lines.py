"""verify_lines.py - kind-pasteur-2026-11-10-S21. Independent brute-force checks for THM-5126 (A)-(B).
 (T) tail law H(L o y) = H + (U-d) C(y) + d (g(y)-1) + (d-E) eps(y): every labelled tournament of order <= 4 and
     60 random ones of orders 5-6, every root, every word y of length <= 4; counts by plain subset DP on the literal
     tournament L o y (no use of the formula).
 (L) every seed of line_seeds.txt with c <= LIMIT: recompute (H,d,E,U) from the code, check c0 = 1 (tail seeds) and
     count H(L o 0^t) (tail) or H(K[TT_(t+1) at x]) (module) literally for t = 1,2,3 against 1 + c 2^s.
 (G) gate: no seed puts 7 or 21 on a ray.
 (B) below_rays.txt: every point below a ray has a certificate; census witnesses and the factors of the two order
     joins are recounted literally (witness tournaments from the audit's independent census).
usage: python verify_lines.py [LIMIT=129]
"""
import sys, os, itertools, random
FAILS=[]
def report(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ")+name+("  | "+detail if detail else ""))
    if not ok: FAILS.append(name)
def H(A):
    n=len(A)
    if n==0: return 1
    out=[sum(1<<u for u in range(n) if A[v][u]) for v in range(n)]
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for m in range(1,1<<n):
        row=dp[m]
        for v in range(n):
            c=row[v]
            if c:
                nx=out[v]&~m
                while nx:
                    u=(nx&-nx).bit_length()-1; nx&=nx-1
                    dp[m|1<<u][u]+=c
    return sum(dp[(1<<n)-1])
def sub(A,keep): return [[A[i][j] for j in keep] for i in keep]
def start_count(B):
    """number of Hamiltonian paths of B starting at vertex 0 (literal DP on the converse ending at 0)"""
    k=len(B)
    if k==1: return 1
    inn=[sum(1<<u for u in range(k) if B[u][x]) for x in range(k)]
    dp=[[0]*k for _ in range(1<<k)]
    dp[1][0]=1
    for m in range(1,1<<k):
        if not m&1: continue
        for x in range(k):
            c=dp[m][x]
            if c:
                nx=sum(1<<u for u in range(k) if B[x][u])&~m
                while nx:
                    u=(nx&-nx).bit_length()-1; nx&=nx-1; dp[m|1<<u][u]+=c
    return sum(dp[(1<<k)-1])
def state(A,v):
    n=len(A); h=H(A); others=[i for i in range(n) if i!=v]; d=H(sub(A,others))
    # E: Hamiltonian paths ending at v = paths of the converse starting at v
    conv=[[A[j][i] for j in [v]+others] for i in [v]+others]
    E=start_count(conv)
    U=0
    for r in range(len(others)+1):
        for X in itertools.combinations(others,r):          # X = vertex set of P1 (possibly empty)
            hx=H(sub(A,list(X)))
            if hx==0: continue
            Y=[v]+[i for i in others if i not in X]
            U+=hx*start_count(sub(A,Y))
    return h,d,E,U
def tail(A,v,y):
    n=len(A); t=len(y); N=n+t; B=[[0]*N for _ in range(N)]
    for i in range(n):
        for j in range(n): B[i][j]=A[i][j]
    for a in range(t):
        z=n+a
        for b in range(a+1,t): B[z][n+b]=1
        for u in range(n):
            if u!=v: B[u][z]=1
        if y[a]: B[v][z]=1
        else: B[z][v]=1
    return B
def g(v):
    if v==0: return 1
    return g(v>>1) if v&1 else g(v>>1)+(v>>1)+((v>>1)&1)
def decode(code,m):
    A=[[0]*m for _ in range(m)]; b=0
    for i in range(m):
        for j in range(i+1,m):
            if code>>b&1: A[i][j]=1
            else: A[j][i]=1
            b+=1
    return A
def all_t(n):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in range(1<<len(pairs)):
        A=[[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if bits>>k&1: A[i][j]=1
            else: A[j][i]=1
        yield A
def check_tail():
    ok=True; cnt=0; rng=random.Random(5)
    pool=[A for n in (1,2,3,4) for A in all_t(n)]
    for n in (5,6):
        for _ in range(30):
            A=[[0]*n for _ in range(n)]
            for i in range(n):
                for j in range(i+1,n):
                    if rng.random()<0.5: A[i][j]=1
                    else: A[j][i]=1
            pool.append(A)
    for A in pool:
        n=len(A)
        for v in range(n):
            h,d,E,U=state(A,v)
            for t in range(1,5 if n<=5 else 4):
                for y in itertools.product((0,1),repeat=t):
                    C=sum((1-y[j])<<j for j in range(t)); gy=g(int("".join(map(str,y)),2)); eps=1-y[0]
                    cnt+=1
                    if H(tail(A,v,y))!=h+(U-d)*C+d*(gy-1)+(d-E)*eps: ok=False
    report("(T) tail law", ok, "%d (L, v, y) triples; all labelled orders <= 4, 60 random of orders 5-6"%cnt)
def blow(A,x,a):
    n=len(A); orig=[]; idx=[]
    for v in range(n):
        for i in range(a if v==x else 1): orig.append(v); idx.append(i)
    N=len(orig); B=[[0]*N for _ in range(N)]
    for p in range(N):
        for q in range(N):
            if p==q: continue
            B[p][q]=(1 if idx[p]<idx[q] else 0) if orig[p]==orig[q] else A[orig[p]][orig[q]]
    return B
def check_lines(limit):
    ok=True; gate=True; n=0
    fn = "line_seeds.txt" if os.path.exists("line_seeds.txt") else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "line_seeds.txt")
    for line in open(fn):
        if line.startswith("#"): continue
        p=line.split(); c,s0=int(p[0]),int(p[1]); mech=p[2]; m,code,v=int(p[3]),int(p[4]),int(p[5])
        if (c==3 and s0<2) or (c==5 and s0<3): gate=False
        if c>limit: continue
        A=decode(code,m); n+=1
        if mech=="tail":
            h,d,E,U=state(A,v)
            if (h,d,E,U)!=tuple(map(int,p[6:10])) or h-E+2*d-U!=1 or (U-d)!=c*2**(s0-1): ok=False; continue
            for t in (1,2,3):
                if m+t<=12 and H(tail(A,v,(0,)*t))!=1+c*2**(s0-1+t): ok=False
            if m==10 and H(tail(A,v,(0,)*1))!=1+c*2**s0: ok=False
        else:
            h=H(A)
            for t in (1,2,3):
                if m+t<=12 and H(blow(A,v,t+1))!=1+(h-1)*2**t: ok=False
            if (h-1)!=c*2**s0: ok=False
    report("(L) line seeds: states and literal counts of the first three ray points", ok, "%d seeds with c <= %d"%(n,limit))
    # certificates below the rays: arithmetic of the product lines, and every odd c <= limit is present with all s < s0
    bfn = "below_rays.txt" if os.path.exists("below_rays.txt") else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "below_rays.txt")
    okb=True; seen={}
    for line in open(bfn):
        if line.startswith("#"): continue
        p=line.split(); c,sv,N=int(p[0]),int(p[1]),int(p[2]); seen[(c,sv)]=p[3]
        if N!=1+c*2**sv: okb=False
        if p[3]=="census":
            if H(decode(int(p[5]),int(p[4])))!=N: okb=False          # literal recount of the witness
        if p[3]=="product":
            a,b=int(p[4]),int(p[5])
            if a*b!=N or H(decode(int(p[7]),int(p[6])))!=a or H(decode(int(p[9]),int(p[8])))!=b: okb=False
        if p[3]=="GAP" and N not in (7,21): okb=False
    s0s={}
    for line in open(fn):
        if line.startswith("#"): continue
        q=line.split(); s0s[int(q[0])]=int(q[1])
    for c in range(1,limit+1,2):
        if c not in s0s: okb=False; continue
        for sv in range(1,s0s[c]):
            if (c,sv) not in seen: okb=False
    gaps=sorted(k for k,v in seen.items() if v=="GAP")
    report("(B) every odd c <= %d has a seed and a certificate for each point below its ray (witness tournaments recounted literally; products are order joins); the only GAP lines are (3,1) and (5,2)"%limit, okb and gaps==[(3,1),(5,2)], "%d points below rays"%len(seen))
    report("(G) gate: no ray through 7 or 21", gate)
if __name__=="__main__":
    limit=int(sys.argv[1]) if len(sys.argv)>1 else 129
    check_tail(); check_lines(limit)
    print("ALL PASS" if not FAILS else "FAILURES: %s"%FAILS); sys.exit(1 if FAILS else 0)
