"""towers.py - the tower (gadget) calculus for portrait states, with brute-force verification.

Pair (B, O): B a tournament, O a set of vertices.  State s(B,O) = (r, M, h, E):
   h = H(B),  r = H(B) - W_B(O),  M = U_B(O),  E = End_B(O)   (Hamiltonian paths of B ending in O).
Closing:  T_t(B,O) = (B => TT_t => v, v -> O, B-O -> v) has H = r + 2^t M  (t >= 1)   [known cut form of the tail law].
Gadget (X, P): X a tournament, P a set of vertices of X.  Level-up:
   B+ = (B => X => v, v -> O, B-O -> v),  O+ = P.
Claim (gadget law):  s(B+, O+) = G(X,P) s(B,O) with the 4x4 integer matrix G built from
   H(X), Pi2(X), r_X(P) = H(X) - W_X(P), kappa_X(P) = Pi2(X) - 2 W_X(P) - 2 W3_X(P),
   Start_X(P), U_X(P), U3_X(P), End_X(P), V_X(P).
"""
import itertools, random, sys
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
    hs=[sum(e) for e in en]; hs[0]=1
    return en,st,hs
def H(A):
    if not A: return 1
    return tables(A)[2][(1<<len(A))-1]
def state(A,O):
    """(r, M, h, E) of the pair (B, O) by direct counting."""
    n=len(A); full=(1<<n)-1; en,st,hs=tables(A); h=hs[full]
    W=0; M=0
    for X in range(1,full):
        Y=full^X
        for e in O:
            if X>>e&1 and en[X][e]:
                for f in O:
                    if Y>>f&1: W+=en[X][e]*st[Y][f]
    for o in O:
        for X in range(1,full+1):
            if X>>o&1: M+=st[X][o]*hs[full^X]
    E=sum(en[full][o] for o in O)
    return (h-W, M, h, E)
def gadget(A,P):
    """constants of the gadget (X, P) and the 4x4 matrix G acting on (r, M, h, E)."""
    n=len(A); full=(1<<n)-1; en,st,hs=tables(A); hX=hs[full]
    pi2=sum(hs[X]*hs[full^X] for X in range(full+1))
    W=0; W3=0; U=0; U3=0
    subs=list(range(full+1))
    for X in range(1,full+1):
        rest=full^X
        # U: designated block X starting at o in P, other block = rest (possibly empty)
        sX=sum(st[X][o] for o in P if X>>o&1)
        if sX:
            U+=sX*hs[rest]
            # U3: two further ordered nonempty paths covering rest
            Y=rest
            while Y:
                Z=rest^Y
                if Z: U3+=sX*hs[Y]*hs[Z]
                Y=(Y-1)&rest
        # W and W3: D1 = X ending at e, D2 = Y starting at f (f != e), third block nonempty for W3
        eX=[(e,en[X][e]) for e in P if X>>e&1 and en[X][e]]
        if eX:
            Y=rest
            while Y:
                Z=rest^Y
                for e,ce in eX:
                    for f in P:
                        if Y>>f&1 and st[Y][f]:
                            if Z==0: W+=ce*st[Y][f]
                            else: W3+=ce*st[Y][f]*hs[Z]
                Y=(Y-1)&rest
    start=sum(st[full][o] for o in P); end=sum(en[full][o] for o in P)
    V=sum(en[X][o]*hs[full^X] for o in P for X in range(1,full+1) if X>>o&1)
    rX=hX-W; kap=pi2-2*W-2*W3
    # state order (r, M, h, E)
    G=[[rX, kap, 0, 0],
       [start+(U-start), 2*start+4*(U-start)+2*U3, start+(U-start), -start],
       [hX, pi2, 0, 0],
       [end, V, -end, 0]]
    return dict(H=hX,pi2=pi2,r=rX,kappa=kap,W=W,W3=W3,start=start,U=U,U3=U3,end=end,V=V,G=G)
def apply(G,s): return tuple(sum(G[i][j]*s[j] for j in range(4)) for i in range(4))
def level_up(A,O,AX):
    """B+ = (B => X => v, v -> O, B-O -> v); vertices: B (0..n-1), X (n..n+k-1), v last."""
    n=len(A); k=len(AX); N=n+k+1; T=[[0]*N for _ in range(N)]
    for i in range(n):
        for j in range(n): T[i][j]=A[i][j]
        for j in range(k): T[i][n+j]=1
    for i in range(k):
        for j in range(k): T[n+i][n+j]=AX[i][j]
        T[n+i][N-1]=1
    for i in range(n):
        if i in O: T[N-1][i]=1
        else: T[i][N-1]=1
    return T
def TT(a): return [[1 if i<j else 0 for j in range(a)] for i in range(a)]
def rand_t(n,rng):
    A=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            if rng.random()<0.5: A[i][j]=1
            else: A[j][i]=1
    return A
if __name__=="__main__":
    rng=random.Random(20261110); bad=0; tests=0
    # 1. closing law H = r + 2^t M
    for _ in range(60):
        n=rng.randint(1,5); A=rand_t(n,rng); O=[i for i in range(n) if rng.random()<0.5]
        r,M,h,E=state(A,O)
        for t in (1,2,3):
            T=level_up(A,O,TT(t)); tests+=1
            if H(T)!=r+2**t*M: bad+=1
    print("closing law H = r + 2^t M: tests",tests,"failures",bad)
    # 2. gadget law on random bases, random gadgets (no unit assumption anywhere)
    bad=0; tests=0
    for _ in range(400):
        n=rng.randint(1,4); A=rand_t(n,rng); O=[i for i in range(n) if rng.random()<0.6]
        k=rng.randint(1,4); AX=rand_t(k,rng); P=[i for i in range(k) if rng.random()<0.6]
        s=state(A,O); g=gadget(AX,P)
        T=level_up(A,O,AX); s2=state(T,[n+p for p in P]); tests+=1
        if s2!=apply(g["G"],s):
            bad+=1
            if bad<=3: print("  FAIL",A,O,AX,P,s,s2,apply(g["G"],s))
    print("gadget law s(B+,O+) = G(X,P) s(B,O): tests",tests,"failures",bad)
    # 3. transitive twins: kappa = 0 iff delta in {1,2}; closed forms
    print("transitive gadget TT_a, pair (i, i+delta), m = a-i-delta after:")
    for delta in (1,2,3,4):
        ks=set(); ok=True
        for i in range(1,4):
            for m in range(0,3):
                a=i+delta+m; g=gadget(TT(a),[i-1,i-1+delta]); ks.add((g["r"],g["kappa"]))
                if delta<=2:
                    pred=[[-1,0,0,0],[(2**delta+1)*2**m, 2*(3**delta+2**delta)*2**(i-1)*3**m, (2**delta+1)*2**m, -(1 if i==1 else 0)],
                          [1,2**a,0,0],[(1 if m==0 else 0),(2**delta+1)*2**(i-1),-(1 if m==0 else 0),0]]
                    if g["G"]!=pred: ok=False
        print("  delta",delta,"(r_X, kappa_X) over i<=3, m<=2:",sorted(ks),("| closed-form matrix matches: %s"%ok if delta<=2 else ""))
