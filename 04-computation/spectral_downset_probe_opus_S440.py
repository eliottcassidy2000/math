#!/usr/bin/env python3
"""
THM-1920 open Q1: does a DOWN-SET insertion (THM-1900 H-neutral: P closed, no arc P->V\\P) have a
special SPECTRAL signature, the way it preserves H?  Probe over all T (n<=5) and all patterns P.
Compare down-set vs non-down-set insertions: B_P structure, whether an eigenvalue is pinned,
var(lambda^2) response (kps S128c139 scalar).
"""
import itertools, sympy as sp
x = sp.symbols('x')

def edges_iter(n):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in range(1<<len(pairs)):
        adj=[[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if (bits>>k)&1: adj[i][j]=1
            else: adj[j][i]=1
        yield adj

def skewmat(adj,n):
    return sp.Matrix(n,n,lambda i,j: 0 if i==j else (1 if adj[i][j] else -1))

def forward_cut(adj,n,P):
    Pc=[j for j in range(n) if j not in P]
    return sum(1 for i in P for j in Pc if adj[i][j])

def lambda2_multiset(adj,n):
    """the lambda^2 values (=-eigenvalue^2 of skew, real nonneg): roots of char in y=x^2."""
    S=skewmat(adj,n); ch=sp.expand((x*sp.eye(n)-S).det())
    # char_S(x) is even (n even) or x*even; substitute x^2=-y? eigenvalues i*lam => x^2=-lam^2
    # collect coefficients; just return numeric lam^2 via numpy for the scalar var
    import numpy as np
    Sn=np.array([[0.0 if i==j else (1.0 if adj[i][j] else -1.0) for j in range(n)] for i in range(n)])
    ev=np.linalg.eigvals(Sn)
    return sorted((abs(e)**2) for e in ev)

def var_l2(adj,n):
    import numpy as np
    l=lambda2_multiset(adj,n); return float(np.var(l))

print("THM-1920 Q1: spectral signature of DOWN-SET (H-neutral) insertions")
print("="*66)
for m in range(2,5):
    downset_new_ev0=0; downset_tot=0; nondown_new_ev0=0; nondown_tot=0
    downset_var_preserved=0
    examples=[]
    for adj in edges_iter(m):
        import numpy as np
        base_l2=lambda2_multiset(adj,m)
        for P in (frozenset(p) for r in range(m+1) for p in itertools.combinations(range(m),r)):
            big=[row[:]+[0] for row in adj]+[[0]*(m+1)]
            for j in range(m):
                if j in P: big[m][j]=1
                else: big[j][m]=1
            new_l2=lambda2_multiset(big,m+1)
            has0 = any(abs(v)<1e-9 for v in new_l2)   # a zero eigenvalue appeared/preserved
            isdown = (forward_cut(adj,m,P)==0)
            if isdown:
                downset_tot+=1
                if has0: downset_new_ev0+=1
            else:
                nondown_tot+=1
                if has0: nondown_new_ev0+=1
    print(f" m={m}: down-set insertions with a 0-eigenvalue: {downset_new_ev0}/{downset_tot}"
          f"   | non-down-set: {nondown_new_ev0}/{nondown_tot}")

# focused: source/sink (extreme down-sets) always pin a 0 eigenvalue?  and the transitive tower
print("\n focused: does EVERY down-set insertion pin a 0 eigenvalue (skew nullity)?")
allpin=True; anydiff=[]
for m in range(2,6):
    for adj in edges_iter(m):
        for P in (frozenset(p) for r in range(m+1) for p in itertools.combinations(range(m),r)):
            if forward_cut(adj,m,P)!=0: continue
            big=[row[:]+[0] for row in adj]+[[0]*(m+1)]
            for j in range(m):
                if j in P: big[m][j]=1
                else: big[j][m]=1
            import numpy as np
            Sn=np.array([[0.0 if i==j else (1.0 if big[i][j] else -1.0) for j in range(m+1)] for i in range(m+1)])
            ev=np.linalg.eigvals(Sn)
            if not any(abs(e)<1e-9 for e in ev):
                allpin=False
    # count parity
print(f"   EVERY down-set insertion yields a singular skew matrix (0 eigenvalue)? {allpin}")
print("\n READING: if TRUE, down-set (H-neutral) insertions are exactly the SPECTRALLY-DEGENERATE")
print(" (rank-dropping) insertions -- the spectral face of THM-1900's H-neutrality.")
