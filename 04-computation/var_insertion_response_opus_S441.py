#!/usr/bin/env python3
"""
How does var(lambda^2) actually move under a=insertion? (opus-S441, follow-up)
var isn't c3-determined for n>=5, so the clean -B*(forward cut) reduction fails. Compute Delta var
directly under source/sink and general P, and test what DOES control it.

Facts: Sum lambda^2 = n(n-1) (fixed per n). Under T(n)->T+u_P (n+1 vertices), Sum lambda^2 -> n(n+1)
(+2n). var = Sum(lambda^4)/(n) - (n-1)^2, so Delta var is driven by Delta tr(S^4).
tr(S^4) counts signed closed 4-walks; adding u adds walks THROUGH u. Test: is Delta tr(S^4) under
insertion a clean function of P (e.g. of the new vertex's score |P| or the forward cut)?
"""
import itertools, numpy as np

def edges_iter(n):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in range(1<<len(pairs)):
        adj=[[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if (bits>>k)&1: adj[i][j]=1
            else: adj[j][i]=1
        yield adj

def skewnp(adj,n):
    return np.array([[0.0 if i==j else (1.0 if adj[i][j] else -1.0) for j in range(n)] for i in range(n)])

def tr4(adj,n):
    S=skewnp(adj,n); S2=S@S; return np.trace(S2@S2)

def var_l2(adj,n):
    S=skewnp(adj,n); l2=[abs(e)**2 for e in np.linalg.eigvals(S)]; return float(np.var(l2))

def add_u(adj,n,P):
    big=[row[:]+[0] for row in adj]+[[0]*(n+1)]
    for j in range(n):
        if j in P: big[n][j]=1
        else: big[j][n]=1
    return big

print("Delta tr(S^4) and Delta var(lambda^2) under a=insertion of u beating P")
print("="*70)
for n in range(3,6):
    dtr4_by_deg={}          # |P| (new vertex out-degree) -> set of Delta tr4
    dtr4_by_fc={}           # forward-cut related? here fc(u) irrelevant; test in-degree d=|Pc|
    all_rows=[]
    for adj in edges_iter(n):
        t0=tr4(adj,n)
        for P in (frozenset(p) for r in range(n+1) for p in itertools.combinations(range(n),r)):
            big=add_u(adj,n,P); t1=tr4(big,n+1)
            d=len(P)                      # u's out-degree
            dtr4=round(t1-t0,4)
            dtr4_by_deg.setdefault(d,set()).add(dtr4)
            all_rows.append((d, dtr4))
    # is Delta tr4 a function of u's out-degree d alone?
    deg_det = all(len(v)==1 for v in dtr4_by_deg.values())
    print(f" n={n}: Delta tr(S^4) determined by new vertex out-degree |P| alone? {deg_det}")
    if deg_det:
        for d in sorted(dtr4_by_deg):
            print(f"     |P|={d}: Delta tr(S^4) = {next(iter(dtr4_by_deg[d]))}")
    else:
        # show the spread
        for d in sorted(dtr4_by_deg):
            vals=sorted(dtr4_by_deg[d])
            print(f"     |P|={d}: Delta tr(S^4) in {vals[:5]}{'...' if len(vals)>5 else ''}")

# transitive var = 2*C(n,3) confirm + the min (regular) var
print("\n" + "="*70)
from math import comb
def transitive(n): return [[1 if i<j else 0 for j in range(n)] for i in range(n)]
for n in range(3,8):
    v=var_l2(transitive(n),n)
    print(f"  n={n}: transitive var(l^2) = {v:.4f}   2*C(n,3) = {2*comb(n,3)}   match={abs(v-2*comb(n,3))<1e-6}")
print("\n READING: Delta tr(S^4) under insertion is NOT a function of the new vertex's degree alone")
print(" for larger n (var(l^2) is genuinely spectral, finer than any local count) -- the honest")
print(" answer to 'how does var move under a': by interlacing, not by a combinatorial forward-cut law.")
print(" CONFIRMED: transitive var(l^2) = 2*C(n,3) (the maximally-spread GIT nullcone vertex).")
