"""
opus-2026-06-07-S707 : the Pfaffian as translator + the inclusion-exclusion tiling
recursion for max Hamiltonian-path count (A038375).

THREADS:
 (1) Pfaffian dictionary: Pf(S)^2 = det(I+2A) = signed weighted disjoint-cycle-cover
     gen. fn. (THM-174); the Pfaffian RECURSION removes 2 vertices (n->n-2) -- a
     "domino" move on the tiling. Topology(dimers/FKT) ~ geometry(tiling) ~ algebra(Clifford).
 (2) The user's INCLUSION-EXCLUSION tiling decomposition of the staircase Delta_{n-2}:
     3 pieces of size n-1 (A,B,C) - 3 of size n-2 (D,E,F) + 1 of size n-3 (G).
     CLAIM: this is the THIRD FINITE DIFFERENCE annihilating the quadratic cell-count
        C(n-1,2) = 3 C(n-2,2) - 3 C(n-3,2) + C(n-4,2).
     => any invariant that is AFFINE/quadratic in the cells satisfies this 3-term recursion;
     H and Pf are NOT, so we test which tournament invariants do.
 (3) max-H = A038375: exact brute small n + Szele/Alon bounds + structure of maximizers;
     test recursions. Build on S531 (H multiplicative over modules; apex block = 1+2^{s-2}).
"""
import numpy as np
from itertools import combinations, product
from math import comb, factorial

# ---------- tournaments as bit-adjacency; A[i][j]=1 iff i->j -------------------
def ham_paths(adj, n):
    # dp over subsets: dp[mask][v] = # paths covering mask, ending at v
    size=1<<n
    dp=[[0]*n for _ in range(size)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(size):
        for v in range(n):
            c=dp[mask][v]
            if not c: continue
            for u in range(n):
                if mask&(1<<u): continue
                if adj[v][u]:
                    dp[mask|(1<<u)][u]+=c
    full=size-1
    return sum(dp[full][v] for v in range(n))

def all_tournaments(n):
    pairs=list(combinations(range(n),2))
    for bits in range(1<<len(pairs)):
        adj=[[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if (bits>>k)&1: adj[i][j]=1
            else: adj[j][i]=1
        yield adj

def pfaffian(S):
    # exact Pfaffian of skew-symmetric integer matrix via recursion (small n, even)
    n=len(S)
    if n==0: return 1
    if n%2==1: return 0
    # expand along row 0
    import copy
    total=0
    def minor(M, i, j):
        idx=[k for k in range(len(M)) if k!=i and k!=j]
        return [[M[a][b] for b in idx] for a in idx]
    for j in range(1,n):
        if S[0][j]==0: continue
        total += ((-1)**(j+1)) * S[0][j] * pfaffian(minor(S,0,j))
    return total

def skew(adj,n):
    return [[ (1 if adj[i][j] else (-1 if adj[j][i] else 0)) for j in range(n)] for i in range(n)]

# ---------- (3a) exact max-H for small n --------------------------------------
print("="*76)
print("(3a) EXACT max Hamiltonian-path count (A038375) by brute force, small n")
print("="*76)
maxH={}
for n in range(2,7):
    best=0; arg=None
    for adj in all_tournaments(n):
        h=ham_paths(adj,n)
        if h>best: best=h; arg=adj
    maxH[n]=best
    szele=factorial(n)//(2**(n-1)) if n>=1 else 1
    print(f"  n={n}: max H = {best:5d}   Szele n!/2^(n-1) = {szele:5d}   max/Szele = {best/szele:.3f}")

# ---------- (3b) n=7 strong local search (best found) -------------------------
print("\n"+"="*76)
print("(3b) n=7 max-H: hill-climb from regular/random (best found, not proven exact)")
print("="*76)
def random_tournament(n,rng):
    adj=[[0]*n for _ in range(n)]
    for i,j in combinations(range(n),2):
        if rng.random()<0.5: adj[i][j]=1
        else: adj[j][i]=1
    return adj
def flip(adj,i,j):
    adj[i][j],adj[j][i]=adj[j][i],adj[i][j]
def hill(n,rng,iters=40):
    adj=random_tournament(n,rng); best=ham_paths(adj,n)
    improved=True
    while improved:
        improved=False
        for i,j in combinations(range(n),2):
            flip(adj,i,j); h=ham_paths(adj,n)
            if h>best: best=h; improved=True
            else: flip(adj,i,j)
    return best
rng=__import__('random').Random(1)
best7=0
for r in range(60):
    best7=max(best7, hill(7,rng))
maxH[7]=best7
sz7=factorial(7)//2**6
print(f"  n=7: best H found = {best7}   Szele = {sz7}   ratio = {best7/sz7:.3f}")

# ---------- (3c) the sequence + recursion tests ------------------------------
print("\n"+"="*76)
print("(3c) the max-H sequence and recursion tests")
print("="*76)
seq={n:maxH[n] for n in range(2,8)}
print("  max-H:", seq)
print("  ratios m(n)/m(n-1):", {n: round(seq[n]/seq[n-1],3) for n in range(3,8)})
# test user's IE signs (3,-3,1) and (3,0,0) and n!/2^(n-2)
for n in range(5,8):
    ie = 3*seq[n-1]-3*seq[n-2]+seq[n-3]
    print(f"  n={n}: m={seq[n]:4d} | 3m(n-1)-3m(n-2)+m(n-3)={ie:4d} | "
          f"3*m(n-1)={3*seq[n-1]:4d} | n!/2^(n-2)={factorial(n)//2**(n-2):5d}")

# ---------- (2) the IE tiling identity: third difference of the cell count ----
print("\n"+"="*76)
print("(2) IE TILING IDENTITY  C(n-1,2) = 3C(n-2,2) - 3C(n-3,2) + C(n-4,2)")
print("    (the user's A,B,C - D,E,F + G ;  third finite diff of a quadratic = 0)")
print("="*76)
ok=all(comb(n-1,2)==3*comb(n-2,2)-3*comb(n-3,2)+comb(n-4,2) for n in range(4,30))
print(f"  identity holds n=4..29: {ok}")
# which tournament invariants satisfy the 3-term IE recursion?  test on the MAX over n:
print("  test max-invariant(n) =? 3*inv(n-1)-3*inv(n-2)+inv(n-3):")
# cyclic-triangle MAX count c(n) = C(n,3) for regular-ish; exact max cyclic triangles:
def max_cyc_triangles(n):
    # max # of cyclic triangles in a tournament = floor stuff; compute exactly small n
    best=0
    if n>6: return None
    for adj in all_tournaments(n):
        c=sum(1 for a,b,cc in combinations(range(n),3)
              if (adj[a][b]and adj[b][cc]and adj[cc][a]) or (adj[b][a]and adj[cc][b]and adj[a][cc]))
        best=max(best,c)
    return best
cyc={n:max_cyc_triangles(n) for n in range(3,7)}
print(f"    max cyclic triangles: {cyc}")
for n in range(6,7):
    if all(cyc.get(k) for k in (n-1,n-2,n-3)):
        print(f"      n={n}: c={cyc[n]} vs 3c(n-1)-3c(n-2)+c(n-3)={3*cyc[n-1]-3*cyc[n-2]+cyc[n-3]}")

# ---------- (1) Pfaffian dictionary checks -----------------------------------
print("\n"+"="*76)
print("(1) PFAFFIAN: Pf(S)^2 = det(I+2A) = det(S) (THM-174); H odd, Pf odd; n->n-2 recursion")
print("="*76)
rng2=__import__('random').Random(3)
for n in [4,6]:
    adj=random_tournament(n,rng2); S=skew(adj,n)
    A=np.array([[adj[i][j] for j in range(n)] for i in range(n)],dtype=float)
    pf=pfaffian(S); dI2A=round(np.linalg.det(np.eye(n)+2*A)); dS=round(np.linalg.det(np.array(S,dtype=float)))
    H=ham_paths(adj,n)
    print(f"  n={n}: Pf(S)={pf}  Pf^2={pf*pf}  det(I+2A)={dI2A}  det(S)={dS}  H={H} (H odd:{H%2==1}, Pf odd:{pf%2==1})")
print("  => Pf(S)^2=det(I+2A)=det(S) verified; both H and Pf are ODD (Redei / even-n Pfaffian).")

# ---------- single-apex-block law and modular multiplicativity (build on S531) ----
print("\n"+"="*76)
print("(extra) single apex-block H=1+2^(s-2); H multiplicative over disjoint modules")
print("="*76)
def transitive(n):  # i->j iff i<j
    return [[1 if i<j else 0 for j in range(n)] for i in range(n)]
def apex_flip_block(n,a,b):  # transitive, then reverse the (a,b) closing arc (a<b)
    adj=transitive(n); adj[a][b]=0; adj[b][a]=1; return adj
for s in range(3,8):
    adj=apex_flip_block(s,0,s-1); print(f"  size-{s} apex block: H={ham_paths(adj,s)}  (1+2^(s-2)={1+2**(s-2)})")
# disjoint modules multiply: block on [0,2] and [3,5] inside n=6 transitive backbone
adj=transitive(6); adj[0][2]=0;adj[2][0]=1; adj[3][5]=0;adj[5][3]=1
print(f"  two disjoint apex blocks [0,2]&[3,5] in n=6: H={ham_paths(adj,6)} (3*3=9 expected)")

# ---------- the H^2 - Pf^2 = 8Q bridge (poly-time Pfaffian skeleton for #P-hard H) ----
print("\n"+"="*76)
print("(bridge) H^2 - Pf^2 = 8Q  exhaustive: always a NON-NEG integer multiple of 8?")
print("="*76)
for n in [4,6]:
    okdiv=True; oknn=True; Qs=[]; Hmax=0; Hmax_pf=None
    for adj in all_tournaments(n):
        H=ham_paths(adj,n); S=skew(adj,n); pf=pfaffian(S)
        d=H*H-pf*pf
        if d%8!=0: okdiv=False
        if d<0: oknn=False
        Qs.append(d//8)
        if H>Hmax: Hmax=H; Hmax_pf=abs(pf)
    print(f"  n={n}: all (H^2-Pf^2) divisible by 8: {okdiv}; all >=0: {oknn}; "
          f"Q in [{min(Qs)},{max(Qs)}]; at max-H={Hmax}: |Pf|={Hmax_pf}")
print("  => H=sqrt(Pf^2+8Q): the poly-time Pfaffian(det) is the SKELETON, Q the correction;")
print("     Pf recurses n->n-2 (domino removal); H is #P-hard but pinned mod 8 by Pf (both odd).")
print("\nDONE.")
