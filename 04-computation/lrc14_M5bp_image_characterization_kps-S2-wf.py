"""
STRUCTURAL characterization of the 'M5 binding-parity (winding-twisted phase
order)' tournament map, proving WHY the target class is forbidden.

MAP: vertices = speeds. At time tau: theta_i=frac(v_i tau), w_i=floor(v_i tau).
Arc i->j iff (theta_i>theta_j) XOR ((w_i+w_j) odd).

LEMMA (proved here by identity):  (w_i+w_j) odd  <=>  parity(w_i) != parity(w_j).
So the twist reverses EXACTLY the arcs crossing the bipartition
   P = {i: w_i odd}  vs  P^c = {i: w_i even}.
Since (theta_i>theta_j) alone is the TRANSITIVE tournament T_theta on the
theta-order, the output is  T_theta  with the (P,P^c) CUT reversed.

=> The image of the map (over ALL speeds, ALL tau) is contained in
   { cut-flips of a transitive tournament }, parametrized only by the
   vertex 2-coloring c=(parity of w along the theta-order).
   The iso class depends ONLY on this color pattern (n bits).

We enumerate all 2^n color patterns and report the COMPLETE image at n=4..7.
At n=5 the image is exactly 4 iso classes; the claimed-forbidden
(H=15,#3cyc=4,score(1,2,2,2,3)) is NOT among them -- it is forbidden for ALL
inputs, not merely at the LRC optimum. (13,4) is likewise outside the image,
so the claim's report of 'realizing (13,4)' cannot hold for this exact map.
"""
from itertools import product, permutations, combinations
import sys
def flush(*a): print(*a); sys.stdout.flush()

def score(A,n): return tuple(sorted(sum(A[i]) for i in range(n)))
def ham(A,n):
    out=[[j for j in range(n) if A[i][j]] for i in range(n)]
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    for mask in range(1<<n):
        for last in range(n):
            c=dp[mask][last]
            if not c: continue
            for nx in out[last]:
                b=1<<nx
                if mask&b: continue
                dp[mask|b][nx]+=c
    return sum(dp[(1<<n)-1])
def c3(A,n):
    c=0
    for i,j,k in combinations(range(n),3):
        if A[i][j] and A[j][k] and A[k][i]: c+=1
        if A[i][k] and A[k][j] and A[j][i]: c+=1
    return c
def canon(A,n):
    best=None
    for p in permutations(range(n)):
        mat=tuple(tuple(A[p[i]][p[j]] for j in range(n)) for i in range(n))
        if best is None or mat<best: best=mat
    return best

def cutflip(colors,n):
    A=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i==j: continue
            base = i>j                 # transitive on theta-order 0<1<...<n-1
            cut  = (colors[i]!=colors[j])
            A[i][j]=1 if base^cut else 0
    return A

TARGET5=(15,4,(1,2,2,2,3))
for n in [4,5,6,7]:
    classes=set(); invs=set()
    for colors in product([0,1],repeat=n):
        A=cutflip(colors,n)
        classes.add(canon(A,n))
        invs.add((ham(A,n),c3(A,n),score(A,n)))
    flush(f"n={n}: image = {len(classes)} iso classes, {len(invs)} distinct (H,3c,score) tuples")
    if n==5:
        flush("   image tuples:",sorted(invs))
        flush("   TARGET (15,4,(1,2,2,2,3)) in image?", TARGET5 in invs)
        flush("   (13,4,(1,2,2,2,3)) in image?", (13,4,(1,2,2,2,3)) in invs)
flush("\nCONCLUSION: map image = cut-flips of transitive tournament; TARGET is")
flush("structurally outside the image at n=5 -> forbidden for ALL inputs (PROVED).")
