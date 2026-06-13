#!/usr/bin/env python3
"""
u(21): (1) corrected elementary cherry/common-neighbour UPPER bound (exact integer),
       (2) the proven literature values (Alexeev-Mixon-Parshall 2024) for comparison,
       (3) an INDEPENDENT numerical relaxation that witnesses unit-distance counts
           beating the triangular-lattice 47 on 21 points (float; marked as such).
monad-explorer-2026-06-06-S710.

NB a bug in s710 v1 cn_upper halved the cherry cap. CORRECT bound:
   #cherries = sum_v C(d_v,2) = sum_{pairs} (#common nbrs) <= 2*C(N,2) = N(N-1).
So sum_v d_v(d_v-1) <= 2*N(N-1), maximised (equal degrees) at
   d <= (1+sqrt(8N-7))/2,   U <= N(1+sqrt(8N-7))/4.
"""
import math, random

# ---- proven / best-known values (Alexeev-Mixon-Parshall 2024, arXiv:2412.11914) ----
# u(n) EXACT for n<=21; bounds for 22..30.
u_exact = {0:0,1:0,2:1,3:3,4:5,5:7,6:9,7:12,8:14,9:18,10:20,11:23,12:27,
           13:30,14:33,15:37,16:41,17:43,18:46,19:50,20:54,21:57}
u_bounds = {22:(60,61),23:(64,66),24:(68,72),25:(72,78),26:(76,84),
            27:(81,90),28:(85,96),29:(89,103),30:(93,110)}

def harborth(N):                       # triangular-lattice (penny) maximum, exact
    return math.floor(3*N - math.sqrt(12*N - 3))

def cn_upper_exact(N):
    """Exact integer max of U=sum d_v/2 over degree seqs with sum d_v(d_v-1) <= 2N(N-1)."""
    cap = 2*N*(N-1)                    # CORRECTED cap on sum d_v(d_v-1)
    deg = [0]*N; used = 0
    while True:
        m = min(deg); idx = deg.index(m); marg = 2*deg[idx]
        if used + marg <= cap:
            deg[idx]+=1; used+=marg
        else:
            break
    return sum(deg)//2

print("="*78)
print("u(N): proven value vs elementary bounds.  [AMP24 = Alexeev-Mixon-Parshall 2024]")
print("="*78)
print(f"{'N':>3} {'tri(Harborth)':>14} {'u(N) PROVEN':>14} {'cherry-upper':>13}  gap(u-tri)  gap(cherry-u)")
for N in range(15, 25):
    H = harborth(N); CU = cn_upper_exact(N)
    if N in u_exact:
        uv = str(u_exact[N]); g1 = u_exact[N]-H; g2 = CU-u_exact[N]
    else:
        lo,hi = u_bounds[N]; uv=f"[{lo},{hi}]"; g1=f"{lo-H}..{hi-H}"; g2=f"{CU-hi}..{CU-lo}"
    mark = "  <== N=21 (SETTLED: u=57)" if N==21 else ""
    print(f"{N:>3} {H:>14} {uv:>14} {CU:>13}  {str(g1):>9}  {str(g2):>11}{mark}")

print()
print("Key: triangular lattice gives 47 at N=21; the PROVEN maximum is 57 (gap 10).")
print("     => triangular/Eisenstein single-norm lattice is NOT optimal at N=21.")
print("     The elementary cherry bound (71) is far above 57: AMP24 needed heavy")
print("     forbidden-subgraph + embedding machinery to close it.")

# ---------------------------------------------------------------------------
# Independent numerical witness: relaxation to maximise near-unit pairs on N pts.
# Soft potential with a narrow attractive well at d=1 plus short-range repulsion.
# float only -> a WITNESS that >47 unit distances are realizable on 21 points.
# ---------------------------------------------------------------------------
def relax_count(N, seed, iters=4000, lr=0.02):
    rng = random.Random(seed)
    P = [[rng.uniform(-3,3), rng.uniform(-3,3)] for _ in range(N)]
    for it in range(iters):
        # anneal well width
        t = it/iters
        sigma = 0.35*(1-t) + 0.02*t          # narrows the unit well over time
        rep   = 0.6*(1-t)                      # repulsion fades
        F = [[0.0,0.0] for _ in range(N)]
        for i in range(N):
            for j in range(i+1,N):
                dx=P[i][0]-P[j][0]; dy=P[i][1]-P[j][1]
                d=math.hypot(dx,dy)+1e-12
                ux,uy=dx/d,dy/d
                # attractive well at d=1: gaussian-derivative pull
                g=math.exp(-((d-1.0)**2)/(2*sigma*sigma))
                fmag=-(d-1.0)/(sigma*sigma)*g          # pulls toward d=1
                # short range repulsion to avoid collapse (d<0.85)
                if d<0.85: fmag += rep*(0.85-d)/d
                F[i][0]+=fmag*ux; F[i][1]+=fmag*uy
                F[j][0]-=fmag*ux; F[j][1]-=fmag*uy
        for i in range(N):
            P[i][0]+=lr*F[i][0]; P[i][1]+=lr*F[i][1]
    # count near-unit pairs (tolerance)
    best_tol=0
    for tol in (0.02,0.01,0.005):
        c=0
        for i in range(N):
            for j in range(i+1,N):
                d=math.hypot(P[i][0]-P[j][0],P[i][1]-P[j][1])
                if abs(d-1.0)<tol: c+=1
        best_tol=max(best_tol,c) if tol==0.02 else best_tol
        if tol==0.005: tight=c
    return c if False else (count_tol(P,N,0.01), count_tol(P,N,0.003))

def count_tol(P,N,tol):
    c=0
    for i in range(N):
        for j in range(i+1,N):
            if abs(math.hypot(P[i][0]-P[j][0],P[i][1]-P[j][1])-1.0)<tol: c+=1
    return c

print()
print("="*78)
print("INDEPENDENT numerical witness (relaxation, FLOAT) on N=21 points:")
print("  best near-unit pair count over random restarts (tol 0.01 / 0.003)")
print("="*78)
best=(0,0);
for s in range(40):
    loose,tight = relax_count(21, s, iters=2500)
    if loose>best[0]: best=(loose,tight,s)
print(f"  best over 40 restarts: {best[0]} pairs at tol 0.01 (tight {best[1]} at 0.003), seed {best[2]}")
print(f"  triangular lattice gives 47; proven max u(21)=57.")
print("  (relaxation is a heuristic witness, NOT an exact construction.)")
