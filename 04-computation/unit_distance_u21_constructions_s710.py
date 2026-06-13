#!/usr/bin/env python3
"""
u(21): the Erdos unit-distance MAXIMUM at N=21 (and neighbours N=17..24).
monad-explorer-2026-06-06-S710.

Distinct from THM-421 (the "first N to beat 3N" floor). Here we want the exact
maximum number of unit distances u(N), pinned by:
  - best explicit construction (exact-integer LOWER bound on u(N)),
  - the cherry/common-neighbour upper bound (exact-integer UPPER bound),
  - comparison with Harborth's triangular-lattice penny number 3N - sqrt(12N-3).

ALL distance counts are EXACT INTEGER. Patch SELECTION uses float ordering only
(ordering does not affect the exactness of the final integer edge count).

Eisenstein coords (a,b) <-> a*(1,0) + b*(1/2, sqrt3/2).
Squared Euclidean distance of difference (da,db): Q(da,db) = da^2 + da*db + db^2.
"unit^2 = D": connect pairs whose difference has Q = D.
"""
import math, itertools

def Q(da, db):
    return da*da + da*db + db*db

# ---- generate Eisenstein integer box ----
def box(R):
    pts = []
    for a in range(-R, R+1):
        for b in range(-R, R+1):
            pts.append((a, b))
    return pts

# float embedding for ORDERING ONLY
def emb(a, b):
    return (a + 0.5*b, (math.sqrt(3)/2)*b)

def count_edges_exact(subset, D):
    """EXACT integer count of pairs in subset with Q(diff)=D."""
    S = set(subset)
    cnt = 0
    for (a, b) in subset:
        # iterate the offset shell of norm D would be faster, but brute O(n^2) is fine & robust
        pass
    # brute exact
    lst = list(subset)
    n = len(lst)
    for i in range(n):
        ai, bi = lst[i]
        for j in range(i+1, n):
            aj, bj = lst[j]
            if Q(ai-aj, bi-bj) == D:
                cnt += 1
    return cnt

def shell(D, Rmax=12):
    return [(da, db) for da in range(-Rmax, Rmax+1) for db in range(-Rmax, Rmax+1)
            if Q(da, db) == D]

def count_edges_shell(subset, D, Rmax=12):
    """Faster exact count using the offset shell."""
    S = set(subset)
    sh = shell(D, Rmax)
    cnt = 0
    for (a, b) in subset:
        for (da, db) in sh:
            nb = (a+da, b+db)
            if nb in S and (a, b) < nb:   # count each edge once
                cnt += 1
    return cnt

# ---- patch: N points closest (Euclidean) to a rational centre (cx,cy) in embedded plane ----
def patch_closest(N, D, centre, R=14):
    cx, cy = centre
    pts = box(R)
    # order by true squared Euclidean distance to centre (float fine for ordering)
    def key(p):
        x, y = emb(*p)
        return (x-cx)**2 + (y-cy)**2
    pts.sort(key=key)
    sub = pts[:N]
    return sub, count_edges_shell(sub, D)

# Harborth triangular penny max
def harborth(N):
    return math.floor(3*N - math.sqrt(12*N - 3))

# CN/cherry exact-integer upper bound: maximize U=sum d_v/2 s.t. sum d_v(d_v-1) <= N(N-1)
def cn_upper(N):
    cap = N*(N-1)               # = 2*C(N,2)
    # greedy: all degrees as equal as possible. Use integer search around d*.
    # We want max sum d_v with sum d_v(d_v-1) <= cap, d_v>=0 integers, N vertices.
    # Marginal cost of raising a vertex from d to d+1 is (d+1)d-(d)(d-1)=2d. Greedy: always raise
    # the vertex with smallest current degree (lowest marginal cost) -> all near-equal.
    deg = [0]*N
    used = 0
    while True:
        # find vertex with min degree
        m = min(deg)
        idx = deg.index(m)
        marg = 2*deg[idx]        # cost to raise from deg to deg+1
        if used + marg <= cap:
            deg[idx] += 1
            used += marg
        else:
            break
    total_d = sum(deg)
    # U = total_d/2 must be integer; if odd, drop one (can't have odd sum of degrees)
    U = total_d // 2
    return U, deg, used, cap

# centres to try (rational, in embedded plane)
def centres():
    cs = []
    # lattice point
    cs.append(("lattice(0,0)", emb(0,0)))
    # edge midpoint of a unit edge
    cs.append(("edgemid", emb(0,0)[0]*0.5+emb(1,0)[0]*0.5, ))  # placeholder
    # build proper rational centres
    out = [("lattice", emb(0,0))]
    # midpoint of (0,0)-(1,0)
    p0=emb(0,0); p1=emb(1,0); out.append(("edge_10", ((p0[0]+p1[0])/2,(p0[1]+p1[1])/2)))
    # centroid of up-triangle (0,0),(1,0),(0,1)
    p2=emb(0,1); out.append(("tri_centroid", ((p0[0]+p1[0]+p2[0])/3,(p0[1]+p1[1]+p2[1])/3)))
    # centre of hexagon = lattice point already; deep hole = centroid
    # midpoint of (0,0)-(1,1) (a norm-3 long edge midpoint)
    p11=emb(1,1); out.append(("mid_11", ((p0[0]+p11[0])/2,(p0[1]+p11[1])/2)))
    return out

print("="*70)
print("u(N) brackets: best Eisenstein/triangular construction vs upper bounds")
print("="*70)
norms = [1, 3, 4, 7, 9, 12, 13, 16, 21, 28]
for N in range(15, 25):
    best = (-1, None, None)
    for D in norms:
        for cname, c in centres():
            sub, e = patch_closest(N, D, c, R=16)
            if e > best[0]:
                best = (e, D, cname)
    H = harborth(N)
    U, deg, used, cap = cn_upper(N)
    star = " <-- N=21" if N==21 else ""
    print(f"N={N:2d}: best_construction={best[0]:3d} (D={best[1]}, centre={best[2]:12s}) "
          f"| Harborth_tri={H:3d} | CN_upper={U:3d}{star}")
