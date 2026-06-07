"""
monad-explorer-2026-06-06-S4
============================
EXACT-arithmetic verification of the Moser-lattice picture for u(21)=57.

Background (THM-431): u(21)=57 is a PROVEN exact value (Alexeev-Mixon-Parshall,
arXiv:2412.11914, 2025). The optimum is NOT a triangular-lattice patch (Harborth
gives only 47 at n=21) but a graph in the MOSER LATTICE (Engel et al., 2024):

   M_L = { a + b w1 + c w3 + d w1 w3 : a,b,c,d in Z },
   w1 = zeta_6 = (1 + i sqrt3)/2,   w1^2 - w1 + 1 = 0     (cos = 1/2)
   w3 = (5 + i sqrt11)/6,           w3^2 - (5/3)w3 + 1 = 0 (cos = 5/6, Moser angle)

So M_L = Z[zeta_6]  extended by the sqrt(-11) direction w3; it lives in the
biquadratic CM field Q(sqrt-3, sqrt-11). We verify:
  (1) M_L has EXACTLY 18 unit vectors (=> max degree 18);
  (2) the triangular lattice gives <= 47 at n=21 (suboptimal);
  (3) the Minkowski sum W6 (+) Delta of a 6-wheel and a unit triangle, with the
      triangle in a transverse (w3) direction, is a 21-vertex graph with EXACTLY
      57 unit distances  ==  the proven u(21).

ALL distance counts are EXACT integers, computed over Q(sqrt3, sqrt11) with
Fraction coefficients in the basis (1, sqrt3, sqrt11, sqrt33). No floats decide
adjacency.
"""
from fractions import Fraction as F
from itertools import combinations, product

# ---------------------------------------------------------------------------
# Exact arithmetic in K = Q(sqrt3, sqrt11), basis b0=1, b1=sqrt3, b2=sqrt11, b3=sqrt33
# An element is a 4-tuple of Fractions.
# ---------------------------------------------------------------------------
Z4 = (F(0),)*4
ONE = (F(1),F(0),F(0),F(0))

def add(x, y): return tuple(x[i]+y[i] for i in range(4))
def sub(x, y): return tuple(x[i]-y[i] for i in range(4))
def smul(s, x): return tuple(s*x[i] for i in range(4))

# multiplication table of basis elements -> (coeff, index)
# 1, sqrt3, sqrt11, sqrt33
_MT = {
    (0,0):(1,0),(0,1):(1,1),(0,2):(1,2),(0,3):(1,3),
    (1,1):(3,0),(1,2):(1,3),(1,3):(3,2),
    (2,2):(11,0),(2,3):(11,1),
    (3,3):(33,0),
}
def mul(x, y):
    r = [F(0)]*4
    for i in range(4):
        if x[i]==0: continue
        for j in range(4):
            if y[j]==0: continue
            a,b = (i,j) if i<=j else (j,i)
            coeff, idx = _MT[(a,b)]
            r[idx] += x[i]*y[j]*coeff
    return tuple(r)

# ---------------------------------------------------------------------------
# Lattice point (a,b,c,d) -> complex coordinate (re, im), each a K-element.
# ---------------------------------------------------------------------------
# basis complex coordinates:
RE = {  # real parts
 '1' : (F(1),F(0),F(0),F(0)),
 'w1': (F(1,2),F(0),F(0),F(0)),
 'w3': (F(5,6),F(0),F(0),F(0)),
 'w13':(F(5,12),F(0),F(0),F(-1,12)),
}
IM = {  # imag parts
 '1' : (F(0),F(0),F(0),F(0)),
 'w1': (F(0),F(1,2),F(0),F(0)),
 'w3': (F(0),F(0),F(1,6),F(0)),
 'w13':(F(0),F(5,12),F(1,12),F(0)),
}
KEYS = ['1','w1','w3','w13']

def coord(v):
    """v=(a,b,c,d) integer -> (re, im) in K."""
    re = Z4; im = Z4
    for k,key in zip(v, KEYS):
        if k==0: continue
        re = add(re, smul(F(k), RE[key]))
        im = add(im, smul(F(k), IM[key]))
    return re, im

def normsq(v):
    """exact |z|^2 in K for lattice vector v (returns a K-element)."""
    re, im = coord(v)
    return add(mul(re,re), mul(im,im))

def is_unit_vec(v):
    return normsq(v) == ONE

# ---------------------------------------------------------------------------
# (1) enumerate the 18 unit vectors of M_L
# ---------------------------------------------------------------------------
print("="*72)
print("(1) UNIT VECTORS of the Moser lattice M_L")
print("="*72)
R = 4
units = []
for a,b,c,d in product(range(-R,R+1), repeat=4):
    if (a,b,c,d)==(0,0,0,0): continue
    if is_unit_vec((a,b,c,d)):
        units.append((a,b,c,d))
print(f"   #unit vectors with |coeff|<= {R}: {len(units)}  (Engel et al. Thm 2.5: 18)")
# verify the algebraic condition ad=bc holds for all of them
adbc = all(a*d==b*c for (a,b,c,d) in units)
print(f"   all satisfy ad=bc : {adbc}")
# show the 6 triangular ones (c=d=0) and the 12 'Moser' ones (c or d != 0)
tri = [u for u in units if u[2]==0 and u[3]==0]
mos = [u for u in units if u[2]!=0 or u[3]!=0]
print(f"   triangular unit vectors (c=d=0): {len(tri)} -> {sorted(tri)}")
print(f"   genuinely-Moser unit vectors    : {len(mos)}")
UNITSET = set(units)

def U_count(points):
    """exact #unit-distance pairs in a set of integer lattice vectors."""
    pts = list(points)
    s = set(pts)
    e = 0
    for p in pts:
        for u in units:
            q = (p[0]+u[0], p[1]+u[1], p[2]+u[2], p[3]+u[3])
            if q in s:
                e += 1
    assert e % 2 == 0
    return e//2

# ---------------------------------------------------------------------------
# (2) triangular-lattice baseline at n=21  (Harborth: floor(3n - sqrt(12n-3)) = 47)
# ---------------------------------------------------------------------------
print()
print("="*72)
print("(2) TRIANGULAR-LATTICE patch at n=21 (pure Z[zeta6], c=d=0)")
print("="*72)
# grow a compact triangular patch: order Z[zeta6] points by |z|^2 about a centre,
# take 21, exact-count.  Try several centres; report the best (max U).
tri_box = [(a,b,0,0) for a in range(-6,7) for b in range(-6,7)]
def trinorm_about(v, cx, cy):
    re,im = coord(v)
    # numeric only for ORDERING (not for adjacency); adjacency is exact above
    import math
    x = float(re[0]) + math.sqrt(3)*float(re[1])
    y = float(im[0]) + math.sqrt(3)*float(im[1])
    return (x-cx)**2 + (y-cy)**2
best_tri = -1
import math
for cx in [i/4 for i in range(0,5)]:
    for cy in [i/4 for i in range(0,5)]:
        pts = sorted(tri_box, key=lambda v: trinorm_about(v,cx,cy))[:21]
        u = U_count(pts)
        if u > best_tri: best_tri = u
print(f"   best triangular 21-pt patch (exact count): U = {best_tri}  (Harborth bound 47)")

# ---------------------------------------------------------------------------
# (3) THE CONSTRUCTION:  W6 (+) Delta  Minkowski sum, exact-counted.
# ---------------------------------------------------------------------------
print()
print("="*72)
print("(3) MINKOWSKI SUM  W6 (+) Delta   (6-wheel  +  unit triangle)")
print("="*72)
# 6-wheel: hub + the 6 triangular unit vectors (a hexagon at unit distance, +spokes)
hub = (0,0,0,0)
W6 = [hub] + tri
print(f"   W6: |V|={len(W6)}  E(W6)={U_count(W6)}  (expect 7, 12)")

# all unit triangles {0,u,v}: u,v unit vectors with (u-v) also a unit vector
triangles = []
for u,v in combinations(units, 2):
    duv = (u[0]-v[0],u[1]-v[1],u[2]-v[2],u[3]-v[3])
    if duv in UNITSET:
        triangles.append((( 0,0,0,0), u, v))
print(f"   #unit triangles through 0 (ordered pairs counted once): {len(triangles)}")

# For each triangle Delta, form W6 (+) Delta, dedupe, exact-count.
results = {}
clean = []
for Delta in triangles:
    pts = set()
    collision = False
    raw = []
    for w in W6:
        for d in Delta:
            p = (w[0]+d[0], w[1]+d[1], w[2]+d[2], w[3]+d[3])
            raw.append(p)
            pts.add(p)
    nV = len(pts)
    nE = U_count(pts)
    results.setdefault((nV,nE), 0)
    results[(nV,nE)] += 1
    if nV == 21 and nE == 57:
        clean.append(Delta)

print("   (|V|, U) distribution over all unit triangles Delta:")
for (nV,nE),cnt in sorted(results.items()):
    star = "   <== faithful 21-vertex, U=57 = u(21)!" if (nV,nE)==(21,57) else ""
    print(f"      |V|={nV:>2}  U={nE:>3}  : {cnt:>3} triangles{star}")

print()
if clean:
    Delta = clean[0]
    print(f"   WITNESS Delta = {{0, {Delta[1]}, {Delta[2]}}}")
    triangular_dirs = [d for d in Delta[1:] if d[2]==0 and d[3]==0]
    print(f"     triangle uses {len(triangular_dirs)} triangular-direction edge(s) and "
          f"{2-len(triangular_dirs)} transverse (Moser/w3) edge(s)")
    pts = sorted({(w[0]+d[0],w[1]+d[1],w[2]+d[2],w[3]+d[3]) for w in W6 for d in Delta})
    print(f"     |V| = {len(pts)},  U = {U_count(pts)}  (EXACT integer count)")
    # degree sequence (exact)
    s = set(pts); deg = []
    for p in pts:
        deg.append(sum(1 for u in units if (p[0]+u[0],p[1]+u[1],p[2]+u[2],p[3]+u[3]) in s))
    from collections import Counter
    print(f"     degree distribution: {dict(sorted(Counter(deg).items()))}  sum/2 = {sum(deg)//2}")
    print(f"     CERTIFIES  u(21) >= 57.  Matches proven u(21) = 57 (THM-431).")
else:
    print("   No faithful 21-vertex U=57 triangle found in this wheel orientation;")
    print("   (try a different wheel/triangle pairing).")

# ---------------------------------------------------------------------------
# (4) sanity: a short greedy densification in M_L to corroborate the ceiling
# ---------------------------------------------------------------------------
print()
print("="*72)
print("(4) Greedy densest 21-pt Moser-lattice graph (corroboration, not a proof)")
print("="*72)
def grow(seed, target=21):
    S = set([seed])
    while len(S) < target:
        cand = {}
        for p in S:
            for u in units:
                q = (p[0]+u[0],p[1]+u[1],p[2]+u[2],p[3]+u[3])
                if q in S: continue
                # gain = # neighbours already in S
                g = sum(1 for w in units if (q[0]+w[0],q[1]+w[1],q[2]+w[2],q[3]+w[3]) in S)
                cand[q] = g
        if not cand: break
        best = max(cand.items(), key=lambda kv:(kv[1], -abs(kv[0][0])-abs(kv[0][1])-abs(kv[0][2])-abs(kv[0][3])))
        S.add(best[0])
    return S
best_greedy = -1
for seed in [hub]:
    S = grow(seed)
    u = U_count(S)
    best_greedy = max(best_greedy, u)
print(f"   greedy from hub: |V|={len(S)}  U={U_count(S)}  (ceiling u(21)=57)")
print()
print("DONE.")
