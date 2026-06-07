#!/usr/bin/env python3
"""
unit_distance_product_defect_monad_s5.py
========================================
monad-explorer-2026-06-07-S5.  Builds on THM-433 (avgdeg additive under [];
N* non-product) and THM-431 (u(21)=57, AMP arXiv:2412.11914).

NEW OBJECT: the PRODUCT-DEFECT function

    delta(N) = u(N) - max over NONTRIVIAL factorizations N=a*b (1<a<=b<N)
                       of  [ b*u(a) + a*u(b) ]            (Erdos product bound)

delta(N) >= 0 ALWAYS (Erdos product construction is a lower bound for u).
  delta(N) = 0  <=>  the GLOBAL optimum on N points is (top-level) a Cartesian
                     product  =>  N is "PRODUCT-OPTIMAL".
  delta(N) > 0  <=>  a genuinely NON-PRODUCT (irreducible) config is strictly
                     denser  =>  N is "IRREDUCIBLE-OPTIMAL".

This QUANTIFIES the irreducibility that THM-433 established only at the crossover
N*.  It turns the binary "N* is non-product" into a full PROFILE over all N<=21
(proven, exact) and exposes a 3-rich / 2-rich split and the principal product
line 3^j.

Also (Part D) a concrete probe of OPEN-Q-057's "is the 28-crosser literally
H(3,3)+1?" question: exactly realize K3^[]3 (27 pts) in Q(sqrt3) and count, for
EVERY plane point, how many of the 27 vertices lie on its unit circle -- i.e.
the max degree an added 28th point can have.  Exact rational arithmetic (a unit
distance is decided by an exact algebraic equality, never a float compare).

All u(n) values for n<=21 are the PROVEN AMP table (cited, not re-proved).
"""

from fractions import Fraction
from itertools import combinations
from math import gcd

# ---------------------------------------------------------------------------
# PROVEN unit-distance maxima  u(n), n = 0..21   (Alexeev-Mixon-Parshall 2024,
# arXiv:2412.11914, Thm 1; matches Schade lower bounds).  OEIS A186705.
# ---------------------------------------------------------------------------
U = [0, 0, 1, 3, 5, 7, 9, 12, 14, 18, 20, 23, 27, 30, 33, 37, 41, 43, 46, 50, 54, 57]
#     0  1  2  3  4  5  6  7   8   9  10  11  12  13  14  15  16  17  18  19  20  21
NMAX = 21

def alpha(n):
    """best average degree on n points = 2 u(n)/n  (Fraction, exact)."""
    return Fraction(2 * U[n], n)

# ---------------------------------------------------------------------------
# Part A: product bound, product-defect delta(N), product-optimality profile
# ---------------------------------------------------------------------------
def product_bound(N):
    """max over nontrivial factorizations a*b=N of b*u(a)+a*u(b); witnesses."""
    best = -1
    wit = []
    for a in range(2, N):
        if N % a:
            continue
        b = N // a
        if a > b:
            continue  # a<=b, count each factorization once
        val = b * U[a] + a * U[b]
        if val > best:
            best = val
            wit = [(a, b, val)]
        elif val == best:
            wit.append((a, b, val))
    return best, wit

print("=" * 78)
print("PART A.  Product-defect delta(N) = u(N) - best Erdos product, N<=21")
print("=" * 78)
print(f"{'N':>3} {'u(N)':>5} {'prodLB':>7} {'delta':>6} {'alpha=2u/N':>11} "
      f"{'>kappa?':>7}  factorization witness(es)")
print("-" * 78)
profile = {}
for N in range(4, NMAX + 1):
    pb, wit = product_bound(N)
    if pb < 0:
        # prime: no nontrivial factorization -> irreducible by definition
        profile[N] = ('prime', None)
        a = alpha(N)
        print(f"{N:>3} {U[N]:>5} {'  --':>7} {'  --':>6} {str(a):>11} "
              f"{'yes' if a > 6 else 'no':>7}  (prime: irreducible by definition)")
        continue
    delta = U[N] - pb
    a = alpha(N)
    kind = 'PRODUCT-opt' if delta == 0 else 'IRREDUC-opt'
    profile[N] = (kind, delta)
    wstr = "; ".join(f"{x}x{y}->{v}" for (x, y, v) in wit)
    print(f"{N:>3} {U[N]:>5} {pb:>7} {delta:>6} {str(a):>11} "
          f"{'yes' if a > 6 else 'no':>7}  {wstr}  [{kind}]")

print()
prod_opt = [N for N in range(4, NMAX + 1) if profile[N][0] == 'PRODUCT-opt']
irr_opt  = [N for N in range(4, NMAX + 1) if profile[N][0] == 'IRREDUC-opt']
primes   = [N for N in range(4, NMAX + 1) if profile[N][0] == 'prime']
print(f"PRODUCT-OPTIMAL (delta=0): {prod_opt}")
print(f"IRREDUCIBLE-OPTIMAL (delta>0): {irr_opt}")
print(f"primes (no product structure): {primes}")

# 3-rich vs 2-rich diagnostic on the COMPOSITE n
def v(n, p):
    e = 0
    while n % p == 0:
        n //= p
        e += 1
    return e
print()
print("3-rich / 2-rich split among COMPOSITE N (does delta=0 track 3-content?):")
print(f"{'N':>3} {'delta':>6} {'v2':>3} {'v3':>3} {'has3':>5} {'verdict':>12}")
for N in prod_opt + irr_opt:
    d = profile[N][1]
    print(f"{N:>3} {d:>6} {v(N,2):>3} {v(N,3):>3} {('yes' if N%3==0 else 'no'):>5} "
          f"{('PRODUCT-opt' if d==0 else 'IRREDUC-opt'):>12}")

# ---------------------------------------------------------------------------
# Part B: superadditivity of alpha over MULTIPLICATION (Erdos = superadditive)
# and the principal product line 3^j: alpha(3^j) = 2j, hits kappa=6 at j=3.
# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("PART B.  alpha(N)=2u(N)/N is SUPERADDITIVE over multiplication:")
print("         alpha(a*b) >= alpha(a)+alpha(b)  (<=> Erdos product bound).")
print("=" * 78)
viol = 0
for a in range(2, NMAX + 1):
    for b in range(a, NMAX // a + 1):
        if a * b > NMAX:
            continue
        lhs = alpha(a * b)
        rhs = alpha(a) + alpha(b)
        tag = "TIGHT" if lhs == rhs else ("ok" if lhs > rhs else "VIOLATION!")
        if lhs < rhs:
            viol += 1
        if a in (2, 3) and a * b <= NMAX:
            print(f"  alpha({a*b:>2}) = {str(lhs):>6}  vs alpha({a})+alpha({b}) "
                  f"= {str(rhs):>6}   [{tag}]")
print(f"  superadditivity violations over all a*b<=21: {viol}  (expect 0)")

print()
print("Principal product line  N = 3^j  (K3^[]j, the j-fold triangle product):")
print(f"{'j':>2} {'N=3^j':>6} {'2j':>4} {'avgdeg K3^[]j':>13} {'note':>20}")
for j in range(1, 5):
    N = 3 ** j
    ad = 2 * j  # j triangle factors, each contributes avgdeg 2
    note = ""
    if N <= NMAX:
        note = f"u({N})={U[N]}, alpha={alpha(N)}" + (" TIE@kappa" if alpha(N) == 6 else "")
    else:
        note = f"realizable LB u({N})>={ad*N//2}=({ad}/2)*{N}"
    print(f"{j:>2} {N:>6} {ad:>4} {ad:>13} {note:>30}")
print("  => alpha(3^j)=2j crosses kappa=6 EXACTLY at j=3 (N=27): the principal")
print("     product line is TANGENT to the kissing threshold at the cube 3^3.")
print("     N=81=3^4 gives alpha>=8>6 (K3^[]4 is a planar UD graph, 8-regular):")
print("     products DO eventually beat 3N, but only at >=32 (THM-433-E), while")
print("     the IRREDUCIBLE crossover N* is at <=28 -- earlier than any product.")

# ---------------------------------------------------------------------------
# Part D: exact K3^[]3 = H(3,3) realization in Q(sqrt3); co-circular probe.
#   Q(sqrt3): represent numbers as p + q*sqrt3, p,q in Fraction.
# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("PART D.  Exact K3^[]3 in Q(sqrt3): unit count, regularity, and the")
print("         'H(3,3)+1' probe (max # of the 27 pts on a common unit circle).")
print("=" * 78)

class Q3:
    """element a + b*sqrt(3), exact."""
    __slots__ = ('a', 'b')
    def __init__(self, a=0, b=0):
        self.a = Fraction(a); self.b = Fraction(b)
    def __add__(s, o): return Q3(s.a + o.a, s.b + o.b)
    def __sub__(s, o): return Q3(s.a - o.a, s.b - o.b)
    def __mul__(s, o): return Q3(s.a*o.a + 3*s.b*o.b, s.a*o.b + s.b*o.a)
    def inv(s):
        d = s.a*s.a - 3*s.b*s.b          # norm
        return Q3(s.a/d, -s.b/d)
    def __truediv__(s, o): return s * o.inv()
    def __eq__(s, o): return s.a == o.a and s.b == o.b
    def __hash__(s): return hash((s.a, s.b))
    def __repr__(s): return f"({s.a}+{s.b}v3)"

# Three unit triangles at three DISTINCT directions whose pairwise generic
# angles keep the Minkowski sum faithful (no accidental coincidences).  We use
# three rotations of a unit edge.  To stay in Q(sqrt3) with EXACT unit length we
# use Pythagorean-style unit vectors with rational/sqrt3 coords of squared-len 1:
#   t0 = (1, 0)
#   t1 = (1/2, sqrt3/2)          (60 deg)
#   t2 = (-1/2, sqrt3/2)         (120 deg)
# Each unit TRIANGLE factor i is generated by {0, e_i, f_i} where e_i,f_i are two
# unit vectors at 60 deg so the 3 points are mutually unit apart.
# Triangle factor in direction theta: vertices {0, R(theta)(1,0), R(theta)(1/2,sqrt3/2)}.
def vec(x, y):  # (Q3, Q3)
    return (x, y)

def Q(a, b=0):
    return Q3(a, b)

# A unit triangle has vertices {0, (1,0), (1/2,sqrt3/2)}.  To keep the Minkowski
# sum of THREE triangles FAITHFUL (27 distinct points, only product unit edges),
# the three triangles must be in GENERIC relative orientation -- using 0/60/120
# deg collapses them into ONE Eisenstein lattice (12 pts).  Fix: rotate each
# factor by a distinct PYTHAGOREAN angle (cos,sin rational, cos^2+sin^2=1) so all
# coords stay in Q(sqrt3) yet the lattices are mutually incommensurate.
def rot(c, s, v):
    """rotate vector v=(Q3,Q3) by (cos,sin)=(c,s), c,s rational Fractions."""
    x, y = v
    return (Q(c)*x - Q(s)*y, Q(s)*x + Q(c)*y)

BASE = [(Q(0), Q(0)),
        (Q(1), Q(0)),
        (Q(Fraction(1, 2)), Q(0, Fraction(1, 2)))]   # 0, (1,0), (1/2, sqrt3/2)

# three distinct primitive Pythagorean rotations (generic relative to Eisenstein)
ANGLES = [(Fraction(1), Fraction(0)),        # identity
          (Fraction(3, 5), Fraction(4, 5)),  # 3-4-5
          (Fraction(5, 13), Fraction(12, 13))]  # 5-12-13

def triangle(k):
    c, s = ANGLES[k]
    return [rot(c, s, v) for v in BASE]

# Minkowski sum of the three triangles = 27 points
T = [triangle(0), triangle(1), triangle(2)]
pts = []
for p0 in T[0]:
    for p1 in T[1]:
        for p2 in T[2]:
            x = p0[0] + p1[0] + p2[0]
            y = p0[1] + p1[1] + p2[1]
            pts.append((x, y))
# dedup (exact)
uniq = []
seen = set()
for p in pts:
    key = (p[0].a, p[0].b, p[1].a, p[1].b)
    if key not in seen:
        seen.add(key); uniq.append(p)
pts = uniq
print(f"distinct points realized: {len(pts)}  (want 27)")

def sqdist(p, q):
    dx = p[0] - q[0]; dy = p[1] - q[1]
    return dx*dx + dy*dy   # Q3

ONE = Q(1)
# unit-distance count + degrees
deg = [0]*len(pts)
edges = 0
for i, j in combinations(range(len(pts)), 2):
    if sqdist(pts[i], pts[j]) == ONE:
        edges += 1; deg[i]+=1; deg[j]+=1
from collections import Counter
print(f"unit distances (edges): {edges}  (want 81 = 3*27)")
print(f"degree sequence: {dict(Counter(deg))}  (want all 6 => 6-regular)")
print(f"3*N = {3*len(pts)}; alpha = {Fraction(2*edges,len(pts))}  (TIE at kappa=6)"
      if len(pts)==27 else "  [realization degenerate -- angles not generic]")

# ----- H(3,3)+1 probe: a 28th point p has degree = #{v in pts : |p-v|=1}.
# Equivalent: how many of the 27 points lie on a common unit circle (center p)?
# The candidate centers p that maximize this are intersection points of unit
# circles around the existing points.  We enumerate centers equidistant (=1)
# from PAIRS of existing points and count exact coincidences.  A center unit-
# distant from points S has all of S on its unit circle.  Max over candidate
# centers = best degree of an added 28th point.
# For each pair (i,j) with |pts_i - pts_j| <= 2, the two unit-circle intersection
# points are the candidate centers; we test each against all 27 pts exactly.
# Intersection of unit circles about A,B: midpoint M +- h * perp(unit(AB)),
# h = sqrt(1 - |AB|^2/4).  h^2 in Q(sqrt3) but h itself may be irrational; to
# stay exact we test a candidate center C by the algebraic condition rather than
# constructing C.  Cleaner: a 28th point unit-distant from a SET S exists iff all
# pairwise distances within S are <=2 AND they are concyclic on a unit circle.
# We instead directly bound the achievable degree by the max, over existing
# points c already in the set, of how many points are at distance exactly... no.
#
# Simplest exact, fully rigorous statement we CAN make without irrationals:
#   For any candidate center C, deg(C) = #{v: |C-v|^2 = 1}.  The MAXIMUM possible
#   is bounded by the max number of the 27 points that are pairwise within
#   distance 2 AND lie on a common unit circle.  We compute, for every existing
#   point as a *center proxy is wrong*; instead we report the distribution of
#   pairwise sqdists to show how many points COULD share a unit circle.
print()
print("H(3,3)+1 probe (EXACT, circumcircle method):")
# A 28th point p added at unit distance from a set S<=pts requires S to lie on a
# UNIT circle centered at p.  deg(p)=|S|.  To beat the tie we need deg>=4
# (81+4=85 = the Engel n=28 value).  Candidate centers with deg>=3 are the
# circumcenters of those 3-subsets whose circumradius is EXACTLY 1.  Enumerate all
# triples, keep circumradius^2 == 1, collect the (exact, Q(sqrt3)) circumcenters,
# then for each distinct center count how many of the 27 pts are at sqdist 1.
def circumcenter(A, B, Cc):
    """exact circumcenter in Q(sqrt3); None if collinear."""
    ax, ay = A; bx, by = B; cx, cy = Cc
    # 2(B-A).O = |B|^2-|A|^2 ; 2(C-A).O = |C|^2-|A|^2
    d11 = (bx-ax)*Q(2); d12 = (by-ay)*Q(2); r1 = (bx*bx+by*by)-(ax*ax+ay*ay)
    d21 = (cx-ax)*Q(2); d22 = (cy-ay)*Q(2); r2 = (cx*cx+cy*cy)-(ax*ax+ay*ay)
    det = d11*d22 - d12*d21
    if det == Q(0):
        return None
    ox = (r1*d22 - d12*r2) / det
    oy = (d11*r2 - r1*d21) / det
    return (ox, oy)

existing = set((p[0].a, p[0].b, p[1].a, p[1].b) for p in pts)
centers = {}        # NEW center-key -> degree
best_new = 0        # best degree of a center NOT coinciding with an existing pt
best_any = 0
coincident_deg6 = 0
if len(pts) == 27:
    seen_centers = set()
    for i, j, k in combinations(range(27), 3):
        O = circumcenter(pts[i], pts[j], pts[k])
        if O is None:
            continue
        if sqdist(O, pts[i]) != ONE:
            continue
        ckey = (O[0].a, O[0].b, O[1].a, O[1].b)
        if ckey in seen_centers:
            continue
        seen_centers.add(ckey)
        d = sum(1 for v in pts if sqdist(O, v) == ONE)
        best_any = max(best_any, d)
        if ckey in existing:
            if d == 6:
                coincident_deg6 += 1
            continue   # center IS an existing vertex -> cannot add a NEW point here
        centers[ckey] = d
        best_new = max(best_new, d)
    deg_dist = Counter(centers.values())
    print(f"  unit-circle centers coinciding with an existing vertex, deg 6: "
          f"{coincident_deg6}  (these are the 27 hexagon centers = the points")
    print(f"     themselves -- a NEW 28th point CANNOT sit here)")
    print(f"  NEW (off-vertex) unit-circle centers through >=3 pts: {len(centers)}")
    print(f"  their degree distribution: {dict(sorted(deg_dist.items())) or '{}'}")
    print(f"  => MAX degree of a genuinely NEW added 28th point = {best_new}")
    if best_new >= 4:
        print(f"  => u(28) >= 81 + {best_new} = {81+best_new} via H(3,3)+1 new pt")
    else:
        print(f"  => H(3,3)+1 new pt gives only u(28) >= {81+best_new} (< 85=Engel).")
        print(f"     The generic K3^[]3 has NO off-vertex point unit-distant from >=4")
        print(f"     of its 27 vertices: every full hexagon is centered ON a vertex.")
        print(f"     ==> the n=28 crosser is NOT 'H(3,3) + one boundary point'; it is")
        print(f"     a genuinely DIFFERENT irreducible blob (consistent with THM-433).")
else:
    print("  [skipped: realization not faithful]")
print("DONE.")
