# opus-2026-07-15-S312 -- THM-856 addendum: the cluster dichotomy at radius 7.
# (1) BEAT LOCALIZATION: mu(D_x cap D_{x+d} cap E) concentrates on
#     Beat_d = {u : ||du|| <= 2/13}; for small d the E-restricted overlap is
#     governed by mu(E cap Beat_d), NOT by (4/169)mu(E).
# (2) MULTI-CLUSTER: if the 7 speeds split into >= 2 clusters separated by
#     >= K0, a spanning tree using ONLY inter-cluster edges has all its
#     overlaps near baseline -> Hunter coercive.
# (3) SINGLE-CLUSTER: 7 speeds within a window: (N+d_i)t are a sub-AP with
#     step t; avoid-measure a(t) > 0 off tiny Farey-1/7 windows; uncovered ->
#     integral_E a(t) dt > 0.
from fractions import Fraction
import itertools, math, random

DELTA = Fraction(1, 13)

def safe_set(P):
    ivs = [(Fraction(0), Fraction(1))]
    for q in P:
        bands = [(Fraction(13*k+1, 13*q), Fraction(13*(k+1)-1, 13*q)) for k in range(q)]
        new = []
        for (a, b) in ivs:
            for (c, d) in bands:
                lo, hi = max(a, c), min(b, d)
                if lo < hi: new.append((lo, hi))
        ivs = sorted(new)
    return ivs

def mu(ivs): return sum(b - a for (a, b) in ivs)

def comb_teeth_in(x, a, b):
    w = Fraction(1, 13*x)
    out = []
    for j in range(math.floor((a - w)*x), math.floor((b + w)*x) + 2):
        lo, hi = max(Fraction(j, x) - w, a), min(Fraction(j, x) + w, b)
        if lo < hi: out.append((lo, hi))
    return out

def pair_measure_on(E, x1, x2):
    tot = Fraction(0)
    for (a, b) in E:
        for (lo, hi) in comb_teeth_in(x1, a, b):
            for (l2, h2) in comb_teeth_in(x2, lo, hi): tot += h2 - l2
    return tot

def subtract_comb(ivs, x):
    out = []
    for (a, b) in ivs:
        cur = a
        for (lo, hi) in sorted(comb_teeth_in(x, a, b)):
            if lo > cur: out.append((cur, min(lo, b)))
            cur = max(cur, hi)
            if cur >= b: break
        if cur < b: out.append((cur, b))
    return [iv for iv in out if iv[0] < iv[1]]

P = [1, 2, 3, 4, 5]
E = safe_set(P)
muE = mu(E)
base = Fraction(4, 169) * muE

# (1) beat localization: overlap vs mu(E cap Beat_d)/mu(Beat_d) * 4/169-scaled
print("(1) BEAT LOCALIZATION -- overlap(x, x+d, E) vs beat-window prediction:")
x = 500
for d in (1, 2, 3, 5, 8, 13, 40, 200):
    ov = pair_measure_on(E, x, x + d)
    # Beat_d = {u: ||du|| <= 2/13}; prediction: (2/13)*mu(E cap Beat_d)... the
    # aligned mass: within the beat window teeth align (joint density ~ 2/13
    # of local measure); outside, near-zero. Predicted ~ (2/13)*mu(E cap Beat_d)
    beat = []
    for k in range(d + 1):
        lo, hi = Fraction(k, d) - Fraction(2, 13*d), Fraction(k, d) + Fraction(2, 13*d)
        beat.append((max(lo, Fraction(0)), min(hi, Fraction(1))))
    # intersect with E
    ebeat = Fraction(0)
    for (a, b) in E:
        for (c, dd) in beat:
            lo, hi = max(a, c), min(b, dd)
            if lo < hi: ebeat += hi - lo
    pred = Fraction(2, 13) * ebeat * Fraction(13, 4)  # triangular alignment avg: (1/2)*(4/13)? calibrate
    print(f"   d={d:4d}: overlap = {float(ov):.6f}  baseline(4/169)muE = {float(base):.6f}"
          f"  mu(E cap Beat_d) = {float(ebeat):.5f}  (beat frac of muE: {float(ebeat/muE):.3f})")

# (2) multi-cluster Hunter with inter-cluster star
print("\n(2) MULTI-CLUSTER packets -- Hunter using only inter-cluster edges:")
def hunter_intercluster(xs, clusters):
    Ei_mass = [mu([iv for (a,b) in E for iv in comb_teeth_in(x, a, b)]) for x in xs]
    # build star-forest: connect every vertex to a hub in the other cluster
    hub_a, hub_b = clusters[0][0], clusters[1][0]
    edges = []
    for i in clusters[0]:
        if i != hub_a: edges.append((i, hub_b))
    for i in clusters[1]:
        if i != hub_b: edges.append((i, hub_a))
    edges.append((hub_a, hub_b))
    # that's len = 7-2+1 = 6 edges but (hub_a,hub_b) is intra?? no: hubs in
    # different clusters -> inter ✓; all listed edges are inter-cluster ✓
    tree = sum(pair_measure_on(E, xs[i], xs[j]) for (i, j) in edges[:6])
    return muE - sum(Ei_mass) + tree, edges[:6]

cases = [
    ([200, 201, 202, 203, 900, 901, 902], [[0,1,2,3],[4,5,6]]),
    ([150, 151, 152, 153, 154, 155, 700], [[0,1,2,3,4,5],[6]]),
    ([100, 101, 350, 351, 352, 800, 801], [[0,1],[2,3,4]]),  # 3 clusters: use first two
]
for xs, cl in cases:
    h, edges = hunter_intercluster(xs, cl)
    U = list(E)
    for x in xs: U = subtract_comb(U, x)
    print(f"   xs={xs}: HUNTER(inter-cluster) = {float(h):+.5f} "
          f"actual uncovered = {float(mu(U)):.5f}  "
          f"{'COERCIVE' if h > 0 else 'not coercive'}")

# (3) single-cluster: a(t) positivity and the uncovered integral
print("\n(3) SINGLE-CLUSTER (consecutive) -- avoid-measure a(t) on E:")
def a_of_t(t, k=7):
    # measure of u avoiding union_{i<k} (W - i t), W = [-1/13, 1/13]
    arcs = []
    for i in range(k):
        c = (-i * t) % 1
        arcs.append(((c - Fraction(1,13)) % 1, (c + Fraction(1,13)) % 1))
    # normalize to sorted arcs on circle; compute complement measure
    segs = []
    for (lo, hi) in arcs:
        if lo <= hi: segs.append((lo, hi))
        else: segs.extend([(lo, Fraction(1)), (Fraction(0), hi)])
    segs.sort()
    cov = Fraction(0); cur = Fraction(0)
    for (lo, hi) in segs:
        if lo > cur: pass
        cov += max(Fraction(0), hi - max(lo, cur))
        cur = max(cur, hi)
    return 1 - cov

# integrate a(t) over E by sampling midpoints of a fine partition (exact-ish scan)
tot = Fraction(0); steps = 0
for (a, b) in E:
    nsteps = max(10, int(300 * float(b - a)))
    for s in range(nsteps):
        t = a + (b - a) * Fraction(2*s+1, 2*nsteps)
        tot += a_of_t(t) * (b - a) / nsteps
        steps += 1
print(f"   integral_E a(t) dt ~ {float(tot):.5f} (positive => uncovered floor "
      f"for ALL large N)  [{steps} sample points]")
print(f"   spot a(t): t=1/2: {float(a_of_t(Fraction(1,2))):.4f}, "
      f"t=1/7+1/5000: {float(a_of_t(Fraction(1,7)+Fraction(1,5000))):.4f}, "
      f"t=2/7+1/300: {float(a_of_t(Fraction(2,7)+Fraction(1,300))):.4f}, "
      f"t=5/13: {float(a_of_t(Fraction(5,13))):.4f}")
