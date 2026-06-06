"""
monad-explorer-2026-06-06-S703b
===============================
ROBUSTNESS of the triangular finite-crossover minimizer (HYP-2267) to PATCH SHAPE.

HYP-2267's honest caveat: the crossover N_cross was measured with disk patches
{Q(g)<=T}. Does triangular still win under DIFFERENT patch shapes?  We re-measure
N_cross (smallest N with U>3N, best radius) under three nested-patch families:
  (disk)  P = {g : Q(g) <= T}                      -- Q-metric ball (S703 default)
  (box)   P = {(i,j) : max(|i|,|j|) <= R}          -- index-box (sheared in plane)
  (eucl)  P = {(i,j) : Euclidean ||embed(g)|| <= R}-- TRUE round disk in the plane

The Euclidean patch is the geometry-fair one: it embeds each form as an actual
2D lattice (Cholesky of the Gram matrix) and takes a genuine round disk, so no
lattice is favoured by an axis-aligned patch.  If triangular wins under ALL
three, the minimizer claim is shape-robust.
"""
import math
from collections import defaultdict

FORMS = {
    "triangular(-3)": (1, 1, 1),
    "square(-4)":     (1, 0, 1),
    "disc-12":        (1, 0, 3),
    "disc-15(d5)":    (1, 1, 4),
    "disc-7":         (1, 1, 2),
    "disc-15(2,1,2)": (2, 1, 2),
}

def rep_counts(a, b, c, Rmax):
    cnt = defaultdict(int)
    B = int(math.isqrt(Rmax // a)) + 3
    By = int(math.isqrt(Rmax // c)) + 3
    for x in range(-B, B + 1):
        for y in range(-By, By + 1):
            R = a*x*x + b*x*y + c*y*y
            if 0 < R <= Rmax:
                cnt[R] += 1
    return cnt

def offsets_norm(a, b, c, D):
    offs = []
    B = int(math.isqrt(D // a)) + 2
    By = int(math.isqrt(D // c)) + 2
    for i in range(-B, B + 1):
        for j in range(-By, By + 1):
            if a*i*i + b*i*j + c*j*j == D:
                offs.append((i, j))
    return offs

def embed(a, b, c):
    """Cholesky: lattice basis vectors e1,e2 with Gram [[2a? no]] ... use Q=a x^2+b xy+c y^2.
    Gram G = [[a, b/2],[b/2, c]] so that (x,y) G (x,y)^T = a x^2 + b xy + c y^2."""
    e1 = (math.sqrt(a), 0.0)
    e2 = (b/2/math.sqrt(a), math.sqrt(c - (b/2)**2/a))
    return e1, e2

def count_unit(pts_set, offsets):
    u = 0
    for (i, j) in pts_set:
        for (di, dj) in offsets:
            if (i+di, j+dj) in pts_set:
                u += 1
    return u // 2

def crossover_disk(a, b, c, D, ratio, Tmax):
    offs = offsets_norm(a, b, c, D)
    for T in sorted(set(rep_counts(a, b, c, Tmax).keys())):
        B = int(math.isqrt(T // a)) + 2; By = int(math.isqrt(T // c)) + 2
        pts = [(i, j) for i in range(-B, B+1) for j in range(-By, By+1)
               if a*i*i + b*i*j + c*j*j <= T]
        if len(pts) < 7:
            continue
        ps = set(pts)
        if count_unit(ps, offs) > ratio*len(pts):
            return len(pts)
    return None

def crossover_box(a, b, c, D, ratio, Rmax):
    offs = offsets_norm(a, b, c, D)
    for R in range(1, Rmax):
        pts = [(i, j) for i in range(-R, R+1) for j in range(-R, R+1)]
        if len(pts) < 7:
            continue
        ps = set(pts)
        if count_unit(ps, offs) > ratio*len(pts):
            return len(pts)
    return None

def crossover_eucl(a, b, c, D, ratio, R2max):
    """genuine round Euclidean disk in the embedded plane."""
    offs = offsets_norm(a, b, c, D)
    e1, e2 = embed(a, b, c)
    def E2(i, j):  # squared Euclidean norm of embedded (i,j) -- equals Q(i,j)!
        x = i*e1[0] + j*e2[0]; y = i*e1[1] + j*e2[1]
        return x*x + y*y
    # Note Q(i,j) IS the squared Euclidean length, so eucl disk == Q disk == crossover_disk.
    # We keep this to confirm the identity numerically.
    prev = None
    for R2 in range(1, R2max):
        B = int(math.isqrt(int(R2 / a))) + 2; By = int(math.isqrt(int(R2 / c))) + 2
        pts = [(i, j) for i in range(-B, B+1) for j in range(-By, By+1)
               if E2(i, j) <= R2 + 1e-9]
        if len(pts) < 7 or len(pts) == prev:
            prev = len(pts); continue
        prev = len(pts)
        ps = set(pts)
        if count_unit(ps, offs) > ratio*len(pts):
            return len(pts)
    return None

print("="*78)
print("PATCH-SHAPE ROBUSTNESS of N_cross (smallest N with U>3N, best radius)")
print("="*78)
print(f"{'form':>16} {'D*':>4} {'dens':>5} {'disk':>6} {'box':>6} {'eucl':>6}")
rows = []
for name, (a, b, c) in FORMS.items():
    cnt = rep_counts(a, b, c, 400)
    pops = [DD for DD in sorted(cnt) if cnt[DD] > 6][:3]
    best = {"disk": None, "box": None, "eucl": None}
    bestD = None
    for DD in pops:
        nd = crossover_disk(a, b, c, DD, 3.0, 50*DD+400)
        nb = crossover_box(a, b, c, DD, 3.0, 40)
        ne = crossover_eucl(a, b, c, DD, 3.0, 50*DD+400)
        if nd and (best["disk"] is None or nd < best["disk"]):
            best["disk"] = nd; bestD = DD
        if nb and (best["box"] is None or nb < best["box"]):
            best["box"] = nb
        if ne and (best["eucl"] is None or ne < best["eucl"]):
            best["eucl"] = ne
    dens = cnt[bestD]/2 if bestD else 0
    rows.append((name, bestD, dens, best))
    print(f"{name:>16} {str(bestD):>4} {dens:>5.0f} "
          f"{str(best['disk']):>6} {str(best['box']):>6} {str(best['eucl']):>6}")

print()
for shape in ["disk", "box", "eucl"]:
    valid = [(r[0], r[3][shape]) for r in rows if r[3][shape] is not None]
    if valid:
        win = min(valid, key=lambda t: t[1])
        print(f"  minimizer under {shape:>5} patch: {win[0]} at N={win[1]}")
print()
print("Note: for a binary quadratic form Q, Q(i,j) IS the squared Euclidean length of")
print("the embedded lattice vector, so the 'eucl' round disk == the 'disk' Q-ball")
print("(they should agree exactly). The 'box' patch is the genuinely different shape.")
