#!/usr/bin/env python3
"""hn_moser_field_tower_s687.py — Hadwiger-Nelson through the fundamental theorem
of algebra: the chromatic obstruction of the plane lives OFF the roots-of-unity
(cyclotomic) locus.

SEED (the user's). A degree-n polynomial has n+1 coefficients and n roots (with
multiplicity) in the complex plane ℂ. ℂ = ℝ² = the Hadwiger-Nelson plane. So the
FTA root-map C^{n+1} → ℂ^n is a map INTO the very plane HN colors. This script
follows that seed and finds where the chromatic obstruction sits on the root
locus.

PART 1 — the χ=3 floor is the cyclotomic locus (z^7 − z).
The minimal Eisenstein unit-distance gadget = the regular unit hexagon + its
center: the 6 vertices are the 6th roots of unity (adjacent chord = 2 sin 30° = 1)
and the center (origin) is at unit distance (radius 1) from all six. Its
unit-distance graph is the WHEEL W6 (hub + 6-cycle), χ = 3. Those 7 points are
EXACTLY the roots of
        z^7 − z = z·(z^6 − 1)            [ 6 roots of unity  +  the origin ].
This is the user's "n roots + the constant": z^6−1 (the hexagon) times z (the
center/origin). And over the field F_7, z^7 − z = ∏_{a∈F_7}(z − a) (Fermat), i.e.
its roots are ALL of F_7 = the 7 colors of the hexagonal 7-coloring. So the SAME
polynomial z^7−z gives, over ℂ, the χ=3 gadget, and over F_7, the 7-coloring
palette. 7 = Φ_3(2) = |PG(2,2)| (already HYP-1043/1104) — here it is the FTA count
6+1.

PART 2 — χ=4 (Moser spindle) requires LEAVING the cyclotomic locus.
The Eisenstein lattice ℤ[ζ6] (triangular lattice) is 3-colorable, so EVERY finite
unit-distance graph inside it has χ ≤ 3 — staying on the roots-of-unity lattice
caps you at the χ=3 floor. The Moser spindle (the minimal χ=4 unit-distance graph,
7 vertices, 11 edges) is two Eisenstein unit-rhombi sharing a vertex, the second
rotated by the spindle angle θ. We build it explicitly, verify χ=4, and identify θ:
        cos θ = 5/6      (NOT in {0,±1/2,±1} ⟹ by Niven's theorem θ is an
                          IRRATIONAL multiple of π),
        e^{iθ} = (5 + i√11)/6  is a root of  3z² − 5z + 3,
        Mahler measure M(3z²−5z+3) = 3 > 1 ⟹ by Kronecker's theorem e^{iθ} is NOT
                          a root of unity (a non-cyclotomic point of the unit
                          circle); leading coefficient = 3 = the prime-3 again.
So the spindle's rhombi live in ℚ(√−3) (Eisenstein) and its connecting rotation in
ℚ(√−11): the spindle is realized in the BIQUADRATIC field ℚ(√−3, √−11). Both −3 and
−11 are HEEGNER numbers. χ≥4 cannot be realized in ℚ(√−3) alone (that lattice is
3-colorable) — the SECOND imaginary-quadratic field is necessary.

THE NEEDLE (new HN statement). On the FTA root locus, the chromatic obstruction is
graded by Mahler measure / number field: Mahler-measure-1 (roots of unity, monic
cyclotomic, ℚ(√−3)) forces only χ≥3; forcing χ≥4 requires a non-cyclotomic
rotation (Mahler measure 3, ℚ(√−11)), i.e. a second imaginary-quadratic field. The
'incommensurate rotation' the repo flagged (S699g, 33.56°) is precisely 'this
algebraic point is not a root of unity', and its field is ℚ(√−11). [CONJECTURE:
χ(ℝ²) is bounded below by how many pairwise-incommensurate imaginary-quadratic
rotations one unit-distance graph can force — de Grey's χ≥5 graph would then be a
tower of several such fields.]

PART 3 — the FTA/spectral bridge (roots of unity = eigenvalues).
The circulant unit-distance graph on the 6m-th roots of unity has eigenvalues =
the connection-set polynomial evaluated at the roots of unity (the FTA locus);
m→∞ recovers the plane's Bessel J0 symbol (the repo's χ≥3.48 floor, S699g). The
hexagon C6(1) eigenvalues are 2cos(πj/3)∈{2,1,−1,−2}, Hoffman χ≥2 (tight).

Session: claude-2026-06-06-S687 (hn-moser-field-tower)."""
import sys; sys.stdout.reconfigure(line_buffering=True)
import math, cmath
from itertools import product

# ---------- generic helpers ----------
def unit_distance_graph(pts, tol=1e-9):
    """edges {i,j} with |pts[i]-pts[j]| == 1."""
    n = len(pts); adj = [[] for _ in range(n)]
    edges = []
    for i in range(n):
        for j in range(i+1, n):
            if abs(abs(pts[i]-pts[j]) - 1.0) < tol:
                adj[i].append(j); adj[j].append(i); edges.append((i, j))
    return adj, edges

def chromatic_number(n, edges, kmax=6):
    """exact chromatic number by backtracking (small graphs)."""
    elist = edges
    for k in range(1, kmax+1):
        color = [-1]*n
        def ok(v, c):
            for (a, b) in elist:
                if a == v and color[b] == c: return False
                if b == v and color[a] == c: return False
            return True
        def bt(v):
            if v == n: return True
            for c in range(k):
                if ok(v, c):
                    color[v] = c
                    if bt(v+1): return True
                    color[v] = -1
            return False
        if bt(0): return k
    return None

# ---------- PART 1: z^7 - z = hexagon + center, and F_7 ----------
print("=== PART 1: the χ=3 floor IS the cyclotomic locus  (z^7 − z) ===")
hexagon = [cmath.exp(2j*math.pi*k/6) for k in range(6)]   # 6th roots of unity
center  = [0+0j]
pts = hexagon + center                                    # 7 points
adj, edges = unit_distance_graph(pts)
chi = chromatic_number(len(pts), edges)
print(f"  hexagon+center: {len(pts)} points (6 sixth-roots of unity + origin),"
      f" {len(edges)} unit edges")
deg = [len(a) for a in adj]
print(f"  degree sequence = {sorted(deg, reverse=True)}  (hub deg 6 + 6-cycle deg 3 = wheel W6)")
print(f"  χ(hexagon+center) = {chi}   [W6: even 6-cycle needs 2, hub forces a 3rd]")
# the 7 points are exactly the roots of z^7 - z = z(z^6-1):
roots_C = [0+0j] + [cmath.exp(2j*math.pi*k/6) for k in range(6)]
maxres = max(abs(r**7 - r) for r in roots_C)
print(f"  these 7 points = roots of z^7 − z = z·(z^6 − 1):  max|r^7 − r| = {maxres:.2e}")
# over F_7: z^7 - z vanishes on all of F_7 (Fermat) -> the 7 colors
f7 = all((pow(a, 7, 7) - a) % 7 == 0 for a in range(7))
print(f"  over F_7:  z^7 − z ≡ ∏_(a∈F_7)(z−a)  (z^7≡z ∀a∈F_7): {f7}  ⟹ roots = the 7 colors")
print(f"  so z^7−z gives the χ=3 gadget over ℂ AND the 7-coloring palette over F_7.")
print(f"  7 = 6+1 (FTA: hexagon z^6−1 times the center z) = Φ_3(2) = |PG(2,2)| (HYP-1043/1104).")

# ---------- PART 2: Moser spindle leaves the cyclotomic locus ----------
print("\n=== PART 2: χ=4 (Moser spindle) requires a NON-cyclotomic rotation ===")
# Eisenstein unit-rhombus A: {0, 1, ζ6, 1+ζ6} in ℤ[ζ6]; far vertex 1+ζ6 has |·|=√3.
z6 = cmath.exp(1j*math.pi/3)              # ζ6 = e^{iπ/3}
A = [0+0j, 1+0j, z6, 1+z6]
# spindle angle: two points at distance √3 from origin, separated by θ, at distance 1:
#   2√3 sin(θ/2) = 1  ⟹  sin(θ/2)=1/(2√3)  ⟹  cosθ = 1 − 2 sin²(θ/2) = 5/6
theta = 2*math.asin(1/(2*math.sqrt(3)))
rot = cmath.exp(1j*theta)
B = [rot*p for p in A]                     # rhombus B = A rotated by θ about 0 (shares vertex 0)
# spindle vertices: 0 shared; the far vertices 1+ζ6 (=A[3]) and rot*(1+ζ6) (=B[3]) get a unit edge
pts2 = [A[0], A[1], A[2], A[3], B[1], B[2], B[3]]   # 7 distinct vertices
adj2, edges2 = unit_distance_graph(pts2)
chi2 = chromatic_number(len(pts2), edges2)
far_gap = abs(pts2[3] - pts2[6])           # |A[3] - B[3]| should be 1 (the spindle edge)
print(f"  Moser spindle: {len(pts2)} vertices, {len(edges2)} unit edges; "
      f"far-vertex gap |(1+ζ6) − e^(iθ)(1+ζ6)| = {far_gap:.6f}")
print(f"  χ(Moser spindle) = {chi2}   (minimal χ=4 unit-distance graph)")
print(f"  spindle angle θ = {theta:.6f} rad = {math.degrees(theta):.3f}°,  cos θ = {math.cos(theta):.6f} (= 5/6)")
# Niven: cosθ rational & θ rational multiple of π ⟹ cosθ ∈ {0,±1/2,±1}. 5/6 ∉ ⟹ θ irrational·π.
print(f"  cos θ = 5/6 ∉ {{0,±1/2,±1}}  ⟹  by Niven θ is an IRRATIONAL multiple of π (θ/π ∉ ℚ).")
# e^{iθ} = (5 + i√11)/6, minimal polynomial 3z² − 5z + 3:
eit = (5 + 1j*math.sqrt(11))/6
print(f"  e^(iθ) = (5 + i√11)/6 = {eit:.6f};  check |e^(iθ)−rot| = {abs(eit-rot):.2e}")
print(f"  3·(e^(iθ))² − 5·e^(iθ) + 3 = {abs(3*eit**2 - 5*eit + 3):.2e}  ⟹ e^(iθ) is a root of 3z²−5z+3")
disc = 25 - 4*3*3
print(f"  3z²−5z+3: discriminant = {disc} = −11; Mahler measure = leading_coeff·∏max(1,|root|) = 3·1·1 = 3 > 1")
print(f"  ⟹ e^(iθ) is NOT a root of unity: rigorous by Niven (cos θ=5/6 rational, θ/π irrational);"
      f"\n    and M=3≠1 ⟹ 3z²−5z+3 is non-cyclotomic (Kronecker/Lehmer), a non-cyclotomic unit-circle point.")
print(f"  rhombi ∈ ℚ(√−3) (Eisenstein), rotation ∈ ℚ(√−11) ⟹ spindle lives in ℚ(√−3, √−11);")
print(f"  both −3, −11 are HEEGNER numbers. χ≥4 needs the 2nd field — ℤ[ζ6] alone is 3-colorable.")

# sanity: confirm rhombus A really lies in the Eisenstein lattice ℤ[ζ6] (a+bζ6, a,b∈ℤ)
def in_eisenstein(z, tol=1e-9):
    # z = a + b ζ6 ; ζ6 = (1 + i√3)/2 ⟹ b = 2 Im z/√3, a = Re z − b/2
    b = 2*z.imag/math.sqrt(3); a = z.real - b/2
    return abs(a-round(a)) < tol and abs(b-round(b)) < tol
print(f"  rhombus A ⊂ ℤ[ζ6]? {all(in_eisenstein(p) for p in A)};  "
      f"rotated far vertex B[3] ∈ ℤ[ζ6]? {in_eisenstein(B[3])}  (False = it left the lattice)")

# ---------- PART 3: roots of unity = eigenvalues (FTA/spectral bridge) ----------
print("\n=== PART 3: FTA/spectral bridge — roots of unity ARE the eigenvalues ===")
# unit-distance graph on the 6th roots of unity = circulant C6(1) (step ±1 = 60° = unit chord)
# eigenvalues of circulant with connection {±1}: λ_j = 2 cos(2π j / 6) = poly z+z^{-1} at ζ_6^j
eig = sorted({round(2*math.cos(2*math.pi*j/6), 6) for j in range(6)})
lam_max, lam_min = max(eig), min(eig)
print(f"  C6(1) (hexagon, no center) eigenvalues = connection poly (z + 1/z) at 6th roots of unity")
print(f"     = {{2cos(πj/3)}} = {eig};  λ_max={lam_max}, λ_min={lam_min}")
print(f"  Hoffman: χ ≥ 1 − λ_max/λ_min = 1 − {lam_max}/{lam_min} = {1 - lam_max/lam_min:.0f} (tight: bipartite 6-cycle)")
print(f"  m→∞ (6m-th roots → the full unit circle): connection poly → Bessel J0 = the plane's")
print(f"  symbol (λ_min≈−0.4028 ⟹ χ(ℝ²)≥3.48, S699g/HYP-2264). The FTA root locus IS the spectrum.")
