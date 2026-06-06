#!/usr/bin/env python3
"""cnp_lee_yang_chromatic_roots_s689.py — the chromatic number of the plane through
the Lee-Yang lens: chromatic roots = Lee-Yang/Fisher zeros of the zero-temperature
antiferromagnetic Potts model on unit-distance graphs.

THE LEE-YANG REFORMULATION OF CNP.
The chromatic polynomial P(G,q) is exactly the zero-temperature antiferromagnetic
Potts partition function Z_G(q, v=−1) = Σ_{S⊆E}(−1)^{|S|} q^{c(S)}. Its complex
roots are the chromatic / Lee-Yang-Fisher zeros. By de Bruijn-Erdős (AC),
        χ(ℝ²) = sup over finite unit-distance graphs G of χ(G).
χ(G) = the least integer q with P(G,q) > 0 (and no larger integer root). So:
        χ(ℝ²) = 1 + max{ integer q that is a chromatic root of some UD graph }.
Hence the {5,6,7} question is a Lee-Yang-zero question:
   • χ ≥ 5  (de Grey 2018):  q = 4 IS a chromatic root of some UD graph.   [KNOWN]
   • χ ≥ 6 :  q = 5 is a chromatic root of some UD graph.                   [OPEN]
   • χ ≥ 7 :  q = 6 is a chromatic root of some UD graph.                   [OPEN]
   • χ ≤ 6  (ELIMINATE 7):  q = 6 is NEVER a chromatic root — i.e. q=6 lies in the
     zero-free region of the entire UD-graph family (P(G,6) > 0 ∀ UD G).   [OPEN]
So 'eliminating 7' = a Lee-Yang ZERO-FREENESS statement at q=6 for unit-distance
graphs; 'eliminating 5' = a zero-FORCING statement (a UD graph with a root at q=5).

THE LEE-YANG CIRCLE OF THE BASIC HN GADGET (verified below).
The minimal Eisenstein gadget is the wheel W6 = unit hexagon + center (χ=3).
P(W_n,q) = q·P(C_n,q−1) = q[(q−2)^n + (−1)^n(q−2)] = q(q−2)[(q−2)^{n−1} + (−1)^{n+1}].
Its 'rim' roots solve (q−2)^{n−1} = ±1 ⟹ |q−2| = 1: the chromatic/Lee-Yang zeros of
EVERY wheel lie on the CIRCLE |q−2|=1 (an antiferromagnetic Lee-Yang circle centred
at q=2). So the canonical HN gadget's Lee-Yang zeros are confined to a circle — and
the chromatic number is read off where that zero-locus crosses the real axis.

This file: builds K3, W6, a triangular patch (χ=3), the Moser spindle (χ=4); gets
their exact chromatic polynomials (Whitney rank sum, exact integers) and chromatic
roots; verifies the W6 Lee-Yang circle; tracks the largest real root and the root
cloud as gadgets force higher χ.

Session: claude-2026-06-06-S689 (cnp-lee-yang-chromatic-roots)."""
import sys; sys.stdout.reconfigure(line_buffering=True)
import math, cmath
import numpy as np

# ---------- geometry / graph builders ----------
def unit_edges(pts, tol=1e-9):
    n = len(pts); E = []
    for i in range(n):
        for j in range(i+1, n):
            if abs(abs(pts[i]-pts[j]) - 1.0) < tol: E.append((i, j))
    return E

def W6_pts():
    return [cmath.exp(2j*math.pi*k/6) for k in range(6)] + [0+0j]

def triangular_patch(R):
    """Eisenstein points a+bζ6 with |a|,|b|≤R — a triangular-lattice patch (χ=3)."""
    z6 = cmath.exp(1j*math.pi/3)
    pts = []
    for a in range(-R, R+1):
        for b in range(-R, R+1):
            pts.append(a + b*z6)
    return pts

def moser_pts():
    z6 = cmath.exp(1j*math.pi/3)
    A = [0+0j, 1+0j, z6, 1+z6]
    theta = 2*math.asin(1/(2*math.sqrt(3))); rot = cmath.exp(1j*theta)
    B = [rot*p for p in A]
    return [A[0], A[1], A[2], A[3], B[1], B[2], B[3]]

# ---------- exact chromatic polynomial via Whitney rank (subset) sum ----------
def chromatic_poly(n, edges):
    """P(G,q) = Σ_{S⊆E}(−1)^{|S|} q^{components(S)}. Returns coeff list low→high."""
    m = len(edges)
    assert m <= 24, "too many edges for subset method"
    coeff = [0]*(n+1)
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for S in range(1 << m):
        for i in range(n): parent[i] = i
        comps = n; bits = S; k = 0; idx = 0
        while bits:
            if bits & 1:
                a, b = edges[idx]; ra, rb = find(a), find(b)
                if ra != rb: parent[ra] = rb; comps -= 1
                k += 1
            bits >>= 1; idx += 1
        coeff[comps] += (-1)**k
    return coeff  # coeff[c] * q^c

def evalP(coeff, q):
    return sum(c*q**i for i, c in enumerate(coeff))

def chromatic_number(coeff, n):
    for q in range(1, n+1):
        if evalP(coeff, q) > 0: return q
    return None

def roots_of(coeff):
    # numpy wants high→low
    hi2lo = coeff[::-1]
    # strip leading zeros
    while len(hi2lo) > 1 and hi2lo[0] == 0: hi2lo = hi2lo[1:]
    return np.roots(hi2lo)

# ---------- run over the gadget tower ----------
def analyze(name, pts):
    n = len(pts); edges = unit_edges(pts)
    coeff = chromatic_poly(n, edges)
    chi = chromatic_number(coeff, n)
    r = roots_of(coeff)
    realroots = sorted([x.real for x in r if abs(x.imag) < 1e-7])
    largest_real = max(realroots) if realroots else float('nan')
    # distances of nonzero roots from q=2 (Lee-Yang circle test)
    nz = [x for x in r if abs(x) > 1e-7 and abs(x-2) > 1e-7]
    d2 = [abs(x-2) for x in nz]
    pq = {q: evalP(coeff, q) for q in (3, 4, 5, 6, 7)}
    print(f"\n--- {name}: {n} vertices, {len(edges)} unit edges ---")
    print(f"  χ = {chi};  P(3)={pq[3]}, P(4)={pq[4]}, P(5)={pq[5]}, P(6)={pq[6]}, P(7)={pq[7]}")
    print(f"  real chromatic roots = {[round(x,4) for x in realroots]};  largest real root = {largest_real:.4f}")
    if d2:
        print(f"  |root − 2| for the {len(d2)} non-degenerate roots: "
              f"min={min(d2):.4f}, max={max(d2):.4f}, mean={sum(d2)/len(d2):.4f}")
    return name, chi, largest_real, d2

print(__doc__.split("Session:")[0])
print("="*78)
res = []
res.append(analyze("K3 (triangle)", [0+0j, 1+0j, cmath.exp(1j*math.pi/3)]))
res.append(analyze("W6 = hexagon+center (basic HN gadget)", W6_pts()))
res.append(analyze("triangular patch R=1 (Eisenstein, χ=3)", triangular_patch(1)))
res.append(analyze("Moser spindle (χ=4, needs √−11)", moser_pts()))

# ---------- verify the W6 Lee-Yang circle exactly ----------
print("\n" + "="*78)
print("LEE-YANG CIRCLE CHECK (W6): rim roots solve (q−2)^5 = −1 ⟹ |q−2| = 1")
n6, e6 = 7, unit_edges(W6_pts()); c6 = chromatic_poly(n6, e6)
r6 = roots_of(c6)
for x in sorted(r6, key=lambda z: (round(z.real,6), round(z.imag,6))):
    on = "ON |q−2|=1" if abs(abs(x-2)-1) < 1e-6 else ("q=0/center" if abs(x)<1e-6 or abs(x-2)<1e-6 else "off")
    print(f"   q = {x.real:+.5f}{x.imag:+.5f}i   |q−2| = {abs(x-2):.5f}   [{on}]")
print("   ⟹ the 5 rim roots lie exactly on the circle |q−2|=1; the real one (q=1) is where")
print("     the Lee-Yang locus meets the real axis below χ=3. The HN gadget's zeros are circular.")

print("\n" + "="*78)
print("THE TREND (Lee-Yang zeros marching toward the integers):")
for (name, chi, lr, d2) in res:
    print(f"   χ={chi}: {name:<42} largest real root = {lr:.4f}")
print("""
 READING: χ(ℝ²) = 1 + (largest integer chromatic root attainable by a UD graph).
 W6/K3/patch (χ=3): largest integer root = 2.  Moser (χ=4): integer root q=3 (P(3)=0).
 de Grey's graph (χ=5, not built here — 1581 vtx): integer root q=4 (P(G,4)=0).
 ELIMINATE 7  ⟺  q=6 is in the UD-family Lee-Yang ZERO-FREE region (P(G,6)>0 ∀ UD G).
 ELIMINATE 5  ⟺  some UD graph FORCES a chromatic root at q=5 (P(G,5)=0).
 Both are open; the Lee-Yang lens says: it is about whether the chromatic-zero locus
 of the UD family pinches the real axis at q=5 and q=6.""")
