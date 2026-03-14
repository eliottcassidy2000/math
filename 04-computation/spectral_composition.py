#!/usr/bin/env python3
"""
spectral_composition.py — opus-2026-03-14-S80

What happens when we compose the tournament spectrum {2,3} with the
Petersen spectrum {3,1,-2}? What algebraic objects emerge?

We explore:
1. The tensor product L(tournament) ⊗ L(Petersen) — Kronecker product of spectra
2. Spectral zeta functions ζ_G(s) = Σ λ_i^{-s}
3. The Ihara zeta function of the Petersen graph
4. Heat kernel traces: Tr(e^{-tA}) and their relation to theta functions
5. Characteristic polynomial composition: what is f(g(z))?
6. Ramanujan property and spectral gaps
7. The "spectral tournament polynomial" — evaluating f at graph eigenvalues
"""

import numpy as np
from math import factorial, comb, pi, e, log, sqrt
from fractions import Fraction

def section(title, n):
    print(f"\n{'='*70}")
    print(f"{n}. {title}")
    print(f"{'='*70}\n")

# ============================================================
section("SPECTRAL DATA RECAP", 1)
# ============================================================

# Tournament polynomial
KEY1, KEY2 = 2, 3
f_roots = [KEY1, KEY2]
print(f"Tournament polynomial: z²-5z+6 = (z-{KEY1})(z-{KEY2})")
print(f"  Roots: {f_roots}")

# Petersen graph adjacency matrix eigenvalues
pet_eigenvalues = [(3, 1), (1, 5), (-2, 4)]  # (eigenvalue, multiplicity)
print(f"\nPetersen adjacency spectrum:")
for ev, mult in pet_eigenvalues:
    print(f"  λ = {ev:+d}, multiplicity = {mult}")

# ============================================================
section("TOURNAMENT POLYNOMIAL AT PETERSEN EIGENVALUES", 2)
# ============================================================

def f_tournament(z):
    return z**2 - 5*z + 6

print("f(z) = z²-5z+6 evaluated at Petersen eigenvalues:")
for ev, mult in pet_eigenvalues:
    val = f_tournament(ev)
    print(f"  f({ev:+d}) = {val:4d}  (mult={mult}, contribution = {mult*val})")

total = sum(mult * f_tournament(ev) for ev, mult in pet_eigenvalues)
print(f"\n  Tr(f(A_P)) = Σ mult·f(λ) = {total}")
print(f"  = f(3)·1 + f(1)·5 + f(-2)·4")
print(f"  = {f_tournament(3)} + {5*f_tournament(1)} + {4*f_tournament(-2)}")
print(f"  = 0 + 10 + 80 = {total}")
print(f"\n  Note: f(3)=0 because 3=KEY₂ is a root!")
print(f"  So the dominant eigenvalue contributes NOTHING.")
print(f"  The trace is entirely from the 'error' eigenvalues.")

# ============================================================
section("PETERSEN MINIMAL POLY AT TOURNAMENT ROOTS", 3)
# ============================================================

def p_petersen(z):
    return z**3 - 2*z**2 - 5*z + 6

print("p(z) = z³-2z²-5z+6 evaluated at tournament polynomial roots:")
print(f"  p(KEY₁) = p(2) = {p_petersen(2)}")
print(f"  p(KEY₂) = p(3) = {p_petersen(3)}")
print(f"  p(KEY₁+KEY₂) = p(5) = {p_petersen(5)}")
print(f"  p(KEY₁·KEY₂) = p(6) = {p_petersen(6)}")
print()
print(f"  p(2) = 8-8-10+6 = {p_petersen(2)}")
print(f"  p(3) = 27-18-15+6 = {p_petersen(3)}")
print(f"  p(2) = {p_petersen(2)} ← NEGATIVE! (2 is between roots 1 and 3)")
print(f"  p(3) = {p_petersen(3)} ← ZERO! (3 IS a root)")
print()
print("The shared root 3 = KEY₂ is the bridge between the two polynomials.")
print(f"  gcd of z²-5z+6 and z³-2z²-5z+6:")
print(f"  z²-5z+6 = (z-2)(z-3)")
print(f"  z³-2z²-5z+6 = (z-3)(z-1)(z+2)")
print(f"  GCD = (z-3)")

# ============================================================
section("COMPOSITION OF POLYNOMIALS", 4)
# ============================================================

print("What is f(p(z))? Tournament polynomial composed with Petersen:")
print("  f(p(z)) = p(z)² - 5·p(z) + 6")
print("  = [p(z) - 2][p(z) - 3]")
print("  = [z³-2z²-5z+4][z³-2z²-5z+3]")
print()

# Compute at integer points
print("  f(p(z)) at small integers:")
for z in range(-3, 6):
    val = f_tournament(p_petersen(z))
    root_str = ""
    if val == 0:
        root_str = " ← ROOT"
    print(f"    z={z:+d}: p(z)={p_petersen(z):6d}, f(p(z))={val:10d}{root_str}")

print("\n  Roots of f(p(z)) = 0: p(z) = 2 or p(z) = 3")
print("  p(z) = 2: z³-2z²-5z+4 = 0")
print("  p(z) = 3: z³-2z²-5z+3 = 0 → (z-3)(z-1)(z+1) = 0? No...")

# Factor p(z) - 3 = z³-2z²-5z+3
# Try z=3: 27-18-15+3 = -3 ≠ 0
# Try z=-1: -1-2+5+3 = 5 ≠ 0

import numpy as np
roots_p2 = np.roots([1, -2, -5, 4])
roots_p3 = np.roots([1, -2, -5, 3])
print(f"\n  Roots of p(z)=2: {np.sort(roots_p2)}")
print(f"  Roots of p(z)=3: {np.sort(roots_p3)}")

print("\n  What about p(f(z))? Petersen composed with tournament:")
print("  p(f(z)) = f(z)³ - 2f(z)² - 5f(z) + 6")
print("  = [f(z)-3][f(z)-1][f(z)+2]")
print("  = [(z-2)(z-3)-3][(z-2)(z-3)-1][(z-2)(z-3)+2]")

for z in range(-3, 8):
    fz = f_tournament(z)
    val = p_petersen(fz)
    root_str = ""
    if val == 0:
        root_str = " ← ROOT"
    print(f"    z={z:+d}: f(z)={fz:4d}, p(f(z))={val:15d}{root_str}")

# ============================================================
section("SPECTRAL ZETA FUNCTIONS", 5)
# ============================================================

print("The spectral zeta function: ζ_G(s) = Σ |λ_i|^{-s}")
print("(only meaningful for non-zero eigenvalues)")
print()

# Petersen: eigenvalues 3 (×1), 1 (×5), -2 (×4)
print("Petersen spectral zeta (with signs):")
for s in range(1, 7):
    zeta = 1 * (3**(-s)) + 5 * (1**(-s)) + 4 * ((-2)**(-s))
    zeta_abs = 1 * (3**(-s)) + 5 * (1**(-s)) + 4 * (2**(-s))
    # Use Fraction for exact
    zeta_f = Fraction(1, 3**s) + Fraction(5, 1) + Fraction(4 * ((-1)**s), 2**s)
    zeta_abs_f = Fraction(1, 3**s) + Fraction(5, 1) + Fraction(4, 2**s)
    print(f"  ζ_P({s}) = 1/3^{s} + 5 + 4·(-2)^{{-{s}}} = {float(zeta_f):.6f}  (= {zeta_f})")

print()
print("Petersen heat kernel trace: Tr(e^{-tA}) = Σ mult·e^{-t·λ}")
print("  = e^{-3t} + 5·e^{-t} + 4·e^{2t}")
for t_val in [0.1, 0.5, 1.0, 2.0]:
    trace = np.exp(-3*t_val) + 5*np.exp(-t_val) + 4*np.exp(2*t_val)
    print(f"  t={t_val:.1f}: Tr = {trace:.6f}")

print()
print("At t=0: Tr = 1+5+4 = 10 = |V(Petersen)| ✓")
print("As t→∞: dominated by 4·e^{2t} (the -2 eigenvalue)")

# ============================================================
section("KRONECKER PRODUCT OF SPECTRA", 6)
# ============================================================

print("Kronecker product T ⊗ P (treat tournament as 2×2 matrix with eigenvalues 2,3):")
print()
print("Eigenvalues of Kronecker product = all products λ_T × λ_P:")
products = []
for lt in [2, 3]:
    for lp, mp in pet_eigenvalues:
        prod = lt * lp
        products.append((prod, mp, lt, lp))
        print(f"  {lt}·{lp:+d} = {prod:+d}  (mult={mp})")

print()
print("Distinct eigenvalues of T⊗P:")
from collections import Counter
eigen_counts = Counter()
for prod, mp, lt, lp in products:
    eigen_counts[prod] += mp
for ev in sorted(eigen_counts.keys()):
    print(f"  λ = {ev:+d}, total multiplicity = {eigen_counts[ev]}")

print()
print(f"Trace of T⊗P = Σ mult·λ = ", end="")
trace = sum(ev * mult for ev, mult in eigen_counts.items())
print(f"{trace}")
print(f"  = Tr(T)·Tr(P) = (2+3)·(3+5-8) = 5·0 = {5*0}... ")
print(f"  Wait: Tr(P) = 3·1 + 1·5 + (-2)·4 = 3+5-8 = 0")
print(f"  So Tr(T⊗P) = Tr(T)·Tr(P) = 5·0 = 0 ✓")

# ============================================================
section("RAMANUJAN PROPERTY", 7)
# ============================================================

print("A k-regular graph is Ramanujan if all non-trivial eigenvalues satisfy")
print("|λ| ≤ 2√(k-1).")
print()
print(f"Petersen: k=3, so bound = 2√2 ≈ {2*sqrt(2):.4f}")
print(f"  Non-trivial eigenvalues: 1, -2")
print(f"  |1| = 1 ≤ {2*sqrt(2):.4f} ✓")
print(f"  |-2| = 2 ≤ {2*sqrt(2):.4f} ✓")
print(f"  Petersen IS Ramanujan!")
print()
print("The Ramanujan bound for k-regular graphs:")
for k in range(2, 12):
    bound = 2*sqrt(k-1)
    print(f"  k={k:2d}: 2√(k-1) = 2√{k-1} = {bound:.4f}")

print()
print("The spectral gap of Petersen: λ₁ - λ₂ = 3 - 1 = 2 = KEY₁")
print("For a 3-regular Ramanujan graph, the BEST possible gap is 3 - 2√2 ≈ 0.17")
print("Petersen has gap 2, far exceeding the Ramanujan minimum.")
print("In fact, gap = KEY₁ and max eigenvalue = KEY₂.")

# ============================================================
section("THE IHARA ZETA FUNCTION", 8)
# ============================================================

print("The Ihara zeta function of a regular graph relates to the adjacency spectrum.")
print()
print("For a q+1 regular graph G with adjacency eigenvalues λ_1,...,λ_n:")
print("  ζ_G(u)^{-1} = (1-u²)^{(q-1)n/2} · Π_i (1 - λ_i u + q u²)")
print()
print("For Petersen (3-regular, so q=2):")
print("  ζ_P(u)^{-1} = (1-u²)^5 · (1-3u+2u²)(1-u+2u²)^5(1+2u+2u²)^4")
print()

# Factor 1 - λu + qu²
print("Individual factors (1 - λu + 2u²):")
for ev, mult in pet_eigenvalues:
    disc = ev**2 - 8
    print(f"  λ={ev:+d}: 1 - ({ev})u + 2u² = 2u² {-ev:+d}u + 1  (disc={disc})")
    if disc >= 0:
        r1 = (ev + sqrt(disc)) / 4
        r2 = (ev - sqrt(disc)) / 4
        print(f"    Roots: u = {r1:.4f}, {r2:.4f}")
    else:
        re = ev / 4
        im = sqrt(-disc) / 4
        print(f"    Roots: u = {re:.4f} ± {im:.4f}i  (complex, |u| = {sqrt(re**2+im**2):.4f})")
        print(f"    |u|² = 1/q = 1/2 ← ON THE RAMANUJAN CIRCLE!")

print()
print("For Ramanujan graphs, all non-trivial zeros lie on |u| = q^{-1/2} = 1/√2")
print(f"  1/√2 = {1/sqrt(2):.6f}")
print(f"  = 1/√KEY₁ = √(1/KEY₁)")
print()
print("This is the graph-theoretic analogue of the Riemann Hypothesis!")
print("The Petersen graph satisfies this GRH analog, as expected for Ramanujan graphs.")

# ============================================================
section("SPECTRAL TOURNAMENT POLYNOMIAL: f(A) WHERE A = ADJACENCY", 9)
# ============================================================

# Build Petersen adjacency matrix
edges = [(0,1),(0,4),(0,5),(1,2),(1,6),(2,3),(2,7),(3,4),(3,8),(4,9),(5,7),(5,8),(6,8),(6,9),(7,9)]
A = np.zeros((10,10), dtype=int)
for i,j in edges:
    A[i][j] = A[j][i] = 1

# Compute f(A) = A² - 5A + 6I
A2 = A @ A
f_A = A2 - 5*A + 6*np.eye(10, dtype=int)

print("f(A_Petersen) = A² - 5A + 6I:")
print()
print("  Eigenvalues of f(A): f(3)=0, f(1)=2, f(-2)=20")
print(f"  So f(A) has eigenvalues 0 (mult 1), 2 (mult 5), 20 (mult 4)")
print()
print(f"  Tr(f(A)) = 0·1 + 2·5 + 20·4 = {0+10+80}")
print(f"  det(f(A)) = 0^1 · 2^5 · 20^4 = 0  (singular!)")
print()

# Check rank
rank_fA = np.linalg.matrix_rank(f_A)
print(f"  rank(f(A)) = {rank_fA}")
print(f"  nullity(f(A)) = {10 - rank_fA}")
print(f"  The null space is the eigenspace of λ=3 (the KEY₂ eigenspace)")
print(f"  = span of the all-ones vector (1,1,...,1)")
print()

# What does f(A) look like?
print("  f(A) matrix (should have eigenvalues 0,2,20):")
for row in f_A:
    print("   ", "".join(f"{int(x):4d}" for x in row))

print()
print(f"  Diagonal entries of f(A): {[int(f_A[i][i]) for i in range(10)]}")
print(f"  All diagonal = {int(f_A[0][0])} = f(0)+d = 6+{int(f_A[0][0]-6)}... ")
print(f"  Actually: (A²)_{{ii}} = degree = 3, A_{{ii}} = 0")
print(f"  f(A)_{{ii}} = 3 - 0 + 6 = 9 = KEY₂² ✓")
print()
print("  Off-diagonal: check (i,j) adjacent vs non-adjacent:")
adj_val = int(f_A[0][1])  # 0-1 are adjacent
non_val = int(f_A[0][2])  # 0-2 are non-adjacent
non_val2 = int(f_A[0][3])  # 0-3 non-adjacent but at distance 2
print(f"  f(A) when i~j: {adj_val}")
print(f"  f(A) when i not~j: {non_val}")
print(f"  (checking another non-adj pair: f(A)_03 = {non_val2})")

# For adjacent i,j: (A²)_{ij} = #common neighbors, A_{ij}=1
# For Petersen: common neighbors of adjacent vertices = 0 (girth 5!)
# So f(A)_{ij} = 0 - 5·1 + 0 = -5 for adjacent
# For non-adjacent: (A²)_{ij} = #common neighbors = 1 (strongly regular)
# f(A)_{ij} = 1 - 5·0 + 0 = 1 for non-adjacent

print()
print("  BEAUTIFUL STRUCTURE:")
print(f"  f(A)_{{i,i}} = {int(f_A[0][0])} = KEY₂²")
print(f"  f(A)_{{i,j}} = {adj_val} = -(KEY₁+KEY₂) when i~j")
print(f"  f(A)_{{i,j}} = {non_val} = 1 when i≁j")
print()
print("  f(A) = 9I - 5A + (J-I-A)... wait")
print(f"  f(A) = A² - 5A + 6I")
print(f"  For Petersen (strongly regular (10,3,0,1)): A² = 3I + 0·A + 1·(J-I-A)")
print(f"  A² = 3I + J - I - A = 2I - A + J")
print(f"  f(A) = (2I - A + J) - 5A + 6I = 8I - 6A + J")
print(f"  = rank(E₈)·I - h(G₂)·A + J")
print(f"  = 8I - 6A + J")
print()
print("  Verify: 8I-6A+J has eigenvalues:")
print(f"  On all-1 vector: 8-6·3+10 = 8-18+10 = 0 ✓ (matches f(3)=0)")
print(f"  On 1-eigenvectors: 8-6·1+0 = 2 ✓ (matches f(1)=2)")
print(f"  On (-2)-eigenvectors: 8-6·(-2)+0 = 20 ✓ (matches f(-2)=20)")
print()
print("  *** f(A_Petersen) = rank(E₈)·I - h(G₂)·A + J ***")
print("  The tournament polynomial of the Petersen adjacency decomposes into")
print("  exactly three Lie-meaningful terms!")

print()
print("="*70)
print("GRAND SYNTHESIS: SPECTRAL COMPOSITION")
print("="*70)
print()
print("1. f(A_P) = 8I - 6A + J = rank(E₈)·I - h(G₂)·A + J")
print("2. The null space of f(A_P) is the all-1 eigenspace (graph regularity)")
print("3. Petersen IS Ramanujan: Ihara zeros on |u|=1/√2=1/√KEY₁")
print("4. The spectral gap of Petersen = KEY₁ = 2")
print("5. Tr(f(A_P)) = 90 = 9·10 = KEY₂²·V(Petersen)")
print("6. The tournament polynomial at Petersen eigenvalues: f(3)=0, f(1)=2, f(-2)=20")
print("7. The shared root KEY₂=3 is the bridge: gcd(f,p) = (z-3)")
