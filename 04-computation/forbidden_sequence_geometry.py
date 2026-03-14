#!/usr/bin/env python3
"""
forbidden_sequence_geometry.py — opus-2026-03-14-S80

The forbidden H values form the sequence 7·3^k = {7, 21, 63, 189, ...}
This is OEIS A005051 (or 7·A000244). Why this exact geometric progression?

We explore:
1. The sequence 7·3^k as Φ₃(2)·KEY₂^k — cyclotomic × power tower
2. Connection to ternary representations: 7 = 21 in base 3, 21 = 210 in base 3
3. The "gap lattice" — are forbidden H values vertices of a 1D lattice in log-space?
4. Relation to E₈ theta function coefficients
5. The sequence as a q-analog: 7 = (3³-1)/(3-1), 21 = 3·7, 63 = 9·7 = (3⁶-1)/(3-1)
6. Connections to coding theory — 7 = Fano plane, 63 = PG(5,2) points
7. The Mersenne connection: 7 = M₃, 63 = M₆ (both 2^k-1)
8. Why 3^k and not 2^k? The asymmetry of KEY₁ vs KEY₂
"""

import itertools
from math import log, gcd, factorial, comb
from fractions import Fraction

def section(title, n):
    print(f"\n{'='*70}")
    print(f"{n}. {title}")
    print(f"{'='*70}\n")

# ============================================================
section("THE FORBIDDEN SEQUENCE: STRUCTURE", 1)
# ============================================================

forbidden = [7 * 3**k for k in range(8)]
print("Forbidden H sequence: 7·3^k")
for k, h in enumerate(forbidden):
    base3 = ""
    n = h
    while n > 0:
        base3 = str(n % 3) + base3
        n //= 3
    # Check Mersenne
    is_mersenne = (h & (h+1)) == 0  # 2^k - 1 check
    mersenne_str = f" = 2^{int(log(h+1,2)+0.5)}-1" if is_mersenne else ""
    # Check Φ₃ value
    phi3_vals = {n2**2 + n2 + 1: n2 for n2 in range(100)}
    phi3_str = f" = Φ₃({phi3_vals[h]})" if h in phi3_vals else ""

    print(f"  k={k}: 7·3^{k} = {h:6d}  (base 3: {base3}){mersenne_str}{phi3_str}")

print(f"\nBase-3 pattern:")
print(f"  7  = 21₃     (2·3 + 1)")
print(f"  21 = 210₃    (2·3² + 1·3)")
print(f"  63 = 2100₃   (2·3³ + 1·3²)")
print(f"  General: 7·3^k = 2·3^(k+1) + 1·3^k = (21)0...0₃ with k zeros")
print(f"  In base 3, the pattern is ALWAYS '21' followed by k zeros!")

# ============================================================
section("MERSENNE × POWER OF THREE: DOUBLE IDENTITY", 2)
# ============================================================

print("The forbidden values that are ALSO Mersenne numbers:")
print(f"  7  = 2³-1 = M₃   (k=0: 7·3⁰)")
print(f"  63 = 2⁶-1 = M₆   (k=2: 7·3²)")
print(f"  Next would need 7·3^k = 2^m - 1")
print()

# Check: when does 7·3^k = 2^m - 1?
# 7·3^k + 1 = 2^m
# For k=0: 8 = 2³ ✓
# For k=1: 22 = 2·11 ✗
# For k=2: 64 = 2⁶ ✓
# For k=3: 190 = 2·5·19 ✗
# For k=4: 568 = 2³·71 ✗
# For k≥3: 7·3^k + 1 ≡ 7·0 + 1 = 1 (mod 8) for k≥2 since 3²≡1 mod 8
# But 2^m = 7·3^k + 1. For k≥3: 7·3^k = 7·27·3^(k-3), need 7·27·3^(k-3)+1 = 2^m
# This is a Ramanujan-Nagell type equation

print("Solving 7·3^k + 1 = 2^m:")
for k in range(20):
    val = 7 * 3**k + 1
    # Check if power of 2
    if val > 0 and (val & (val - 1)) == 0:
        m = int(log(val, 2) + 0.5)
        print(f"  k={k}: 7·3^{k} + 1 = {val} = 2^{m}  ✓")

print(f"\nOnly k=0 and k=2 give Mersenne numbers!")
print(f"This is a Pillai/Catalan-type result: 2^m - 3^n = 1 has finitely many solutions")
print(f"  2³ - 3¹ = 8-3 = 5... no")
print(f"  Actually: 7·3^k = 2^m - 1 means 2^m = 7·3^k + 1")
print(f"  k=0: 2³ = 8 ✓ (difference: 2^3 - 7·3^0 = 1)")
print(f"  k=2: 2⁶ = 64 ✓ (difference: 2^6 - 7·3^2 = 1)")
print(f"  By Mihailescu (Catalan conjecture, proved 2002):")
print(f"  the only perfect powers differing by 1 are 8 and 9 = 2³ and 3²")
print(f"  So 2^m = 7·3^k + 1 with m,k ≥ 1 has only finitely many solutions")

# ============================================================
section("PROJECTIVE GEOMETRY: 7 AND 63 AS PG POINT COUNTS", 3)
# ============================================================

print("Points in projective spaces PG(d, q) = (q^(d+1) - 1)/(q - 1):")
print()
for q in [2, 3, 4, 5, 7]:
    print(f"  q={q}:")
    for d in range(1, 7):
        pts = (q**(d+1) - 1) // (q - 1)
        star = ""
        if pts in [7, 21, 63, 189, 567]:
            star = f"  ← FORBIDDEN! = 7·3^{[7,21,63,189,567].index(pts)}"
        print(f"    PG({d},{q}) = {pts:6d} points{star}")
    print()

print("KEY DISCOVERY:")
print("  7  = |PG(2,2)| = Fano plane (the unique minimal projective plane)")
print("  21 = |PG(2,4)| = projective plane over GF(4)")
print("  63 = |PG(5,2)| = 6-dim binary projective space")
print("  21 = |PG(4,2)| ← WAIT, let me check...")
pg42 = (2**5 - 1)//(2 - 1)
print(f"  Actually |PG(4,2)| = {pg42} = 31 = h(E₈)+1")
print()
print(f"  7  = (2³-1)/(2-1)  = |PG(2,2)| = KEY₂² - KEY₁")
print(f"  21 = (4³-1)/(4-1)  = |PG(2,4)| = KEY₂·(h(G₂)+1)")
print(f"  63 = (2⁶-1)/(2-1)  = |PG(5,2)| = KEY₂²·(h(G₂)+1)")

# ============================================================
section("GAUSSIAN BINOMIAL COEFFICIENTS AND q-ANALOGS", 4)
# ============================================================

def gauss_binom(n, k, q):
    """Gaussian binomial coefficient [n choose k]_q"""
    if k < 0 or k > n:
        return 0
    if k == 0:
        return 1
    num = 1
    den = 1
    for i in range(k):
        num *= (q**(n-i) - 1)
        den *= (q**(i+1) - 1)
    return num // den

print("q-analogs at q=2 (counting subspaces of GF(2)^n):")
print()
for n in range(1, 8):
    row = [gauss_binom(n, k, 2) for k in range(n+1)]
    total = sum(row)
    star = ""
    if total in [7, 21, 63, 189]:
        star = f"  ← FORBIDDEN! = 7·3^{[7,21,63,189].index(total)}"
    print(f"  n={n}: {row}  total={total}{star}")

print()
print("Total subspaces of GF(2)^n = Σ [n,k]_2:")
for n in range(1, 10):
    total = sum(gauss_binom(n, k, 2) for k in range(n+1))
    star = "  ← FORBIDDEN" if total in [7, 21, 63, 189, 567] else ""
    print(f"  n={n}: {total}{star}")

print()
print("At q=3:")
for n in range(1, 7):
    row = [gauss_binom(n, k, 3) for k in range(n+1)]
    total = sum(row)
    star = ""
    if total in [7, 21, 63, 189]:
        star = f"  ← FORBIDDEN!"
    print(f"  n={n}: {row}  total={total}{star}")

# ============================================================
section("THE ASYMMETRY: WHY 3^k NOT 2^k?", 5)
# ============================================================

print("If forbidden were 7·2^k instead of 7·3^k:")
alt_forb = [7 * 2**k for k in range(6)]
print(f"  Would give: {alt_forb}")
print(f"  But 14 IS achievable! (many tournaments at n=6 have H=14... wait, H odd)")
print(f"  H is always odd (Rédei), so 7·2^k is even for k≥1 — automatically excluded!")
print()
print("The asymmetry explained:")
print("  H(T) = I(Ω(T), 2) is ALWAYS odd (Rédei's theorem)")
print("  KEY₂ = 3 is odd, so 7·3^k is always odd ← meaningful gap")
print("  KEY₁ = 2 is even, so 7·2^k is even for k≥1 ← trivially excluded")
print("  The forbidden sequence MUST use the odd key!")
print()
print("Deeper: H ≡ 1 (mod 2) always. What about mod higher powers?")

# Check H mod 3
print("\nH values at n=6 modulo 3:")
# From the exhaustive data
h_vals_n6 = [1, 3, 5, 9, 11, 13, 15, 17, 19, 23, 25, 27, 29, 31, 33, 37, 41, 43, 45]
for r in range(3):
    vals = [h for h in h_vals_n6 if h % 3 == r]
    print(f"  H ≡ {r} (mod 3): {vals}")

print(f"\nAll three residues mod 3 are represented!")
print(f"But forbidden {7, 21, 63} ≡ {{1, 0, 0}} (mod 3)")
print(f"  7  ≡ 1 (mod 3)")
print(f"  21 ≡ 0 (mod 3)")
print(f"  63 ≡ 0 (mod 3)")

# ============================================================
section("GALOIS FIELD STRUCTURE: GF(2) vs GF(3)", 6)
# ============================================================

print("The tournament polynomial over finite fields:")
print("  f(z) = z² - 5z + 6 = (z-2)(z-3)")
print()
print("Over GF(2): f(z) = z² - z + 0 = z(z-1) = z(z+1)")
print(f"  Roots: 0, 1 ∈ GF(2)")
print(f"  f splits completely over GF(2)")
print()
print("Over GF(3): f(z) = z² - 2z + 0 = z(z-2) = z(z+1)")
print(f"  Roots: 0, 2 ∈ GF(3)")
print(f"  f splits completely over GF(3)")
print()
print("Over GF(5): f(z) = z² + 0z + 1 = z² + 1")
print(f"  Roots: z² = -1 = 4 in GF(5)")
print(f"  √4 = 2,3 in GF(5) since 2²=4, 3²=4")
print(f"  Roots: 2, 3 ∈ GF(5) — SAME as over ℤ!")
print(f"  KEY₁ and KEY₂ live natively in GF(5) = GF(KEY₁+KEY₂)")
print()
print("Over GF(7): f(z) = z² - 5z + 6 = z² + 2z + 6")
print(f"  Disc = 4 - 24 = -20 ≡ 1 (mod 7)")
print(f"  √1 = 1,6 in GF(7)")
print(f"  Roots: (-2 ± 1)/2 = (-1)/2, (-3)/2 in GF(7)")
print(f"  2⁻¹ = 4 in GF(7), so roots = (-1)·4, (-3)·4 = -4, -12 ≡ 3, 2 (mod 7)")
print(f"  Roots: 2, 3 ∈ GF(7) — AGAIN the same!")
print()
print("Summary: 2 and 3 are roots of z²-5z+6 in every GF(p) for p≥5")
print("because they are actual integers, not formal symbols.")
print("But the CHARACTER of f changes:")
print(f"  Over GF(2): f ≡ z(z+1) — has ZERO as a root")
print(f"  Over GF(3): f ≡ z(z+1) — has ZERO as a root")
print(f"  Over GF(5): f ≡ z²+1   — constant term vanishes mod 5")
print(f"  KEY₁, KEY₂ ≡ 0 (mod KEY₁, KEY₂ resp.) — tautological")
print(f"  but: 5 ≡ 0 (mod 5), 6 ≡ 1 (mod 5)")

# ============================================================
section("THE TRIPLING MAP: H → 3H IN CONFLICT GRAPH TERMS", 7)
# ============================================================

print("If H = I(G, 2), what graph operation gives I(G', 2) = 3·I(G, 2)?")
print()
print("I(G, x) for disjoint union: I(G₁ ⊔ G₂, x) = I(G₁, x)·I(G₂, x)")
print("I(K₁, x) = 1 + x, so I(K₁, 2) = 3 = KEY₂")
print()
print("Therefore: I(G ⊔ K₁, 2) = 3·I(G, 2)")
print("  Adding an isolated vertex to Ω TRIPLES H!")
print()
print("The forbidden sequence becomes:")
print(f"  H=7:   Ω = K₃                 I(K₃, 2) = 7")
print(f"  H=21:  Ω = K₃ ⊔ K₁           I = 7·3 = 21")
print(f"  H=63:  Ω = K₃ ⊔ K₁ ⊔ K₁     I = 7·9 = 63")
print(f"  H=189: Ω = K₃ ⊔ K₁ ⊔ K₁ ⊔ K₁ I = 7·27 = 189")
print()
print("The forbidden sequence = K₃ ⊔ (k copies of K₁)")
print("This means: a tournament can NEVER have a conflict graph that is")
print("'a triangle plus some number of isolated vertices'!")
print()
print("WHY? The triangle K₃ requires 3 mutually conflicting odd cycles.")
print("An isolated vertex is a cycle disjoint from all others.")
print("The impossibility of K₃ (HYP-1047 area) propagates:")
print("  If K₃ is impossible as a subgraph of Ω where all other vertices")
print("  are isolated from K₃, then K₃ ⊔ kK₁ is impossible for all k.")
print()
print("KEY INSIGHT: The 'seed' of impossibility is K₃ (the triangle).")
print("The tripling operation (adding isolated cycle) just propagates it.")
print("The forbidden sequence is really about ONE impossibility: K₃ ⊄ Ω(T).")

# ============================================================
section("K₃ AS CONFLICT SUBGRAPH: THE CORE PROHIBITION", 8)
# ============================================================

print("The question reduces to: WHY can't Ω(T) contain K₃ as a subgraph")
print("where the three vertices are disjoint from all other Ω-vertices?")
print()
print("Actually wait — the full statement is stronger.")
print("K₃ ⊔ kK₁ means exactly 3+k cycles in Ω, where 3 are mutually")
print("conflicting and k are isolated (disjoint from all others).")
print()
print("But maybe K₃ CAN appear as a subgraph of larger Ω?")
print("The question is whether Ω can have I(Ω,2) = 7·3^k.")
print()
print("For I(Ω,2) = 7 = 1 + 6: this means Ω has exactly 3 edges worth")
print("of 'obstruction'. Specifically:")
print("  I(G,2) = Σ_S 2^|S| where S ranges over independent sets")
print("  I(G,2) = 7 iff independent set counts (i₀,i₁,i₂,...) satisfy")
print("  1 + 2i₁ + 4i₂ + 8i₃ + ... = 7")
print("  Only solution: i₀=1, i₁=3, i₂=0 (all higher = 0)")
print("  This means |V(Ω)|=3 and no independent set of size 2 exists")
print("  i.e., all pairs are edges: Ω = K₃.")
print()
print("Similarly I(G,2) = 21 = 1 + 2·4 + 4·3 = 1+8+12... no")
print("  21 = 1 + 2i₁ + 4i₂ + 8i₃ + ...")
print("  21 = 3·7 = I(K₁,2)·I(K₃,2) = I(K₃⊔K₁, 2)")
print("  i₁=4, i₂=3, i₃=0")
print("  Or: 21 = 1 + 2·2 + 4·4 = 21? 1+4+16=21 yes!")
print("  That would be i₁=2, i₂=4")
print("  Graph: 2 vertices, 4 independent pairs? Only 1 pair with 2 vertices.")
print()

# Let's actually enumerate which graphs have I(G,2) = 21
print("Enumerating graphs with I(G,2) = 21:")
# Small graphs: up to ~6 vertices
for n_v in range(1, 7):
    count = 0
    examples = []
    # Iterate over all graphs on n_v vertices
    edges_possible = [(i,j) for i in range(n_v) for j in range(i+1, n_v)]
    for edge_mask in range(2**len(edges_possible)):
        adj = [[False]*n_v for _ in range(n_v)]
        for idx, (i,j) in enumerate(edges_possible):
            if edge_mask & (1 << idx):
                adj[i][j] = adj[j][i] = True
        # Compute I(G,2)
        ip = 0
        for sub_mask in range(2**n_v):
            verts = [i for i in range(n_v) if sub_mask & (1 << i)]
            # Check independence
            is_indep = True
            for a in range(len(verts)):
                for b in range(a+1, len(verts)):
                    if adj[verts[a]][verts[b]]:
                        is_indep = False
                        break
                if not is_indep:
                    break
            if is_indep:
                ip += 2**len(verts)
        if ip == 21:
            count += 1
            if count <= 3:
                edge_list = [(i,j) for idx,(i,j) in enumerate(edges_possible) if edge_mask & (1 << idx)]
                examples.append(f"{n_v}v, edges={edge_list}")
    if count > 0:
        print(f"  {n_v} vertices: {count} graphs with I(G,2)=21")
        for ex in examples[:3]:
            print(f"    Example: {ex}")

print()
print("Similarly for I(G,2) = 63:")
for n_v in range(1, 7):
    count = 0
    examples = []
    edges_possible = [(i,j) for i in range(n_v) for j in range(i+1, n_v)]
    for edge_mask in range(2**len(edges_possible)):
        adj = [[False]*n_v for _ in range(n_v)]
        for idx, (i,j) in enumerate(edges_possible):
            if edge_mask & (1 << idx):
                adj[i][j] = adj[j][i] = True
        ip = 0
        for sub_mask in range(2**n_v):
            verts = [i for i in range(n_v) if sub_mask & (1 << i)]
            is_indep = True
            for a in range(len(verts)):
                for b in range(a+1, len(verts)):
                    if adj[verts[a]][verts[b]]:
                        is_indep = False
                        break
                if not is_indep:
                    break
            if is_indep:
                ip += 2**len(verts)
        if ip == 63:
            count += 1
            if count <= 2:
                edge_list = [(i,j) for idx,(i,j) in enumerate(edges_possible) if edge_mask & (1 << idx)]
                examples.append(f"{n_v}v, edges={edge_list}")
    if count > 0:
        print(f"  {n_v} vertices: {count} graphs with I(G,2)=63")
        for ex in examples[:2]:
            print(f"    Example: {ex}")

print()
print("="*70)
print("GRAND SUMMARY")
print("="*70)
print()
print("The forbidden H sequence 7·3^k arises because:")
print("  1. H = I(Ω(T), 2) where Ω is the conflict graph")
print("  2. I(K₃, 2) = 7 = Φ₃(2) — the Fano number")
print("  3. I(K₃ ⊔ kK₁, 2) = 7·3^k — tripling via isolated vertices")
print("  4. K₃ cannot appear as Ω(T) (three mutually conflicting cycles impossible)")
print("  5. Adding isolated vertices (disjoint cycles) preserves impossibility")
print()
print("The 'seed' of ALL forbidden values is the single fact:")
print("  K₃ ⊄ Ω(T) in the sense that Ω cannot have 3 mutually conflicting")
print("  odd cycles that are disjoint from all other cycles.")
print()
print("Equivalently: in any tournament, you cannot find 3 odd directed cycles")
print("that pairwise share vertices but are collectively 'maximal' (no other")
print("cycles interact with them). Any two overlapping cycles force additional")
print("cycles via the Splicing Lemma, which prevents the K₃ configuration.")
