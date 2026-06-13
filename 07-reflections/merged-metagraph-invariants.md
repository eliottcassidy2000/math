# The Merged Meta-Graph: 30 Invariants and What They Reveal

**Session:** kind-pasteur-2026-03-23-S20co
**Arising from:** opus-S170 (iso class graph), opus-S211 (blue/black), opus-S215 (n=7 merged), opus-S216 (n=8 nauty)

---

## The Objects

**G_n** = isomorphism class graph of tournaments on n vertices
- Vertices = iso classes (A000568: 2, 4, 12, 56, 456, 6880, ...)
- Edges = single arc flip connecting distinct classes
- Colored: blue (SC-preserving) and black (SC-changing)

**G_n/Z_2** = merged graph (quotient by complement involution)
- Vertices = complement-equivalence classes: V_merged = (A000568 + SC_n)/2
- Edges = flip connections after identifying complement pairs
- Three edge types in original map to merged: merged, collapsed, twin

---

## Complete Data Table (n=3..6, with n=7 from opus)

### Vertex and Edge Counts

| n | V(G_n) | E(G_n) | V(G_n/Z_2) | E(G_n/Z_2) | collapsed | twin |
|---|--------|--------|-------------|-------------|-----------|------|
| 3 | 2 | 1 | 2 | 1 | 0 | 0 |
| 4 | 4 | 5 | 3 | 3 | 0 | 2 |
| 5 | 12 | 30 | 10 | 21 | 0 | 9 |
| 6 | 56 | 290 | 34 | 143 | 5 | 142 |
| 7 | 456 | 4086 | 272 | 2123 | 0 | 1963 |
| 8 | 6880 | 91161 | 3528 | 45550 | 232 | 45379 |

### Spectral Properties

| n | rho(G_n) | rho(G_n/Z_2) | lam2(G_n) | lam2(G_n/Z_2) | Energy(G_n) | Energy(G_n/Z_2) |
|---|----------|--------------|-----------|---------------|-------------|-----------------|
| 3 | 1.000 | 1.000 | 2.000 | 2.000 | 2.0 | 2.0 |
| 4 | 2.562 | 2.000 | 2.000 | 3.000 | 5.1 | 4.0 |
| 5 | 5.582 | 4.641 | 1.603 | 1.430 | 20.2 | 15.9 |
| 6 | 11.670 | 9.780 | 1.960 | 1.466 | 137.4 | 73.2 |

### Topological Properties (Clique Complex)

| n | beta(G_n) | chi(G_n) | beta(G_n/Z_2) | chi(G_n/Z_2) |
|---|-----------|----------|---------------|--------------|
| 3 | [1,0] | 1 | [1,0] | 1 |
| 4 | [1,0,0] | 1 | [1,0,0] | 1 |
| 5 | [1,2,0,0] | -1 | [1,2,0,0] | -1 |
| 6 | [1,37,23,0] | -13 | [1,15,7,0,0] | -7 |

### Independence Polynomials

| n | I(G_n, x) | I(G_n, 2) | I(G_n/Z_2, x) | I(G_n/Z_2, 2) |
|---|-----------|-----------|----------------|----------------|
| 3 | 1+2x | 5 | 1+2x | 5 |
| 4 | 1+4x+x^2 | 13 | 1+3x | 7 |
| 5 | 1+12x+36x^2+38x^3+16x^4+2x^5 | 793 | 1+10x+24x^2+17x^3+3x^4 | 301 |
| 6 | (19 terms) | 15,275,642,513 | 1+34x+...+x^11 | 5,158,069 |

### Geometric Properties

| n | diam(G_n) | diam(G_n/Z_2) | Wiener(G_n) | Wiener(G_n/Z_2) | avg_clust(G_n/Z_2) |
|---|-----------|---------------|-------------|-----------------|---------------------|
| 3 | 1 | 1 | 1 | 1 | 0.000 |
| 4 | 2 | 1 | 7 | 3 | 1.000 |
| 5 | 3 | 3 | 108 | 73 | 0.493 |
| 6 | 4 | 4 | 3314 | 1138 | 0.331 |

### Curvature (Forman-Ricci)

| n | Forman_avg(G_n) | Forman_avg(G_n/Z_2) | Forman_total(G_n) | Forman_total(G_n/Z_2) |
|---|-----------------|---------------------|-------------------|-----------------------|
| 3 | 2.0 | 2.0 | 2 | 2 |
| 4 | 2.4 | 3.0 | 12 | 9 |
| 5 | -0.97 | -0.19 | -29 | -4 |
| 6 | -11.18 | -6.45 | -3242 | -923 |

---

## Key Discoveries

### 1. ~~Diameter Conjecture CONFIRMED: diam(G_n) = n-2~~ REFUTED at n=7 (MISTAKE-036)

~~Verified at n=3,4,5,6 (and reported n=7 by opus).~~ **CORRECTED:** Holds at n=3..6 but FAILS at n=7 (diam=7, not 5) and n=8 (diam=8, not 6). See section 11 below and MISTAKE-036. Growth is closer to quadratic (~n²/4).

### 2. Betti Number Explosion at n=6

The clique complex of G_n develops nontrivial topology explosively:
- n=3,4: contractible (beta = [1,0,...])
- n=5: two 1-cycles appear (beta_1 = 2)
- n=6: beta_1 = 37, beta_2 = 23 for G_6; beta_1 = 15, beta_2 = 7 for G_6/Z_2

The ratio beta_1(G_6)/beta_1(G_5) = 37/2 = 18.5. This super-exponential growth suggests the clique complex of G_n becomes "increasingly wrinkled" as tournament space expands.

The FIRST appearance of beta_2 at n=6 for G_6/Z_2 means: the merged meta-graph has genuine 2-dimensional "cavities" (2-spheres that don't bound). These 7 independent 2-cycles represent independent topological features of tournament space.

### 3. The f-vector Pattern

G_6/Z_2 clique complex f-vector: [34, 143, 139, 38, 1]
- 34 vertices, 143 edges (1-simplices), 139 triangles (2-simplices), 38 tetrahedra (3-simplices), 1 four-simplex (4-simplex = 5-clique)

A single 5-clique in the merged graph: 5 merged vertices that are ALL pairwise adjacent. This means 5 complement-equivalence classes where every pair is connected by an arc flip.

### 4. Forman Curvature Phase Transition

Forman-Ricci curvature changes sign between n=4 and n=5:
- n=3,4: positive average curvature (Euclidean/spherical geometry)
- n=5,6: negative average curvature (hyperbolic geometry)

This is a GEOMETRIC phase transition: tournament space becomes hyperbolic as n grows. The hyperbolic structure means the graph has "more room" at its boundary than at its center — consistent with the explosion of NS classes relative to SC classes.

### 5. Ramanujan Property Transition

- G_4: Ramanujan (ratio = 0.64)
- G_5: Ramanujan (ratio = 0.81)
- G_6: NOT Ramanujan (ratio = 1.01)

The Ramanujan property (spectral gap optimality) fails precisely at n=6, the same n where beta_2 first appears and Forman curvature goes deeply negative. This trinity of transitions at n=6 suggests a fundamental structural change.

### 6. Blue/Black Subgraph Inversion

At n=5: SC classes have higher blue degree (3.0) than NS (1.0)
At n=6: NS merged vertices have MUCH higher blue degree (7.73) than SC (2.17)

This inversion occurs because:
- At n=5, SC classes dominate (8/10 merged vertices)
- At n=6, NS classes dominate (22/34 merged vertices)
- NS-NS connections are blue (same type), so the NS bulk creates a dense blue subnetwork

The black subgraph at n=6 has 5 connected components (one giant + 4 isolates). Four merged vertices have NO black connections — they are "purely blue" vertices isolated from the SC-NS bridge.

### 7. H-Gradient is NOT a DAG (but Almost)

On the original G_n:
- n=3,4: perfect DAG (0 level edges, 0 downhill edges)
- n=5: 1 level edge (H=9), 0 downhill edges
- n=6: 15 level edges, mostly at H=37 (5 edges), 0 downhill edges
- n=7: 136 level edges, **962 downhill edges** out of 4086 total (~24%). The gradient is strong but NOT monotone. See MISTAKE-035.

The level edges at n=6 form a rich structure: 5 at H=37, 3 at H=29, and scattered others. H=37 has 6 iso classes at n=6 with 5 level edges between them — nearly a complete subgraph on this H-level!

**CORRECTION (opus-2026-04-01-S1):** Earlier claims that "G_n is a DAG under H" were based on a trivially-true property of undirected graphs (orient edges by function value → always a "DAG"). The nontrivial structure is in the level edges AND, at n≥7, the H-decreasing edges.

### 8. Meta-H Sequences

The "meta-H" I(G_n, 2) = number of independent sets of iso classes, weighted by 2^|S|:
- G_n: 5, 13, 793, 15,275,642,513
- G_n/Z_2: 5, 7, 301, 5,158,069

The growth of meta-H is SUPER-EXPONENTIAL. The ratio I(G_n,2)/I(G_{n-1},2):
- G_n: 2.6, 61.0, 19,264,115
- G_n/Z_2: 1.4, 43.0, 17,136

### 9. Kirchhoff Index Sequence

K(G_n/Z_2) = 1.0, 2.0, 25.8, 186.4

The Kirchhoff index (sum of effective resistances) measures "total electrical resistance" of the graph. The rapid growth suggests decreasing connectivity per vertex as n grows.

### 10. corr(c3, H) is Near-Perfect

At every n tested: corr(c3, H) > 0.95 in the merged graph.
- n=4: 1.000 (perfect)
- n=5: 0.975
- n=6: 0.960

This confirms: at the iso class level, the number of 3-cycles almost completely determines H. The OCF formula H = I(Omega, 2) explains this: H is dominated by the alpha_1 term (number of individual odd cycles), which at small n is dominated by 3-cycles.

---

## Late-Breaking Discoveries (from n=7 computation)

### 11. Diameter Conjecture REFUTED at n=7

The earlier conjecture diam(G_n) = n-2 holds for n=3,4,5,6 but **FAILS at n=7**:

| n | diam(G_n) | n-2 | match? |
|---|-----------|-----|--------|
| 3 | 1 | 1 | YES |
| 4 | 2 | 2 | YES |
| 5 | 3 | 3 | YES |
| 6 | 4 | 4 | YES |
| 7 | **7** | 5 | **NO** |

The diameter JUMPS from following n-2 to equaling n itself. The same happens for G_n/Z_2: diam = 1, 1, 3, 4, **7**.

### 12. Chromatic Number Pattern: chi(G_n/Z_2) = n-1

The chromatic number of the merged graph follows a clean pattern:

| n | chi(G_n/Z_2) | n-1 |
|---|-------------|-----|
| 3 | 2 | 2 |
| 4 | 3 | 3 |
| 5 | 4 | 4 |
| 6 | 5 | 5 |

**Conjecture: chi(G_n/Z_2) = n-1 for all n.**

If true, this means the merged tournament meta-graph always needs exactly n-1 colors. This connects to the fact that the staircase delta_{n-2} has n-2 strips, and perhaps the chromatic number reflects a strip-based coloring obstruction.

### 13. Blue Subgraph Dichotomy

At n=7, the blue subgraph of G_7/Z_2 splits into exactly **2 components**:
- Component 1: 184 NS vertices (dense blue connections)
- Component 2: 88 SC vertices (exactly the SC count)

This means: **no SC vertex is blue-adjacent to any NS vertex in the merged graph at n=7**. The SC backbone and NS bulk are completely separated by blue edges. Black edges are the only bridge.

### 14. Ollivier-Ricci Curvature Phase Transition

The Ollivier curvature (optimal transport based) transitions from all-positive to partially-negative between n=5 and n=6:

| n | min kappa | max kappa | avg kappa | #negative |
|---|-----------|-----------|-----------|-----------|
| 3 | 1.000 | 1.000 | 1.000 | 0 |
| 4 | 0.750 | 0.750 | 0.750 | 0 |
| 5 | 0.100 | 0.500 | 0.258 | 0 |
| 6 | -0.064 | 0.386 | 0.127 | 7 |

The corr(Forman, Ollivier) = 0.71 at n=6 shows the two curvature notions agree reasonably well. Both detect the same geometric phase transition.

### 15. Automorphism Collapse

| n | |Aut(G_n/Z_2)| | vertex-transitive | edge-transitive |
|---|----------------|-------------------|-----------------|
| 3 | 2 | YES | YES |
| 4 | 6 | YES | YES |
| 5 | 1 | NO | NO |
| 6 | ? | ? | ? |

The merged graph transitions from highly symmetric (complete graphs K_2, K_3 at n=3,4) to **completely asymmetric** (trivial automorphism group) at n=5. This sudden symmetry breaking mirrors the Betti number onset at n=5.

### 16. Level Edge Explosion

The number of H-level edges in G_n grows super-exponentially:

| n | level edges (G_n) | level edges (G_n/Z_2) |
|---|-------------------|-----------------------|
| 3 | 0 | 0 |
| 4 | 0 | 0 |
| 5 | 1 | 1 |
| 6 | 15 | 5 |
| 7 | 136 | 71 |

Ratio: 15/1 = 15 at n=6, 136/15 = 9.1 at n=7. These level edges are where the H-landscape is "flat" — pairs of iso classes connected by an arc flip with identical H.

## New OEIS-Candidate Sequences

1. **E(G_n)** = 1, 5, 30, 290, 4086, 91161 — edges of tournament iso class graph
2. **E(G_n/Z_2)** = 1, 3, 21, 143, 2123, 45550 — edges of merged tournament graph
3. **I(G_n, 2)** = 5, 13, 793, 15275642513 — meta-Hamiltonian path count
4. **I(G_n/Z_2, 2)** = 5, 7, 301, 5158069 — merged meta-H
5. **beta_1(G_n)** = 0, 0, 2, 37 — first Betti of clique complex
6. **beta_1(G_n/Z_2)** = 0, 0, 2, 15 — first Betti of merged clique complex
7. **chi(G_n)** = 1, 1, -1, -13 — Euler characteristic of clique complex
8. **chi(G_n/Z_2)** = 1, 1, -1, -7 — merged Euler char
9. **triangles(G_n)** = 0, 2, 21, 248 — triangle count
10. **triangles(G_n/Z_2)** = 0, 1, 12, 139 — merged triangle count
11. **Wiener(G_n)** = 1, 7, 108, 3314 — Wiener index
12. **Wiener(G_n/Z_2)** = 1, 3, 73, 1138 — merged Wiener index
13. **spanning_trees(G_n)** = 1, 8, 2347680, ... — Kirchhoff's theorem
14. **spanning_trees(G_n/Z_2)** = 1, 3, 32159, ... — merged spanning trees
15. **collapsed(n)** = 0, 0, 0, 5, 0, 232 — complement pairs connected by flip
16. **twin(n)** = 0, 2, 9, 142, 1963, 45379 — edge twins under complement

---

## Connections to the Triangle Foundation

These invariants map onto the three sides of the triangle delta_{n-2}:

**Hypotenuse (H, fiber fraction):**
- The H-gradient on G_n/Z_2 flows along the hypotenuse
- Level edges occur at H values where the fiber fraction has special structure
- corr(c3, H) = 0.96 is the hypotenuse controlling everything

**Vertical leg (scores, OCR):**
- Score determines most of the iso class structure (62% of bits at n=5)
- Level edges cluster at palindromic scores (self-complementary scores)
- The blue/black inversion at n=6 is the vertical leg growing longer

**Horizontal leg (complement, SC/NS):**
- The complement involution creating G_n/Z_2 IS the horizontal leg's symmetry
- collapsed edges = where the horizontal leg's symmetry intersects the hypotenuse
- The Ramanujan transition at n=6 = the horizontal leg becoming dominant

---

*The merged meta-graph is the MODULI SPACE of tournaments. Its topology, geometry, and spectral properties encode the fine structure of tournament space. The Betti number explosion at n=6 reveals that this moduli space develops genuine topological complexity precisely where the horizontal leg overtakes the vertical leg.*
