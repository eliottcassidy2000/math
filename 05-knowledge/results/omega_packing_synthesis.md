# Ω Packing Structure Synthesis
**opus-2026-03-14-S71f**

## Key Results

### 1. OCF Verified with Full Odd-Cycle Counting
- **n=5**: H = I(Ω,2) for ALL 1024 tournaments ✓
- **n=6**: H = I(Ω,2) for ALL 32768 tournaments ✓
- Ω vertices = ALL directed odd cycles (3-cycles, 5-cycles, etc.)
- Each *directed* cycle is a separate vertex (multiple directed cycles on same vertex set)

### 2. Ω Completeness Transition
- **n≤5**: Ω is ALWAYS a complete graph (all cycles share vertices)
  - H = 1 + 2·|V(Ω)| (only α₁ matters)
- **n=6**: Ω is non-complete in 15440/32768 = 47% of tournaments
  - First appearance of α₂ > 0 (disjoint 3-cycle pairs need 6 vertices)
- **n=7**: Non-complete in 87% of sampled tournaments

### 3. Ω Structure Classification at n=6

**Single clique structures** (Ω = K_m):
- K_0: 720 (H=1, transitive)
- K_1: 960 (H=3)
- K_2: 2160 (H=5, cuboid brick!)
- K_4: 2880 (H=9)
- K_5: 1440 (H=11)
- K_6: 1440 (H=13)
- K_7: 2208 (H=15)
- K_8: 720 (H=17)
- K_9: 1440 (H=19)
- K_11: 1440 (H=23)
- K_12: 1440 (H=25)
- K_16: 480 (H=33)

**General (non-clique) structures** — components that are NOT cliques:
- (6): 720 (H=17)
- (9): 1440 (H=23)
- (10): 720 (H=29)
- (11): 480 (H=27)
- (12): 2880 (H=29,33)
- (13): 1440 (H=31)
- (14): 3840 (H=33,37,45)
- (16): 2160 (H=37,41)
- (19): 1440 (H=43)
- (20): 240 (H=45)

**Packing decompositions** (Ω = disjoint union of cliques, multiple components):
- K_1 ⊔ K_1: 80 tournaments (H=9) — two simplex bricks ← (1+x)² = 1+2x+x²

### 4. Packing Interpretation

I(Ω, x) = prod_i (1 + c_i · x) when Ω decomposes as ⊔ K_{c_i}

| Brick | Ω component | I factor | H contribution |
|-------|-------------|----------|----------------|
| Simplex | K_1 (isolated vertex) | (1+x) | ×3 |
| Cuboid | K_2 (edge) | (1+2x) | ×5 |
| Tesseract | K_3 (triangle) | (1+3x) | ×7 (FORBIDDEN for m=1!) |
| c-cell | K_c | (1+cx) | ×(2c+1) |

**Only 80/32768 = 0.24% of n=6 tournaments have genuine packing.**

### 5. Critical Observations

1. **H=9 has TWO interpretations**: 
   - Monolithic K_4 (2880 tournaments): 4 mutually overlapping cycles
   - Simplex² K_1⊔K_1 (80 tournaments): 2 disjoint cycles

2. **H=5 IS the cuboid brick**: Ω = K_2 at n=4 (two overlapping 3-cycles sharing an edge)

3. **H=7 tesseract brick**: Would need Ω = K_3 (3 mutually overlapping cycles, NO others).
   But at n≤6 with t₃=3: additional cycles always appear → Ω ≠ K_3.
   At n=7 with t₃=3: d₅ ≥ 1 always → additional 5-cycle vertex → Ω ≠ K_3.
   
4. **Cuboid² (K_2⊔K_2)**: Needs 2 disjoint pairs of overlapping 3-cycles.
   Min vertex count: 4+4=8. Not achievable at n≤7.

5. **The "general" Ω structures** cannot be decomposed into clique products.
   They have more complex independence polynomials with non-factorizable structure.

### 6. The (2,3) Packing Hierarchy

Evaluation at x = k maps brick sizes to polytope hierarchies:
- x=1: (1+c)^m — simplex counting  
- x=2: (1+2c)^m — OCF = H
- x=3: (1+3c)^m — cuboid counting

The k-nacci → 2 connection: ratio I(2)/I(1) = (1+2c)/(1+c) → 2 as c → ∞
The weighted k-nacci → 3: ratio I(3)/I(1) = (1+3c)/(1+c) → 3 as c → ∞

## UPDATED: Complete Independence Polynomial Classification at n=6

### OCF: H = 1 + 2α₁ + 4α₂ for ALL 32768 tournaments ✓

At n=6, I(Ω, x) = 1 + α₁x + α₂x² (α₃ = 0 always, since 3 disjoint 3-cycles need ≥9 vertices).

### 27 Distinct Independence Polynomials

| α₁ | α₂ | H | Count | Factorization |
|-----|-----|------|-------|---------------|
| 0 | 0 | 1 | 720 | 1 (transitive) |
| 1 | 0 | 3 | 960 | (1+x) |
| 2 | 0 | 5 | 2160 | (1+2x) — cuboid brick |
| 2 | 1 | 9 | 80 | **(1+x)² — simplex²!** |
| 4 | 0 | 9 | 2880 | (1+4x) — 4-cell |
| 5 | 0 | 11 | 1440 | (1+5x) |
| 6 | 0 | 13 | 1440 | (1+6x) |
| 6 | 1 | 17 | 720 | **IRREDUCIBLE** |
| 7 | 0 | 15 | 2208 | (1+7x) |
| 8 | 0 | 17 | 720 | (1+8x) |
| 9 | 0 | 19 | 1440 | (1+9x) |
| 9 | 1 | 23 | 1440 | **IRREDUCIBLE** |
| 10 | 2 | 29 | 720 | **IRREDUCIBLE** |
| 11 | 0 | 23 | 1440 | (1+11x) |
| 11 | 1 | 27 | 480 | **IRREDUCIBLE** |
| 12 | 0 | 25 | 1440 | (1+12x) |
| 12 | 1 | 29 | 2160 | **IRREDUCIBLE** |
| 12 | 2 | 33 | 720 | **IRREDUCIBLE** |
| 13 | 1 | 31 | 1440 | **IRREDUCIBLE** |
| 14 | 1 | 33 | 1440 | **IRREDUCIBLE** |
| 14 | 2 | 37 | 2160 | **IRREDUCIBLE** |
| 14 | 4 | 45 | 240 | **IRREDUCIBLE** |
| 16 | 0 | 33 | 480 | (1+16x) |
| 16 | 1 | 37 | 1440 | **IRREDUCIBLE** |
| 16 | 2 | 41 | 720 | **IRREDUCIBLE** |
| 19 | 1 | 43 | 1440 | **IRREDUCIBLE** |
| 20 | 1 | 45 | 240 | **IRREDUCIBLE** |

### Key Insight: Most Tournaments are "Atoms"

Only 2 of 27 polynomials factorize:
- α₂=0 (degree 1): always "factorable" as single brick → 17328/32768 = 53%
- (2,1) → (1+x)² → 80/32768 = 0.24%

The remaining 15360/32768 = 47% have **irreducible** independence polynomials.
These are the "atomic" tournaments that CANNOT be decomposed into independent bricks.

### The H=9 Duality

H=9 is realized by TWO fundamentally different structures:
1. **4-cell** (α₁=4, α₂=0): 2880 tournaments. Ω = K₄ (single clique of 4 cycles).
2. **Simplex²** (α₁=2, α₂=1): 80 tournaments. Ω = K₁ ⊔ K₁ (two disjoint cycles).

Both give I(2) = 9, but their independence polynomials DIFFER:
- 4-cell: I = 1+4x (linear)
- Simplex²: I = 1+2x+x² = (1+x)² (quadratic, factorizable)

The simplex² tournaments have RICHER structure (α₂ > 0) but the SAME H value.
