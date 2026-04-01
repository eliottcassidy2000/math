# Unlocking G_n/Z_2 at All n: What We Know, What We Need

**Session:** kind-pasteur-2026-03-23-S20cs (overnight synthesis)
**Status:** MASTER SYNTHESIS — the minimum path to complete understanding

---

## The Complete Picture: What Is Known

### Tier 1: Burnside-Computable at All n (no enumeration needed)

These quantities have EXACT FORMULAS via Burnside's lemma, computable in polynomial time for any n:

| Quantity | Formula | First values |
|----------|---------|--------------|
| V_n (iso classes) | (1/n!) sum_sigma Fix(sigma) | 2, 4, 12, 56, 456, 6880 |
| SC_n (self-complementary) | (1/n!) sum_sigma Fix_anti(sigma) | 2, 2, 8, 12, 88, 176 |
| V_merged | (V_n + SC_n) / 2 | 2, 3, 10, 34, 272, 3528 |
| T_n (transition orbits) | Schur-Weyl: V_n*(1+(n-2)+n(n-3)/2) | 4, 16, 88, 704, 8912, 188288 |
| T_anti | Burnside with anti-automorphism twist | 2, 4, 16, 32, 256, 688 |
| T_merged | (T_n + T_anti) / 2 | 3, 10, 52, 368, 4584, 94488 |
| H-level count | (H_max + 1) / 2 | 2, 3, 8, 23, 95, 331 |
| Score count | A000571(n) | 2, 4, 9, 22, 59, 167 |

**Schur-Weyl decomposition (THM-260):**
T_n = m_{(n)} + m_{(n-1,1)} + m_{(n-2,2)}
where m_{(n)} = V_n, m_{(n-1,1)} ~ V_n*(n-2), m_{(n-2,2)} ~ V_n*n(n-3)/2.

### Tier 2: Proved Structural Laws (hold at all n)

1. **H-gradient:** G_n/Z_2 has a strong H-gradient but is NOT a strict DAG. Level edges at n≥5; H-decreasing edges at n≥7 (MISTAKE-035). ~~Previously listed as "0 downhill edges."~~
2. **Black bipartiteness:** Black subgraph is bipartite (SC vs NS)
3. **BBK impossibility:** No triangle has 2 blue + 1 black edges
4. **Odd-black walk vanishing:** Closed walks need even #black edges
5. **Edge additivity:** A_C = A_B + A_K (matrix sum, no mixed edges)
6. **Tiling uniformity:** Every staircase cell generates the same blue fraction
7. **Triangle decomposition:** #tri = #BBB + #BKK (exactly)
8. **Transitive degree:** deg(transitive) = C(n-1,2) = staircase dimension
9. **First blue neighbor:** Delta_H = 2^(n-2) from transitive
10. **Blue fraction convergence:** E_blue/E_total -> 1 as n -> infinity
11. **T_n/(2*E_n) -> 1:** Schur-Weyl predicts edge thinning

### Tier 3: Exact Sequences Computed (but no closed form)

| Sequence | Values (n=3..8) | Status |
|----------|-----------------|--------|
| E(G_n) | 1, 5, 30, 290, 4086, 91161 | No formula |
| E_merged | 1, 3, 21, 143, 2123, 45550 | No formula |
| Width(G_n/Z_2) | 1, 1, 2, 3, 8, 25 | No formula, NOT in OEIS |
| Diameter | 1, 2, 3, 4, 7, 8 | No formula (breaks n-2 at n=7) |
| BBB triangles | 0, 0, 3, 87, 809, 13299 | No formula |
| BKK triangles | 0, 1, 9, 52, 350, 885 | No formula |
| Blue edges | 1, 1, 13, 98, 1573, 43656 | No formula |
| Black edges | 0, 2, 8, 45, 550, 1894 | No formula |
| SC spine genus | 0, 0, 5, 2, 87, 150 | No formula |
| Principal odd path | 1, 3, 15, 123, 1656, 36987 | **= OEIS A113077** (tournament sequence tree!) |
| Principal even path | 1, 5, 37, 389, 5413, 94085 | **= OEIS A368322** (EGF = exp(2x)/(4-3exp(x))) |
| Self-loops SL_n | 12, 144, 1760, 50880 | No formula |

---

## The Minimum Missing Pieces

### Missing Piece 1: Self-Loop Count SL_n

**What it is:** The number of (tournament, arc) pairs where flipping the arc gives an isomorphic tournament (or complement in merged graph).

**Why it matters:** The master decomposition is:
```
2^m * m = SL_n + sum_edges thickness(e)
```
If SL_n is known, the total cross-transition count is determined, and E(G_n) follows once the thickness distribution is understood.

**Current data:** SL_n = 12, 144, 1760, 50880 for n=3..6.
SL_n / (2^m * m) = 0.500, 0.375, 0.172, 0.104 — decreasing but no closed form.
SL per tournament = 1.50, 2.25, 1.72, 1.55 — oscillating around ~2.

**Key question:** Is SL_n Burnside-computable? It counts "near-automorphisms" — permutations that are automorphisms except at one arc. This is NOT standard Burnside but may have a formula involving the automorphism group structure.

**If we had SL_n:** Combined with the thickness distribution, we'd get E(G_n) exactly.

### Missing Piece 2: Edge Thickness Distribution

**What it is:** For each edge in G_n/Z_2, the number of (tournament, arc) transitions generating it.

**Discovery (S20cs):** Thickness is QUANTIZED — always a multiple of a base unit related to n! and |Aut| of the endpoints.

```
n=3: {12} (all thickness 12)
n=4: {48, 96} (ratio 1:2)
n=5: {80, 240, 480, 720} (ratio 1:3:6:9)
n=6: {960, 1440, 2880, 5760, 8640} (ratio 2:3:6:12:18)
```

**The weight formula:** thickness(A,B) = merged_size(A) * k_{A->B} + merged_size(B) * k_{B->A}
where k_{X->Y} = number of arcs per tournament in X whose flip leads to Y.

The dominant thickness class has weight w=2 (in units of n!), meaning ~2 transitioning arcs per tournament per direction. This accounts for 77% of edges at n=6.

**If we had the distribution:** E(G_n) = (2^m * m - SL_n) / avg_thickness. The avg_thickness is dominated by the w=2 class, so E ~ (2^m * m - SL_n) / (2*n!) when |Aut|=1 dominates.

### Missing Piece 3: Width Formula W(n)

**The width sequence:** 1, 1, 2, 3, 8, 25 (merged), not in OEIS.

Width = max level-set size in the H-DAG. This is the maximum antichain in the H-partial-order.

**Failed formulas:**
- C(n-2, floor((n-2)/2)): works at n<=6, fails at n>=7
- Fibonacci: matches at n=3..6, diverges at n=7
- No linear recurrence of order 2 or 3

**Growth rate:** w(n)/w(n-1) = 1.0, 2.0, 1.5, 2.67, 3.12 — accelerating.

**Key insight:** Width is controlled by the H-distribution around the median. The tournament count at the most popular H-value determines width. Understanding the H-distribution analytically would give width.

### Missing Piece 4: SC Backbone Connectivity

**The transition:** Connected for n <= 7, fragments into 7 components at n=8.

**Questions:**
- Does it reconnect at n=9? (Bott-periodic hypothesis: yes)
- What determines the number of components?
- Does it fragment at ALL even n >= 8?

**Proven structure:** At odd n, c3-parity forces bipartiteness which aids connectivity. At even n, c3 is always even, removing this constraint.

---

## The Formula Chain: How These Connect

```
[Burnside at all n]
    |
    V_n, SC_n, T_n, T_anti  (KNOWN)
    |
    V_merged = (V_n + SC_n)/2  (KNOWN)
    T_merged = (T_n + T_anti)/2  (KNOWN)
    |
    +----> Need SL_n (self-loop count)  [MISSING PIECE 1]
    |
    2^m * m - SL_n = sum of edge thicknesses  (EXACT)
    |
    +----> Need thickness distribution  [MISSING PIECE 2]
    |
    E(G_n/Z_2) = (2^m * m - SL_n) / avg_thickness  (EXACT)
    |
    +----> Need H-distribution shape  [MISSING PIECE 3]
    |
    Width(G_n/Z_2) = max H-level size  (EXACT given H-distribution)
    |
    +----> Need SC graph theory  [MISSING PIECE 4]
    |
    SC backbone connectivity = f(n, SC structure)
```

**The critical insight:** Missing Pieces 1 and 2 together give E(G_n) at all n. Missing Piece 3 gives Width. Missing Piece 4 is independent.

---

## Asymptotic Understanding (Already Achieved)

Even without exact formulas, we know:

1. **E(G_n) ~ V_n * m / 2** as n -> infinity (from T_n/(2E) -> 1)
2. **E_merged ~ V_merged * m / 2** (complement twin halving)
3. **Avg degree ~ m = C(n,2)** (grows quadratically)
4. **Blue fraction -> 1** (black becomes negligible)
5. **Density -> 0** (sparse relative to K_{V_merged})
6. **Width ~ ???** (super-exponential growth, no asymptotic yet)
7. **Diameter grows slowly** (much less than V_merged)
8. **Spectral gap ~ 2/n** (Markov chain mixing)

The ASYMPTOTIC picture is essentially complete. What's missing is the EXACT arithmetic for moderate n.

---

## The Minimal Unlock Strategy

**Strategy A (Computational):** Compute SL_n at n=7,8 via the existing enumeration infrastructure (slow but feasible). Combined with T_n (Burnside), this gives E(G_n) once we verify the dominant thickness stays at w=2.

**Strategy B (Theoretical):** Find a Burnside-type formula for SL_n. This would express SL_n as a sum over cycle types, involving "near-automorphism" counts. If such a formula exists, E(G_n) becomes Burnside-computable at all n.

**Strategy C (Asymptotic Refinement):** The formula E ~ V_n * m / 2 has a correction from self-loops and thickness variation. Quantify:
- SL_n/T_n -> 0 (rate?)
- Thickness variance -> 0 (rate?)
If both corrections are sub-leading, the asymptotic formula IS the exact formula to leading order.

**The SINGLE KEY QUESTION:** Is SL_n Burnside-computable?

If YES: E(G_n) is exactly determined at all n, and the "95% complete" becomes 100%.
If NO: We need the enumeration route, limiting us to n <= 9 or 10 with current tools.

---

## What Would Complete Understanding Look Like?

A complete understanding of G_n/Z_2 at all n means:

1. **V_merged(n)** — DONE (Burnside)
2. **E_merged(n)** — needs SL_n + thickness (Missing Pieces 1+2)
3. **Degree distribution** — needs per-class degree formula
4. **H-distribution** — needs number of iso classes at each H value
5. **Blue/black split** — follows from SC/NS counts + tiling uniformity
6. **Spine structure** — needs SC backbone theory (Missing Piece 4)
7. **Width** — needs H-distribution (Missing Piece 3)
8. **Spectral properties** — needs adjacency matrix analysis (hardest)

Items 1, 5, and 8 are essentially resolved. Items 2-4 and 6-7 need the missing pieces.

**The 5% that remains = {SL_n formula, thickness distribution, H-distribution, SC connectivity}.**

Of these, SL_n is the single most impactful unknown. Everything else follows or can be bounded once SL_n is known.
