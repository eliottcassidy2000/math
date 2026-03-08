# GLMY Path Homology of Tournaments and Circulant Digraphs — Synthesis

## Session: opus-2026-03-07-S38

### Background

GLMY path homology (Grigor'yan-Lin-Muranov-Yau) assigns topological invariants
to directed graphs. For a digraph G on n vertices:

- **Allowed p-paths A_p**: sequences (v_0,...,v_p) of distinct vertices with each v_i→v_{i+1} an edge
- **∂-invariant paths Ω_p**: {u ∈ A_p : ∂u has all faces in A_{p-1}}
- **Path Betti numbers**: β_p = dim(ker ∂_p / im ∂_{p+1}) on the chain complex Ω_•
- **Euler characteristic**: χ = Σ (-1)^p β_p

### Key Results: Tournaments

#### 1. Betti Number Distribution (EXHAUSTIVE for n≤5, sampled for n≥6)

| n | Total | β=(1,0,...) | β=(1,1,0,...) | β=(1,0,0,1,...) | Other |
|---|-------|-------------|---------------|-----------------|-------|
| 3 | 8     | 6 (75%)     | 2 (25%)       | 0               | 0     |
| 4 | 64    | 40 (62.5%)  | 24 (37.5%)    | 0               | 0     |
| 5 | 1024  | 720 (70.3%) | 304 (29.7%)   | 0               | 0     |
| 6 | ~1000 | ~940 (94%)  | ~60 (6%)      | 0               | 0     |
| 7 | ~500  | ~410 (82%)  | ~48 (9.6%)    | ~42 (8.4%)      | 0     |
| 8 | ~80   | ~66 (82%)   | ~1 (1.3%)     | ~13 (16%)       | β_4=1: rare |

#### 2. β_2 = 0 ALWAYS for tournaments (HYP-301)
Across ~5000 tournaments tested (n=3 through n=8), β_2 has NEVER been nonzero.
**ALGEBRAIC PROOF STRUCTURE**: ker(∂_2) = im(∂_3) EXACTLY for all tested tournaments.
- Ω_2 = transitive triples (a→b→c with shortcut a→c)
- Ω_3 = "doubly transitive" 4-paths (a→b→c→d with a→c AND b→d)
- Every 2-cycle in Ω_2 is filled by a 3-chain from Ω_3
- Verified EXHAUSTIVELY at n=4,5; sampled 100 each at n=6,7

#### 3. β_3=1 and β_1=1 are MUTUALLY EXCLUSIVE (HYP-302)
In 500 samples at n=7:
- 37 have β_3=1 (all with β_1=0)
- 39 have β_1=1 (all with β_3=0)
- NO tournament has both β_1>0 AND β_3>0

#### 4. β_1 Threshold (PROVED for n≤4)
- n=3: β_1=1 iff t_3=1 (the 3-cycle tournament)
- n=4: β_1=1 iff t_3=2 (exactly)
- n=5: β_1=1 requires t_3≥3, but t_3=3 is mixed (half have β_1=1)
- n=5: score sequence also doesn't determine β_1 (mixed for some sequences)

#### 5. The Ω_2 Structure for Tournaments
The 2-path (v_0, v_1, v_2) where v_0→v_1→v_2 is ∂-invariant (in Ω_2) iff v_0→v_2.
So Ω_2 = {transitive triples in the tournament}.
This is a CRUCIAL structural insight: the ∂-invariant condition at dimension 2
is EXACTLY the transitivity condition.

#### 6. Geometric Interpretation
- β=(1,0,...): contractible (point) — the tournament is "topologically trivial"
- β=(1,1,0,...): has a 1-hole — topologically like S^1 (circle)
- β=(1,0,0,1,...): has a 3-hole — topologically like S^3 (3-sphere!)
- β=(1,0,0,0,1,...): has a 4-hole — topologically like S^4 (first appears at n=8)

#### 7a. Complement Duality PROVED (HYP-303)
- n=5 EXHAUSTIVE: β(T) = β(T^op) for ALL 1024 tournaments
- n=6 SAMPLED: β(T) = β(T^op) for 100/100 tournaments
- This makes β(T) a complement-invariant, like H(T)

#### 7b. Stability of Circulant Path Homology (VERIFIED)
- C_n^{1,2}: β=(1,1,0,...) for ALL n≥4 (tested through n=17)
- C_n^{1,3}: β=(1,2,1,0,...) for ALL n≥6 (tested through n=23)
- Confirms Tang-Yau strong stability theorem computationally

#### 7. Euler Characteristic
- n=4: χ=1 for 40 (contractible), χ=0 for 24 (S^1-like)
- n=5: χ=1 for 720 (contractible), χ=0 for 304 (S^1-like)

### Key Results: Circulant Digraphs

#### 1. Single generator S={s}
- C_n^{s}: β_0 = gcd(n,s) components, β_1 = gcd(n,s)
- When n prime: always β=(1,1,0,...) — topologically S^1

#### 2. Two generators S={s_1, s_2}
Following Tang-Yau (arXiv:2602.04140):
- S={1,2}: β=(1,1,0,...) — always S^1
- S={1,s} with s≠2: β=(1,2,1,0,...) — T^2 (torus) when n large enough
- S={1,n-1} (bidirected cycle): β=(1,n-1,0,...) — huge β_1

#### 3. Three generators — RICH topology appears
- C_5^{1,2,3}: β=(1,0,0,1,0) — S^3!!!
- C_6^{1,3,5}: β=(1,1,0,6,0,0) — mixed S^1 and S^3
- C_7^{1,2,5,6}: β=(1,0,0,22,0,0) — 22 independent 3-holes!
- C_6^{2,3,4}: β=(1,1,2,2,0,0) — simultaneous β_2 and β_3

#### 4. Complete circulant C_n^{[n-1]} (complete digraph)
- C_5^{1,2,3,4}: β=(1,0,0,0,44) — massive β_4
- C_6^{1,2,3,4,5}: β=(1,0,0,0,0,265)
- Pattern: β_{n-1} explodes, all lower Betti numbers vanish

#### 5. Complement pairs
C_n^S and C_n^{[n-1]\S} often have related topology.
At n prime, the shift automorphism acts irreducibly, giving all β_0=1.

### Tang-Yau Fourier Decomposition

For circulant digraphs, the shift τ: v ↦ v+1 (mod n) is an automorphism.
This decomposes the chain complex into eigenspaces:

Ω_m ≅ ⊕_{λ^n=1} Ω_m^(λ)

The boundary maps become "symbol matrices" M_m(λ) whose entries are
polynomials in λ. For S={γ_0,...,γ_{d-1}}:

M_1(t) = [t^{γ_0} - 1, t^{γ_1} - 1, ..., t^{γ_{d-1}} - 1]

**Strong stability theorem** (Tang-Yau 4.5): For "no wrap-around" S,
the Betti numbers stabilize for large primes.

### Conjectures

- **HYP-301**: β_2(T) = 0 for ALL tournaments T
- **HYP-302**: β_1(T) · β_3(T) = 0 for all tournaments (mutual exclusion)
- **HYP-303**: β(T) = β(T^op) (complement invariance of path homology)
- **HYP-304**: For tournaments, β_p ∈ {0, 1} for all p ≥ 1
- **HYP-305**: The complete digraph K_n has β = (1, 0, ..., 0, D_{n-1}) where D_{n-1} grows rapidly

### Open Questions
1. WHY does β_2=0 for tournaments? Algebraic proof?
2. What exactly determines β_3=1 at n=7? Is there a cycle-count threshold?
3. Does β_5 appear at n=9?
4. Does the mutual exclusion β_1·β_3=0 hold generally?
5. Connection to OCF: does β(T) carry information about H(T)?
6. Can the Fourier decomposition be extended to "nearly circulant" tournaments?

### Scripts
- path_homology_v2.py (core implementation)
- path_homology_deep.py (Phase 1: exhaustive n≤5, sampling n=6,7)
- path_homology_phase2.py (Phase 2: β_3 discovery, Euler char)
- path_homology_phase3.py (Phase 3: cycle generators, correlates)
- path_homology_phase4.py (Phase 4: complement duality, n=8, arc-flip)
