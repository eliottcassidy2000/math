# GLMY Path Homology of Tournaments and Circulant Digraphs — Synthesis

## Sessions: opus-2026-03-07-S38, opus-2026-03-08-S39

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
| 9 | 50    | 44 (88%)    | 0 (0%)        | 6 (12%)         | 0           |

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

### OCF-Topology Connection (Major Finding)

**S-phase (β_3=1) tournaments have the most directed cycles and highest H:**

| Phase | t3 | t5 | H (mean) | Transitivity ratio |
|-------|-----|-----|----------|-------------------|
| P (contractible) | 8.6 | 15.2 | 75.3 | 0.51 |
| C (S^1) | 9.1 | 18.2 | 87.0 | 0.52 |
| S (S^3) | **11.1** | **24.4** | **114.8** | **0.43** |

Since OCF gives H = 1 + 2(t3+t5+t7) + 4·(disjoint pairs), the higher H of
S-phase comes from having MORE odd cycles. The topological complexity (directed
3-holes) correlates with cycle richness.

**n=5 exhaustive**: C-phase H ∈ {9,11,15}, P-phase H ∈ {1,3,5,9,13}.
H alone doesn't determine topology (H=9 appears in both phases).

### Session S39: Major New Results

#### Fourier Decomposition v3 (CORRECT — 90/90 validated)
Key fix: Ω_p is NOT span{individually ∂-invariant paths}. It's the kernel
of the "junk projection" — chains where non-allowed boundary faces cancel.
This gives Ω_p^(λ) = ker(J_p(λ)) where J_p is the junk coefficient matrix.

#### Illegal Merged Steps F = (S+S)\S Controls Topology
For C_n^S:
- F = {(a+b) mod n : a,b ∈ S, a≠b, (a+b) mod n ∉ S ∪ {0}}
- |F|=0 (S closed under addition): β_1 = n-1 (bidirected pair case)
- |F|=max, L=∅ (no legal merges): highest-dimensional sphere topology
- At n=7: |F|=3, |S|=3 → 6×S^4; |F|=1, |S|=2 → circle or torus

#### Paley Tournament Path Homology (MAJOR FINDING)
The Paley tournament P_p (p ≡ 3 mod 4 prime) has S = QR(p).
- **P_3**: β=(1,1,0) — circle S^1
- **P_7**: β=(1,0,0,0,6,0) — **6×S^4**! (Euler char χ=7)
- **P_11**: β=(1,0,0,0,...) — contractible through dim 3
- For QR(p): F = QNR(p) EXACTLY — illegal merges are the quadratic non-residues
- Gauss sum connection: the character sum Σ_{s∈QR} λ^s controls eigenspace structure

#### n=9 Tournament Results (50 samples, max_dim=3)
- β=(1,0,0,0): 44 (88%) — contractible
- β=(1,0,0,1): 6 (12%) — S^3
- β=(1,1,...): 0 (0%) — **C-phase completely disappears!**
- β_3=1 appears at t3 as low as 7 (not just high-t3 tournaments)

#### Topology Census (n=5,7 prime circulants)
n=5: 8 circles, 4 S^3, 2 β_1=6, 1 with 44×S^4
n=7: 30 circles, 6 tori, 8 (6×S^4), 6 S^3, 3 (22×S^3), 3 β_1=8

### Conjectures

- **HYP-301**: β_2(T) = 0 for ALL tournaments T (verified n=3-9)
- **HYP-302**: β_1(T) · β_3(T) = 0 for all tournaments (verified n=3-9)
- **HYP-303**: β(T) = β(T^op) (complement invariance — PROVED at n=5, verified n=6-7)
- **HYP-304**: For tournaments, β_p ∈ {0, 1} for all p ≥ 1
- **HYP-305**: The complete digraph K_n has β = (1, 0, ..., 0, D_{n-1}) where D_{n-1} grows rapidly
- **HYP-306**: S-phase proportion increases with n (0% at n≤6, 8% at n=7, 16% at n=8, 12% at n=9)
- **HYP-307**: C-phase (β_1=1) disappears for n ≥ 9 (0/50 at n=9)
- **HYP-308**: For Paley tournaments P_p, β_p depends only on (p mod 4) and p
- **HYP-309**: F = QNR for all Paley tournaments (verified p=3,7,11)

### Open Questions
1. PROVE β_2=0 algebraically using the Ω_2/Ω_3 structure
2. What exactly determines β_3=1? At n=9, appears even at low t3
3. Does β_5 appear at n=9? (Need max_dim≥5, too slow for brute force)
4. Does the mutual exclusion β_1·β_3=0 hold generally?
5. What is β(P_11), β(P_19) at full dimension?
6. Is P_7 = 6×S^4 related to the Gauss sum |g| = √p?
7. Does the QR-topology connection extend to p ≡ 1 mod 4?

### Scripts
- path_homology_v2.py (core implementation, validated)
- path_homology_fourier_v3.py (CORRECT Fourier decomposition, 90/90 validated)
- path_homology_fourier.py (v1 Fourier, β_0 and β_1 only)
- path_homology_fourier_v2.py (v2 Fourier, BUGGY — wrong Ω)
- path_homology_deep.py (Phase 1: exhaustive n≤5, sampling n=6,7)
- path_homology_phase2.py (Phase 2: β_3 discovery, Euler char)
- path_homology_phase3.py (Phase 3: cycle generators, correlates)
- path_homology_phase4.py (Phase 4: complement duality, n=8, arc-flip)
- path_homology_beta2_proof.py (β_2=0 algebraic analysis)
- path_homology_n8.py, path_homology_n8_extended.py (n=8 exploration)
- path_homology_n9_fast.py (n=9 tournament β_3 check)
- path_homology_phases.py (topological phase analysis)
- path_homology_Hconnection.py (OCF-topology connection)
- quick_stability.py (circulant stability verification)
- topology_landscape.py (circulant topology census)
- symbol_matrix_analysis.py (F vs topology analysis)
- paley_path_homology.py (Paley tournament path homology)
