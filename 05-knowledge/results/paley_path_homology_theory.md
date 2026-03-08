# Paley Tournament Path Homology — Theoretical Summary

## Sessions: S38, S39, S40

### Setup
For p ≡ 3 mod 4 prime, the Paley tournament P_p is the circulant tournament
on Z_p with connection set S = QR(p) (quadratic residues).

Properties:
- |S| = (p-1)/2
- S ∩ (-S) = ∅ (since -1 ∈ QNR when p ≡ 3 mod 4)
- P_p is self-complementary (complement has S' = QNR)
- P_p is vertex-transitive (circulant shift τ: v → v+1)

### Computational Results

| p | QR | β | Geometric type | χ |
|---|-----|---|---------------|---|
| 3 | {1} | (1,1,0) | S^1 (circle) | 0 |
| 7 | {1,2,4} | (1,0,0,0,6,0) | ∨^6 S^4 | 7 |
| 11 | {1,3,4,5,9} | Computing... | ? | predicted 11 |

### Per-Eigenspace Decomposition

The cyclic shift τ decomposes the chain complex into p eigenspaces
indexed by k = 0,...,p-1 with eigenvalue λ_k = ω^k (ω = e^{2πi/p}).

**P_7 (COMPLETE):**
- k=0: Ω dims = [1,3,6,9,9,6,3], β = (1,0,0,0,0,0,0)
- k=1,...,6: ALL have Ω dims = [1,3,6,9,9,6,3], β = (0,0,0,0,1,0,0)
- Total: β = (1,0,0,0,6,0,0), χ = 7

**P_11 (PARTIAL — dims 0-6 complete):**
- All eigenspaces k=1,...,10 have the SAME Ω dims (proved theoretically)
- Ω dims (k≠0) through dim 6: [1, 5, 20, 70, 205, 460, 700, ...]
- IMPORTANT: Inner sequence [5,20,70,205,460,700,...] is NOT palindromic
  (460 ≠ 700), unlike P_7's [3,6,9,9,6,3]

**P_3:**
- Special case: β_1 = 1 comes from k=0 eigenspace (not k≠0)
- This is because |S|=1 gives only 1 step sequence at each dim

### Why All Eigenspaces Have Same Ω Dims

**Theorem** (proved theoretically, verified computationally):
For Paley tournament P_p, all non-trivial eigenspaces (k=1,...,p-1) have
identical Ω_m dimensions.

**Proof sketch:**
1. For k ∈ QR: multiplication by k permutes QR, so the map
   φ_k: step sequence s → ks permutes step sequences. The junk matrix
   J_m(ω^k) = P^T J_m(ω) P for some permutation P. Hence same nullity.

2. For k ∈ QNR: multiplication by k sends QR → QNR, giving a bijection
   between QR-step-sequences and QNR-step-sequences. By self-complementarity
   of P_p, the complement tournament (with S = QNR) has the same structure.

3. The character sums confirm: Σ_{s∈QR} ω^{ks} has exactly ONE distinct
   value for all k ∈ QR, and one for all k ∈ QNR. Verified for
   p = 3, 7, 11, 19, 23, 31, 43, 47.

### Euler Characteristic

For P_p with p ≥ 7:
- Each non-trivial eigenspace has alternating sum of Ω dims = 1
- So per-eigenspace χ^(k) = 1 for k ≠ 0
- k=0 contributes χ^(0) = 1 (from β_0 = 1)
- Total χ = 1·1 + (p-1)·1 = p

This is CONSISTENT with β being concentrated at a single even dimension d:
χ = 1 + (p-1)·(-1)^d = p requires (-1)^d = 1, i.e., d even.
For p ≡ 3 mod 4, p-3 ≡ 0 mod 4, so d = p-3 is even. ✓

### Illegal Merged Steps F = (S+S)\S

For Paley: F = QNR exactly (verified p = 3, 7, 11).
- |F| = (p-1)/2 (maximal), L = ∅ (no legal merges when p > 3)
- This gives the "maximal F" case from the F-topology classification
- Maximal F ↔ highest-dimensional sphere topology

### Circulant Tournament Comparison (p = 7)

8 total circulant tournaments on Z_7:
- 6 have β = (1,1,0,...,0), χ = 0 → S^1
- 2 have β = (1,0,0,0,6,0), χ = 7 → ∨^6 S^4
  (these are QR and QNR, both Paley up to complement)

Only Paley has χ = p among all circulant tournaments at p = 7.

### Gauss Sum Connection

The Gauss sum g = Σ χ(a)ω^a satisfies:
- |g| = √p
- g² = -p (since p ≡ 3 mod 4)
- Σ_{s∈QR} ω^{ks} = (χ(k)·g - 1)/2

The character sum controls the junk matrix entries and hence the Ω dims.
The uniform eigenspace structure follows from g being a single algebraic
quantity that determines all character sums via the Legendre symbol.

### β_2 = 0 for Tournaments — Structural Analysis (S40)

**Theorem** (verified exhaustively through n=6, randomly at n=7):
For any tournament T, H_2^{path}(T) = 0.

**Key structural discoveries:**

1. **Ω_2 ≠ span(transitive triples).** Ω_2 also includes "cancellation chains":
   linear combinations of non-transitive 2-paths where non-allowed face terms
   cancel. Example: (a,b₁,c) - (a,b₂,c) with c→a, where the -(a,c) terms cancel.
   The gap dim(Ω_2) - |TT| ranges from 0 to 5 at n=5.

2. **Cancellation chains never individually in ker(∂_2)** (0 out of 2649 at n=5).
   But **mixed elements** (TT + cancellation parts) CAN be 2-cycles.
   So ker(∂_2|Ω_2) > ker(∂_2|TT) for many tournaments.

3. **DT 4-paths**: (a,b,c,d) with a→b→c→d AND a→c AND b→d.
   - 4-clique DT (also a→d): always individually in Ω_3
   - Non-clique DT (d→a): also individually in Ω_3 (all faces in A_2)
   - DT ⊋ 4-clique paths; the extra non-clique DT paths are needed

4. **Simplicial H_2 ≠ 0**: The simplicial complex of transitive subsets has
   nontrivial H_2 for 40/1024 tournaments at n=5. But GLMY H_2 = 0 because
   Ω_3 is larger than the simplicial 3-chain space.

5. **β_2 = 0 verified with correct Ω_2**: Using compute_omega_basis for the
   TRUE Ω_2 and Ω_3, β_2 = 0 confirmed at n=4,5,6 (all), n=7 (sample).

### β_2 = 0 Proof — Deeper Analysis (S41)

6. **DT = {p ∈ A_3 : all faces in A_2}** (proved for n≤5).
   This is the set of allowed 3-paths where every face is also allowed.
   Equivalent to: (a,b,c,d) with a→b→c→d, a→c, b→d.

7. **dim(Ω_2) = |TT| + Σ(mult-1)** where mult = multiplicity of each non-allowed
   face among non-TT 2-paths. EXACT formula, verified 100% at n=5,6.

8. **DT sufficiency FAILS at n=6**: For 960/32768 tournaments,
   im(∂_3|DT) < ker(∂_2|Ω_2). The gap is 1 or 2.
   - Gap=1 cases: 720, all have score (1,2,2,3,3,4), t3=6
   - Gap=2 cases: 240, all have score (2,2,2,3,3,3), t3=8
   Full Ω_3 (including cancellation 3-chains) STILL gives β_2=0.

9. **Cancellation 3-chains are essential**: At n=6, Ω_3 \ DT contains
   multi-term linear combinations (6-14 paths) where non-A_2 face terms
   cancel. Their boundaries fill the remaining ker(∂_2) gap.

10. **Cone construction**: Coning from a source vertex (out-deg n-1) kills ALL
    2-cycles at n=5. But for non-source vertices, the residual is nonzero.
    Multi-vertex cone needed for general tournaments.

### Circulant Census at n=11 (S41)

- 32 circulant tournaments on Z_11 with |S|=5
- Paley S=[1,3,4,5,9]: β=[1,0,0,0,0] through dim 4, χ=1
- NO other circulant has χ=11 (or even χ>1) through dim 4
- At p=7: only Paley (and complement) had χ=p among circulants
- Paley's sphere topology is unique among circulants

### Open Questions
1. What are the full Ω dims for P_11? (Ω_7=690 known; Ω_8 computing)
2. At which dimension d does β_d = p-1 for P_11?
3. Is d = p-3 always, or does it vary?
4. Why is P_7 palindromic but P_11 is not?
5. Can the Gauss sum theory predict d directly?
6. β_2 = 0 algebraic proof for general tournaments — WHAT APPROACH?
   - NOT simplicial acyclicity (simplicial H_2 can be nonzero)
   - NOT DT sufficiency alone (fails at n=6)
   - Must use FULL Ω_3 including cancellation 3-chains
   - Need homological argument working with chain complex directly
   - Possible: spectral sequence, Mayer-Vietoris, or filtration argument
