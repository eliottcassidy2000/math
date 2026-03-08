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

### Open Questions
1. What are the full Ω dims for P_11? (dims 7-10 computing)
2. At which dimension d does β_d = p-1 for P_11?
3. Is d = p-3 always, or does it vary?
4. Why is P_7 palindromic but P_11 is not?
5. Can the Gauss sum theory predict d directly?
6. β_2 = 0 algebraic proof for general tournaments
