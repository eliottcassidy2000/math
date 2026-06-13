# THM-447 — Skew-Sylvester Doubling of Tournaments: three copies and one negated copy

**Status:** PROVED (block algebra below) — computational verification in progress (skew_doubling_core_kps1.py).
**Source:** kind-pasteur-2026-06-09-S1. User directive: "n*2 tournament doubling recursion … three
copies of the original, and one negated copy … skew Hadamard matrices, normalized with all elements
in the first row equal to 1, analogous to the tiling model with the Hamiltonian path fixed."
**Related:** OPEN-Q-045 (T[K_2] blowup — all-positive doubling), OPEN-Q-044 (SC blowup — transposed
cross blocks), THM-067 (Mersenne vanishing), THM-250..253 (spectral bridge / rapidity), T767,
HYP-2332..2337.

## Setup

For a tournament T on n vertices let A be its 0/1 adjacency matrix and
**M = A − A^T** its dominance matrix: skew-symmetric, entries ±1 off the diagonal, 0 on the
diagonal. Set **S = M + I** (the ±1 "skew-type" matrix: S + S^T = 2I). T is recovered from M.
A skew-Hadamard matrix of order m is a ±1 matrix S with S + S^T = 2I and S S^T = mI.

## Definition (the skew-Sylvester double)

**D(T)** is the tournament on 2n vertices {1..n} ∪ {1'..n'} with dominance matrix

```
M' = [ M     M+I ]
     [ M−I   −M  ]
```

Concretely: all three vertex-pair types (i,j), (i,j'), (i',j) with i≠j are oriented as in T
(**three copies of the arc set**); pairs (i',j') are oriented as in T^op (**one negated copy**);
and each vertex beats its own twin: i → i'.

This is the exact skew analog of Sylvester/Walsh doubling H_{2m} = [[H, H], [H, −H]]: the four
blocks are M, M, M, −M up to the forced diagonal corrections ±I (which encode the twin arcs).

## Claims and proofs

**(1) M' is a valid dominance matrix (D(T) is a tournament).**
M'^T = [[M^T, (M−I)^T], [(M+I)^T, −M^T]] = [[−M, −M−I], [−M+I, M]] = −M'. So M' is
skew-symmetric; its diagonal blocks have zero diagonal (M has zero diagonal), and every
off-diagonal entry is ±1: within blocks M, ±1 off-diagonal; in block M+I the diagonal is +1
(twin arcs i → i') and off-diagonal entries are M_ij = ±1. ∎

**(2) S' = M' + I_{2n} satisfies S' + S'^T = 2I; if S = M + I is skew-Hadamard of order n, then
S' is skew-Hadamard of order 2n.**
In block form S' = [[S, S], [S−2I, 2I−S]]. Skew-type: (S) + (S)^T = 2I on diagonal blocks
(2I−S has (2I−S)+(2I−S)^T = 4I−2I = 2I), and off-diagonal: S + (S−2I)^T = S + S^T − 2I = 0. ✓
Orthogonality (assuming S S^T = nI):
- (1,1): S S^T + S S^T = 2nI ✓
- (1,2): S(S−2I)^T + S(2I−S)^T = S(S^T−2I) + S(2I−S^T) = 0 ✓
- (2,2): (S−2I)(S−2I)^T + (2I−S)(2I−S)^T = 2(S S^T − 2(S+S^T) + 4I) = 2(nI − 4I + 4I) = 2nI ✓
So S' S'^T = 2n·I_{2n}. ∎
(The base case S = [1] of order 1 is skew-Hadamard; iterating D gives a **skew-Walsh tower** of
orders 2^k. Compare: plain Sylvester doubling does NOT preserve skewness — [[H,H],[H,−H]] of a
skew H has (H'+H'^T) = 2I ⊗ J ≠ 2I.)

**(3) Spectral law — Chebyshev T_2.**
```
M'^2 = [ M² + (M+I)(M−I)         M(M+I) + (M+I)(−M)  ]   = [ 2M² − I     0      ]
       [ (M−I)M + (−M)(M−I)      (M−I)(M+I) + M²     ]     [ 0           2M² − I ]
```
i.e. **M'² = I_2 ⊗ (2M² − I)**. Hence every eigenvalue μ of M (pure imaginary, μ = iλ) yields
eigenvalues μ' of M' with **μ'² = 2μ² − 1 = T_2(μ)** — the degree-2 Chebyshev polynomial — each
with doubled multiplicity. In skew form: **λ → √(2λ² + 1)**, so the quantity

```
λ² + 1   EXACTLY DOUBLES under D.
```

For a DRT-type spectrum λ² + 1 = (Hadamard order), this is conservation of the Hadamard order
under doubling: (2λ²+1) + 1 = 2(λ²+1). Iterating k times from λ:
λ_k² = 2^k (λ² + 1) − 1. From the 1-vertex tournament (λ=0): λ_k = √(2^k − 1) — the Mersenne
spectral radii. Asymptotic stretch factor √2 — the triangle's hypotenuse/leg ratio.

**(4) Scores.** score_{D(T)}(i) = s_i (copy 1) + s_i (cross, beats those it beat) + 1 (twin)
= **2s_i + 1**; score_{D(T)}(i') = (n−1−s_i) (copy 2, T^op) + s_i (cross) = **n − 1**,
independent of T. The second copy is "flattened" to constant score n−1 regardless of T.
If T is regular on n = 2k+1 vertices, D(T) is near-regular on 2n with scores n−1, n each n times.

**(5) Canonical Hamiltonian path (the Walsh/Gray reflection).** If P: n → n−1 → ⋯ → 1 is a
Hamiltonian path of T, then

```
n → n−1 → ⋯ → 1 → 1' → 2' → ⋯ → n'
```

is a Hamiltonian path of D(T): the twin arc 1 → 1' exists by construction, and k' → (k+1)' in
copy 2 = T^op iff k+1 → k in T, which is the base-path arc. The doubled base path = (path,
twin-arc, REVERSED path) — the reflect-and-append of Gray/Walsh ordering.

**(5-CORRECTED, in-session; see MISTAKE-065 and THM-452.)** In this frame the tiling of D(T)
splits as: copy-1 block = tiles x of T (identity copy); copy-2 block = **grid transpose of x
WITHOUT negation** — the T^op negation is exactly absorbed by traversing copy 2 along the
reversed path (the original claim "grid-transposed negated" was wrong); cross block = T's full
ordered arc matrix A[i][j], whose σ-partner pairs (i,j)↔(j,i) carry **complementary** bits —
this is where the single negated Sylvester copy lives — plus an all-ones twin anti-diagonal
(the hypotenuse). Net: three positive copies + one negated copy, exactly [[H,H],[H,−H]] at tile
level, with the negation in the cross block, not the copy-2 block. Verified 100% (1098 framed +
1096 labeled cases, n=3..6). Full statement: THM-452 (σ-coset law).

## Verification record (this session)

- C1–C4 (tournament+scores, M'² law, skew-type, canonical path): exhaustive n=3,4,5 — ALL PASS
  (`skew_doubling_core_kps1.out`).
- Skew-Hadamard tower clean orders 2..64 (`drt_mersenne_tower_kps1.out`).
- Eigenvalue map λ → √(2λ²+1) numerically exact (sample + tower spectra: T15 has
  spec {0, ±i√15 ×7}, M² = J − 15I exact).
- D(border(C_3)) normalized core ≅ Paley T_7, H = 189 (explicit isomorphism found).

## The normalization dictionary (skew-Hadamard ↔ tiling model)

A skew-Hadamard matrix of order m, normalized so the first row is all +1 (always possible by
column negations + skew-restoring row negations), has: first column −1 below the corner, diagonal
+1 — the remaining (m−1)×(m−1) core is S(T) of a doubly regular tournament T on m−1 vertices.
The border row of 1s = an added dominating source vertex. Fixing additionally the superdiagonal
of the core (the base Hamiltonian path arcs, n−1 of them) leaves exactly C(n−1,2) free ±1 bits =
the **tiles**. So:

```
normalized first row        ↔  dominating source / score hierarchy anchor (cut space)
fixed base-path arcs        ↔  the staircase frame (the Hamiltonian path of the tiling model)
free C(n−1,2) core bits     ↔  the tiles (cycle space, the hypercube Q_m)
```

The user's analogy is exact: **normalizing a skew-Hadamard matrix = fixing the Hamiltonian path
in the tiling model.** The doubling D acts on this frame as the Sylvester step.

## Verification

- skew_doubling_core_kps1.py: exhaustive n=3..5 checks of (1)–(5) — IN PROGRESS, see
  05-knowledge/results/skew_doubling_core_kps1.out.
- THM-448 (stub): the DRT/Mersenne tower C_3 → Paley T_7 → DRT_15 → ⋯ via bordered doubling.
