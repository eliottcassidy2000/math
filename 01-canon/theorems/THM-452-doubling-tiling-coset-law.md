# THM-452 — The Tiling Coset Law: D(T)'s tiling has CONSTANT grid-symmetry defect; grid transpose implements op; the tower's path is the Gray code

**Status:** PROVED (block-transform derivation) + VERIFIED 100% (1098 framed + 1096 labeled
cases n=3..6; tower orders 2..16 exact). The original "blue conjecture" version is REFUTED
(see MISTAKE-065). Adversarially re-verified (`verify_C_tiling_kps1.out`).
**Source:** kind-pasteur-2026-06-09-S1 (branch C + verifier). Resolves HYP-2335.
**Related:** THM-447(5), blue/black strict definitions (CLAUDE.md), grid-sym fraction formula,
THM-024 (SC ⟹ involutive anti-aut).

## Setup

Frame D(T) by its canonical Hamiltonian path (THM-447(5)): relabel u_i ↦ n+i, v_j' ↦ n+1−j
(path position t gets label 2n−t). Tiles of the 2n-staircase split into: copy-1 block (labels
n+1..2n), copy-2 block (labels 1..n), cross block (n²−1 tiles, containing the n−1 twin tiles
on the anti-diagonal X+Y = 2n+1). Let σ = grid-symmetry involution (X,Y) ↦ (2n+1−Y, 2n+1−X),
x = tiling vector of T in its base frame.

## (1) The block transform (verified 100%)

```
copy-1 block:  t(n+i, n+j)   = x(i,j)            identity copy
copy-2 block:  t(X, Y)       = x(σ_n(X,Y))       grid transpose, NO negation
cross block:   t(n+i, n+1−j) = A[i][j]  (i≠j)    T's full ordered arc matrix
twin tiles:    ≡ 1                                all-ones anti-diagonal (the hypotenuse)
```
The T^op negation in copy 2 is EXACTLY absorbed by traversing the reversed path (MISTAKE-065:
"grid-transposed negated" was wrong). The negated Sylvester copy lives in the CROSS block:
σ-partner tiles (i,j) ↔ (j,i) carry the complementary bits A[i][j] vs A[j][i]. Net tile-level
schema = three positive copies + one negated = [[H,H],[H,−H]]. ✓

## (2) The coset law (the repaired blue conjecture)

The canonical tiling of D(T) is grid-symmetric for NO tournament (0/2194 tested — the blue
conjecture fails uniformly). Instead, the defect is CONSTANT, independent of T:
```
σ(t(D(T))) = t(D(T)) XOR c_n,     c_n = indicator of the n(n−1) off-twin cross tiles
```
Diagonal blocks are σ-symmetric (σ swaps copy-1 ↔ copy-2 bits, equal by (1)); the cross block
is σ-ANTI-symmetric. So x ↦ t(D(x)) is an injective GF(2)-affine map of rank C(n−1,2) whose
image is a section of the SINGLE fixed coset (blue subspace) + c_n. Dim(blue subspace at 2n) =
n(n−1) = |support(c_n)| (a numerical coincidence worth watching). Relative density of the
D-image in its coset: 2^{−(n−1)(n+2)/2}.

## (3) Grid transpose = the op-functor; op acts on the D-image as translation by c_n

σ(t(D(T))) is precisely the canonical tiling of copy-swap-relabeled D(T)^op
(dominance [[M, −M+I], [−M−I, −M]] = P(−M')P). Verified 74/74. It is never itself a D-image
(that would force A[i][j] = A[j][i]).

## (4) Blue-framability is (almost) Paley-exclusive

Over ALL Hamiltonian-path frames of D(T): blue (grid-symmetric) frames exist ONLY for T = C₃
(9/45 paths blue) and T = regular T₅ (185/15505) among the 18 iso classes n ≤ 5; all other
doubles are PURE BLACK. Provable necessary condition (verifier strengthening): blue-framability
of D(T) requires n odd AND T regular (score-multiset reversal argument). And the tower:
**T15 has ZERO anti-automorphisms of any order** (complete search) ⟹ its class is pure black,
NOT self-converse — unlike Paley DRTs. Tower doubles W₄, W₈, W₁₆ likewise have 0
anti-involutions.

## (5) The tower's canonical path is the binary reflected GRAY CODE

The recursive (path, twin, reversed path) construction = reflect-and-append = Gray code; twin
arcs at doubling scale d are exactly the Gray bit-flips at that scale. Closed-form arc law of
the tower (exact, orders 2..16): a → b iff b_L XOR parity(popcount((a AND b) >> (L+1))) = 1,
L = lowest differing bit — the "skew-Walsh function" (Sylvester character truncated at the twin
scale). Pure Walsh parity(popcount(a AND b)) describes only 52/105 of W16's tiling bits —
REFUTED as descriptor; the truncated form gets 105/105.

## (6) The fractal hypotenuse and the tiling weight recursion

In the tower member W_{2^k}'s tiling, EVERY scale's twin anti-diagonal appears as an all-ones
fractal hypotenuse. Tiling Hamming weight obeys w_{2n} = 2w_n + (n−1) + C(n,2):
w = 0, 0, 2, 13, 61 at orders 2,4,8,16 (61/105 for W16). T15 insertion-frame weight 48/91.

## Scripts

`skew_doubling_tiling_walsh_kps1.py`, `verify_C_tiling_kps1.py` (+ .out; ASCII tiling grids in
the .out show the fractal hypotenuse).
