# THM-452 — The Tiling Coset Law: D(T)'s tiling has CONSTANT grid-symmetry defect; grid transpose implements op; the tower's path is the Gray code

**Status:** PROVED (block-transform derivation) + VERIFIED 100% (1098 framed + 1096 labeled
cases n=3..6; tower orders 2..16 exact). The original "blue conjecture" version is REFUTED
(see MISTAKE-065). Adversarially re-verified (`verify_C_tiling_kps1.out`).
**Source:** kind-pasteur-2026-06-09-S1 (branch C + verifier). Resolves HYP-2335.
**Related:** THM-447(5), blue/black strict definitions (CLAUDE.md), grid-sym fraction formula,
THM-024 (SC ⟹ involutive anti-aut).

**2026-11-07 correction (Session222, independently audited):** sections1 and6 use
FORWARD-bit1, the global complement of the usual BACKWARD-bit1 gauge. The constant
reflection syndrome is unchanged by that switch, but is not itself a representative
of its nonzero affine fiber. Nonblue/disjoint-image assertions require input n>=2.
The five tower weights belong to orders1,2,4,8,16. See the
[audit and exact controls](../../05-knowledge/results/session222_tiling_recursions_20261107/audit/REPORT.md#3-thm-452-correction-and-safe-patch).

## Setup

In this file's block formulas let bit1 mean a FORWARD nonconsecutive arc along the
path, and let x use that same gauge. This differs from THM-280/554's backward-bit1 convention.

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

For the usual BACKWARD-bit1 convention complement both input and output tilings:
the two copy formulas remain identity/mirror copies, off-twin cross entries become
1−A[i][j], and twin tiles become0. Every twin arc is forward; one is a fixed-path
edge and the other n−1 are free tiles.

## (2) The coset law (the repaired blue conjecture)

For input n>=2 the canonical tiling of D(T) is grid-symmetric for NO tournament
(the original tested universe had0/2194). For n=1 there are no free tiles and the
tiling is blue. At every order the defect is CONSTANT, independent of T:
```
σ(t(D(T))) = t(D(T)) XOR c_n,     c_n = indicator of the n(n−1) off-twin cross tiles
```
Diagonal blocks are σ-symmetric (σ swaps copy-1 ↔ copy-2 bits, equal by (1)); the cross block
is σ-ANTI-symmetric. So x ↦ t(D(x)) is an injective GF(2)-affine map of rank C(n−1,2) whose
image lies in the SINGLE affine fiber `{t : (1+σ)t=c_n}=t0+B_(2n)`, where
`B_(2n)=ker(1+σ)` and t0 is any one D-image, for example the canonical double of
the transitive input in the declared gauge. Since σ(c_n)=c_n, c_n is itself blue;
the former expression `B_(2n)+c_n` equals B_(2n), not the desired nonzero fiber.
Dim(blue subspace at 2n) =
n(n−1) = |support(c_n)| (a numerical coincidence worth watching). Relative density of the
D-image in its coset: 2^{−(n−1)(n+2)/2}.

## (3) Grid transpose = the op-functor; op acts on the D-image as translation by c_n

σ(t(D(T))) is precisely the canonical tiling of copy-swap-relabeled D(T)^op
(dominance [[M, −M+I], [−M−I, −M]] = P(−M')P). Verified 74/74. For n>=2 it is never itself a D-image
(that would force A[i][j] = A[j][i]). At n=1 the unique empty tiling is a D-image.

## (4) Blue-framability is (almost) Paley-exclusive

Over ALL Hamiltonian-path frames of D(T): blue (grid-symmetric) frames exist ONLY for T = C₃
(9/45 paths blue) and T = regular T₅ (185/15505) among the 18 iso classes with 3 ≤ n ≤ 5;
all other doubles in that tested universe are PURE BLACK. The singleton input n=1
has a blue empty tiling. Provable necessary condition (verifier strengthening): blue-framability
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

In the FORWARD-bit1 gauge, every scale's twin anti-diagonal is all ones. Weight obeys
w_f(2n)=2w_f(n)+(n−1)+C(n,2), with w_f=0,0,2,13,61 at orders1,2,4,8,16
(61/105 for W16). The historical T15 insertion-frame weight48/91 is in this gauge.

In the standard BACKWARD-bit1 gauge the same twins are zeros, and
w_b(2n)=2w_b(n)+C(n,2), with w_b=0,0,1,8,44 at orders1,2,4,8,16.
Indeed w_b(n)+w_f(n)=C(n−1,2). Each copy contributes its old weight; the off-twin
cross tiles form C(n,2) complementary pairs; the n−1 free twins contribute only
in the forward gauge. The reflection-syndrome law is unchanged because σ(1)=1.
These are tiling Hamming weights, not Hamiltonian-path counts.

## Scripts

`skew_doubling_tiling_walsh_kps1.py`, `verify_C_tiling_kps1.py` (+ .out; ASCII tiling grids in
the .out show the fractal hypotenuse).
