---
id: THM-639
title: The Hamiltonian-path (step-sequence) frame for runner families — (A) gap laws are step-reversal invariant, so unique minimizers are palindromes; (B) the AP is the unique wall-count minimizer (the "transitive tournament" of the runner world); (C) the gap-process law is Haar on the annihilator loop and the balanced lattice L_0(E) is a complete invariant (= rational-affine class = primitive step sequence up to reversal)
status: PROVED (A, B, C below — elementary/classical, proofs included). REFUTED-AS-ATTEMPTED (D): the single-parameter balanced-girth deficit bound fails on constants (rank-11 packing count; and 13-point integer sets are pigeonhole-forced to girth 4 in all bounded-range regimes) — documented as a dead end with the successive-minima repair direction.
source: mac-mini-2026-07-07-S43 (owner directive: apply the tournament tiling-model frame — viewing the object from one of its Hamiltonian paths — to lonely runners; HYP-4887)
depends_on:
  - THM-530   # mu_{1/7}, m_P, the (P,E) architecture
related:
  - HYP-4837  # PZ-on-U (S41): the balanced lattice and "pairs never contribute"
  - HYP-4797  # kps-S59: diameter monotonicity (the frame's subset/refinement move)
  - HYP-4817  # monad-S2: tail-diameter theorem
  - HYP-4791  # klein-S155: successor-digraph/tournament frame at the k=8 criticality
external: Kronecker–Weyl (closed subgroups of tori); three-distance theorem; the project's tournament tiling model (CLAUDE.md, "everything is the triangle").
---

# THM-639 — The Hamiltonian-path frame for runner families

**The dictionary being formalized.** In the tournament tiling model, an n-player
tournament is viewed from one of its Hamiltonian paths: the path's n−1 arcs are the
spanning tree (cut space, the score hierarchy), and the remaining C(n−1,2) arcs are the
tiles (cycle space), whose binary pattern determines everything interesting; the
transitive tournament is the all-aligned configuration. Here the same split is applied
to a runner family `E = {e_1 < … < e_k} ⊂ ℤ` (k = 13 for the LRC(14) pure-cluster leg):

| tournament world | runner world |
|---|---|
| players | runners (speeds), incl. the stationary runner |
| the binary relation (who beats whom) | the pairwise difference `e_j − e_i` (the relative runner) |
| base Hamiltonian path (n−1 arcs) | the sorted order: steps `s_i = e_{i+1} − e_i` (k−1 steps) |
| tiles = cycle space | the **balanced lattice** `L_0(E) = {m : Σ m_i e_i = 0, Σ m_i = 0}` (rank k−2) |
| transitive tournament (all tiles aligned) | the AP (all steps equal) — Theorem B |
| grid/transpose symmetry of tilings | step reversal `E ↦ −E` — Theorem A |
| iso classes (tilings / fibration) | lattice classes = rational-affine classes — Theorem C |
| tile flip / metagraph edge | order-cell wall crossing in the x-sweep (a collision `(e_j−e_i)x ∈ ℤ`) |

Throughout, `g(E,x)` is the multiset of circular gaps of the configuration
`C_E(x) = {frac(e·x) : e ∈ E}`, and a *gap functional* is any quantity of the form
`x ↦ F(g(E,x))` or its law under `x ~ U[0,1)` (so `μ_θ`, `E[U_θ]`, `E[maxgap]`, all
gap moments).

## Theorem A (step reversal; palindromic minimizers)

**(i)** For a.e. `x`, `g(−E, x) = g(E, x)` as multisets. Hence every gap functional has
identical law for `E` and `−E`.
**(ii)** In sorted order, `−E` (affine-normalized) has the REVERSED step sequence
`(s_{k−1}, …, s_1)`.
**(iii)** Consequently the minimizer set of any gap functional is closed under step
reversal, and **any minimizer that is unique as an affine class has a palindromic step
sequence**.

*Proof.* (i) For `e·x ∉ ℤ`, `frac(−e·x) = 1 − frac(e·x)`; so `C_{−E}(x)` is the image of
`C_E(x)` under the circle isometry `y ↦ −y`, which preserves the gap multiset. (ii) Sort
`−E`: `−e_k < ⋯ < −e_1`, consecutive differences `e_{k−j+1} − e_{k−j}` read the `s_i`
backwards; translation by `e_k` normalizes to partial sums of the reversed sequence.
(iii) Immediate from (i)+(ii): if `F(E) = min` then `F(−E) = min`; uniqueness forces
`E ∼ −E` (affine), i.e. steps = reversed steps. ∎

*Census (S43 script).* All the burst's mu/E[maxgap]-record families are palindromic:
AP `(1^12)`, death-star `(2,2,2,2,2,1,1,2,2,2,2,2)`, monad's record
`(2,2,2,2,1,1,1,1,2,2,2,2)`, the S41 PZ-minimizer `(2,2,1,1,1,1,1,1,1,1,2,2)`, the deep
interlace. The S41 E[U]-minimizer is NOT palindromic — so by (iii) it is not a unique
minimizer (it comes as a mirror pair), and indeed the palindrome-constrained floor sits
strictly higher (0.0988 vs 0.0938): **the E[U] landscape has a spontaneously broken
reversal symmetry**, unlike μ (whose unique minimizer, the AP, is a palindrome).
Practical use: palindromicity halves adversarial search spaces for symmetric objectives
and is a structural check on any claimed unique extremal.

## Theorem B (wall count; the AP is the "transitive tournament")

For sorted `E` with steps `s_1, …, s_{k−1}`:
> `Σ_{i<j} (e_j − e_i) = Σ_{r=1}^{k−1} s_r · r(k−r)`,
and over families of k distinct integers this is uniquely minimized (per affine class) by
the AP, with value `Σ_r r(k−r) = C(k+1, 3)` (= 364 at k = 13).

The quantity counts, with multiplicity, the collision walls `{x ∈ (0,1] : (e_j−e_i)x ∈ ℤ}`
of the x-sweep — the total number of "tile flips" the configuration's cyclic order
undergoes. So the AP is the family with the COARSEST order-cell complex: the minimal-
complexity frame, exactly as the transitive tournament is the minimal class of the
tiling model. (This is why the exact Farey-cell engines are tractable on the AP, why its
witness structure is three-distance, and — heuristically — why it is the master extremal
object of the whole program.)

*Proof.* `e_j − e_i = Σ_{r=i}^{j−1} s_r`; the step `s_r` is counted by the pairs
`i ≤ r < j`, of which there are `r(k−r)`. Each coefficient `r(k−r) > 0` and each
`s_r ≥ 1` (distinct integers), so the sum is `≥ Σ_r r(k−r)` with equality iff all
`s_r = 1`. Pair `(i,j)` collides exactly `e_j − e_i` times in `(0,1]`. ∎

## Theorem C (the gap law is Haar on the annihilator loop; complete invariant)

Let `E` be primitive (gcd of differences 1, wlog after affine normalization) and
`L(E) = {m ∈ ℤ^k : Σ m_i e_i = 0}`, `L_0(E) = L(E) ∩ {Σ m_i = 0}`.

**(i)** The phase curve `x ↦ (e_1 x, …, e_k x) mod 1` is a closed 1-dimensional subgroup
of `T^k`, equal to the annihilator `Ann(L(E))`, and the law of the phase vector under
`x ~ U[0,1)` is **Haar measure on that subgroup**.
**(ii)** The law of the gap multiset process is the pushforward of Haar on the image
subgroup in `T^k/T_diag`, whose annihilator (inside the quotient's character group
`{Σ m = 0}`) is `L_0(E)`. Hence **every gap functional's law is determined by `L_0(E)`**
(as a sublattice of `ℤ^k`, up to index permutation).
**(iii)** Conversely `L_0(E)` determines `span_ℚ{e, 𝟙}`, hence `E` up to rational-affine
maps and index permutation. So: **lattice class = rational-affine class = primitive step
sequence up to reversal** — the isomorphism classes of runner families, forced (as in the
tiling model) by the invariances of the pairwise relation itself: translation (`Σm = 0`)
and dilation (`Σ m_i e_i = 0` homogeneous).

*Proof.* (i) The map `φ(x) = x·e mod 1` is a continuous homomorphism `ℝ → T^k` with
`φ(1) = 0`, so it factors through `T¹` and its image is a compact (closed) subgroup.
By Pontryagin duality a closed subgroup equals the annihilator of its annihilator; the
annihilator of the image is `{m : m·(xe) ∈ ℤ ∀x} = {m : m·e = 0} = L(E)`. Primitivity
makes `φ` injective on `[0,1)`, and `φ` pushes Lebesgue to the normalized arc-length
(constant speed) on the image = Haar on a 1-dimensional compact group. (ii) Gap
multisets are invariant under the diagonal rotation, so the gap map factors through
`T^k/T_diag`; the pushforward of Haar under a quotient homomorphism is Haar on the
image, and duality of the quotient identifies its annihilator as `L(E) ∩ {Σm = 0}`.
(iii) `L_0` has rank `k−2` (two independent linear constraints on a non-constant
integer vector), so its rational orthogonal complement is the 2-plane `span{e, 𝟙}`;
recovering the plane recovers `e` up to `e ↦ ae + b·𝟙`, `a,b ∈ ℚ`; normalizing to a
primitive sorted integer vector leaves exactly the sign choice = step reversal. ∎

*Remarks.* (1) Theorem A(i) is the special case `L_0(−E) = L_0(E)` of C(ii). (2) The
project's piecemeal invariances (dilation, translation, boxeph's speed-vs-co-offset
reflection, "pairs never contribute": a support-2 balanced relation needs
`m(e_i − e_j) = 0`, impossible) are all instances of C. (3) C is the formal foundation
under the deficit program: the deviation of ANY gap law from the iid law is a function
of `L_0` alone, graded by the norms of its vectors — which is what S41 (mac-mini) and
S154/S155 (klein) compute layer by layer.

## D (attempted and refuted as stated): the balanced-girth tail bound

The natural quantitative form of C — "if `L_0(E)` has ℓ1-girth `W(E) ≥ W*` then
`|E[U_{1/7}] − (6/7)^13| ≤ T(W*) ≤ 0.086`, closing the sparse lane" — FAILS on
constants with the naive ingredients: the rank-11 packing count `(1+2n/W)^{11}`
contributes a `3^{11} ≈ 1.8·10^5` factor that no geometric decay `u^{n/2}` recovers
(computed table: `T(22) ≈ 323`). Moreover the hypothesis is nearly vacuous: 13 bounded
integers always carry additive quadruples (`e_a + e_b = e_c + e_d`, norm-4 balanced
relations), so `W(E) = 4` for every family in the program's range (census: all bank
families girth 4; a 13-element Sidon set only reaches girth 6). **Do not pursue the
single-parameter girth route.** The repair direction (open): the full successive-minima
profile `λ_1 ≤ … ≤ λ_{11}` of `L_0` with a covolume-aware count `Π(1 + 2n/λ_i)`, i.e.
"how many short cycles at each scale" — consistent with the S41 finding that the
deficit is a graded, signed, cross-weight object.

## Scripts
`04-computation/lrc14_hampath_frame_macmini_S43.py` (+ `.out` in 05-knowledge/results):
palindrome/wall/girth census, palindrome-constrained descent, T(W) table.
