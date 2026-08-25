---
id: THM-1950
title: "H ≥ disc REDUCED TO THE STRONGLY-CONNECTED BASE, via the skew-determinant SCC-composition law — the exact structural parallel of THM-1860 (c₃ ≤ H) for death-star's HYP-8636. For every tournament T, the Rédei Hamiltonian-path count H(T) is at least the poly-time skew-determinant disc(T)=|det(I+K)|/2^{n−1} (K=A−Aᵀ), with equality iff T is transitive. This is PROVED modulo the strongly-connected base H(C) ≥ max(1,s(C))·disc(C), where s(C)=𝟙ᵀ(I+K)⁻¹𝟙=‖(I+K)⁻¹𝟙‖² is the total inverse-response. The reduction rests on three PROVED facts: (i) the SCC-composition law disc(C₁⇒⋯⇒C_r)=∏disc(Cᵢ)·[∏(1+sᵢ)+∏(1−sᵢ)]/2^r (Schur complement; disc is SUPER-multiplicative, unlike H which is exactly multiplicative, and unlike char_A which factors — THM-1925); (ii) the velocity-addition law s(C₁⇒C₂)=(s₁+s₂)/(1+s₁s₂) (a Möbius/SL₂ action on the SCC composition); (iii) the elementary kernel inequality max(1,x)·max(1,y) ≥ max(1+xy,x+y)/2 for x,y≥0. Peeling the top strong component with the invariant P(T)=max(1,s(T))·disc(T) closes the induction: H(T)≥P(T)≥disc(T)."
status: >
  (0) VERIFIED: H ≥ disc exhaustive n ≤ 6 (33864 tournaments, 0 violations; reconfirms
  death-star HYP-8636) and sampled n = 7,8; equality locus = transitive only.
  (1) PROVED (this session), all machine-verified:
    - SCC-composition law for disc: exact to 1.8e-15 over all 10976 reducible tournaments n≤6.
    - s velocity-addition s(C₁⇒C₂)=(s₁+s₂)/(1+s₁s₂): derived by block solve, verified.
    - kernel inequality max(1,x)max(1,y) ≥ max(1+xy,x+y)/2: proved by 3 cases, 500k-sample slack ≥ 0.
    - the invariant statement H(T) ≥ P(T)=max(1,s(T))·disc(T): exhaustive n≤6 (all 2^{C(n,2)}).
  (2) THE RESIDUAL BASE (open, the whole remaining content): for strongly connected C,
    H(C) ≥ max(1,s(C))·disc(C). Verified exhaustively through n=7: n≤6 on labelled
    tournaments and n=7 on all 353 strong isomorphism classes (THM-4051), with 0 violations.
    Tight only at C₃ (H=3, s=3, disc=1 ⟹ max(1,s)disc=3). The former sampled n=7
    minimum ratio 4.22 and monotone-growth inference are SUPERSEDED: the exact minimum is
    27/8 at the Paley tournament, so no monotonic margin claim is licensed. Note s(C)≥1 FAILS for strong C at
    n≥7 (min s = 0.667,0.556 at n=7,8), which is exactly why the base must carry the max(1,·) and
    the kernel is the two-sided max-inequality, not the naive ∏sᵢ bound.
  (3) So this is "H ≥ disc reduced to strongly-connected", NOT "H ≥ disc proved" — the exact
    status of THM-1860 (c₃ ≤ H) for its sibling conjecture; the composition machinery is fully
    proved and the base is the residual.
source: klein-2026-07-21-S400 (owner: synthesise repo progress; work high-leverage open math to the fundamental; small improvements toward proofs)
depends_on:
  - THM-1860    # c₃ ≤ H reduced to strong — the structural template
related: [THM-1858, THM-1925, THM-1810, THM-1885, THM-474]
answers: HYP-8636   # death-star-S78's H ≥ disc conjecture — now reduced to the strong base
script: 04-computation/h_ge_disc_reduction_to_strong_klein_S400.py, disc_composition_law_klein_S400.py, h_ge_disc_reduction_klein_S400.py (+ .out)
---

# THM-1950 — H ≥ disc, reduced to the strongly-connected base

death-star-S78 (HYP-8636) conjectured `H(T) ≥ disc(T)` with equality iff transitive, where
`H(T)` = # directed Hamiltonian paths (Rédei, #P-hard, always odd) and
`disc(T) = |det(I+K)|/2^{n−1} = ∏_j(1+μ_j²)/2^{n−1}` (`K = A−Aᵀ`, `±iμ_j` its eigenvalues) is the
poly-time skew-determinant (THM-474). It is the WOWII "easy invariant bounds hard invariant" shape.
Here it is reduced to a clean strongly-connected base, with the entire SCC-composition machinery
proved — the exact structural parallel of THM-1860.

## The engine: the total inverse-response s

For a tournament `T` set `x = (I+K)⁻¹𝟙` and `s(T) = 𝟙ᵀx`. Because `K` is skew,
`s = xᵀ(I+K)x = xᵀx = ‖x‖² ≥ 0`, and Cauchy–Schwarz gives `s ≤ n`; the singular values of `I+K`
are `√(1+σ²)≥1`, so `‖(I+K)⁻¹‖≤1`. A singleton has `s=1`.

## The three proved facts

1. **SCC-composition law (Schur complement).** For the ordered sum `T = T₁ ⇒ T₂` (all of `T₁`
   beats all of `T₂`), writing `I+K = [[I+K₁, J],[−Jᵀ, I+K₂]]` with `J` all-ones, the Schur
   complement on the `(1,1)` block gives `det(I+K)=det(I+K₁)·det(I+K₂+s₁J)` and then (matrix
   determinant lemma) `det(I+K)=det(I+K₁)det(I+K₂)(1+s₁s₂)`. Dividing by `2^{n−1}`:
   > `disc(T₁ ⇒ T₂) = disc(T₁)·disc(T₂)·(1+s₁s₂)/2.`
   Iterating over the transitive SCC order `C₁ ⋗ ⋯ ⋗ C_r`:
   > `disc(T) = (∏ᵢ disc(Cᵢ))·[∏ᵢ(1+sᵢ) + ∏ᵢ(1−sᵢ)]/2^r.`
   Since `sᵢ ≥ 0`, the correction `≥ 1`: **disc is super-multiplicative over strong components**
   (contrast: `H` is exactly multiplicative — THM-1860; `char_A` factors — THM-1925; `disc` does
   *not*). Verified to `1.8e-15` over all 10976 reducible tournaments `n ≤ 6`.

2. **Velocity-addition law for s.** Solving the block system gives
   > `s(T₁ ⇒ T₂) = (s₁ + s₂)/(1 + s₁ s₂),`
   the relativistic velocity-addition / Möbius law (`s = tanh` of an additive "rapidity"). So the
   SCC composition acts on `s` by `SL₂` Möbius maps — the disc/`s` pair sits inside the a/b `SL₂`
   frame (THM-1810/1885). (This is why reducible `s` can drop below 1 even when every component
   has `s ≥ 1`.)

3. **Kernel inequality (elementary, 3 cases; FORMALIZED in Lean).** For all `x, y ≥ 0`,
   > `max(1,x)·max(1,y) ≥ max(1+xy, x+y)/2.`
   Cases `x,y≥1` (⟺ `xy≥1`), `x,y≤1` (⟺ `xy≤1`), `x≥1≥y` (⟺ `x≥y`) each reduce to a triviality.
   Machine-checked, `sorry`-free, kernel-pure `[propext, Classical.choice, Quot.sound]`:
   `HgeDiscKernel.kernel_ineq` and the peel form `HgeDiscKernel.peel_step`
   (`04-computation/lean/TournamentH7/TournamentH7/HgeDiscKernel.lean`, in the root manifest) —
   the exact analogue of THM-1860's `SumLeProd.lean` arithmetic kernel.

## `s` is a regularity coordinate: `s(T)=n ⟺ T regular`

Because a regular tournament has `K𝟙 = 0` (row sums `= out−in = 0`), `x=(I+K)⁻¹𝟙=𝟙` and `s=n`; and
`s=n` forces `x∥𝟙` (equality in Cauchy–Schwarz `s=𝟙ᵀx≤√n‖x‖=√n√s`) hence `K𝟙=0`, i.e. regular.
Verified exhaustively n=3,5,7 (0 mismatches). So `s ∈ [0,n]` with the top value pinned to regularity;
regular tournaments (Paley, rotational) uniquely maximize `s`, and there the base reads `H ≥ n·disc`
(e.g. Paley-7: `H=189 ≥ 7·8=56`). Under `⇒` the velocity law drives `s` into `(−1,1)` (reducible
`s<1`), so `s` measures how far the top strong component is from regular, in a Möbius-additive
coordinate.

## The reduction

Let `P(T) = max(1, s(T))·disc(T)`. **Claim: `H(T) ≥ P(T)` for every tournament, given the base.**

- **Base (strong `C`):** `H(C) ≥ max(1,s(C))·disc(C) = P(C)`. *(the residual, verified n≤7.)*
- **Step:** peel the top strong component, `T = C₁ ⇒ T'`. Then `H(T)=H(C₁)H(T')` (Rédei/SCC) and
  `disc(T)=disc(C₁)disc(T')(1+s₁s')/2` with `s'=s(T')`. Using the base and the inductive
  hypothesis `H(T')≥P(T')=max(1,s')disc(T')`,
  `H(T) ≥ max(1,s₁)max(1,s')·disc(C₁)disc(T')`
  `     ≥ [max(1+s₁s', s₁+s')/2]·disc(C₁)disc(T')`  (kernel inequality)
  `     = max(1,s(T))·disc(C₁)disc(T')(1+s₁s')/2`  (since `s(T)(1+s₁s')=s₁+s'`)
  `     = max(1,s(T))·disc(T) = P(T).`
- Hence `H(T) ≥ P(T) ≥ disc(T)`. ∎ (modulo the base)

The peel treats `T'` as a single (possibly reducible) block; the Schur law needs only that `C₁`
beats all of `T'`, which the top strong component does. Verified in invariant form `H ≥ P(T)`
exhaustively for all `2^{C(n,2)}` tournaments `n ≤ 6`.

## The residual base and its texture

All that remains is **`H(C) ≥ max(1,s(C))·disc(C)` for strongly connected `C`.** When `s(C) ≥ 1`
this is the genuinely stronger `H ≥ s·disc`; when `s(C) < 1` (which happens for strong `C` at
`n ≥ 7`) it is just `H ≥ disc`. Tight only at `C₃`. The base is the sibling of THM-1860's open
strong base `c₃ ≤ H`: a strongly-connected tournament is vertex-pancyclic (Moon) with many
Hamiltonian paths, so `H` is large where `disc` (bounded by `∏(1+μ_j²)/2^{n−1}` with the fixed
energy `Σμ_j² = C(n,2)`) is small. A likely route: bound `s·disc` by a Hamiltonian-path injection,
or an eigenvalue-product argument on the strong spectrum (Perron + the isotropic pairs, THM-1858).

**Exact order-seven correction (THM-4051).** The complete 456-class quotient
universe contains 353 strong classes, all satisfying the base strictly. Its
minimum additive slack is `21`, but its minimum ratio is `27/8`, attained by
the regular Paley tournament with `(H,disc,s·disc)=(189,8,56)`. The old
sampled value near `4.22` missed this hostile symmetry class; in particular,
the finite minima do not justify a monotone-growth claim.

## Why this matters

- It puts death-star's HYP-8636 on the same rigorous footing as `c₃ ≤ H`: a proved reduction to a
  strongly-connected base plus a proved algebraic kernel, honest about the open base.
- The **SCC-composition law for disc** and the **velocity-addition law for s** are new structural
  theorems in their own right — they explain *why* disc is not SCC-multiplicative (the fleet's
  multiplicativity results are all for `A`/`char_A`, not `K`/disc) and tie the skew-determinant to
  the `SL₂`/Möbius (a/b, THM-1885) frame.
