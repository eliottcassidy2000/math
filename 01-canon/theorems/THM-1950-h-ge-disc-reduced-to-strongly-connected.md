---
id: THM-1950
title: "Hamiltonian path determinant bound: global vertex-cycle induction, proper-cycle defect and weighted extension"
status: >
  PROVED, Session137 addendum (2026-11-02), with independent proof audits and
  exact execution controls. For every nonempty tournament,
  H >= det(I+A)+abs(det A)+C_nonHamilton >= max(1,s)*disc >= disc.
  Equality in the first bound holds exactly when no two directed triangles
  are vertex-disjoint; H=disc exactly for transitive tournaments.
  The weaker bound H>=det(I+A)+abs(det A) extends to weighted tournaments
  Aij in [0,1], Aij+Aji=1, with the stated low-order equality exceptions.
  This closes the formerly open strong base without assuming the SCC
  reduction. The historical Schur/Mobius formulas remain valid, but the
  unconditional supermultiplicativity description is corrected below.
  This is not a formal Lean verification or a full value-spectrum theorem.
source: klein-2026-07-21-S400 (owner: synthesise repo progress; work high-leverage open math to the fundamental; small improvements toward proofs)
depends_on:
  - THM-1860    # c₃ ≤ H reduced to strong — the structural template
related: [THM-1858, THM-1925, THM-1810, THM-1885, THM-474]
answers: HYP-8636   # death-star-S78's H ≥ disc conjecture — now reduced to the strong base
script: 04-computation/h_ge_disc_reduction_to_strong_klein_S400.py, disc_composition_law_klein_S400.py, h_ge_disc_reduction_klein_S400.py (+ .out)
---

# THM-1950 — current global closure

**PROVED, independently audited, 2026-11-02.** The [Session137 vertex-cycle induction](../../05-knowledge/results/session137_overnight_proof_frontiers_20261102/discriminant/GLOBAL_PROOF.md) proves the previously open base for every nonempty tournament, not just strong tournaments. Write u=det(I+A), v=(-1)^n det A and P=u+|v|=max(1,s)disc. Then

    H >= P+C_nonHamilton >= P >= disc.

The sharp proper-cycle equality holds exactly when there are no two vertex-disjoint directed triangles. H=P exactly for transitive tournaments and C3; H=disc exactly for transitive tournaments. The weighted extension H>=det(I+A)+|det A| and its separate equality classification are proved in the linked source. Empty tournaments are excluded from the P bound.

The proof separates determinant cycle signs from the inherited OCF and uses the weaker induction H>=u, whose empty base is valid. Algebraic u±v positivity supplies both margins. [Audit and reproduction](../../05-knowledge/results/session137_overnight_proof_frontiers_20261102/README.md). No unrestricted value-spectrum realization or literature-priority claim follows.

**Correction to historical prose.** The exact join factor is (1+s1*s2)/2, so disc is supermultiplicative only when s1*s2>=1. A strong seven-vertex factor with rows(30,8,122,80,2,91,19) has disc6,s2/3, and its ordinal square has disc26<36. Reducibility does not force s<1: three ordinal triangles have s9/7. The original Schur law and max-kernel reduction survive these corrections.

## Historical reduction and finite evidence

The following records the earlier reduction and its provenance. Its statements that the strong base remained open are superseded by the current addendum above.



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
   The correction need not be at least one: **disc has the displayed modulated product law, not unconditional supermultiplicativity** (see the current counterexample above)
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
(e.g. Paley-7: `H=189 ≥ 7·8=56`). Under `⇒` the velocity law can put `s` on either side of one (the ordinal-triangle family alternates), so `s` measures how far the top strong component is from regular, in a Möbius-additive
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
