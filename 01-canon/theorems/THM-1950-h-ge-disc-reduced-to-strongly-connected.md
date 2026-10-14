---
id: THM-1950
title: "H >= disc reduced to the strongly-connected base via exact skew-determinant and inverse-response composition"
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
    In the audited strong orders 3 through 7, tight only at C₃ (H=3, s=3, disc=1 ⟹ max(1,s)disc=3). The former sampled n=7
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

**2026-10-14 response-domain correction.** The exact Schur laws and the
conditional reduction below are unchanged. The former universal
super-multiplicativity, real `SL_2`/rapidity, and reducible `s<1` claims were
false. Canonical order-seven representative `299` has `(H,disc,s)=(31,6,2/3)`;
adjoining a sink lowers `disc` to `5`. See the exact witness and minimality
certificate in
[`tournament_response_correction_20261014.md`](../../05-knowledge/results/tournament_response_correction_20261014.md).

death-star-S78 (HYP-8636) conjectured `H(T) ≥ disc(T)` with equality iff transitive, where
`H(T)` = # directed Hamiltonian paths (Rédei, #P-hard, always odd) and
`disc(T) = |det(I+K)|/2^{n−1} = ∏_j(1+μ_j²)/2^{n−1}` (`K = A−Aᵀ`, `±iμ_j` its eigenvalues) is the
poly-time skew-determinant (THM-474). It is the WOWII "easy invariant bounds hard invariant" shape.
Here the inequality is reduced to a clean strongly-connected base, with the entire SCC-composition machinery
proved — the exact structural parallel of THM-1860. The conjectured all-order equality characterization
is not established by the weak base inequality alone.

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
   The two-block correction is at least one **iff `s₁s₂>=1`**; mere
   nonnegativity is insufficient. In particular `disc` is not universally
   super-multiplicative over strong components. If every component has
   `sᵢ>=1`, the displayed many-component factor is at least one: expand its
   numerator as twice the sum of the `2^(r-1)` even-cardinality products.
   The formula itself was verified numerically to `1.8e-15` over the recorded
   reducible tournaments through order six. The correction certificate is exact.

2. **Velocity-addition law for s.** Solving the block system gives
   > `s(T₁ ⇒ T₂) = (s₁ + s₂)/(1 + s₁ s₂),`
   an exact rational composition law on nonnegative responses. For fixed
   `s₁`, the representing matrix `[[1,s₁],[s₁,1]]` has determinant `1-s₁²`:
   the map is constant at `s₁=1`, and for `s₁>1` cannot be rescaled into
   real `SL_2`. Real additive rapidity `arctanh(s)` applies only in
   `0<=s<1`. The globally defined coordinate
   `q(T)=(1-s(T))/(1+s(T))` instead satisfies `q(T₁⇒T₂)=q(T₁)q(T₂)`.
   It retains the sign distinguishing responses above and below one.

3. **Kernel inequality (elementary, 3 cases; FORMALIZED in Lean).** For all `x, y ≥ 0`,
   > `max(1,x)·max(1,y) ≥ max(1+xy, x+y)/2.`
   The cases `x,y>=1`, `x,y<=1`, and `x>=1>=y` (or its symmetric case) each reduce to a trivial inequality.
   Machine-checked, `sorry`-free, kernel-pure `[propext, Classical.choice, Quot.sound]`:
   `HgeDiscKernel.kernel_ineq` and the peel form `HgeDiscKernel.peel_step`
   (`04-computation/lean/TournamentH7/TournamentH7/HgeDiscKernel.lean`, in the root manifest) —
   the exact analogue of THM-1860's `SumLeProd.lean` arithmetic kernel.

## `s` is a regularity coordinate: `s(T)=n ⟺ T regular`

Because a regular tournament has `K𝟙 = 0` (row sums `= out−in = 0`), `x=(I+K)⁻¹𝟙=𝟙` and `s=n`; and
`s=n` forces `x∥𝟙` (equality in Cauchy–Schwarz `s=𝟙ᵀx≤√n‖x‖=√n√s`) hence `K𝟙=0`, i.e. regular.
Verified exhaustively n=3,5,7 (0 mismatches). So `s ∈ [0,n]` with the top value pinned to regularity;
regular tournaments (Paley, rotational) uniquely maximize `s`, and there the base reads `H ≥ n·disc`
(e.g. Paley-7: `H=189 ≥ 7·8=56`). Reducibility does not force `s<1`:
a singleton factor gives `s=1`, while joining `C₃` to the response-`2/3`
strong witness gives `s=11/9>1`. The exact sign rule is
`s(T₁⇒T₂)-1=-(s₁-1)(s₂-1)/(1+s₁s₂)`.

## A globally typed multiplicative pair

The same Schur formulas, with no additional hypothesis, give

```text
D_+(T)=disc(T)(1+s(T))/2 > 0,
D_-(T)=disc(T)(1-s(T))/2,
D_±(T₁⇒T₂)=D_±(T₁)D_±(T₂).
```

This is a reformulation of the displayed composition law, not an additional
conjecture or a novelty claim. It recovers `disc=D_++D_-` and
`disc*s=D_+-D_-`; the ratio `q=D_-/D_+` gives the global signed coordinate.
The second coordinate may be negative and must not be discarded.

In particular `P(T)=max(1,s(T))*disc(T)=D_+(T)+|D_-(T)|`. Writing the
coordinates of two factors as `(a,b)` and `(c,d)`, with `a,c>0`, gives

```text
P(T₁⇒T₂)=ac+|bd| <= (a+|b|)(c+|d|)=P(T₁)P(T₂).
```

Thus the invariant used by the reduction is submultiplicative even though
`disc` itself is neither universally submultiplicative nor universally
super-multiplicative. For example `C₃⇒C₃` has correction factor `5`, while
the order-seven witness followed by a singleton has correction `5/6`.

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
`n ≥ 7`) it is just `H ≥ disc`. Among the audited strong orders 3 through 7,
equality occurs only at `C₃`. The base is the sibling of THM-1860's open
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
