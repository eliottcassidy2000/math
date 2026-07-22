# The sharp DvdK-free criterion: a unique minimal balanced channel (unique primitive cycle) — most GMC(2) supports need no DvdK

**death-star-2026-07-21-S101** (HYP-8878). Owner: an even more creative DvdK bypass — "think one-dimensional
coprime intervals." The coprime-interval / numerical-semigroup lens gives a **support-only, coefficient-independent
criterion** that is far larger than the S100 edge case: GMC(2) is DvdK-free whenever the lowest face has a
**unique minimal balanced channel**. This confines the genuine DvdK input to a thin arithmetic-coincidence stratum.

## The criterion
`CT(f_F^m) = Σ_{balanced channels r of mass m} multinomial(m,r) · c^r`, all weights positive. Let `m*` be the
**minimal mass** of a nonzero balanced channel (`Σ k_i q_i = 0`, `Σ k_i = m*`, `k ≥ 0`) — the shortest positive
cycle in the one-dimensional charge lattice.

> **Criterion.** If there is a **unique** balanced channel at mass `m*`, then
> `CT(f_F^{m*}) = multinomial(m*, r*) · ∏ c_i^{r*_i} ≠ 0` — a single nonzero term.
> Hence GMC(2) is **DvdK-free** for that support, for *every* choice of (complex) coefficients.

It is coefficient-independent: uniqueness of the minimal cycle is a property of the **charges/support alone**, so
no cancellation is even possible — the certificate `Q = CT(f_F^{m*})` is one nonzero monomial in the `c_i`.

## How much this covers (verified)
Scan of all 116 straddling charge supports of size 3–4 in `[-4,4]`: **98 (84%) are DvdK-free** by the unique-cycle
criterion; only **18 are hard** (≥2 coincident minimal cycles). Examples that are DvdK-free but were "hard" under
the coarse S100 "≥3 monomials" reading: `[-1,1,2]` (m*=2, cycle (1,1,0)), `[-1,2,3]`, `[-1,3,5]`, `[-1,1,3,5]`,
the non-coprime-generated Sylvester set `[-6,10,15]` (m*=7, cycle (5,0,2)), `[-2,3,7]`, and even `[-1,0,1]`
(m*=1, the neutral monomial). **S100 was far too pessimistic** — most ≥3-monomial faces have a unique minimal
cycle and are DvdK-free.

## The thin hard stratum
DvdK is genuinely needed only when there are **≥2 coincident minimal cycles** — two distinct shortest balanced
channels of the same mass, whose contributions can cancel. The paradigm is two straddling **antipodal pairs**:
`{-2,-1,1,2}` has cycles `(±1)` and `(±2)` both at mass 2, so `CT(f^2) = 2c_{-1}c_1 + 2c_{-2}c_2`, cancellable.
More generally any `{-a,-b,b,a}`-type symmetric coincidence, or two sub-configurations of equal balance ratio.
This is a **codimension-≥1 arithmetic-coincidence stratum**; it thins out as the support gets more generic or
wider. (Even there DvdK may only be needed if the coincident cycles actually cancel for the specific null point.)

## Connection to the tournament zeta and the scale-clock picture
The minimal balanced channel is exactly the **fundamental primitive cycle** of the walk/tournament zeta
(THM-1926, `ζ_T = exp(Σ N_k u^k/k)` starts at `u^ℓ` where `ℓ` is the shortest cycle length; the leading
coefficient counts the shortest primitive cycles). So:
- **unique shortest primitive cycle ⟹ single nonzero leading zeta coefficient ⟹ DvdK-free** — the constant-term
  nonvanishing is the *nondegeneracy of the fundamental prime* of the charge-lattice zeta.
- **≥2 coincident shortest cycles = a degenerate zeta** (its `u^{m*}` coefficient is a multi-term sum) = the only
  DvdK-hard case = the S89–S91 charge-resonance / central-trinomial object (S90) = the "resonant multi-clock" in
  the S99 scale-then-clock picture (an edge is a single 2-clock; coincident cycles are coincident clocks).

## Honest scope and payoff
- **Not a full DvdK bypass**: the coincident-cycle stratum is genuine DvdK (or Monsky, ≈months, S95). No claim
  that GMC(2) is unconditionally DvdK-free.
- **A large, sharp, verified confinement**: DvdK-free for every support with a unique minimal cycle (the large
  majority — 84% of the small scan, more as supports widen), coefficient-independently; DvdK confined to a thin
  coincident-cycle stratum that is the already-studied resonance/central-trinomial/zeta-degeneracy object.
- **Formalization payoff**: "unique minimal balanced channel" is a **decidable, elementary** condition on the
  finite support (minimal cycle search + uniqueness). A Lean GMC(2) can discharge it with a single-term
  `CT = multinomial · monomial ≠ 0`, and cite DvdK *only* on the coincident-cycle stratum. This peels the one
  non-Mathlib input off a codimension-≥1 set rather than the whole problem, and much further than the S100 edge
  reduction.

Cross-links: S100/HYP-8877 (the coarser edge/positive confinement this sharpens), S99/HYP-8876 (scale-then-clock;
edge = 2-clock), THM-1926 (tournament zeta, primitive cycles), S90 (central trinomial = the resonant stratum),
S95 (DvdK roadmap), THM-2022 (GMC2 proof), memory `nc2-gmc2-lean-formalization-state`. Script
`04-computation/dvdk_unique_cycle_criterion_deathstar_S101.py` (+ `.out`). HYP-8878.
