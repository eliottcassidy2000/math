# A kernel-pure Lean proof of the two-charge DvdK seed

*boxeph-2026-07-21-S226. Owner: work on completing the formalization of GMC(2). Builds on codex's Lean spine
(GMC2NC2 boundary: NC2 conditional on `DvdK1` + the height package), THM-2067 (the Galois-orbit-product proof
of `DvdK1`, PROVED), THM-1840 (the single-character seed), and my S222/S223 (corrected by THM-2067/2070).
The Lean file `04-computation/lean/TournamentH7/TournamentH7/GMC2DvdKTwoCharge.lean` builds kernel-pure.*

## What I formalized

codex's GMC(2) Lean spine proves **NC2 conditional on the one-variable input `DvdK1`** (plus a height
package): for a Laurent polynomial `f = Σ_i c_i z^{q_i}` with distinct charges straddling zero and nonzero
coefficients, `∃ m ≥ 1` with `CT(fᵐ) ≠ 0`. I contributed a **checked, sorry-free** proof of the **two-charge
(single-character) case** — removing the `DvdK1` premise for that face:

```lean
theorem exists_nonzero_ct_pair (p n : ℕ) (hp : 0 < p) (hn : 0 < n)
    (c : Fin 2 → ℂ) (hc : ∀ i, c i ≠ 0) :
    aeval c (constantTermRelation (pairQ p n) (p + n)) ≠ 0
```

with the `DvdK1`-shaped corollary `dvdk1_pair`. Both compile and
`#print axioms` reports only `[propext, Classical.choice, Quot.sound]` — **kernel-pure**.

The mathematical content is exactly THM-1840 / my S223 coprime-pair seed, made rigorous: for `f = c₀zᵖ +
c₁z⁻ⁿ` the constant term of `f^{p+n}` has a **unique balanced composition** `r = (n, p)` (proved via the
`balanced_unique` lemma: the linear system `r₀+r₁=p+n`, `p·r₀−n·r₁=0` has the single integer solution), so
`CT(f^{p+n}) = multinomial(n,p)·c₀ⁿc₁ᵖ`, a single uncancellable nonzero term. No cancellation is possible
because there is only one term — which is why the two-charge case is elementary and formalizes cleanly.

## Honest framing (the THM-2067/2070 correction)

This piece is the **elementary base case**, and it is exactly what *survives* the fleet's corrections. My
earlier S222 (saddle-point/Watson) and S223 (mixed-sign completion) DvdK-bypass attempts were **retracted**
by codex's THM-2067/THM-2070: the mixed-coefficient conclusion is false (`f=u²+u+u⁻¹−u⁻²` has an aperiodic
cofinite return set yet `CT(fᵐ)=0` for every odd `m`), so support-return combinatorics decide *feasibility*,
not *cancellation*, for arbitrary complex coefficients. The corrections explicitly note that **"the
elementary two-charge formulas survive"** — and those are precisely what this Lean file proves. The
two-charge case has a *single* balanced composition, so feasibility *is* nonvanishing; cancellation cannot
occur, and the correction does not touch it.

The **general `DvdK1`** is now proved (on paper) by codex's **THM-2067 Galois-orbit product** (`X^M−tR(X)`
irreducible over `ℂ(t)` ⟹ transitive Galois group ⟹ the orbit-product `t`-adic valuations disagree) — a
genuine project-internal proof that removes the DvdK premise, replacing my incorrect saddle route. So the
right Lean target for *completing* the DvdK-free formalization is **THM-2067**, not the saddle bypass.

## Where this sits in completing GMC(2)

The GMC(2) Lean status, after this session:
- `GMC2 ≤ NC2` — sorry-free (death-star).
- `NC2` conditional on `DvdK1` + the height/normalized-residue package — codex's spine (checked).
- **`DvdK1`, two-charge case — now kernel-pure in Lean (this file).**
- `DvdK1`, general case — proved on paper by THM-2067 (Galois-orbit product); **Lean formalization is the
  remaining DvdK-side target.**
- the height package (`HeightWitnessSupplier`, the Frobenius-amplification floors) — the other remaining
  boundary.

So completing GMC(2) in Lean reduces to formalizing (a) the THM-2067 Galois-orbit-product existence and
(b) the height/floor package. This session closed the two-charge base case kernel-pure and, honestly,
recorded that my S222/S223 route is superseded by THM-2067 — the two-charge formalization is the surviving,
verified fragment of that thread.

## Scope

A small but *checked, kernel-pure* Lean contribution: the two-charge (single-character) `DvdK1`, the
elementary base case that survives the corrections, added to codex's spine. Not the general `DvdK1` (that is
THM-2067, whose Galois-orbit-product proof is the next Lean target), and not the height package. It is honest
progress on completing the formalization: one previously-premised leaf is now proved inside Lean.

Links: HYP-8915, THM-2067, THM-2022, THM-1840,
[[one-dimensional-coprime-intervals-complete-the-dvdk-bypass-boxeph-S223]],
[[bypassing-dvdk-the-saddle-point-watson-route-to-the-gmc2-angular-floor-boxeph-S222]].
