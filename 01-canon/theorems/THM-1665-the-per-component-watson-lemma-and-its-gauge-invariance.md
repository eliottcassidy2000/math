---
id: THM-1665
title: "THE PER-COMPONENT WATSON LEMMA, PROVED IN ONE LINE, AND WHY GAUGE INVARIANCE MAKES IT WELL-POSED. The Laplace-layer nullcone (HYP-8350: L(p^m)=0 ∀m ⟹ p=0, L(v^k)=k!) was reduced by THM-1610(E) to Watson–Nevanlinna, with the Gevrey-1 bound the one missing input (HYP-8445, flagged 'not verified'). IT IS ELEMENTARY: |p(v)| ≤ C₀(1+v)^D with C₀ = sup|p|/(1+v)^D finite (→|a_D| at ∞), so |L(p^m)| ≤ C₀^m ∫(1+v)^{Dm}e^{−v}dv = eC₀^m∫₁^∞ w^{Dm}e^{−w}dw ≤ e·C₀^m·(Dm)! — exactly Gevrey-1 in τ = t^{1/D}. The SAME estimate bounds the analytic remainder of the resolvent Ψ(t)=∫[1/(1−tp)−1]e^{−v}dv (its tail is ∫(tp)^N/(1−tp)·e^{−v}dv, and THM-1610(E)'s rotated contour keeps |1−tp| off 0), so BOTH Watson–Nevanlinna hypotheses now hold in the same τ: the sector opening π(1+D) > π AND Gevrey-1. THE GAUGE: the nullcone carries a scaling gauge p ↦ cp under which the Watson DATA (degree D, sector π(1+D), Gevrey type 1) is INVARIANT — a function of D alone — while C₀ is COVARIANT (↦ |c|C₀). So the lemma descends to the gauge quotient (degree-D polynomials mod scaling) and 'per-component = per-degree' is well-defined — the same abstract move as the repo's cut/switching gauge. CONSEQUENCE: HYP-8350's reduction is now COMPLETE — L(p^m)=0 ∀m ⟹ Ψ ≡ 0 (Phragmén–Lindelöf on the sub-exponential decay e^{−c/|τ|}) — leaving the single residual Ψ≡0 ⟹ p≡0, which is DvdK's own Liouville step and is IMMEDIATE on the sign-definite locus (sign-definite p has Re(e^{iα}L(p))>0, so it fails the nullcone already at m=1)."
status: >
  The Gevrey-1 coefficient bound |L(p^m)| ≤ e C₀^m (Dm)! is PROVED (one-line estimate) and
  VERIFIED for several p, m = 2..8. The gauge-covariance of C₀ (scales by |c|) and
  gauge-invariance of the sector/Gevrey data are PROVED and verified. The remainder bound
  for the analytic resolvent follows from the same estimate GIVEN THM-1610(E)'s contour
  keeps |1−tp| bounded below — stated, with that contour input cited not re-proved. The
  sign-definite sub-case is PROVED. HYP-8445 is CLOSED; HYP-8350 is REDUCED to the single
  residual Ψ≡0 ⟹ p≡0 (DvdK's Liouville step), NOT closed.
source: mac-mini-2026-07-20-S145 (owner: "think gauge invariance abstractly, prove the
  per-component watson lemma via the standing route")
depends_on:
  - THM-1610  # (E) the sector opening pi(1+D) via the rotated 'standing' contour
  - THM-1600  # the Laplace layer at degree 1 (derangements); HYP-8350
related:
  - THM-1645  # GMC(2)'s radial layer, of which this is the tool
  - THM-474   # the cut/switching gauge -- the same abstract move
closes: HYP-8445
---

# THM-1665 — the per-component Watson lemma, and its gauge invariance

## The gauge (abstract)

The Laplace-layer nullcone problem — `L(p^m) = 0` for all `m ⟹ p = 0`, `L(f) = ∫₀^∞ f e^{−v}dv`,
`L(v^k) = k!` — carries a **scaling gauge** `p ↦ cp` (`c ≠ 0`): the condition `L((cp)^m) = c^m
L(p^m) = 0` is invariant. Under it:

- the **Watson data** — degree `D`, sector opening `π(1+D)`, Gevrey type `1` in `τ = t^{1/D}` —
  depends on `D` **alone**, hence is **gauge-invariant**;
- the **size constant** `C₀ = sup_{v≥0} |p(v)|/(1+v)^D` is **gauge-covariant**, `C₀ ↦ |c|·C₀`
  (verified: scales by exactly `|c| = 3`).

> So the per-component Watson lemma is a statement on the **gauge quotient** — degree-`D`
> polynomials modulo scaling — and "per-component = per-degree" is well-defined. This is the
> same abstract move as the repo's cut/switching gauge (THM-474, THM-1420): name the group,
> show the invariant data descends to the quotient, prove one representative per orbit.

## The lemma (the "standing route")

`THM-1610(E)` rotates the **standing contour** `v = ρe^{iφ}`, `|φ| < π/2`, to continue
`Ψ(t) = Σ_{m≥1} L(p^m) t^m` to a sector of opening `π(1+D)` in `τ = t^{1/D}` — past Watson's
threshold `π`. The one missing Watson–Nevanlinna hypothesis was the Gevrey-1 bound (HYP-8445).
It is elementary:

> **Lemma.** `C₀ := sup_{v≥0}|p(v)|/(1+v)^D` is finite (the ratio `→ |a_D|` as `v → ∞`), and
> `|L(p^m)| ≤ C₀^m ∫₀^∞(1+v)^{Dm}e^{−v}dv = e·C₀^m∫₁^∞ w^{Dm}e^{−w}dw ≤ **e·C₀^m·(Dm)!**`.

That is exactly **Gevrey-1 in `τ = t^{1/D}`** (the coefficient of `τ^{Dm}` is `L(p^m) ≤
e(C₀^{1/D})^{Dm}(Dm)!`). Verified `m = 2..8` for `v−1`, `v²+1`, `v³−3v+2`, with `|L(p^m)|/(Dm)!`
staying `O(1)`.

The **remainder** of the analytic resolvent `Ψ(t) = ∫₀^∞[1/(1−tp(v)) − 1]e^{−v}dv` after `N`
terms is `∫₀^∞ (tp)^N/(1−tp)·e^{−v}dv`, and the same `(1+v)^D` estimate — together with
THM-1610(E)'s rotated contour keeping `|1−tp|` bounded below — gives
`|Ψ − Σ_{m<N}| ≤ C·|τ|^{DN}·C₀^N·(DN)!`, the Gevrey-1 **remainder** bound. So **both**
Watson–Nevanlinna hypotheses hold in the same `τ`.

## Consequence: the reduction is complete

With the zero series, optimal truncation gives `|Ψ(τ)| ≤ inf_N e C₀^N N! |τ|^N ~ e^{−c/|τ|}`
(verified numerically), and sub-exponential decay in a sector of opening `> π` forces, by
Phragmén–Lindelöf / Watson uniqueness,

> **`L(p^m) = 0` for all `m ≥ 1` ⟹ `Ψ ≡ 0`.**

## The residual, and where it is immediate

`Ψ ≡ 0` means the pushforward `μ = p_*(e^{−v}dv)` has `∫ w^m dμ = 0` for all `m ≥ 1` with
`∫dμ = 1`. Concluding `μ = δ₀` (hence `p ≡ 0`) is **DvdK's own Liouville/monodromy step** and
is the single remaining piece — it is not automatic because `μ` can be non-compact /
indeterminate.

> **It is immediate on the sign-definite locus.** If `p([0,∞))` lies in a half-plane
> `{Re(e^{iα}w) > 0}`, then `Re(e^{iα}L(p)) = ∫ Re(e^{iα}p) e^{−v}dv > 0`, so `L(p) ≠ 0` —
> a sign-definite `p` fails the nullcone already at `m = 1`. (Verified: `v+1`, `v²+1` give
> `L(p) = 2, 3`; a sign-changing `v²−3v+1` has `L(p) = 0`.)

So the nullcone can contain only **sign-changing** `p`, and for those the residual is the
determinacy step.

## Honest scope

- **HYP-8445 (the Gevrey-1 bound) is CLOSED** — that is the concrete deliverable, and it is a
  one-line estimate.
- **HYP-8350 is NOT closed.** It is now reduced, with both Watson hypotheses met, to the single
  statement `Ψ ≡ 0 ⟹ p ≡ 0` — the analytic resolvent vanishing forces the polynomial to vanish.
  That is DvdK's Liouville step (a real theorem, not reproduced here) and is trivial only on the
  sign-definite locus.
- The remainder bound uses THM-1610(E)'s contour keeping `|1−tp|` off `0`; I cite that, I do not
  re-derive the contour placement (it is the DvdK residue/branch avoidance).
- The gauge here (scaling) is **not** the cut/switching gauge of the tiling side; the analogy is
  structural (invariant data descends to a quotient), not an identity of groups.
- Nothing here bears on GMC(2)'s *cross-shell* coupling (HYP-8470) — that is the other radial
  piece and is untouched. GMC(2) still needs both.

*Artifacts:* `04-computation/per_component_watson_macmini_S145.py` (+out).
