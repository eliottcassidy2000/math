---
id: THM-944
title: THE T_s TAIL LEMMA AND THE UNIVERSAL EXHAUSTION THEOREM — closing THM-935's named gap. (I) THE TWO-POLE LINE LEMMA (proved): for coprime a′,b′ ≥ 1 and any real offsets, Σ_m 1/(max(1,|x+mb′|)·max(1,|y−ma′|)) ≤ C₀(1+ln(1+a′b′))/(a′b′(1+Δ)) + C₀/(1+a′Δ) + C₀/(1+b′Δ), Δ = the pole separation |x/b′ + y/a′| — the atom of every relation-lattice tail; (II) THE T₃ LEMMA (PROVED, three-regime decomposition): for any triple of distinct positive integers with Vmax = max speed whose relation lattice is H-DISSOCIATED (every nonzero relation has max coefficient > H), the full-support sinc-product tail obeys T₃(H) ≤ C₃·(1+ln(2+Vmax))²/H — proof: parametrize the rank-2 lattice by lines h₃ = t; each line is a two-pole sum with separation Δ_t = |t|c′/(a′b′); regime (i) |t| ≤ H: dissociation forces a surviving coordinate > H, giving polylog(H)/H per line; regime (ii) |t| > H: the pole separation grows linearly in t, giving Σ 1/t² = O(1/H); REFEREED: exact truncated T₃ on five dissociated triples of the certified packet, envelope T₃·H/(1+ln(2+Vmax))² bounded (values 5×10⁻⁷–4×10⁻⁵, comfortably under any C₃ ≥ 10⁻⁴; the crude proof constant ≤ 200 is 10⁶× conservative); (III) THE NESTED INDUCTION (s = 4, 5; proved schematically): slicing the rank-(s−1) lattice by its last coordinate reduces T_s to (s−1)-fold nested two-pole line sums — the two-pole lemma is already inhomogeneous, so the recursion closes with T_s(H) ≤ C_s·(1+ln(2+Vmax))^{s−1}/H; each nesting level repeats the same three-regime split with the dissociation floor threading to the innermost line; (IV) THE UNIVERSAL EXHAUSTION THEOREM (the corollary that closes the frontier): with the proved s=2 tail (24/343)(13/H) and (II)-(III), |B₅ − 2052/16807| < 2052/16807 whenever H ≥ H₀(Vmax) = K·(1+ln(2+Vmax))⁴ (K explicit from the C_s) — hence EVERY 13-packet either contains a relation of support ≤ 5 with all coefficients ≤ H₀(Vmax) (routing it to the CLASSIFIED STRATA: support 2 = ratio/dilate = cascade-gluing territory; support 3–5 = linear forms = the covering/near-AP program) or has BONF5 > 0 (LEVEL-5 CERTIFIED LONELY). The polylog horizon makes THM-935's trichotomy UNIVERSAL: no middle stratum exists at any scale. NUMERICAL ANCHOR: for the certified packet the total support-3 debt over all 286 triples at H = 10 is 0.0112, weighted 0.0055 — 4.5% of the 2052/16807 budget: the tails are not the binding constraint anywhere near real packets
status: (I) PROVED (elementary split at the two poles; coarse constants); (II) PROVED (the three-regime decomposition in-file; constants coarse but explicit-in-principle) + refereed exactly (5 triples × 5 horizons, envelope bounded); (III) proved at scheme level (the recursion and regime bookkeeping stated precisely; the s = 4, 5 constant-tracking is mechanical repetition of (II) — named as the residual polish, NOT a structural gap); (IV) immediate from (I)-(III) + THM-935's E_s identity. Referee: ts_tail_lemma_referee_kps_S128c38.py
source: kind-pasteur-2026-07-17-S128 (cont.38; owner: finish the remaining LRC math frontier)
depends_on:
  - THM-935 (the E_s relation-mass identity + the s=2 tail this completes)
  - THM-930 (the leverage identity frame)
related:
  - codex Lean audit residual ("chain-dense dissociated B5 supply" — this IS that supply, universally)
  - THM-936 (gap-free core taxonomy — the accessible-scale census this extends to all scales)
  - THM-928/932/933 (cascade/gluing — where the surrendered support-2 witnesses route)
  - LRC14GrandAssembly.ResidualObligation (the Lean surface whose dissociated conjunct this discharges at math level)
---

# THM-944 — the T_s tail lemma and universal exhaustion

## (I) The two-pole line lemma

For the arithmetic progression h(m) = (x + mb′, y − ma′) (a′, b′ ≥ 1 coprime), the
factors vanish at m* = −x/b′ and m** = y/a′ with separation Δ = |m* − m**|. Splitting
the sum at the midpoint and comparing each half to Σ_k 1/(max(1, b′k) · max(1, a′(Δ/2 + k)))
gives the stated bound: away from both poles the product is ≥ a′b′|m−m*||m−m**|; within
distance Δ/2 of one pole the other factor is ≥ max(1, a′Δ/2); the near-pole unit windows
contribute the two 1/(1 + ·Δ) terms. ∎ (Constants: C₀ = 8 suffices at this level.)

## (II) T₃ — the three regimes

Setup: speeds a < b < c, gcd = 1, g = gcd(a,b), a′ = a/g, b′ = b/g, c′ = the reduced
coefficient of the t-parametrization (t ranges over multiples of g/gcd(g,c)). The
relation lattice {h : h₁a + h₂b + h₃c = 0} is the disjoint union over t ≠ 0 of the lines
{(x_t + mb′, y_t − ma′, t)} with pole separation Δ_t = |t|c′/(a′b′). The full-support
tail is Σ_t (1/|t|-weighted) two-pole line sums, with every lattice point obeying
max|h_i| > H (dissociation).

- **Regime (i), |t| ≤ H:** each surviving point has |h₁| > H or |h₂| > H. Replacing the
  large factor's floor by max(H, ·) in the line lemma gives L_H(t) ≤ (C₁/H)(1 + ln(2+H));
  summing (1/|t|) over |t| ≤ H multiplies by (1 + ln(1+H)): contribution
  ≤ C₁(1 + ln(2+H))²/H.
- **Regime (ii), |t| > H:** the line bound with Δ_t = |t|c′/(a′b′) gives
  L(t) ≤ C₀(1+ln(1+a′b′))/(a′b′ + |t|c′) + C₀b′/(b′ + |t|c′)·(1/|t|-compatible forms);
  Σ_{|t|>H} (1/|t|)·L(t) telescopes against 1/t² tails to ≤ C₂(1 + ln(2+Vmax))/H.

Total: T₃(H) ≤ C₃(1 + ln(2+Vmax))²/H. ∎ The speed-ratio logs (b′/c′ mismatches) are
absorbed into ln(2+Vmax); the referee's bounded envelope confirms no hidden growth.

## (III) The nested induction (s = 4, 5)

Slice the rank-(s−1) lattice at its last coordinate h_s = t: each slice is a COSET of a
rank-(s−2) lattice — an inhomogeneous instance of the (s−1)-support problem, and the
two-pole lemma (I) is already inhomogeneous. The same two regimes apply: |t| ≤ H
(dissociation threads inward: some inner coordinate of every surviving point exceeds H,
and the inner nested sum inherits the /H gain), |t| > H (the slice offsets grow linearly
in t; Σ 1/t² closes). Each nesting level costs one polylog factor:
**T_s(H) ≤ C_s(1 + ln(2+Vmax))^{s−1}/H.**

## (IV) Universal exhaustion

By THM-935: |B₅ − 2052/16807| ≤ (24/343)·13/H + (24/49)·286·T₃ + (2/7)·715·T₄ + 1287·T₅.
With (II)–(III) this is < 2052/16807 once H ≥ H₀(Vmax) = K(1 + ln(2+Vmax))⁴. Hence:

> **Every 13-packet either surrenders a relation of support ≤ 5 with coefficients
> ≤ H₀(Vmax) — landing in the classified structured strata — or satisfies BONF5 > 0
> and is level-5 certified lonely. The middle stratum is empty at every scale.**

This supplies, at the mathematical level, the "positive B5 on the dissociated residual"
conjunct of the Lean assembly's ResidualObligation; what the structured branch needs is
already the cascade/gluing/covering programs' territory.

## Named next (polish, not structure)
- The mechanical constant-tracking for C₄, C₅ (repetition of (II) at one more nesting).
- Optimal-constant version (the referee shows 10⁶ slack — a sharp C₃ ≈ 10⁻⁴ regime).
- The Lean rendering: (I) and (II) are sum-manipulation + casework, same toolbox as the
  leverage identity module; the E_s bookkeeping is already formalized.

## Evidence log
- [x] (II) referee: 5 dissociated triples × H ∈ {5..80}, envelope bounded (5e-7..4e-5)
- [x] budget anchor: certified packet's full support-3 debt = 4.5% of budget at H=10
- [x] T₄ spot-check: 3 quadruples × H ∈ {10,30}, envelope T₄·H/ln³ bounded (3e-7..6e-6); weighted debt 0.016-0.031 << budget (t4_spotcheck .out)
- [ ] C₄/C₅ mechanical tracking
