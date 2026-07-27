---
id: THM-2547
title: "External cut-character/root-control pairing and the delta_2 cancellation boundary"
status: >
  FINITE-EXACT EXTERNAL CONTROL + REFUTED PHYSICAL INTERPRETATION.
  The displayed unrelated coefficient arrays pair nontrivially in all
  432 placements (both leg conventions), exactly in Q(zeta_91), but this
  is not a common-ancestry or physical current.  A Boolean delta_2 hostile
  proves exact non-universality: full support of both factors does not
  prevent 108/432 pairing zeros.  No target/owner retention, typed bridge,
  live-row transplant, or LRC(14) branch exclusion follows.
source: opus-2026-07-27; hostile correction gain-slice-bridge-2026-07-27
depends_on:
  - THM-2508-affine-cut-bundle-covariance-and-carry-permutation
  - THM-2506-punctured-stalk-primitive-module-saturation-and-thirteen-primary-pushforward-no-go
  - THM-2471-owner-first-collision-weighted-root-service-and-temporal-atom-boundary
related:
  - THM-2478 (typed future-collision/old-component equality; not a cut bridge)
  - THM-2507 (Section 7 proposed the test)
  - HYP-9032 (live-row transplant trichotomy and owner-clock re-homing)
  - MISTAKE-281 (failed physical-current promotion)
script: 04-computation/lrc14_cut_bundle_ancestry_pairing_opus_20260727.py
output: 05-knowledge/results/lrc14_cut_bundle_ancestry_pairing_opus_20260727.out
script_sha256: c1d18785e4282006fa262b11d8165832a49d39f54e0eaeb8bef8c42d222ccd05
output_sha256: 5fb872fd3c834de26233bb41f48d221fff0ab064dab3e67a7cc0bbf726069d96
hash_basis: working-tree bytes (LF)
---

# THM-2547 -- external pairing and its cancellation boundary

## Exact external-control theorem

Let

```text
C = (6/25)(delta_3 + delta_6 + delta_10 + delta_12)
```

be the synthetic weighted root-correlation control reconstructed from
THM-2471's exact companion, and let `d` be THM-2506's canonical two-row
defect

```text
d(0,5)=1, d(0,3)=-1, d(1,5)=1, d(1,4)=-1.
```

Write `Chat` for the `C_13` transform of `C` and `Psi` for THM-2508's
cut-bundle transform of `d`. After making the external numerical label
alignment `alpha=alpha`, define

```text
P^-_{tau,a}(beta)
  = sum_{alpha in F_13^*} Chat(-alpha) Psi_{tau,a}(alpha,beta),
P^+_{tau,a}(beta)
  = sum_{alpha in F_13^*} Chat(+alpha) Psi_{tau,a}(alpha,beta).
```

Then `P^-` and `P^+` are nonzero in `Q(zeta_91)` for every

```text
(tau,a,beta) in F_13^* x F_7^* x F_7^*.
```

Thus each convention has exactly `432/432` nonzero placements. An exact
dual witness is

```text
100 P^-_{1,1}(1) = 312 zeta_91^39 - 312 zeta_91^52.
```

This is an existence theorem for two explicitly displayed finite arrays.
It is not a theorem about a physical LRC row.

## Exact non-universality theorem

There is a valid Boolean first-collision root control with one predecessor
packet at root `0` and one disjoint arrival packet at root `2`. Its
correlation is `C=delta_2`, so `C(0)=0`, its service is positive, and all
twelve nontrivial root Fourier colours are nonzero. The same canonical
defect has all `5,184` primitive cut coefficients nonzero. Nevertheless:

```text
#{(tau,a,beta): P^-_{tau,a}(beta)=0} = 108,
#{(tau,a,beta): P^+_{tau,a}(beta)=0} = 108.
```

For the dual convention the zero `tau` values are `{4,6,12}`; for the
untwisted convention they are `{1,3,7}`. For each listed `tau`, all 36
choices of `(a,beta)` vanish. In particular, at dual `tau=4`,

```text
R_{4,a,c}(2)=0                    for every a in F_7^*, c in F_7,

sum_{alpha!=0} zeta_13^(2 alpha) Psi_{4,a}(alpha,beta)
  = 13 sum_c R_{4,a,c}(2) xi_7^(-beta c)
  = 0.
```

The identity uses `sum_v R_{tau,a,c}(v)=0`. This isolates the failed
implication: nonzero individual factors do not prevent cancellation in the
sum over `alpha`.

## Type and retention audit

The original physical interpretation is refuted for three independent
reasons.

1. The root control and cut defect do not live on an exhibited common
   Boolean base. The first is a synthetic companion control; the second is a
   canonical packet from the already-empty high-septimal branch. No atomwise
   product or common ancestry morphism is constructed.
2. THM-2478 (21a)-(21b) identifies a printed future collision-colour label
   with one fixed old target-component label. It does not identify the
   temporal root character with THM-2506's static cut-row character. The
   `alpha` alignment above is external bookkeeping.
3. `K_{0,beta}=0` excludes the trivial Fourier character `alpha=0`. It is
   not the physical root-coordinate condition `C(0)=0`. Likewise, owner and
   deep margins are computed separately and never enter the pairing.

The covariance and kappa/cut torsor-collapse calculations remain exact for
the external algebraic pairing. They show that this array is equivariant and
that invariant averages die; they do not supply physical sidecars.

## Independent controls

The `Q(zeta_91)` computation reduces every placement modulo `Phi_91`, checks
the full `5,184` primitive-mode factorization, `864` beta-zero vanishings,
`504` invariant readouts, all 91 translations, and both torsor collapses. A
separate direct Radon/DFT implementation evaluates the 864 displayed
nonzero placements at each split prime

```text
547, 911, 1093, 2003, 2549
```

and finds no modular zero; it also reproduces all `144` beta-zero controls
per prime. Normal and `python -O` runs are required to be byte-identical.

## Exact remaining branch-transplant obligation

For a physical exclusion of the live low-septimal bank, the needed
quantifier is at least the following. For every live row `R` among the 165
rows, construct one positive-measure Boolean ancestry base `Omega_R` on
which collision, owner, word, deep, and cut data all derive from the same
atoms, together with a lawful typed intertwiner from static cut data to the
chosen collision/target component. Then choose a row-dependent but
`omega`-independent `(tau_R,a_R,beta_R)` and prove

```text
integral_{Omega_R} sum_{alpha!=0}
  Gamma_R(omega,-alpha)
  Psi_{d_R(omega),tau_R,a_R}(alpha,beta_R) dmu(omega) != 0.
```

The product must be formed before marginalization, Fourier transform, and
integration, while retaining the same owner, target, word, and deep/source
factors. Pointwise nonvanishing of the two factors is insufficient because
both the `alpha` sum and the later `omega` integral can cancel.

HYP-9032 refines the likely shape of the typed intertwiner: a live row cannot
recover the old vertical `C_7` column, and `13`-adic clock dilation transports
only the parity subgroup.  Its only native full `F_7` candidate is the owner
clock.  Even if that trichotomy supplies a live array satisfying the row-zero,
dichotomy, and magnitude laws, the common-base/integrated noncancellation
quantifier above remains the semantic step that this external control did not
prove.

## Correction lineage

The filename preserves the historical positive-promotion slug. The former
claim that this was the "first physical cut-character current" is refuted by
MISTAKE-281 and replaced by the two exact statements above. LRC(14) remains
open.
