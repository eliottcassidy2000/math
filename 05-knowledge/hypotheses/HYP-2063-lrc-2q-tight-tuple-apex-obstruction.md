---
id: HYP-2063
status: PARTIAL PROOF — reduction + apex obstruction + tight-tuple repair PROVED; residual characterised (iff); full 2q corrector FALSE
source: opus-2026-06-02-S559
related:
  - HYP-2056
  - HYP-2059
  - HYP-2061
  - THM-369
---

# HYP-2063: the k+1=2q polynomial-method obstruction is the apex zero-divisor q

The literature's accelerator for LRC(n≥8) is the polynomial method (Sungkawichai-
Trakulthongchai, arXiv:2604.23906, Prop 4.1): for `k+1` an odd prime, the tight
tuple `(1,…,k)` is a universal corrector over the field `ℤ_{k+1}`. It FAILS at
`k+1=2q` (n=14=2·7), the published wall. This pinpoints why.

## Proved

- **Reduction Theorem.** `ℤ_{2q}≅ℤ_2×ℤ_q`; every unit mod 2q is odd ⇒ the mod-2
  component of `s·v+r·(1,…,2q-1)` is FORCED (`≡v_i+i`), mod-q is free over `ℤ_q^×`.
  So the corrector into `{1,…,2q-2}^{2q-1}` holds **iff** `∃ s',r'∈ℤ_q^×:
  s'w_i+r'c_i ≠ f_i` ∀i, where `w_i=v_i mod q`, `c_i=i mod q`, `f_i=0` if `v_i+i`
  even else `q-1`. (0 mismatches, q=3,5,7,11.)
- **Apex Obstruction.** The unique speed `≡0 mod q` in `{1,…,2q-1}` is the apex
  `q=(k+1)/2=n/2`. Its constraint is parameter-free in `r'`; for the tight tuple
  (`v_q=q≡0`, `f_q=0`) it is `0≠0`, unsatisfiable ⇒ the full corrector is FALSE for
  the tight tuple (all q). The apex IS the zero-divisor where `ℤ_{2q}` stops being
  a field — the localized cause of Prop 4.1's failure, = the repo's apex/co-observer.
- **Tight-tuple repair.** The apex is `1/2`-safe at the base `t=s/(2q)`; excluding
  it, `s=r=1` corrects all non-apex runners (land at `2i mod 2q ∈{2,…,2q-2}`).
  Verified `excl-apex=True`, all q. ⇒ the hardest tuple (the paper's `c=k+1` lift
  bottleneck) is handled analytically once the apex is set aside.

## Characterised (iff, verified)

For parity-matched `v` (all `f_i=0`), apex excluded: un-correctable **iff** the
ratios `{-c_i/w_i}` cover `ℤ_q^×`. The tight tuple collapses ratios to `{-1}`
(easiest); the residual hard inputs are the ratio-spread tuples.

## Not claimed / open

The full 2q Prop-4.1 is FALSE (ratio-spread residual). Open: (a) clear the
residual using the `r/p` freedom of `t=s/(2q)+r/p`; (b) hand the apex/zero-divisor
to the pinch/shield route (HYP-2061, THM-396) — apex = the (q,q) shield; (c) the
story transfers to all `k+1=2·prime` (n=10,14,22,…).

**See:** `07-reflections/lrc-2q-tight-tuple-analogue-the-apex-is-the-field-failure-s559.md`,
`04-computation/lrc_2q_tight_tuple_analogue_s559.py` (+.out),
`05-knowledge/results/lrc_2q_{apex_repair,residual_characterization}_s559.out`;
HYP-2056 (even-fold), HYP-2059/2061 (pinch/shield), THM-369; arXiv:2604.23906.
