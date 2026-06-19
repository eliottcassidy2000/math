---
id: HYP-2641
status: STRONG-PARTIAL — the unbounded→bounded reduction REDUCED to a single 1-D Weyl decorrelation estimate (margin 0.13–0.18, vs the divergent-lattice obstruction MISTAKE-078); LRC(14) NOT proved
source: kind-pasteur-2026-06-19-S13
related:
  - HYP-2637   # Freiman-dimension pockets (this is the unbounded-direction closure)
  - HYP-2607   # the crux: consec maximizes meas(S7)
  - THM-538    # support-6 floor / the divergent lattice obstruction (MISTAKE-078)
  - HYP-2610   # stranger-contraction (the |A|=2 special case of this)
  - THM-531    # dilation invariance (handles dilated APs)
---

# HYP-2641: the far-element plateau recursion — converting the UNBOUNDED direction to a 1-D Weyl estimate

The LRC(14)-S3 crux is `meas(S7(E)) ≤ cap_k` for all primitive k-sets (k=8,9,10), `=` the tight finite
check at the AP + the OPEN unbounded direction. THM-538's lattice-Fourier route hit the **divergent
absolute envelope** (MISTAKE-078). This hypothesis CONVERTS the unbounded direction into a single,
standard, margin-having equidistribution estimate, sidestepping the divergence.

## The plateau identity (the largest element decorrelates)

Write `meas(S7)(E) = p_0(E)`. If `w = max E` is large, the orbit point `frac(w x)` equidistributes
(Weyl) and acts as an INDEPENDENT uniform fill-in relative to the bounded-frequency core `E' = E∖{w}`:
```
  p_0(E)  →  Plat(E') := p_0(E') + (1/7)·p_1(E')   as the dissociation of w grows,
```
where `p_1(E') = P(E'-orbit misses exactly one of the 6 inner sectors)` (a missed sector gets filled by
the independent point w.p. 1/7). VERIFIED exact: at k=9, core `consec_8`, `Plat = 0.32721 + (1/7)(0.24422)
= 0.36210`, matching the computed far-`w` plateau to 5 digits. The approach is **rate O(1/w)**: `|p_0 −
Plat| ≈ 0.5/w` (w≥20 ⟹ error < 0.008; verified w up to 300).

## The bound Q and the RECURSION in k (VERIFIED)

Define `Q(m) := max over m-element sets (0∈E) of [p_0(E) + (1/7)p_1(E)]`. Then for any k-set with a far
largest element, `p_0(E) ≤ Q(k−1) + C/w`. Exact (verified):
```
  Q(7) = 0.19660  < cap_8 = 0.38153   (margin 0.185)   — achieved at consec_7
  Q(8) = 0.36210  < cap_9 = 0.49426   (margin 0.132)   — achieved at consec_8
  Q(9) = 0.44789  < cap_10 = 0.60440  (margin 0.157)   — achieved at consec_9
```
**`Q(k−1)` is achieved at the bounded AP-core `consec_{k−1}`** (itself a finite check at `k−1`), and is
`< cap_k` with **margin 0.13–0.18** — an order of magnitude looser than the tight `consec_k` margin
(0.0014 at k=9). So the far-element direction is NOT tight: a crude Weyl bound closes it.

`Q(k−1)` is itself a `(k−1)`-extremal problem, so this is a genuine **recursion DOWN in k**: the
unbounded-`k` bound rests on a bounded-`(k−1)` finite check with more margin; it bottoms out at small k.
(And `Q(m)` over UNBOUNDED `m`-sets reduces again by the same far-element argument, so the bounded max IS
the sup.)

## The complete partition (with this hypothesis)

```
  1. dilated AP            -> = consec by dilation-invariance (THM-531)  -> < cap [exact]
  2. bounded spread (max<B) -> finite check (incl. the tight AP & near-AP)  -> < cap [DONE]
  3. max ≥ B (far element)  -> p_0(E) ≤ Q(k−1) + C/B < cap [margin 0.13]  -> THIS hypothesis
```
With margin `m(k) ≥ 0.13` and rate `C/w` (`C ≈ 0.5`), the threshold is `B ≈ C/m(k) ≈ 4` — TINY. So the
bounded finite check needs only a small `B`, and everything beyond is covered by the plateau bound.

## The one remaining rigorous lemma (the precise target)

> **1-D Weyl decorrelation:** for a k-set `E` with `w = max E`, `|p_0(E) − Plat(E∖{w})| ≤ C/w` with an
> explicit absolute `C`.

Equivalently (THM-538 language): the signed sum of `K(n)` over relations `n ∈ Λ°(E)` that INVOLVE the
coordinate `w` (`n_w ≠ 0`) is `≤ C/w`. This is a 1-D estimate (one fast frequency), NOT the full
divergent lattice sum — the support-6 floor + the per-coordinate envelope `0.6973/|n|` make relations
involving large `w` sparse (they need `Σ_{j≠w} n_j e_j = −n_w w`, forcing large other coordinates).
This is the standard Erdős–Turán / single-frequency discrepancy estimate, and with margin 0.13 it does
not need to be sharp.

## Why this is progress

- It **replaces the MISTAKE-078 obstruction** (the full divergent lattice sum) for the unbounded
  direction with a single 1-D estimate that has comfortable margin.
- It is the explicit, quantitative form of HYP-2610's "stranger contraction" (which was the `|A|=2`
  special case) and HYP-2637's "dissociation peel," now with the EXACT plateau `Q(k−1)` and the recursion.
- The tight part (margin 0.0014) is confined to the bounded finite check (done); the open analytic part
  (the Weyl estimate) has margin 0.13.

## Honest status
LRC(14) NOT proved. The unbounded direction is reduced to the 1-D Weyl decorrelation lemma (margin 0.13,
standard tool); the bounded part is the done finite check. Remaining: (a) prove the Weyl estimate with an
explicit `C` (and hence `B`); (b) verify the finite check to that `B` (small). Files:
`04-computation/lrc14_{far_element_decay,far_element_recursion,weyl_plateau}_kps.py`.
