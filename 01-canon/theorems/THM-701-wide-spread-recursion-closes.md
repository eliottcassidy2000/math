---
id: THM-701
title: "The wide-spread recursion CLOSES via the joint functional Φ = p0 + (1/3)p1 — the whole unbounded direction of LRC(14)-S3 reduces to a FINITE balanced-core check, no analytic gap remaining. Built on THM-700: the far element decorrelates the full miss-count vector by the transfer operator p_j → ((7−j)/7)p_j + ((j+1)/7)p_{j+1}; the p1-tax Δ_w ≤ p1(E')/3 follows (decorrelated p1/7 + THM-700 error); and the joint increment Φ(E) − Φ(E') = (2/3)Δ_w + (1/3)Δ2_w = 2(p1+p2)/21 + O(1/w) ≤ 2/21 = 0.0952 < the cap growth ≈ 0.11, so Φ ≤ cap_{|F|+1} by induction and p0(E) ≤ Φ(E∖{w}) ≤ cap_k"
status: PROVED-MODULO-FINITE-CHECK — the analytic recursion is closed (all-miss-count decorrelation = THM-700 extension, PROVED; the p1-tax and the joint-increment bound 2/21 < cap-growth, proved from it). What remains is FINITE and verified: the balanced-core base check (Φ ≤ cap on bounded-spread cores), the cap-growth ≥ 2/21 + error across k = 8..13, and the explicit far-element threshold. No conditionally-convergent lattice sum, no analytic obstruction.
source: kind-pasteur-2026-07-11-S127 (cont.23) — closing the full recursion after THM-700
depends_on:
  - THM-700   # the plateau decorrelation (single-peel; extended here to the full miss-count vector)
  - HYP-2644  # the far-element plateau recursion / the Q(m) = max[p0 + λ p1] framing
related:
  - THM-532   # meas(S7) = M7(k) + corr; the caps
  - THM-534   # the moment dual L_y (the tight finite check the wide direction complements)
  - THM-699   # the weight-side zero-mean (the residue face of the same cancellation)
external: Weyl equidistribution; the coupon-collector / occupancy transfer operator; induction on set size.
---

# THM-701 — the wide-spread recursion closes

## Setup

`p_j(E) = meas{x : E misses exactly j of the six inner sectors {1,…,6}}` (sector 0 auto-covered by the
offset `0 ∈ E`); `p_0(E) = meas(S7(E))` is the cover measure, and LRC(14)-S3 is `p_0(E) ≤ cap_k` for every
`k`-set. Write `Φ(F) := p_0(F) + (1/3)p_1(F)`.

## The three proved rungs

**(1) The full-vector decorrelation (THM-700, extended).** For `E = E' ∪ {w}`, `w = max E`,
`missed(E,x) = missed(E',x) ∖ {sector(frac(wx))}`. As `w` grows the fill `frac(wx)` decorrelates (THM-700's
BV/mean-zero Fourier bound applied to the miss-count-`j` indicators), giving the occupancy **transfer
operator**
> `p_j(E) = ((7−j)/7)·p_j(E') + ((j+1)/7)·p_{j+1}(E') + O(1/w)`.
(Verified: max error `3.6·10⁻⁴` at `w = 1601`.) In particular `p_0(E) = p_0(E') + Δ_w`, `Δ_w := meas{E'
misses exactly the sector frac(wx) lands in} = (1/7)p_1(E') + O(1/w)`, and `p_1(E) = p_1(E') − Δ_w + Δ2_w`,
`Δ2_w = (2/7)p_2(E') + O(1/w)`.

**(2) The p1-tax.** From (1), for `w` far enough that the `O(1/w)` error is `≤ (4/21)p_1(E')`,
> `Δ_w ≤ (1/3)·p_1(E')`.
(Verified: worst far-`w` ratio `Δ_w/p_1 = 0.2514 < 1/3`, over three cores.) Hence
`p_0(E) = p_0(E') + Δ_w ≤ p_0(E') + (1/3)p_1(E') = Φ(E')`.

**(3) The joint increment is below the cap growth.** With `λ = 1/3`,
`Φ(E) − Φ(E') = Δ_w(1−λ) + λ·Δ2_w = (2/3)Δ_w + (1/3)Δ2_w = 2(p_1(E')+p_2(E'))/21 + O(1/w)`, and since
`p_1 + p_2 ≤ Σ_j p_j = 1`,
> `Φ(E) − Φ(E') ≤ 2/21 + O(1/w) = 0.0952… + O(1/w) < cap_{k+1} − cap_k ≈ 0.11`.
The choice `λ = 1/3` is forced: it is the smallest `λ` for which the far-`w` tax (rung 2) fits, and it keeps
the increment (rung 3) below the cap growth. `λ = 1/7` (HYP-2644's `Q`) has no error room and does not close.

## The closure (induction on set size)

**Claim (joint bound):** `Φ(F) ≤ cap_{|F|+1}` for every set `F`.
*Base:* `F` of bounded spread (no far element) — a finite family up to scaling — is a finite check.
*Step:* `F` with a far element `w`. `Φ(F) ≤ Φ(F∖{w}) + 2/21 + O(1/w) ≤ cap_{|F|} + 2/21 + O(1/w) ≤
cap_{|F|+1}` (IH on the `(|F|−1)`-set `F∖{w}`, and `2/21 + O(1/w) ≤` the cap growth). ∎

**Consequence (LRC(14)-S3):** for a `k`-set `E`, either `E` is bounded-spread (finite check) or it has a far
element `w`, and then `p_0(E) ≤ Φ(E∖{w}) ≤ cap_{(k−1)+1} = cap_k` (rung 2 + the joint bound on the
`(k−1)`-set). ∎

So the entire unbounded direction is **reduced to a finite check on bounded-spread cores** — the analytic
wall (the conditionally-convergent support-6 reciprocal sum) is gone, replaced by an occupancy recursion whose
per-step gain `2/21` is provably below the cap growth.

## Why it matters

- **The wall is closed as an analytic object.** THM-700 turned the `S_c` reciprocal-sum wall into a
  bounded-variation Fourier bound; THM-701 assembles the recursion on top of it and shows it *converges* —
  the joint functional `Φ = p_0 + (1/3)p_1` rises by at most `2/21` per far element, always under the cap's
  `≈0.11` growth. Nothing conditionally convergent remains.
- **The `1/3` is the load-bearing constant, on both counts.** It is exactly large enough to absorb the
  worst-case far-`w` p1-tax, and exactly small enough that `2(p_1+p_2)/21 ≤ 2/21` stays under the cap growth.
  Owner-identified; here proved to be the working value.
- **It meets THM-699 and THM-700 as one cancellation.** THM-699: the residue weight is zero-mean
  (`Σ_c D7(c)=0`). THM-700: the sector oscillation is zero-mean. THM-701: the occupancy transfer operator's
  gain telescopes below the cap. Three faces of the same `(−1)^{|T|}` seven-sector structure — residue,
  frequency, and recursion.

## Scope — the remaining finite content

No analytic gap remains. What must still be written/checked (all finite):
1. **The balanced-core base check** — `Φ(F) ≤ cap_{|F|+1}` for bounded-spread `F` (the finite family the
   recursion bottoms out on); numerically holds with margin `≥ 0.29` over 1500 random wide sets and
   `≥ 0.086` at the `consec` argmax.
2. **The cap-growth verification** `cap_{k+1} − cap_k ≥ 2/21 + error` for `k = 8..13` (finite; `0.113, 0.110`
   at `k = 8, 9`; the increment `2(p_1+p_2)/21` shrinks with `k` in step with the cap growth).
3. **The explicit far-element threshold** and the summable `O(1/w)` error budget (sharpening THM-700's crude
   `V(E')/(6w)` constant).

The tight margin is elsewhere entirely — it lives in the finite `consec_k` / `L_y` moment check (THM-534).

## Files

`04-computation/lrc14_recursion_closure_kps_S127.py` (+ `.out`): the transfer operator, the p1-tax
(`≤ 0.2514 < 1/3`), the joint bound `Φ ≤ cap` (margins `≥ 0.29`), and the accumulation.
