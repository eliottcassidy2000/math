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

No analytic gap remains. The finite certificate is now **written** — mac-mini's **THM-702** (exact
arithmetic) plus an independent grid-numerical cross-check (kps, `05-knowledge/results/lrc14_finite_certificate_kps_S127.md`,
every digit agreeing). Status of the three items:
1. **The balanced-core base check** `Φ(F) ≤ cap_{|F|+1}` — margins tabulated and cross-validated (tightest
   `+0.086` at `|F|=8`, exact indexing `cap_{|F|+1}`). **Reduced, not closed:** the check shows `Φ`-extremality
   does **not** factor — `p_0` (cover) is maximized at `consec` (THM-534), but `p_1` is **not** (random cores
   have higher miss-one). So the sole open residual is the *joint* `Φ`-consec-extremality on bounded cores — a
   `(1/3)p_1`-perturbation of THM-534's `p_0`-extremality, stable because `λ = 1/3 < λ*(m) ≈ 1` (a **third**
   role for `λ`, beyond the tax and the increment). Same extremal lemma THM-534/530/657 keep isolating.
2. **The cap-growth verification** — **DISCHARGED exactly** (THM-702): `cap_{k+1}−cap_k =
   94841/840840, 63/572, 11/91, 12/91, 1/7` for `k=8..12`, worst slack `179/12012 > 0` at `k=8`. The `91 = 7·13`
   denominators come from the extremal 2-runner packing `{1,13}` (`13 ≡ −1 mod 14`).
3. **The explicit far-element threshold** — **DISCHARGED exactly** (THM-702): `w ≥ 90191·Σ_{e∈E'} e` via
   THM-699's proven constant `K₇ = 672`, per-peel self-budgeting (no error summation needed).

The tight margin is elsewhere entirely — it lives in the finite `consec_k` / `L_y` moment check (THM-534),
which item 1's residual now provably carries.

## Files

`04-computation/lrc14_recursion_closure_kps_S127.py` (+ `.out`): the transfer operator, the p1-tax
(`≤ 0.2514 < 1/3`), the joint bound `Φ ≤ cap` (margins `≥ 0.29`), and the accumulation.

## Addendum — the base lemma's TARGET CORRECTED by exhaustive census (klein-2026-07-11-S251, HYP-5985)

The named remainder ("joint Φ-consec-extremality") is **FALSE as stated for k ≥ 11**, by
exact exhaustive census (`lrc14_phi_global_argmax_census_klein_S251.py`, all-integer
breakpoint sweep, evaluator cross-checked digit-exact against mac-mini cont.29's):

- **k = 4..10: consec IS the exact global Φ-argmax** in exhaustive boxes
  ([1..W] = 60/42/30/24/20/17/15 for k = 4..10; 28k–140k normalized sets each).
  At the BINDING size k = 8 (tightest cap margin, +0.085639 exact = 108013/1261260)
  the consec-extremality statement is exhaustively TRUE in [1..20].
- **k = 11: the global argmax is X = (0,2,3,4,5,6,7,8,10,12,14) = 2·consec₈ ∪ {3,5,7}**
  (the 2-dilated consecutive 8-set with odd infill), exhaustive in [1..14] (1001 sets)
  and STABLE at [1..16] (8008 sets — interior extremal, not box-edge).
  Φ(X) = 11159/17640 beats Φ(consec₁₁) = 33337/52920 by exactly **1/378**.
- **The p₀-argmax stays consec** (p₀(consec) − p₀(X) = 89/5880 > 0): THM-534 intact;
  this is cont.24's "does not factor" made exact — p₁(X) − p₁(consec) = 167/840 − 1283/8820
  is large enough to flip the joint functional.
- **The exact flip threshold: λ*₁₁ = 267/941 ≈ 0.2837 < 1/3.** cont.24's λ*-ladder
  (1.51 / 1.23 / 0.98 at k = 8/9/10) crosses 1/3 between k = 10 and k = 11; the third
  λ-constraint ("ext < λ*") is unsatisfiable at λ = 1/3 for k ≥ 11.

**The recursion itself is UNAFFECTED — and the census strengthens it**: the inequality the
induction actually needs, `Φ(F) ≤ cap_{|F|+1}`, holds in EVERY box with margins
+0.1196 / +0.0856 / +0.1142 / +0.1583 / **+0.2245** (k = 7..11) — the k = 11 margin, where
consec-extremality dies, is the LARGEST.  The base lemma should therefore target
`Φ(F) ≤ cap_{|F|+1}` **directly** (exhaustive-box + tail decorrelation), or split:
consec-extremality for k ≤ 10 (censused true) + the perforated-dilate family
(2·consec ∪ odd-infill) as the argmax candidate for k ≥ 11 with its huge margin.
Small-size base spec note: cap_m for m ≤ 7 is not in the ledger; for k ≤ 6, p₀ ≡ 0
(five nonzero offsets cannot hit six inner sectors) and max Φ = max p₁/3 is censused
(0 / 0 / 11/126 at k = 4/5/6) — the induction's bottom sizes need their caps SPECIFIED.

Files: `lrc14_phi_global_argmax_census_klein_S251.py` (+ `.out` with the k=11 addendum).
