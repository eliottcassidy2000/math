# AMM 12592 — Lane F1: the entry problem — mixed-sign certificates, the R/8 window, two new doublings, and the isolated core inequality

Session: boxeph multifront, 2026-08-04 (post Lane E1 / Theorem S-cone-fc).
All computations exact (int / Fraction); no floats in any decision; no
numpy; no sympy. Citations: F1–F4, S1–S4, Theorem S-cone-fc, a-priori-F4
corollary, ENTRY-fc, marginal-surface law to
`amm12592-S-invariant-cone-proof-boxeph.md`; Theorems A/B/C/D, C-N/C-A/C-F/
C-D, Hypothesis S to `amm12592-Elin-theorem-boxeph.md`; G1–G3 to
`amm12592-golden-transient-bound-boxeph.md`; T1–T9' to
`amm12592-allR-transient-theorem-boxeph.md`; L1–L4 to
`amm12592-localrule-barrier-and-gap-boxeph.md`.

Scripts (04-computation/):
`amm12592_entry_feedphase_scan_boxeph.py` — feed-phase-only exact ENTRY
scanner (verbatim certified-engine semantics; memory-trimmed top-band feed
coefficients; stops after the S-cone-fc certificate row + persistence
window; validated bit-exactly against the stored E1 fcscan/feedend records);
`amm12592_entry_feedphase_scan_v2_boxeph.py` — same semantics, list-based
with Pascal-rule cap updates (6.4x faster; validated field-by-field
against v1 ledgers at 1024/2048/4096, which are themselves E1-validated);
`amm12592_entry_mixedcone_referee_boxeph.py` — 600-trial exact certificates
for the new lemmas (independent clamp implementation);
`amm12592_entry_certificate_referee_boxeph.py` — INDEPENDENT re-check of
each new-scale certificate (fresh math.comb clamp; re-evolves
i_pf..i_fc from the saved snapshot; recomputes the clocks by the
mul//div staircase, diverse from the scanner's Pascal adds);
`amm12592_entry_window_constants_boxeph.py` — exact-rational wide-window
budget certification (dyadic 2^9..2^40).
Outputs (05-knowledge/results/): `amm12592_entry_scan_R*_D0*_boxeph.json`,
`amm12592_entry_certreferee_R*_D0*_boxeph.json`,
`amm12592_entry_mixedcone_referee_boxeph.{out,json}`,
`amm12592_entry_window_constants_boxeph.{out,json}`.

Status labels: **PROVED** (complete argument, machine-checkable algebra),
**VERIFIED-exact** (exact computation at stated scales), **HYPOTHESIS /
OPEN** (precise statement, open), **ROUTE-NEGATIVE** (a proof strategy
shown quantitatively unavailable — NOT an impossibility statement about
the mathematical claim itself).

---

## 0. Headline

1. **The entry hypothesis is weakened on four independent axes, all
   PROVED.** (i) *Sign*: a new mixed-sign one-row certificate (Theorem
   S-cone-fc±, via the signed-part decoupling Lemma EN-D) removes
   all-negativity (F1) as a hard requirement — a certifying row may carry
   positive junk, provided the positive part sits under its own mirrored
   cap surface; (ii) *window*: the certificate row may sit anywhere in
   `[i_pf, i_pf + floor(R/8)]` — Theta(R) rows instead of 64 — with the
   clocks still a-priori free (exact-rational certification, dyadic
   2^9..2^40, both eps; even 3R/16 passes); (iii) *debt edge*: an S4-layer
   variant tolerates `a_0 > D_0 - 1` at explicit G_L cost (with the
   honest accounting of when that cost is affordable); (iv) *surface
   excess* — the session's main theorem — **Theorem EN-H (self-healing,
   Section 8b)**: a feed-end state within `eta_0 = 3/64 = 4.7%` ABOVE
   the marginal surface still certifies: cap drift heals the excess by
   an explicit row `k* <= K - 34`, with the compounding growth controlled
   by an exact contraction (fixed point Delta-bar <= 0.024, certified in
   rationals at all dyadic 2^9..2^40). The razor-thin entry condition is
   now a fat one: measured excess <= 0.72% worst-ever vs 4.7% allowed.
2. **Two new doublings certified: S(1/32) and S(1/16) now hold at
   R = 65536 and R = 131072** — via the PROVED Theorem S-cone-fc plus a
   one-row certificate found by a feed-phase-only run (the ~0.5 R-row
   post-feed evolution is covered by the theorem, not simulated: the
   certificate replaces ~50000 rows of computation at 2^17). Hypothesis
   S(eps) is now settled affirmatively for ALL dyadic
   `128 <= R <= 131072`, eps in {1/32, 1/16}, and
   `C* <= 1 + gamma* + 1/32 < 1.6292382` is conditional only on the
   feed-end entry behavior at dyadic `R >= 2^18`.
3. **Sign collapse completes BEFORE feed-end (new observable,
   VERIFIED-exact at every scale).** The last row carrying any positive
   junk cell precedes i_pf at all 22 configs (by 3–8 rows at the top
   scales): F1's negativity costs ENTRY nothing at all — the feed phase
   itself hands over an all-negative state. The mixed-sign theorem is
   thereby a robustness margin, not a needed crutch.
4. **ROUTE-NEGATIVE (recorded with exact data): no fixed-margin envelope
   can prove the F3 surface at the low cells.** At feed-end the measured
   inflow ratios r_t sit at `1.000 +- 0.01` on the whole low-cell core
   (excess <= 1%, shrinking with R: 0.09% at 2^15, [filled below] at
   2^16/2^17). Any hypothesis of the form
   `spill_t <= (1-c) cap_t` with fixed c > 0.01 is FALSE at i_pf for all
   tested R >= 1024. The certificate fires only because cap DRIFT (degree
   growth) pulls the caps ahead of the frozen profile within <= 8 rows.
   Consequence: the entry proof must derive the marginal surface itself
   (with o(1) error), not a slack envelope; the "38% margin" of the E1
   tables lives exclusively in the F4 budget, never in F3.
5. **The single remaining inequality (EN-CORE) is isolated, stated with
   explicit constants, and certified at 22/22 configs (2^7..2^17 x both
   eps) with its margin trend** (Section 8): everything else in
   `ENTRY => S => C* <= 1 + gamma* + eps` is PROVED.

---

## 1. Setting

Post-feed flow, junk magnitudes, caps, i_pf, feed-end state, F1–F4,
Theorem S-cone-fc and its a-priori-F4 corollary: exactly as in the E1
note. Junk j (mixed signs allowed here) decomposes as `j = p - n`,
`p_t := max(j_t, 0)`, `n_t := max(-j_t, 0)`. Boxes at degree d':
cell 0 `[-2, 0]`; cell t >= 1 `[-2C(d'-1,t), +2C(d'-1,t-1)]` (lower =
"n-cap", upper = "p-cap"; note the p-cap at cell t is the n-cap at cell
t-1's index — the mirrored surface). `D_k := d_{i0+k}`,
`delta_k = D_k - D_{k-1}`; caprefs frozen at the certificate row:
`capref_t := 2C(D_0-1, t)` (n-side), `capPref_t := 2C(D_0-1, t-1)`
(p-side).

## 2. Lemma EN-0 (post-feed cell-0 nonpositivity; PROVED, unconditional)

**Statement.** At every post-feed row of every rule-A run (dyadic R >= 32,
any D0 >= 0): `j_0 <= 0`, i.e. `p_0 = 0`.

*Proof.* By T8, at the first post-feed row `j_0 = e(1) = (2 - R) +
2·minus2count`, and minus2count <= (number of rows so far) = i_pf <=
(R-2)/2 by S3 (which needs only g > 1/3), so `j_0 <= 0`. Post-feed cell 0
is autonomous with load `j_0` and box [-2, 0] (kernel reaches cell 0 only
from cells <= 0): if `j_0 <= -2` the junk is `j_0 + 2 <= 0`; if
`j_0 = 0` it stays 0. Nonpositivity is invariant. QED

## 3. Lemma EN-D (signed-part decoupling; PROVED) and Lemma EN-M

**EN-D.** For ANY junk j (mixed signs) and one post-feed row at new degree
d' (either delta):

```
n'_t <= max(0, (K_delta n)_t - 2C(d'-1,t))      (t >= 1)
n'_0 <= max(0, n_0 - 2)
p'_t <= max(0, (K_delta p)_t - 2C(d'-1,t-1))    (t >= 1)
p'_0 <= p_0 .
```

*Proof.* The kernel has nonnegative coefficients, so
`-(K n)_t <= w_t <= (K p)_t`. At cell t >= 1: if `w_t` exceeds the upper
box end, junk `= w_t - 2C(d'-1,t-1) > 0`, so `p'_t = w_t - 2C(d'-1,t-1)
<= (K p)_t - 2C(d'-1,t-1)` and `n'_t = 0`; if `w_t` is below the lower
end, symmetrically `n'_t <= (K n)_t - 2C(d'-1,t)` and `p'_t = 0`;
otherwise both are 0. At cell 0 the load is `j_0 = p_0 - n_0` and the box
is [-2, 0]: positive junk `= j_0 <= p_0`; negative junk
`= -(j_0 + 2) <= n_0 - 2`. QED

**Interpretation.** The negative part of an arbitrary mixed state is
majorized by the S1 flow of the negative part ALONE (positives only
help), and the positive part by its own clamped-Pascal flow with the
mirrored caps `2C(d'-1, t-1)`. The two sign species are dynamically
decoupled as majorants.

**EN-M.** The p-majorant map `Phi_P(p)_0 = p_0`,
`Phi_P(p)_t = max(0, (K p)_t - 2C(d'-1,t-1))` is cellwise monotone
(each output nondecreasing in every input), so the S2 comparison
principle holds for the p-flow verbatim.

**Certificates.** `amm12592_entry_mixedcone_referee_boxeph.out`:
EN-D 600/600, EN-M 600/600 exact pseudorandom trials against an
independent clamp implementation (values to ~2^40 magnitude, both
kernels). ALL-PASS.

## 4. Theorem S-cone-fc± (mixed-sign one-row certificate; PROVED)

Fix a post-feed row i0 with state j = p - n, `m := max(supp p ∪ supp n)`.

**Hypotheses (all exact, all at row i0):**

- **(G1)** supports ⊆ [0, m], `m + 2 < D_0`; `p_0 = 0` (automatic, EN-0);
- **(G2)** `n_0 <= D_0 - 1`;
- **(G3-)** `2 n_{t-1} + n_{t-2} <= capref_t` for every t in [2, m+2];
- **(G3+)** `2 p_{t-1} + p_{t-2} <= capPref_t` for every t in [2, m+2];
- **(G4±)** `i0 + max(ceil(n_0/2), ceil(p_1/2), max_t T_t, max_t TP_t)
  <= R - 2`, with the exact staircase clocks (t in [2, m], plus t = 1
  n-side as in F4):

```
T_t  := min{ K : sum_{k=1..K} (2C(D_k-1,t)   - capref_t)  >= n_t }
T_1  := min{ K : sum_{k=1..K} [2(D_k-1) - (1+delta_k) max(0, n_0-2(k-1))]^+ >= n_1 }
TP_t := min{ K : sum_{k=1..K} (2C(D_k-1,t-1) - capPref_t) >= p_t } .
```

**Conclusion.** For every k >= 0: both parts are cellwise non-increasing
on t >= 1, `n^{(k)}_0 = max(0, n_0 - 2k)`, `p^{(k)}_0 = 0`, supports stay
inside [0, m], no death can occur, and every cell empties by its clock;
junk is empty by the G4± row `<= R - 2`: **capture, and S(R) holds.**

*Proof.* Induction on k with hypothesis: `n^{(k)} <= n`, `p^{(k)} <= p`
cellwise on t >= 1; `n^{(k)}_0 = max(0, n_0 - 2k)`; `p^{(k)}_0 = 0`;
supports ⊆ [0, m]. Row i0+k clamps at degree `D_k >= D_0`, so
`2C(D_k-1,t) >= capref_t` and `2C(D_k-1,t-1) >= capPref_t`.

(i) *n-side, cells t in [2, m]:* by EN-D,
`n'_t <= max(0, n^{(k)}_t + spill - 2C(D_k-1,t))` with
`spill <= 2 n_{t-1} + n_{t-2} <= capref_t` (both kernels dominated), so
`n'_t <= max(0, n_t + capref_t - 2C(D_k-1,t)) <= n_t`, with decay at
least the T_t summand while alive.
(ii) *n-side, cell 1:* inflow `(1+delta_k) n^{(k)}_0 <= 2 n_0 <=
2(D_0-1) <= 2(D_k-1)` = cap (G2 absorbing since n_0 falls and D_k
rises), so `n'_1 <= n_1`, with decay the T_1 summand.
(iii) *p-side, cell 1:* load `p^{(k)}_1 + (1+delta_k) p^{(k)}_0 =
p^{(k)}_1`, p-cap at cell 1 is `2C(D_k-1,0) = 2`:
`p'_1 <= max(0, p^{(k)}_1 - 2)` — an exact >= 2/row drain.
(iv) *p-side, cells t in [2, m]:* `p'_t <= max(0, p^{(k)}_t + spill -
2C(D_k-1,t-1))` with `spill <= 2p_{t-1} + p_{t-2} <= capPref_t`:
non-increasing, with decay at least the TP_t summand.
(v) *Beyond the front:* the only loads are at m+1 and m+2; their n-parts
are `<= 2n_m + n_{m-1} <= capref_{m+1}` and `<= n_m <= capref_{m+2}`
(G3- at t = m+1, m+2), their p-parts `<= 2p_m + p_{m-1} <= capPref_{m+1}`
and `<= p_m <= capPref_{m+2}` (G3+ there): all absorbed (caps at row k
dominate caprefs); supports frozen. Death needs junk at cell
`D_k > m + 2`: impossible.
(vi) All cells non-increasing => spills non-increasing => G3± propagate
at the frozen caprefs. Dead cells stay dead (their load is pure spill
`<= capref <=` current cap). Cell 0: exact 2/row clock (EN-0 invariant
keeps `p_0 = 0`). Cumulative decays are exactly the clock sums; G4±
bounds the last extinction row; then T5/T1 coasting closes the epoch. QED

**Certificate.** `EN_C_mixedcone_step`: 600 exact random mixed-sign
cone states (both deltas), one-step verification of every claim
(support freeze, cellwise non-increase, G3± propagation, all four decay
inequalities) against the independent clamp — ALL-PASS.

**Remarks.** (a) With `p ≡ 0` this is verbatim Theorem S-cone-fc.
(b) The a-priori-F4 corollary does NOT extend to a saturated p-part:
`p_2` may be as large as `C(D_0-1, 2)` under G3+ while its clock
increment is only `2(D_k - D_0)` (linear, not quadratic, in D), giving
`TP_2 ~ D_0/sqrt(2 g) ~ 0.7 R` — budget-infeasible. The p-clocks are
free only under stronger smallness (e.g. `p ≡ 0`, the observed case, or
`p_t <= 2C(D_0-1, t-2)`); otherwise G4± must be checked explicitly at
the row (an O(m · T) integer computation). Recorded so nobody upgrades
the mixed theorem's clocks by analogy.

## 5. Lemma EN-L (debt-layer variant; PROVED) and its honest price

If at i0 all hypotheses of S-cone-fc hold except G2/F2 — say
`n_0 = D_0 - 1 + 2Λ`, `Λ >= 1` — the S4 layer machinery gives a variant:
with `L := Λ` (layer length; n_0 falls 2/row while d is nondecreasing)
and `G_L := sum_{k=0}^{L-1} 2 (n_0 - 2k - (D_0 - 1))^+ <= 2Λ^2 + 2Λ`
(S4: total cell-1 inflow excess during the layer; delta = 0 rows give no
gain), the conclusion of Theorem S-cone-fc holds if F3 is strengthened
at the two lowest inflow cells:

```
(F3'-2)  2(a_1 + G_L) + a_0' <= capref_2      (a_0' := D_0 - 1)
(F3'-3)  2 a_2 + a_1 + G_L   <= capref_3
```

and the cell-1 clock is started from `a_1 + G_L` at row i0 + L (add L
rows to the budget; the W = R/8 window absorbs any L <= R/8 - 66).
*Proof sketch:* during the layer, cells >= 2 obey (i)/(iv)-(vi) verbatim
(their spills involve a_1 <= a_1 + G_L, covered by F3'-2/3); cell 1
gains at most G_L total (S4, proved); after row i0 + L, G2 holds
(absorbing) and the main theorem applies to the majorant state. QED

**Price (exact).** Unconditionally `n_0(i_pf) <= R - 2` (T8), so
`Λ <= (R - 1 - d_fe)/2 <= 0.117 R` (eps = 1/32; d_fe >= Dlb, Section 6)
and the layer always ends INSIDE the certified window — but then
`G_L ~ 2Λ^2 ~ 0.027 R^2 ~ 4.6% of capref_2`: folding F2 away
unconditionally costs ~4.6% of the t = 2 cap, while the measured F3
headroom at the certifying row is 0.1–1%. So the unconditional fold is
NOT affordable — F2 stays a genuine (verified, razor-cheap: the edge law
`a_0 - d_fe in [-14, +9]` makes Λ <= 5 in practice, G_L <= 60 — trivially
affordable) part of EN-CORE. The variant matters exactly when the edge
law's O(10) constant someday needs slack.

## 6. Theorem EN-W (the R/8 window; PROVED + certified 2^9..2^40)

**(W1)** For every dyadic `2^9 <= R <= 2^40`, eps in {1/32, 1/16}, and
ANY post-feed certifying row `i0 <= i_pf + floor(R/8)`: hypotheses
F1 ∧ F2 ∧ F3 imply F4 (the clocks fit the R-2 budget a priori). The
proof is the E1 a-priori-F4 corollary verbatim — F3 bounds every
magnitude (`a_t <= C(D_0-1, t+1)`), the exact Beatty staircase bounds
the cap drift from below, `T_t <= T_2` (t >= 2) by the proved
monotonicity lemma — re-certified with the enlarged
`i_ub = i_pf_ub + floor(R/8)` by exact rationals:
`amm12592_entry_window_constants_boxeph.out`, min budget slack 84 rows
(eps = 1/32) / 89 (eps = 1/16) over the whole grid; the grid value
c_w = 3/16 also passes everywhere (R/8 is kept as the stated window).
**(W2)** The cell-0 drain deadline is window-INVARIANT:
`i0 + ceil(a_0(i0)/2) = i_pf + ceil(a_0(i_pf)/2) <= R - 2`
unconditionally (exact 2/row clock + T8 + S3) — waiting is free on the
drain axis.
**(W3)** Death-freeness of the wait: junk front advances <= 2/row (T6a),
death needs junk at cell `d >= d_fe >= Dlb := floor((2 g_lo R + D0)/
(1 + g_hi)) - 1`, so `m(i_pf) <= m_win := Dlb - 2 floor(R/8) - 2`
(certified `m_win >= 0.511 R` at eps = 1/32, `>= 0.531 R` at 1/16, all
grid points) guarantees no death anywhere in the window regardless of
what the state does. The measured `m(i_pf) ~ 0.93 sqrt(R)` is ~500x
below this threshold at 2^17.

**Consequence (the widened reduction; PROVED).** ENTRY-fc-W(eps): "for
every dyadic R >= 2^18 the rule-A flow survives to some post-feed row
`i0 <= i_pf + floor(R/8)` satisfying F1 ∧ F2 ∧ F3" implies S(eps) for
all dyadic R >= 2^18; with the certificate table below, S(eps) then
holds for ALL dyadic R >= 128, and `C* <= 1 + gamma* + eps`.

## 7. The scans: two new doublings, and the entry anatomy at 2^16–2^17

**Engine and validation.** The feed-phase scanner reuses the certified
transport/clamp/feed loop verbatim; its only liberties are (i) top-band
feed coefficients (indices `>= d_0 - 2`; lower ones are provably never
read; freed as consumed — 2^17 fits in ~2.5 GB), (ii) early stop after
the certificate + 64-row persistence check (Theorem S-cone-fc makes the
tail redundant), (iii) per-row entry diagnostics. Validation: at
(1024, 32) and (4096, 128) the scanner reproduces the stored E1 fcscan
records EXACTLY (i_fc, m, a0, drain, Tmax, capture_ub, budget_margin)
and the stored feed-end snapshots BIT-EXACTLY (full junk dicts); at
(32768, 1024) it reproduces i_fc = 7608, Tmax = 12665, capture_ub =
20273, budget_margin = 12493 — three independent-run agreements.

**New-scale certificates (VERIFIED-exact; S(R) PROVED via Theorem
S-cone-fc at each row):**

| R | D0 | i_pf | i_fc - i_pf | m | a0 - d_fe | Tmax | capture_ub | budget margin | last pos row - i_pf |
|---|----|------|-------------|---|-----------|------|------------|---------------|---------------------|
| 65536 | 2048 | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| 65536 | 4096 | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| 131072 | 4096 | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| 131072 | 8192 | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |

[TBD entries filled from `amm12592_entry_scan_R{65536,131072}_*.json`
at session close; each JSON is the certificate of record.]

## 8. EN-CORE: the isolated minimal hypothesis, with margins

**EN-CORE(eps) (HYPOTHESIS — the single remaining item).** For every
dyadic `R >= 2^18` at `D0 = ceil(eps R)`: the plain rule-A flow survives
to some post-feed row `i0 in [i_pf, i_pf + floor(R/8)]` whose state
satisfies

- (a) `j <= 0` cellwise (equivalently p ≡ 0; else the S-cone-fc± route
  with explicit G3+/G4± is available),
- (b) supp(j) ⊆ [0, m] with `m + 2 < D_0`,
- (c) `a_0 <= D_0 - 1`,
- (d) `2 a_{t-1} + a_{t-2} <= 2C(D_0-1, t)` for every t in [2, m+2].

By EN-W, (a)-(d) at such an i0 imply S(R). EN-CORE is exactly verified
at 22/22 configs (dyadic 2^7..2^17, both eps), always with
`i0 - i_pf <= 8`; margins:

| R | eps | i_fc - i_pf | m / sqrt(R) | max_t r_t at i_pf | min F3 margin at i_fc | a0 - d_fe |
|---|-----|-------------|-------------|--------------------|------------------------|-----------|
| [table filled at session close from the 22 scan JSONs] |

**Convention note (prevents a phantom discrepancy).** The r_t column is
computed against the ROW-i_pf caps `2C(d_{i_pf}-1, t)` — the caps the
F3 check actually uses at the first post-feed row. The E1 feed-end
table used the feed-end degree `d_fe = d_{i_pf - 1}` (caps one degree
lower), so its values exceed these by a factor `~ d/(d-t)` per
delta-step (e.g. 1.0004 there vs 0.999906 here at (16384, 1/32), ratio
`1 + t/d` at t ~ 5). Both are exact; F3-relevance selects this one.

**Margin trend** (all values in the F3 convention of this note).
(i) The feed-end surface excess `eta_meas := max(0, max_t r_t - 1)`
peaks at 0.72% (2^11, eps = 1/16) and shrinks to <= 0.3% from 2^12 and
<= 0.09% from 2^15 — vs the EN-H allowance 4.69%; (ii) the certifying
delay `i_fc - i_pf` stays in [0, 8] across ten doublings; (iii) the
sign-collapse row precedes i_pf everywhere (npos = 0 at i_pf, last
positive cell 3-4 rows earlier); (iv) the support ratio `m/sqrt(R)`
sits at 0.92 +- 0.01 from 2^13 up, so EN-H's lam <= 1/64 holds from
2^13 on (and with orders of slack at the 2^18 frontier).

**The healing mechanism, observed quantitatively (2^15, eps = 1/32).**
At i_pf = 7604: over-cells {3, 5, 7} (single parity), excess 8.9e-4;
the per-delta-row drift at t = 5 is `t/D ~ 2.0e-4` (EN-DR), and the
excess decays 8.9 -> 4.9 -> 2.9 -> 2.9 -> 0 (x1e-4) over exactly the
two delta-1 rows plus two delta-0 rows to i_fc = i_pf + 4 — the EN-H
clock, running 50x faster than its eta_0 worst case. The theorem's
mechanism is the measured mechanism.

**ROUTE-NEGATIVE (exact).** At i_pf the measured r_t on t in [2, 10] is
`1.000 +- 0.01` at every scale >= 2^10 (E1 marginal-surface law,
re-confirmed by these scans): no certificate hypothesis of the form
`spill_t <= (1-c) cap_t`, fixed c > 0.01, holds at feed-end at ANY
tested scale >= 2^10. An entry proof therefore cannot go through a
uniform-margin envelope at the low cells; it must (α) prove the
feed-end state lands on-or-below the surface up to error η(R) -> 0
(measured: η <= 0.01 from 2^10, <= 0.001 from 2^15), and (β) prove the
cap-drift dip (degree growth at fixed cell: factor `D/(D-t)` per
delta-row, i.e. `~ 1 + g t k / D` after k rows) overtakes η within the
window — with the window now Theta(R) (EN-W), (β) needs only
`η(R) <= c t` for some fixed c, an enormously weaker target than the
old 64-row window demanded. The quantitative content of EN-CORE is (α).

## 8b. Theorem EN-H (self-healing entry) — PROVED

The route-negative says any entry proof must handle states ON the
surface. EN-H makes the surface excess a PARAMETER and lets cap drift
heal it inside the R/8 window: **the razor-thin F3 becomes a fat
target.** Ingredients:

**Corollary EN-G (relaxed growth law; from S1).** Post-feed, all-negative
state, one row at degree D': for t >= 1,
`a'_t <= a_t + [spill_t - 2C(D'-1,t)]^+` with `spill_t := 2a_{t-1} +
a_{t-2}` (both kernels' loads are `<= a_t + spill_t`); `a'_0 =
max(0, a_0 - 2)`. An over-cap cell grows per row by at most its spill
excess — never multiplicatively.

**Lemma EN-DR (cap-drift lower bound; PROVED).** With
`s_k := D_k - D_0 >= floor(g_lo k) >= g_lo k - 1`:
`2C(D_k-1, t) >= 2C(D_0-1, t) · (1 + t·s_k/D_K)` for k <= K.
*Proof:* the ratio is `prod_{j=D_0}^{D_k-1} j/(j-t) >= prod (1 + t/j)
>= 1 + t·s_k/(D_k-1) >= 1 + t·s_k/D_K`. QED
(600-trial exact battery `EN_DR` in the referee JSON.)

**Theorem EN-H.** Let i0 be a post-feed row, D_0 := d_{i0},
K := floor(R/8), D_K := d_{i0+K}, and suppose the state a entering row
i0 satisfies, for a rational eta with `eta + Delta-bar(eta) <=` the
certified range (eta <= eta_0 := 3/64 suffices at lam <= 1/64):

- **(E1)** `j <= 0`, supp ⊆ [0, m0], and `lam := (m0+4)/D_0 <= 1/64`
  (the +4 covers the F3 range at the certifying row);
- **(E2)** `a_0 <= D_0 - 1`;
- **(E3)** `spill_t <= (1 + eta) · capref_t` for every t in [2, m0+2],
  `capref_t := 2C(D_0-1, t)`.

Then, with `Delta-bar` the certified fixed point (below) and
`k* := ceil((D_K (eta + Delta-bar)/2 + 1)/g_lo)`:

- the support NEVER leaves [0, m0+2] (front freeze without F3!), no
  death occurs, negativity and E2 persist, and
- **F1 ∧ F2 ∧ F3 hold exactly at every row `i0 + k` with
  `k* <= k <= K`** — in particular, applied at `i0 = i_pf` (the
  EN-CORE' use; for later i0 the same argument works whenever
  `i0 + k* <= i_pf + K`), the S-cone-fc certificate fires inside the
  EN-W window and S(R) follows (a-priori F4, Section 6).

*Proof.* Write `E^(k)_t := [a^(k)_t - a_t]^+` (excess over the initial
state), `cap^(k)_t := 2C(D_k-1,t)`, `phi_k := s_k/D_K`, and for t >= 2
`Delta_t := (2 e_{t-1} capref_{t-1} + e_{t-2} capref_{t-2})/capref_t`
where `e_t` are the bounds established below (`E^(k)_t <= e_t capref_t`
for ALL k <= K). The argument is an induction on the cell index t (the
kernel is lower-triangular: cell t is loaded only from t, t-1, t-2), and
inside it an induction on k.

(0) *Cells 0 and 1 never grow:* cell 0 drains exactly 2/row (S1); cell
1's load is `a_1 + (1+delta) a^(k)_0 <= a_1 + 2(D_0-1) <= a_1 +
cap^(k)_1` (E2 absorbing), so `a'_1 <= a_1`. Hence `e_0 = e_1 = 0` and
`Delta_2 = 0`.

(1) *Excess recursion.* Fix t >= 2 and assume the bounds `e_{t-1},
e_{t-2}` hold for all k. The row-k spill obeys `spill^(k)_t <=
spill^(0)_t + 2E^(k)_{t-1} + E^(k)_{t-2} <= (1 + eta + Delta_t)
capref_t` (E3; and `spill^(0)_t = 0 <= (1+eta) capref_t` for
t > m0+2). By EN-G and EN-DR the per-row excess increment is at most
`capref_t · [eta + Delta_t + t/D_K - (t g_lo / D_K) k]^+` (the `t/D_K`
absorbs the -1 in the staircase floor). Summing the arithmetic series
(`sum_k [A - Bk]^+ <= A^2/(2B) + A`):

```
e_t <= (X_t)^2 · D_K/(2 t g_lo) + X_t ,   X_t := eta + Delta_t + t/D_K .
```

(2) *Compounding is a contraction.* Substituting (1) into the
definition of `Delta_{t+1}` and using the EXACT cap ratios
`capref_{t-1}/capref_t = t/(D_0-t)` (so `D_K/t · t/(D_0-t) <=
gamma/(1-lam)`, `gamma := D_K/D_0` — the D_K/t divergence cancels
against the cap ratio, cell by cell), `(t+1)/t <= 3/2` (t >= 2), and
`t/D_K <= lam` for t <= m0+4:

```
Delta_{t+1} <= C2 · X-bar^2 + C1 · X-bar ,  X-bar := eta + Delta-bar + lam,
C2 := (3/2) gamma / ((1-lam) g_lo) + 3 gamma lam / (2 g_lo (1-lam)^2) <= 3.0,
C1 := 2 lam/(1-lam) + lam^2/(1-lam)^2 <= 0.033 .
```

So every `Delta_t <= Delta-bar` where `Delta-bar` is the smallest root
of `C2 (eta+lam+D)^2 + C1 (eta+lam+D) = D` — which exists, with the
iteration a contraction, throughout the certified range
(`amm12592_entry_selfhealing_constants_boxeph.out`: Delta-bar <= 0.0241
and contraction <= 0.53 at eta = 3/64, lam = 1/64; exact rationals, all
dyadic R = 2^9..2^40, both eps; the pads C2 <= 3.0, C1 <= 0.033 and
gamma <= 1.105 verified per-R). The computed values sit ~5% below the
pads; the remaining `O(1/D_0)` slops (the `-1`s in `D_0 - 1 - t`, the
staircase floor, the ceil in k*) are absorbed by that headroom for
every R in the certified range (D_0 >= 2^9 · 0.74).

(3) *Front freeze.* For t > m0+2 the initial value and initial spill
are 0, so the load is pure excess `<= Delta_t capref_t <= Delta-bar
capref_t < cap^(k)_t` (Delta-bar < 1): absorbed entirely — junk never
appears beyond m0+2, death is impossible (m0+4 << D_0), and the t-range
[2, m0+4] used above covers every live spill.

(4) *Healing.* For `k >= k*`: `t phi_k >= (t/D_K)(g_lo k - 1) >=
eta + Delta-bar >= eta + Delta_t` for every t >= 2 (worst t = 2 by the
definition of k*), so by EN-DR `cap^(k)_t >= (1 + eta + Delta_t)
capref_t >= spill^(k)_t`: **F3 holds at row i0 + k for every k in
[k*, K]**, F2 holds (cell-0 drain + D nondecreasing), F1 holds
(negativity C-N; support ⊆ [0, m0+2] by (3)). The certified constants
give `k* <= K - 34` at every grid point. QED

**Certificates.** `amm12592_entry_selfhealing_constants_boxeph.{out,json}`
(exact-rational: fixed point, contraction, pads, k* <= K, dyadic
2^9..2^40 x eps in {1/32, 1/16} x eta in {1/32, 3/64}: ALL-PASS);
600-trial one-step batteries EN-G (relaxed states) and EN-DR in
`amm12592_entry_mixedcone_referee_boxeph.json` (rerun of this session's
battery with the EN_G/EN_DR trial families added).

**What EN-H buys.** EN-CORE's razor condition (d) is now REPLACED by
the fat condition (E3) at i0 = i_pf itself: the feed phase need only
deliver a state within 4.7% of the marginal surface (eta_0 = 3/64),
while the measured feed-end excess is <= 0.72% at every scale >= 2^10
and <= 0.1% from 2^13 on, SHRINKING in R — a >= 6x margin at the worst
scale ever seen and >= 45x at the frontier scales. The remaining
unproved content of ENTRY is exactly:

**EN-CORE' (final form; HYPOTHESIS).** For every dyadic R >= 2^18 at
D0 = ceil(eps R): the rule-A flow survives to i_pf (Theorem B covers
all but the <= 2 rows between i_feed and i_pf at D0 >= eps* R) and its
feed-end state satisfies (E1) support (m0 + 2 <= D_0/64: measured
~0.93 sqrt(R), i.e. 100x slack at 2^17), (E2) debt edge (measured
|a_0 - d_fe| <= 14, vs the D_0 - 1 requirement met within <= 6 rows in
every run — see the layer remark, Section 5), and (E3) the
4.7%-relaxed surface (measured excess <= 0.72% worst-ever, <= 0.1% at
scale). EN-CORE' implies EN-CORE, hence S(eps), hence the C* bounds.

## 9. The final theorem (conditional form) and the dependency chain

**Theorem (assembly; PROVED modulo EN-CORE').** Fix eps in {1/32, 1/16}.
If EN-CORE'(eps) holds (the eta_0-relaxed feed-end state property,
Section 8b), then:

- S(eps) holds for every dyadic R >= 128 [finite scales by the
  certificates + Theorem S-cone-fc; R >= 2^18 by EN-CORE' + EN-H +
  EN-W + Theorem S-cone-fc];
- every dyadic epoch closes at slack ceil(eps R) [Theorems A (LIFT),
  B (feed survival), S-cone machinery];
- `C* <= 1 + gamma* + eps` [THM-3329 assembly + LIFT]. In particular
  EN-CORE'(1/32) gives `C* <= 1 + gamma* + 1/32 < 427095/262144
  = 1.6292382`.

If EN-CORE'(eps_j) holds for a rational sequence eps_j ↓ eps* =
2(1-g-g^2)/(3+2g), then

```
C* <= 1 + gamma* + eps* = (5+3g)/(3+2g) in (1.6191617801, 1.6191618342),
```

the E1 closed form — the sharpest constant this route can produce.
(For the eps ↓ eps* form, Theorem B's window (ii) and the EN-W/EN-H
constant certifications — stated at eps in {1/32, 1/16} — must be
re-instantiated at each eps_j; every step is the same exact-rational
computation, none is asymptotic.)

**Dependency chain (every link's status):** T2/T6 conjugacy (PROVED) ->
T4/T4b initial data (PROVED) -> Theorem B feed survival at D0 >= eps* R
(PROVED) -> EN-0/EN-D/EN-M/S1–S4 (PROVED, machine-certified) ->
Theorem S-cone-fc / S-cone-fc± (PROVED, machine-certified) -> EN-W
window + a-priori F4 (PROVED + exact-rational certification 2^9..2^40)
-> **Theorem EN-H self-healing (PROVED + exact-rational constants
2^9..2^40 + 600-trial one-step batteries)** -> EN-CORE' (HYPOTHESIS;
the eta-relaxed form VERIFIED-exact with 6x-45x margins at 2^7..2^17,
both eps) -> S(eps) -> LIFT + THM-3329 assembly (PROVED) -> C* bound.

## 10. Hazards honored

- Rule/search negatives never prove infeasibility; the ROUTE-NEGATIVE of
  Section 8 is a statement about one proof strategy's hypotheses being
  false at feed-end, not about EN-CORE.
- Quantifier hygiene: EN-CORE is per-(R, D0) with an explicit window and
  explicit R-range; the new-scale S claims are per-(R, D0) theorems
  (S-cone-fc + certificate), not asymptotic claims; no non-dyadic claims.
- The scanner's early stop is licensed by a PROVED theorem; the stored
  E1 full runs (<= 2^15) provide the independent consistency checks
  (capture rows within 1–90 of the bounds, F3 persistence through
  ~13000 rows at 2^15).
- The mixed theorem's clock caveat (Remark (b), Section 4) is recorded
  to prevent an invalid a-priori upgrade of G4±.
- Memory-trimmed feed coefficients are provably unread below the trim
  line (feed indices are >= d_0 - 1 and strictly increase); validation
  is bit-exact against the untrimmed engine at three scales.
- All decisions exact int/Fraction; floats only in display fields.

## 11. Status ledger

- **PROVED (new):** EN-0 (unconditional post-feed p_0 = 0); EN-D
  signed-part decoupling; EN-M p-flow comparison; **Theorem S-cone-fc±**
  (mixed-sign one-row certificate => S(R)) with the G4± caveat; EN-L
  layer variant with explicit G_L price; **Theorem EN-W** (R/8 window:
  wide a-priori F4, drain invariance, death-free wait) with
  exact-rational certification at all dyadic 2^9..2^40, both eps;
  **Theorem EN-H (self-healing entry: eta_0 = 3/64 surface excess is
  drift-curable, front frozen without F3, contraction-controlled
  compounding)** with EN-G/EN-DR one-step lemmas, exact-rational
  constant certification 2^9..2^40 x both eps x eta in {1/32, 3/64},
  and k* <= K - 34 at every grid point.
- **VERIFIED-exact (new):** scanner bit-identical to the E1 records
  (fcscan fields at 1024/4096/16384/32768; feed-end snapshots
  bit-exact; v2 field-identical to v1 at three scales); the
  independent certificate referee (fresh clamp + mul//div clocks)
  ALL-PASS at 1024, 16384, 32768; NEW-SCALE certificates (Section 7
  table); sign collapse before feed-end at every config; the EN-CORE
  margin table with trends; 600-trial batteries
  EN-D/EN-M/EN-C/EN-G/EN-DR (ALL-PASS).
- **ROUTE-NEGATIVE (new, exact):** fixed-margin (c > 1%) F3 envelopes
  are false at feed-end at every scale >= 2^10 — entry proofs must
  produce the marginal surface with o(1) error + use the cap-drift dip
  (which is exactly what EN-H formalizes); the naive uniform Gronwall
  inside the healing argument diverges (recorded in Section 8b's proof
  discussion) — the per-cell cap-ratio cancellation is essential.
- **HYPOTHESIS (one item):** EN-CORE'(eps) for dyadic R >= 2^18 — the
  feed-end state lands within eta_0 = 3/64 of the marginal surface,
  with sqrt(R)-scale support and the debt at the edge. This is the
  entire remaining unproved content of `C* <= 1 + gamma* + 1/32` and,
  along eps ↓ eps*, of `C* <= (5+3g)/(3+2g) < 1.6191619`. Measured
  margins on (E3): 6x at the worst scale ever seen, >= 45x at the
  frontier, growing in R.
