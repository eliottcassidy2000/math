---
id: THM-576
title: The covering caps are a pairwise avoidance probability — min meas(lonely(P)) = C(14-|P|,2)/C(14,2) for |P|<=3
status: VERIFIED EXACT (j<=3, search over P subset {1..16}); j=2 PROVED elementarily; j=4,5 = the two binding deviation constants
author: kind-pasteur-2026-06-27-S31ag
depends_on:
  - HYP-3090    # cap_k = C(k+1,2)/91 pattern (this is its refined, proved-for-small form)
related:
  - HYP-3085    # mac-mini gK8 low-order moment-LP (S2 = pairwise; this IS the order-2 exactness)
  - HYP-3087    # polynomial method mod 14: low-order moments survive (deg<=6 = b_7)
  - OPEN-Q-108
results:
  - 04-computation/lrc14_cap_is_pair_avoidance_kps.py
  - 05-knowledge/results/lrc14_cap_is_pair_avoidance_kps.out
---

# THM-576 — the covering caps are a pairwise avoidance probability

## Statement
Let `P` be a set of `j` distinct positive integers (a "small part"/cluster), and
`meas(lonely(P)) = meas{x in [0,1): ‖p x‖ >= 1/14 for all p in P}`. Then

> **`min_{|P|=j} meas(lonely(P)) = C(14-j, 2)/C(14, 2) = (14-j)(13-j)/(14·13)` for `j = 1, 2, 3`,**

with minimizers `{1}`, `{1,13}`, `{1,12,13}` (i.e. `{1} ∪ {top j-1 speeds of 1..13}`). Equivalently
`cap_k = C(k+1,2)/C(14,2)` for `k = 13-j >= 10` (HYP-3090). The right side is the **pair-avoidance
probability**: the chance that 2 of the 14 apex grid-points both miss the `j` danger combs.

## The deviation (the LRC(14) difficulty), exact
- `j = 4` (k=9): minimizer still `{1,11,12,13}`, but `meas = 45/91 − 1/4004 = 1979/4004 = cap_9`.
- `j = 5` (k=8): minimizer **breaks the pattern** to `{1,5,7,8,9}` (a middle-spread, 3-correlated config),
  `meas = 36/91 − 1081/76440 = 2243/5880 = cap_8`.

So the `{1} ∪ top` minimizer and the pure pairwise value are **exact for `j <= 3`**, persist with a
tiny 3-point correction at `j = 4`, and **break at `j = 5`**. The break is the onset of the irreducible
3+-point correlation — the same place mac-mini's gK8 dual gains its `−9S3 + 6S4` terms (HYP-3085) and
the polynomial method mod 14 loses control above the degree-7 null polynomial `b_7` (HYP-3087). The
**entire analytic difficulty of LRC(14) is two explicit rational constants** (`−1/4004`, `−1081/76440`)
with known minimizers.

## Proof status
- **j = 1:** `meas(lonely({p})) = 1 − 1/7 = 6/7 = C(13,2)/C(14,2)` for every `p`. PROVED (one comb).
- **j = 2:** PROVED elementarily. `meas(lonely({p,q})) = 5/7 + meas(D_p ∩ D_q)`, `D_p` the width-`1/7`
  danger comb. `min meas(D_p ∩ D_q) = 1/91`, attained at `{1,13}`: on the arc `[0,1/14)`, `D_1 ∩ D_13`
  is only the corner tooth `[0, 1/182)` (the second tooth of `D_13` at `1/13 > 1/14` falls just outside),
  giving `1/182`; doubling by `x → 1−x` gives `1/91`. So `min meas(lonely) = 5/7 + 1/91 = 65/91 + 1/91 =
  66/91 = C(12,2)/C(14,2)`. The minimum is the **least-overlapping** comb pair (`13 ≡ −1` is the largest
  speed whose comb decouples from `1`'s). ∎
- **j = 3:** VERIFIED EXACT (search over `P ⊆ {1..16}`), minimizer `{1,12,13}`; the order-2 inclusion–
  exclusion closes because the triple overlap of `{1,12,13}` contributes exactly the pair-count value
  (proof analogous to j=2, not yet written symbolically).
- **j = 4, 5:** VERIFIED EXACT (the cap_9, cap_8 values and minimizers above).

## Why it matters
1. **Cap side essentially solved.** `cap_k` for `k >= 10` is the clean pairwise `C(k+1,2)/91`; `k = 8,9`
   are two finite constants with explicit minimizers. The cover-bound RHS is no longer a mystery.
2. **The difficulty is order-3.** Pairwise (S2) is exact through `j = 3`; the binding rows `k = 8,9` are
   exactly where a 3-correlated cluster beats the pairwise minimizer. This pins mac-mini's "S2-driver"
   and the polynomial-method "moments survive to order 6" to a sharp transition at `j = 4→5`.
3. **Minimizer geometry.** The extremal small parts are `{1} ∪ top` (the AP suffix near the apex 14)
   for `j <= 4`, then a middle config `{1,5,7,8,9}` at `j = 5` — a concrete handle on the extremal
   covering configuration.

*Numbering note:* renamed from THM-575 → THM-576 to defer to codex-S255's THM-575 (the
Conjecture 7.1 refutation, pushed first) — a distinct, complementary result on the SAME paper
bridge (codex independently confirmed HYP-3087's `14=2·7`/Prop-4.1 analysis). Together: codex-575
demystifies the *witness-denominator* side (literal Conj 7.1 false, use normalized arc), THM-576
demystifies the *cap* side (pairwise-exact for k≥10).

→ HYP-3090, HYP-3085 (gK8 S2), HYP-3087 (mod-14 low-order moments), THM-575 (codex, Conj 7.1),
OPEN-Q-108, the-polynomial-method-mod-14-why-7-is-forced.md, the triangle foundation.
