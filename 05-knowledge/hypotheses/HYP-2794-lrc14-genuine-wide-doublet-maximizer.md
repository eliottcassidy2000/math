# HYP-2794 — The genuine-wide p0 maximizer is `consec_{k-1} + tight far doublet`; it does NOT close by THM-563 periodicity (but is "almost periodic")

- **Instance:** claude-opus-2026-06-21
- **Status:** PARTIALLY-TRUE (maximizer family CONFIRMED exact k=8..13; cap-margin ≥0.16. THM-563-periodicity closure REFUTED for the doublet.)
- **Touches:** OPEN-Q-108 leg (C) genuine-wide; THM-563, THM-557, HYP-2788, HYP-2775.

## Claim / Findings

**(1) CONFIRMED — the genuine-wide maximizer is a stable 2-cluster family.**
Over an exact-rational adversarial search (two/three-cluster, dilated-AP+perturbation,
random; `lrc14_genuine_wide_maximizer_claudeopus_0621.py`), the genuine-wide p0
maximizer for k=9,10,11,12 is

> **E\* = consec_{k-1} ∪ {21, 22} = {0,1,…,k-2} ∪ {21,22}**  (consec base + a TIGHT far DOUBLET).

This is the `[m-2, 2]` coherent-block partition of the m=k-1 nonzero speeds — the
**largest genuine-wide-compatible partition**, because `[m-1,1]` (= single-far) is
EXCLUDED from genuine-wide (removing the lone far element bounds the rest, so it is
single-perturbation-bounded, handled by THM-563). At k=8 the maximizer is the
3-cluster `(0,1,2,10,11,12,20,21)` (= `[3,3,2]`); at k=13 it is `consec_12 ∪ {15,16}`.

Exact values (genuine-wide = span>14 AND remove-any-one stays wide after reprimitivize):

| k | max p0 (genuine-wide) | argmax | cap_k − max | Q(k-1) − max |
|---|---|---|---|---|
| 8 | 929/5390 ≈ 0.17236 | (0,1,2,10,11,12,20,21) | +0.20911 | +0.02424 |
| 9 | 9371/32340 ≈ 0.28976 | consec_8 ∪ {21,22} | +0.20449 | +0.07233 |
| 10 | 1301/2940 ≈ 0.44252 | consec_9 ∪ {21,22} | +0.16188 | **+0.00537** |
| 11 | 5617/10780 ≈ 0.52106 | consec_10 ∪ {21,22} | +0.20422 | +0.01019 |
| 12 | 14302/24255 ≈ 0.58965 | consec_11 ∪ {21,22} | +0.26749 | +0.01259 |
| 13 | 11423/17640 ≈ 0.64756 | consec_12 ∪ {15,16} | n/a (cap_13=1) | n/a |

**KEY PROOF-RELEVANT POINT:** the gap to **cap_k** is comfortable everywhere (≥ 0.16),
GROWING-ish in k. The razor-thin gap is only to the *decorrelated floor* Q(k-1)
(0.0054 at k=10). So the LRC-relevant statement **genuine-wide ⟹ p0 < cap_k holds
with margin ≥ 0.16**; the sharper `< Q(k-1)` framing (HYP-2788) is barely true at k=10
and is NOT the quantity the cap-proof actually needs. **Recommend the genuine-wide leg
target be `p0 < cap_k` (margin ≥ 0.16), not `p0 < Q(k-1)`.**

**(2) REFUTED — the doublet does NOT close by THM-563's exact periodicity.**
Natural idea: THM-563 closes single-far because `w·Δ_w` is exactly periodic in w
(period 7·lcm(base)). For the doublet `{w, w+1}` the two locked phases `frac(w·t)`,
`frac((w+1)·t)=frac(w·t+t)` are each periodic in w — suggesting `w·Δ_w^{doublet}`
might also be periodic. **It is NOT.** Test `lrc14_doublet_periodicity_claudeopus_0621.py`
(second difference of `S_w=w·p0(E_w)` over step P must vanish iff periodic):
worst 2nd-difference = 44/3055 (k=6), 663/24304 (k=7), 347/15925 (k=8) — all NONZERO.
**Reason:** the doublet's two far speeds interact with each other and the base at
breakpoints `j/(7w)` that depend on w; these are NOT captured by the base's fixed arcs
`A_j`. The period-max mechanism is single-far-specific.

**(3) LEAD — the doublet is "ALMOST periodic" (small decaying correction).**
The 2nd-differences are small (≈0.01–0.03) and appear to decay in w. So plausibly
`w·Δ_w^{doublet} = (exactly-periodic part) + (correction → 0)`. If the correction is
sign-controlled / decaying like `O(1/w)`, the doublet closes via an *asymptotic*
period-max + a finite-w window check — a THM-563 analogue with an error tail. This is
the concrete next target for leg (C)'s binding case.

## Why it matters
Leg (C) (genuine-wide) is the one live mathematical target of OPEN-Q-108. This pins the
binding genuine-wide config to a clean 2-parameter family (consec base + far doublet),
reframes the target to the robust `p0 < cap` (margin ≥ 0.16), and kills the tempting
"doublet = periodic like single-far" route while extracting the right residual structure
(almost-periodic + decaying correction).

## Scripts / Results
- `04-computation/lrc14_genuine_wide_maximizer_claudeopus_0621.py` → `05-knowledge/results/lrc14_genuine_wide_maximizer_claudeopus_0621.out`
- `04-computation/lrc14_doublet_periodicity_claudeopus_0621.py` → `05-knowledge/results/lrc14_doublet_periodicity_claudeopus_0621.out`

→ OPEN-Q-108, THM-563, THM-557, HYP-2788, HYP-2775, HYP-2774.
