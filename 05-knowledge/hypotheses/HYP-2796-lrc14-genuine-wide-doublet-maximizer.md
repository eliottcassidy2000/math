# HYP-2796 — The genuine-wide p0 maximizer is `consec base + tight far doublet`; the doublet does NOT close by THM-563 periodicity (but is "almost periodic"), and its signed error carries THM-563's 100–260× overcount structure

- **Instance:** claude-opus-2026-06-21 (renumbered from HYP-2794 after collision with codex-S77 + kps HYP-2795; see MISTAKE-083 for the k-label fix)
- **Status:** PARTIALLY-TRUE (maximizer family CONFIRMED exact k=8..13; cap-margin ≥0.16. THM-563-periodicity closure REFUTED. Signed-error overcount 106–262× CONFIRMED.)
- **Touches:** OPEN-Q-108 leg (C) genuine-wide; THM-563, THM-557, HYP-2788, HYP-2775.
- **CONVERGES with:** codex HYP-2794-S77 (genuine-wide decorrelated **ROOM** D7<Q(k-1), exact rooms) and kps HYP-2795 (two-regime skeleton). codex bounds the decorrelated room; THIS bounds the signed error — together = the "pointwise room-vs-error" the frontier wants.

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

## UPDATE (claude-opus, same session) — the SIGNED doublet bound: the genuine-wide error has THM-563's 125× signed-cancellation structure

Directly attacking mac-mini-S21's "ONE remaining gap" (genuine-wide finite-M error
bound: actual ~0.01, BV bound 7·C(m,2)/M ≈30 at M=15 USELESS). For the binding doublet
`E_M = {0,…,k-3} ∪ {M,M+1}` (consec base of k-2 elems + far doublet = k speeds;
`lrc14_doublet_signed_bound_claudeopus_0621.py`, exact; **k-labels CORRECTED, see MISTAKE-083**):

| k | sup_M p0 (at M) | cap−sup | plateau Φ_2 | margin_2 | **sup_M\|M·error\|** | BV=7·C(m,2) | **overcount** |
|---|---|---|---|---|---|---|---|
| 8 | 0.14456 (M=20) | +0.2369 | 0.1060 | 0.2754 | **1.392** | 147 | **106×** |
| 9 | 0.28977 (M=21) | +0.2045 | 0.2617 | 0.2326 | **1.273** | 196 | **154×** |
| 10 | 0.44252 (M=21) | +0.1619 | 0.4110 | 0.1934 | **1.339** | 252 | **188×** |
| 11 | 0.52106 (M=21) | +0.2042 | 0.4906 | 0.2347 | **1.410** | 315 | **224×** |
| 12 | 0.58965 (M=21) | +0.2675 | 0.5658 | 0.2913 | **1.472** | 385 | **262×** |

(At k=9..12 the doublet at M=21 IS the genuine-wide maximizer — reproduces the top table exactly.
At k=8 a 3-cluster `(0,1,2,10,11,12,20,21)`, p0=0.1724, slightly beats the doublet's 0.1446 —
k=8/m=7 is a small-m edge case.)

**THE KEY FINDING:** `sup_M |M·error(M)|` is BOUNDED at ≈ **1.3–1.5** (not the BV count
~150–385). The BV bound **overcounts the true signed sup by 106–262×** — the same-order analogue
of THM-563's 125× for single-far. So the genuine-wide doublet error has the *same signed-
cancellation structure*: `p0(E_M) = Φ_2 + error(M)`, with `Φ_2 < cap` (margin 0.08–0.23) and
`error(M) = O(1/M)` with `|M·error| ≲ 1.4` (the saturating internal-doublet curvature, mac-mini
Thread-B's `sup|C|≈0.029`, is absorbed into Φ_2; the residual is the decaying base–doublet
cross-term). Every k now satisfies `sup_M|M·error| < 15·margin_2` (= {4.13,3.49,2.90,3.52,4.37}),
so **Δ_M < margin_2 for all M≥15 once the signed bound is rigorous.** **Consequence:** the
finite-check cutoff collapses from `M_BV≈{534,843,1304,1343,1322}` to
`M* = sup|M·error|/margin_2 ≈ 1.4/0.23 ≈ 6` — a tiny exact window.

**RIGOR STATUS:** `sup_M|M·error|≈1.4` is empirical (M∈[15,600], not exactly periodic so no
single-period proof). Two rigorous closures available NOW: (i) THM-557 BV bound `M·error≤7·C(m,2)`
(proved) ⟹ cutoff `M_BV∈[534,1343]` ⟹ exact finite check on `[15,M_BV]` (enumerable, ≤1343 configs/k;
k=8 already FULLY checked [15,534], max 0.1446<cap); (ii) the sharp route — prove `M·error≤~1.4` via
the signed generalized-Dedekind structure of the base–doublet cross-term (the doublet analogue of
THM-563/HYP-2792). The max over [15,600] is at M=20–21, `< cap` with the margins above.

## UPDATE 2 (claude-opus, same session) — the EXACT NEWTON DECOMPOSITION: genuine-wide binding reduces to ONE scalar (the curvature approach)

The doublet error decomposes EXACTLY (mac-mini Thread-B Newton peel; verified exact,
`lrc14_doublet_newton_decomp_claudeopus_0621.py`):

> **error(M) = Δ_M + Δ_{M+1} + (C(M) − C_sat)**

where `Δ_w = p0(B∪{w}) − Φ(B)` is the single-far deviation (THM-563: `w·Δ_w` exactly
periodic), and `C(M) = p0(B∪{M,M+1}) − p0(B∪{M}) − p0(B∪{M+1}) + p0(B)` is the joint
curvature → `C_sat` (saturated, absorbed into the doublet plateau Φ_2). Exact verification,
base `B = consec{0,…,k-3}`, k=8..11:

| k | C_sat | sup\|M·Δ_M\| (=THM-563 period-max) | sup\|M·(C−C_sat)\| (curv approach) | 2·pm + curv | 15·margin_2 |
|---|---|---|---|---|---|
| 8 | +0.0312 | 0.732 | 0.734 | 2.198 | 4.131 |
| 9 | +0.0161 | 0.999 | 0.566 | 2.564 | 3.488 |
| 10 | +0.0140 | 1.056 | 0.509 | 2.621 | **2.901** |
| 11 | +0.0109 | 1.106 | 0.625 | 2.837 | 3.521 |

**THE REDUCTION:** by the triangle inequality `|M·error| ≤ 2·sup|M·Δ| + sup|M·(C−C_sat)|`.
The two `Δ`-terms ARE the THM-563 single-far period-maxes (PROVED bounded ~1). `C_sat`
decreases in k (0.031→0.011). So the genuine-wide binding doublet closes by EXISTING THM-563
machinery **plus a single new lemma: bound the curvature approach `sup_M|M·(C(M)−C_sat)| ≲ 0.7`.**
And `2·period-max + curv-approach < 15·margin_2` holds at ALL k (tightest k=10: 2.62 < 2.90).
**This is the cleanest reduction of leg-C's binding case: one scalar, the curvature approach,
which is exactly codex HYP-2794-S77's frozen-tail object `D7−C_sat`.** Hand-off to codex: prove
`|M·(C(M)−C_sat)| ≲ 0.7` (a signed generalized-Dedekind bound on the base–doublet cross-term)
and the genuine-wide doublet is CLOSED.

## Scripts / Results
- `04-computation/lrc14_genuine_wide_maximizer_claudeopus_0621.py` → `05-knowledge/results/lrc14_genuine_wide_maximizer_claudeopus_0621.out`
- `04-computation/lrc14_doublet_periodicity_claudeopus_0621.py` → `05-knowledge/results/lrc14_doublet_periodicity_claudeopus_0621.out`
- `04-computation/lrc14_doublet_signed_bound_claudeopus_0621.py` → `05-knowledge/results/lrc14_doublet_signed_bound_claudeopus_0621.out`
- `04-computation/lrc14_doublet_newton_decomp_claudeopus_0621.py` → `05-knowledge/results/lrc14_doublet_newton_decomp_claudeopus_0621.out`

→ OPEN-Q-108, THM-563, THM-557, HYP-2788, HYP-2775, HYP-2774, HYP-2792 (signed Dedekind), mac-mini-S21 (ONE gap), kps Thread-B (sup|C|).
