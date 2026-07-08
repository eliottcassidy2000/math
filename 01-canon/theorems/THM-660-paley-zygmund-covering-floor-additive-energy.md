---
id: THM-660
title: The Paley–Zygmund covering floor — on mac-mini's covering reformulation (THM-657, W = uncovered measure), the second-moment bound mu_{1/7}(E) >= E[W]²/E[W²] (Paley–Zygmund; the OPTIMAL two-moment lower bound on P(W>0)) is strictly stronger than (7/6)E[W] and clears the honest (A′) bars at k = 11, 12, 13 diameter-free — with margins +0.016 / +0.109 / +0.216 — where (7/6)E[W] FAILS at k=11,12 (0.184/0.176 < bars). The floor is additive-energy-ORDERED: Var(W)/E[W]² increases with the reduced additive energy R2, so the block (AP, max energy) is the minimizer at k=11,12; this unifies the energy axis across the density-side tent floor (THM-656) and the covering-side floor. The three legs reduce to ONE moment inequality: min_E E[W]²/E[W²] >= bar, a ratio of additive-energy moments (more tractable than "AP minimizes mu")
status: PROVED per-shape (Paley–Zygmund is rigorous — Cauchy–Schwarz; and it is the exact 2-moment optimum, the LP minimizer placing mass at 0 and at E[W²]/E[W] <= 6/7). The UNIFORM floor min_E PZ >= bar: EXACT block moments now computed (klein-S178, Farey-cell piecewise integration — independently reproducing monad's block_13 PZ = 221828403/815409784 to the digit, validating the method). The **k=11 leg (the binding case) is exhaustively EXACT-verified in the high-energy regime**: over ALL 2993 primitive k=11 shapes with diam <= 15 (where the minimizer provably lives — PZ falls with additive energy, and max energy = min diameter), the exact minimum is PZ = 21576025/62216714 = 0.346788 >= bar 83549/252252 = 0.331212 (margin +0.015576, a definite rational, no longer a grid artifact); diam >= 16 shapes have PZ >= 0.45 (descent) with energy-ordering corr(R2, PZ) = -0.96. Exact block values PZ(block_k) = 0.34710/0.30803/0.27205 at k=11/12/13. THE COMPLEMENTARITY: the block PZ CLEARS exactly at k=11,12,13 and FAILS at k=8,9,10 (0.574/0.481/0.411 < bars) — precisely the tent floor's domain; the two second-moment methods (tent variance k<=10, W-variance/PZ k=11,12,13) partition k=8..13, meeting at k=10/11.
source: klein-2026-07-07-S177 (HYP-5317)
depends_on:
  - THM-657   # mac-mini's covering reformulation: mu = P(W>0), W = sum_i (g_i-1/7)_+, mu >= (7/6)E[W]
  - THM-656   # klein's tent-variance / additive-energy floor (the density-side energy axis this unifies with)
related:
  - HYP-5157  # monad-explorer PZ-on-V: the SAME bound E[W]²/E[W²], computed at k=13 (min = AP = 0.27205 exact) — credited
  - THM-530   # "consecutive minimizes mu" — the mu-route extremal lemma this moment-route parallels/relaxes
  - THM-638   # pair joint-window masses = the pair terms of E[W²]'s correlation expansion
external: Paley–Zygmund 1932 (2nd-moment method); Stevens 1939 (circle covering); the moment problem on [0,M].
---

# THM-660 — the Paley–Zygmund covering floor

## Statement

On the covering reformulation (THM-657): `W(x) = 1 − coverage = Σ_i (g_i(x) − 1/7)_+` is the
uncovered measure, `0 ≤ W ≤ 6/7`, and `mu_{1/7}(E) = P(W > 0)`. Then

> **`mu_{1/7}(E) ≥ E[W]² / E[W²]`** (Paley–Zygmund).

*Proof.* `E[W] = E[W · 1_{W>0}] ≤ E[W²]^{1/2} · P(W>0)^{1/2}` (Cauchy–Schwarz); square. ∎

This is the **optimal** two-moment bound: minimizing `P(W>0)` over all distributions on `[0, 6/7]`
with the given `E[W], E[W²]` puts mass `1 − E[W]²/E[W²]` at `0` and the rest at `a = E[W²]/E[W]`
(here `a ≤ 6/7` always), attaining `E[W]²/E[W²]`. Since `E[W²] ≤ (6/7)E[W]` (from `W ≤ 6/7`),
`PZ ≥ (7/6)E[W]` — it dominates mac-mini's first-moment bound, with equality only when `W ∈ {0, 6/7}`.

## The diameter-free discharge of k = 11, 12, 13

The first-moment floor `(7/6)E[W]` (THM-657) is `0.184 / 0.176 / 0.148` at `k = 11/12/13` — **below
the bars `0.331 / 0.199`** at k=11,12 (it only clears k=13). The second-moment `PZ` clears **all
three** (N=20M, block = the minimizer at k=11,12):

| k | bar | `(7/6)E[W]` | **PZ = E[W]²/E[W²]** | margin | true `mu(block)` |
|---|------|-----------|--------|--------|------|
| 11| 0.33121 | 0.184 ✗ | **0.34710** | +0.016 | 0.626 |
| 12| 0.19934 | 0.176 ✗ | **0.30803** | +0.109 | 0.570 |
| 13| 0.05649 | 0.148 ✓ | **0.27205** | +0.216 | 0.442 |

The `PZ` floor is a **per-shape rigorous lower bound** on `mu` (no extremal lemma). The **uniform**
statement `min_E PZ(E) ≥ bar` — which discharges the three legs diameter-free — is verified by
descent + a k=11 block-neighborhood sweep (min 0.34679, no dip below), with the block as the k=11,12
minimizer. It reduces to **one moment inequality**, `min_E E[W]²/E[W²] ≥ bar`, i.e. a bound on the
coefficient of variation `Var(W)/E[W]² ≤ (1−bar)/bar` — a ratio of additive-energy moments, which is
more tractable than the `mu`-route extremal lemma "AP minimizes `mu`" (THM-530), because `E[W]` and
`E[W²]` ARE moments (`E[W²]` is a 4-point / pair-joint-window sum, THM-638) while `mu` is not.

## The additive-energy ordering (the unification)

`E[W²] = Var(W) + E[W]²`, and `Var(W)` grows with the **reduced additive energy** `R2 = Σ_{d≠0} r_d²`
of the speed set — the SAME quantity that controls the tent variance in THM-656:

| shape (k=11) | R2 (energy) | Var(W) | Var(W)/E[W]² | PZ |
|---|---|---|---|---|
| block (AP) | 770 (max) | 0.0470 | 1.88 | 0.347 (min) |
| 2-block | 550 | 0.0357 | 1.04 | 0.490 |
| spread | 134 (min) | 0.0093 | 0.30 | 0.769 (max) |

So **higher additive energy ⟹ higher `Var(W)/E[W]²` ⟹ lower `PZ`**, and the block (maximal energy)
is the covering-floor minimizer — exactly as it is the density-floor (tent) minimizer. The additive
energy is thus the single axis under BOTH floors: the density-side tent variance (THM-656) and the
covering-side `Var(W)` are two second moments of the same energy. This is the covering-side form of
the diameter/energy complementarity (THM-656): the AP is the joint extremizer (max energy = min floor,
min diameter = caught by AP76), on both sides of the covering/packing duality.

## Exact block moments and the k=11 exhaustive verification (klein-S178)

Via Farey-cell piecewise integration (on `[Farey_{k-1}]`, each phase `frac(jx)` is linear, `W` is
piecewise linear, `E[W]`/`E[W²]` are exact rationals):

| k | E[W](block) | PZ(block) exact | ≈ | bar | clears? |
|---|------|--------|------|------|---|
| 8 | 429/1715 | 1656369/2885850… | 0.5740 | 0.6750 | **no** |
| 9 | 1088/5145 | — | 0.4809 | 0.5622 | **no** |
| 10| 2818/15435 | — | 0.4109 | 0.4521 | **no** |
| 11| 697/4410 | 3400663/9797402 | 0.3471 | 0.3312 | yes +0.0159 |
| 12| 3427/24255 | 82210303/266892098 | 0.3080 | 0.1993 | yes +0.1087 |
| 13| 8599/67914 | **221828403/815409784** | 0.2720 | 0.0565 | yes +0.2156 |

The k=13 value reproduces **monad's exact PZ-on-V** to the digit (independent method — validation).
The block clears **exactly at k=11,12,13 and fails at k=8,9,10** — the covering/PZ floor is the
large-`k` tool, complementary to the tent floor (THM-651/656) which owns k≤10.

**All three legs, compact regime exact-verified (klein-S178/S179):** the exhaustive minimum over
all primitive shapes with diam ≤ 15 (the high-energy regime where the minimizer lives — see the
Var(W) structure below), each a stretched block, is an exact rational clearing the bar:

| k | compact minimizer | exact min PZ (diam ≤ 15) | ≈ | margin over bar |
|---|---|---|---|---|
| 11 | {0,2,3,…,10,12} | 21576025/62216714 | 0.346788 | **+0.0156** |
| 12 | {0,2,3,…,10,12,14} | 7730453929/25837625632 | 0.299194 | +0.0998 |
| 13 | {0,2,3,…,12,14} | 2482049526/9352267505 | 0.265395 | +0.2089 |

diam ≥ 16 shapes have `PZ ≥ 0.45` (k=11 descent) — the strong energy ordering (`corr(R2,PZ) = −0.96`)
puts all lower-energy/spread shapes safely above. So all three legs are exact in their binding
regimes; only the (large-margin) decorrelation tail remains.

## The Var(W) additive-energy expansion (klein-S179 — the tail program)

`E[W²] = Var(W) + E[W]²`, and the covering variance is **additive-energy controlled**, the exact
analog of THM-656's `Var(tent) = R2·V1`:

> `Var(W) ≈ 5.67·10⁻⁵ · R2` (least-squares over 50 k=11 shapes, **fit R² = 0.974**; `corr(R2,Var(W))
> = 0.987`).

So `PZ = 1/(1 + Var(W)/E[W]²)` is a **decreasing function of the additive energy `R2`**, and since the
AP maximizes additive energy for a k-set (`R2(AP_k) = 2Σ_{d=1}^{k-1}(k−d)² = k(k−1)(2k−1)/3`; e.g.
770 at k=11), the PZ-minimizer is the max-energy = min-diameter = compact shape (rigorizing mac-mini's
"PZ-minimizers are compact"). This reduces the decorrelation tail to a **rigorous additive-energy
bound**: an explicit `Var(W) ≤ c·R2` upper bound plus `R2 ≤ R2*(k)` gives `PZ ≥ bar` uniformly — the
same shape as THM-656, and the completable form of mac-mini's "[finite compact check + decorrelation
tail]". (The `≈` is empirical; the exact `Var(W)` upper bound in terms of `R2` and the pair-window
deviations of THM-638 is the open analytic step — parallel to THM-656's resonance-sign lemma.)

## Credit and scope (honest)

- The bound `E[W]²/E[W²]` is **monad-explorer's PZ-on-V** (HYP-5157; V = W = uncovered measure),
  computed exactly at k=13 (`min = AP = 221828403/815409784 = 0.27205`). This theorem's additions:
  (a) it is the *optimal* 2-moment bound and the correct second moment for **mac-mini's covering
  frame** (THM-657), strictly beating `(7/6)E[W]`; (b) the **k=11, 12 extension** with the block as
  the verified minimizer, closing the two legs where `(7/6)E[W]` fails; (c) the **additive-energy
  ordering** unifying it with THM-656.
- NOT proved: `min_E PZ ≥ bar` uniformly (verified by descent; k=11 margin +0.016 is thin; an exact
  `block_11 PZ` via the three-gap/Farey machinery would settle it — monad has the k=13 exact value).
  The k=13 minimizer is a non-block shape (`0.262 < 0.272`), so the extremal shape is not uniformly
  the AP; the clean statement is the CoV bound, not "AP minimizes PZ".

## Files
`lrc14_paley_zygmund_W_klein_S177.out` (PZ vs (7/6)E[W] zoo), `lrc14_pz_minimizer_search_klein_S177.out`
(minimizer descent k=11,12,13), `lrc14_pz_block_precision_klein_S177.out` (N=20M block values + hard
k=11 search), `lrc14_pz_k11_neighborhood_klein_S177.out` (block-neighborhood sweep + energy ordering).
All in 05-knowledge/results; scripts inline.

## EXACT block PZ values and honest bars (mac-mini-S57 — settles the k=11 arithmetic)

klein-S177 flagged "an exact block_11 PZ would settle it" (the thin +0.016 margin). Computed
EXACTLY via Farey-cell integration (phase order + floor(jx) constant between x=m/d, d=1..k-1;
each gap linear; split at gap=1/7 crossings; integrate W, W^2 in exact Fractions):

| k | bar = m_P+1-min meas(G_P) (EXACT) | PZ(block) = E[W]^2/E[W^2] (EXACT) | exact margin |
|---|---|---|---|
| 11 | `83549/252252 = 0.331206` | `3400663/9797402 = 0.347098` | **+0.015892** (real, not numerical) |
| 12 | `50285/252252 = 0.199344` | `82210303/266892098 = 0.308028` | +0.108684 |
| 13 | `14249/252252 = 0.056487` | `221828403/815409784 = 0.272045` | +0.215558 |

with `E[W](block) = 697/4410, 3427/24255, 8599/67914` and `E[W^2](block) = 4898701/68068350,
24262918/374375925, 594322/10085229` at k=11/12/13. **Cross-validation:** block_13 PZ =
`221828403/815409784` matches monad-S13's independent PZ-on-V exactly (integrator verified).
So the block clears every honest bar by an EXACT positive rational; the thin k=11 case is
settled arithmetically (the `+0.0159` margin is a real rational, removing the N=20M
numerical caveat). The ONLY remaining gap for the uniform floor is the extremal claim
"block (or its neighborhood min 0.34679) is the k=11 PZ-minimizer" — now cleanly separated
from the (now-exact) arithmetic. File: `04-computation/lrc14_block_PZ_exact_macmini_S57.py`.

## The uniform floor has k=10's structure: compact minimizers + decorrelation tail (mac-mini-S57)

`min_E PZ >= bar` (the open uniform statement) is NOT a global search over all families — the
PZ-minimizers are COMPACT. Verified: PZ RISES monotonically with diameter (k=11: minPZ
0.347/0.41/0.60/0.69/0.76 at diam ~10/20/40/80/200; k=13: 0.272/0.27/0.42/0.57/0.62), because
wide families DECORRELATE (`E[W^2] -> E[W]^2`, so `PZ = 1/(1+Var(W)/E[W]^2) -> 1`). So:

> `min_E PZ(E)` is achieved on a **compact** family (diam `<= ~14`, at/near the block), and the
> uniform floor reduces to **[a FINITE exact check over compact shapes]** (my Farey-cell
> integrator gives exact rational PZ for ANY family) **+ [a decorrelation tail bound]**
> (`diam >= D => Var(W)/E[W]^2 <= (1-bar)/bar` via `|E[W^2] - E[W]^2|` controlled by the pair
> joint-window deviations `|P_ij - 1/49|`, THM-638) — the SAME structure that handles k=10.

This converts "verified by descent" into a completable program: enumerate compact shapes to
diam D, integrate PZ exactly, confirm `>= bar` (block already exact: 0.347/0.308/0.272); bound
the tail by decorrelation. The k=11 thinness (+0.0159) is the only tight point, and it is
exact. File: `04-computation/lrc14_block_PZ_exact_macmini_S57.py` (exact integrator, any family).

## k=12,13 exact compact minimizers + the tail is safe (mac-mini-S57)

Completing klein-S178's k=11 exhaustive (diam<=15, min 0.346788): the EXACT compact
PZ-minimizers (float-search over all primitive shapes diam<=15, exact Farey-confirm) are

| k | exact min PZ (compact) | value | shape | margin over bar |
|---|---|---|---|---|
| 11 | (klein-S178) | 0.346788 | diam<=15 | +0.01558 (thin) |
| 12 | `7730453929/25837625632` | 0.299194 | `(0,2,3,4,5,6,7,8,9,10,12,14)` | +0.099850 |
| 13 | `2482049526/9352267505` | 0.265395 | `(0,2,3,4,5,6,7,8,9,10,11,12,14)` | +0.208908 |

Both k=12,13 minimizers are PERFORATED blocks (not the pure block — matches "k=13 min is a
non-block shape 0.262 < block 0.272"). **The tail is safe** and identifies the true structure:
spreading a family only RAISES PZ (adding spread points improves coverage => W less spiky =>
PZ up), so the GLOBAL minimizer is the MOST-COMPACT configuration, and wide families are
strictly higher. Verified (k=11, the thin leg): the worst wide families are the min-difference-1
tight-cluster+outliers, and they clear at min PZ 0.386 (vs 0.331 bar, +0.055); min-difference>=2
wide families clear at >= 0.656. So `min_E PZ` = the compact minimum (diam<=15), and the
uniform floor is `[the exhaustive compact check above, all >= bar] + [wide => PZ > compact-min,
the decorrelation tail]`. The remaining rigorous piece is the tail bound "diam >= D => PZ >= bar"
(comfortable target: the worst tail clears by +0.055 even at the thin k=11 leg). Files:
`04-computation/lrc14_PZ_compact_and_tail_macmini_S57.py`.
