# The CV gatekeeper is miscalibrated on covering families — the deep well passes, the real floor-deficit family passes too, and speed 7 is where the floor bottoms

*opus-2026-07-03-S58. The owner asked me to lower-bound the tight-family measure floor. The fleet's
crisp open target is THM-579's covering-floor CV criterion, with klein's HYP-3554 leaving one sharp
open sub-question (its next-step #2). I answered it — exactly — and in doing so found the CV
gatekeeper is not merely non-uniform (known) but **miscalibrated**: it fails on families whose actual
floor is fine and passes on the family whose actual floor is worst. Honest, and it repositions the
route, not closes it.*

## The target, precisely

THM-579 (mac-mini, PROVED) separates the covering floor: for covering `S = R ∪ 14Q`,
`R'(S) = meas(lonely S)/(m_R m_Q) ≥ 1 − CV(N_R)·√((1−m_Q)/m_Q)`, so the floor is positive whenever the
**gatekeeper** `CV(N_R)² < m_Q/(1−m_Q)` holds, where `N_R(t) = #{a∈0..13 : t+a/14 is R-safe}` is the
14-sheet count. klein (HYP-3554) scanned 1828 sets: `CV(N_R)²` is **unbounded** (sup 8.74 at
`R={1..13}\{12}`, `m_R→0`, amplified by the speed-7 resonance `7·a/14 = a/2`), so a *uniform* CV bound
is false — but the actual `R'` stayed positive everywhere. klein's **next-step #2**: *do the
gatekeeper-failing dense `R` actually occur as the `R`-part of a genuine 13-speed covering family?*

## The r=1 simplification, and the deep well

For `r=1` (a single far runner `14q`), `m_Q = meas(lonely{q}) = 6/7` **always**, so the threshold is
exactly `m_Q/(1−m_Q) = 6`. The gatekeeper becomes simply **`CV(N_R)² < 6`**. And the covering-min
extremal — the **deep well `{1..12,182}`** — *is* an `r=1` config: `R={1..12}`, `Q={13}`,
`14·13 = 182`. mac-mini's verification table ran `r=2..6` and skipped it. Filling it:
`CV²(R={1..12}) = 1.095 < 6` — **the deep well PASSES the gatekeeper**, with CS bound `+0.573` and
actual `R' = 0.818`. The tightest known covering family is *not* where the CV criterion struggles.
(Exact `CV²` via `Var(N_R) = 14·Σ_{d=0}^{13} A(d/14) − (14 m_R)²`, `A` = the R-safe autocorrelation.)

## klein #2, answered: YES — and the CV proxy is miscalibrated

Testing the dense `R = {1..13}\{j}` and completing each to a covering 13-family `R ∪ {14q}`:

| `R` | `CV²` | `<6`? | covering `S` | `meas(S)` | actual `R'` |
|---|---|---|---|---|---|
| `{1..12}` (deep well) | 1.095 | ✓ | `{1..12,182}` | 0.0239 | 0.818 |
| `{1..13}\{7}` | 3.63 | ✓ | `{1..6,8..13,14}` | 0.0227 | **0.315** |
| `{1..13}\{6}` | 7.76 | ✗ | `{1..5,7..13,14}` | 0.0082 | 1.167 |
| `{1..13}\{12}` | 8.74 | ✗ | `{1..11,13,84}` | 0.0054 | 0.514 |

Two facts fall out:

1. **YES to klein #2.** `R={1..13}\{6}` and `R={1..13}\{12}` — both gatekeeper-*failing* (`CV²>6`) —
   **do** complete to genuine 13-speed covering families (`{1..5,7..13,14}`, `{1..11,13,84}`). So the CV
   gatekeeper is not uniform *even on the real covering family*; those families need the exact SPEC /
   `Γ₀(N)` route (HYP-3553), not the CV proxy. And their actual `R'` stays `>0` (1.167, 0.514) — klein's
   robustness, now confirmed in exact arithmetic.

2. **The gatekeeper is miscalibrated.** The covering family with the **lowest actual floor**,
   `R' = 0.315`, is `{1..13}\{7}` (drop speed 7) — and it **passes** the gatekeeper (`CV²=3.63<6`).
   Meanwhile the gatekeeper-failing families have *higher* `R'` (0.514, 1.167). So `CV²>6` does **not**
   track the true floor deficit: it flags the wrong families. The exact-SPEC / `Γ₀(N)` route is not
   merely "more uniform" than CV — it tracks a **different, correct** quantity. A future uniform floor
   proof cannot be a repair of the CV bound; it must replace the proxy.

## Speed 7 is where the floor bottoms

The lowest floor drops when speed **7** is removed. This is the `14 = 2·7` signature from the other
side: the band `1/14 = 1/(2·7)` has its Fourier zeros at `7ℤ` (`âhat(k)=−sin(πk/7)/(πk)`, my S57
7-Fourier-zeros / THM-579's `ahat`), speed 7 folds the 14-sheets onto `a/2` (klein's `CV`-amplifier),
and dropping 7 from the covering family is what minimizes `R'`. The three roles of 7 — Fourier-zero of
the band, sheet-folder in `Var(N_R)`, floor-minimizer in `R'` — are one phenomenon. The apex-7 face of
14 is exactly where the covering floor is thinnest.

## Three axes, and where the crux actually is

On covering families, three "tightness" notions that were tacitly conflated are **distinct**:
(i) small safe-measure `meas(S)` (deep well 0.024; `{1..11,13,84}` 0.005 — *smaller*);
(ii) CV-gatekeeper failure `CV²>6` (`\{6},\{12}`);
(iii) low actual floor `R'` (`\{7}`, 0.315).
They do not coincide. "The tight-family measure floor" is not a single family — it forks by which axis
you minimize. The genuinely binding object is **(iii)**, the actual `R'`, and its minimizer (drop-7) is
invisible to both the CV proxy (ii) and the naive tightness (i).

## Status

- **Verified (exact):** deep well `r=1` passes (`CV²=1.09<6`, `R'=0.818`); two gatekeeper-failing dense
  `R` complete to covering families with `R'>0`; drop-7 gives the min floor `R'=0.315` yet passes the
  gatekeeper. (`lrc14_floor_CV_r1_deepwell_opus_S58.py`, `lrc14_floor_dense_R_covering_crux_opus_S58.py`.)
- **Answered:** klein HYP-3554 next-step #2 (YES, gatekeeper-failing `R` occur on the real family) — the
  CV route cannot be made uniform even restricted to coverings; and it is *miscalibrated*, so the
  `Γ₀(N)`/exact-SPEC route (HYP-3553, THM-580) is the necessary one.
- **Not a proof / open (= the tight-family measure floor):** a *uniform* `R'>0` over covering families.
  This session maps where it is hard (drop-7, the apex-7 face) and rules out the CV proxy as the tool;
  it does not lower-bound `R'` uniformly.

Note the alignment with kps HYP-4060 (concurrent): the deep well `{1..12,182}` is the **unique
primitive covering-min extremizer** (`M=14/183`) — and it is exactly the `r=1` family that *passes* the
CV floor criterion here. So the `M`-extremal family is benign for the floor; the floor-deficit minimizer
(drop-7) is a *different* family. kps's covering-min `M` is a fourth axis, and its extremizer does not
coincide with the floor-deficit one either. (HYP-4060 also settles the primitivity subtlety — imprimitive
`14·{1..13}` is covering+tight but gcd-reduces to non-covering — so all of this is a gcd=1 statement.)

Related: THM-579 (the CV criterion), HYP-3554/klein (CV unbounded + its #2, here answered), HYP-3553 +
THM-580 (the `Γ₀(N)`/congruence route, now confirmed as the necessary one), THM-576 (the `m_Q` cap),
opus HYP-4058/S57 (the 7-Fourier-zeros of the band = the same `ahat`), OPEN-Q-108 (near-tight cores),
MISTAKE-078 (the absolute-bound divergence). Scripts: `lrc14_floor_CV_r1_deepwell_opus_S58.py`,
`lrc14_floor_dense_R_covering_crux_opus_S58.py`. HYP-4061.
