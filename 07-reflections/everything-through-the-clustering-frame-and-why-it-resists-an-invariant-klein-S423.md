---
source: klein-2026-07-24-S423
status: RE-READING of the project's results through one frame (clustering / spread-vs-structured), as a
  catalogue of instances, PLUS a new negative that explains the frame's limits: no simple numeric invariant
  can capture tightness, because the required invariant must be dilation-INVARIANT and translation-VARIANT,
  and each natural candidate fails on a specific family. Quantitative test included.
tags: [lrc14, clustering-frame, meta, additive-energy, dilation-invariance, re-reading, honest-negative]
---

# Everything through the clustering frame — and why the frame resists an invariant

**klein-2026-07-24-S423.** One frame keeps reproducing across every method the project has tried.

## The frame
> For a family of periodic obstacle sets `D_v` (each of density `2h`), whether they can COVER splits on the
> structure of `{v}`: **spread / dissociated ⇒ the obstacles behave independently ⇒ uncovered measure ≈
> ∏(1−2h) > 0 ⇒ no covering (loose). Clustered / structured ⇒ obstacles resonate ⇒ efficient covering ⇒
> tightness is possible.** In short: **tightness requires clustering.**

## Past work, re-read as instances of the one frame
| result | spread side | clustered side |
|---|---|---|
| **THM-518** (Riesz two-route diagnosis) | Riesz product certifies loose sets | **stalls at 1.0096 on AP-cores** |
| **Bedert's `dim₂²/n³` gain** | large additive dimension ⇒ gain | dies at `dim₂ ≈ 2–3` (AP) |
| **THM-503** almost-Sidon class | `Σg²/(3v_av_b) < 36/49` ⇒ **proved loose** | AP-cores carry abundant `|T|≥3` relations |
| **kps-S138** stranger decoupling | generic strangers decouple to **k=24** | APs/dilates **fail at k=16** |
| **klein-S416** counting-vs-decoupling | uncovered tracks `(1−2h)^13` for generic far parts | counting needed for resonant ones |
| **klein-S422** clustering law | geometrically spread far speeds ⇒ **contradiction** | far part must cluster |
| **THM-504** conditional-convergence wall | within-level cancellation fails | the `|T|≥3` relations are the AP's |
| **tight locus `{AP, GW}`** | everything else is loose | the two most clustered configurations |

Four independent methods hit the *same* wall on the *same* sets. That convergence is the real evidence that
additive/multiplicative structure is the invariant of the hard locus, not an artifact of any one technique.

## The new content: the frame cannot be reduced to a numeric invariant
Measured across 22 configurations (tight, translates, dilates, dissociated, random), correlation with `gap`:

| candidate invariant | corr with gap | fails because |
|---|---|---|
| additive energy `E(S)` | −0.336 | **translation-invariant** — `{20..32}` has the *same* `E = 0.669` as `{1..13}` yet `gap` is **5.38×** tight |
| `Σ 1/v` (origin-sensitive) | **−0.610** (best) | **not dilation-invariant** — `3·{1..13}` is tight but `Σ1/v = 1.06` |
| divisibility pairs | −0.132 | geometric set `{1,2,4,…,4096}` is maximally divisible (6.00) yet loose (4.62×) |
| residue coverage mod 14 | −0.353 | `{20..32}` scores 0.929, same as the AP, yet is the loosest tested |

**The obstruction, stated cleanly.** `gap` is **dilation-invariant** (`gap(cS) = gap(S)`) and **translation-variant**
(`{1..13}` tight, `{20..32}` loose). So any invariant that predicts tightness must be *both*. Additive energy has
the first property but not the second; `Σ1/v` has the second but not the first; divisibility and residue coverage
each break on an explicit family. **No natural candidate has both.** That is precisely why the frame is an
excellent heuristic and a bad theorem — and why OPEN-Q-108 is genuinely Diophantine rather than a
counting statement waiting to be found.

## Reversals worth remembering
- **Φ-consec-extremality is FALSE at k ≥ 11** (repo census): the most-clustered (consecutive) configuration is
  extremal only for `k ≤ 10`; a perforated-dilate family overtakes it after. So "most clustered = most extremal"
  is *not* a law even inside the project's own data.
- **`{20..32}`**: maximal additive clustering, maximal looseness. The single cleanest counterexample to the naive
  frame; worth keeping as the standard test case for any proposed invariant.
- **Geometric `{1,2,4,…}`**: maximal divisibility clustering, loose. Kills the divisibility reading.

## What the frame is good for anyway
1. **Predicting where a method will fail** — any counting/averaging method will stall on AP-cores. Confirmed four
   times; stop re-discovering it.
2. **Choosing the split for a proof** — dissociated part by decoupling, structured part by parameter families
   (klein-S420), which is the only route that covers all `k` at once.
3. **Recognising artifact ceilings** — a `k < 1/(ch)` bound is bookkeeping; a structural obstruction is a family.

→ klein-S415/S416/S418/S419/S422 (lemma, dichotomy, unification, tight-threshold theorems, clustering law),
kps-S138, THM-503/504/518, OPEN-Q-108.
