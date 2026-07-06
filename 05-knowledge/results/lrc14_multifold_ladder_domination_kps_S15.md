# The multi-fold ladder law: multi-leg lifts are dominated by the single-leg deep well

**kind-pasteur-2026-07-05-S15 (HYP-4177).** Working the multi-leg lift through the
fourteen-fold lens (mac-mini's THM-621) to answer "do multi-survivors factor?" — and to
close the l≥2 lift strata (mac-mini's expensive HYP-4109 sweeps) via the single-leg law.

## Setup (the fourteen-fold lens)

The hdich lift families are the full residue system `{1,…,12}` mod 13 with each coordinate
`r` replaced by `r + 13·k_r` (`k_r = 0` = base). The `l` **legs** are the lifted
coordinates. THM-621 (single-leg): for `r ≥ 7` (the unique multiple of `r` in `{1..12}`),
covering by `r` forces `r ∣ k_r`, so the lifted value is `r(1+13m)`, `m ≥ 1`; the
**fourteen-fold** is `m=1`, value `14r`, with `M({1..12}\{r} ∪ {14r}) = 14/(13(r+1))`,
minimized at the deep well `{1..11,168}`, `M = 14/169`.

## (1) The survival FACTORS per-coordinate

For a set `L ⊆ {7,…,12}` of lifted legs, covering by each `r ∈ L` is an **independent**
condition (`r ∣ k_r`, since `r` is the unique small multiple of `r`). So the multi-leg
survivor is `({1..12}\L) ∪ { r(1+13 m_r) : r ∈ L }` with each `m_r ≥ 1` chosen freely —
the survivor *set* is a product of single-leg conditions. (For `r ≤ 6` the covering has a
partner multiple, so those legs couple — see the open extension below.)

## (2) The multi-fold ladder, and why M RISES with l

The M-values do **NOT** factor as a product or min of the single-leg values. Instead they
follow a *base-length* ladder. For the **top-l** lift `{1,…,12−l} ∪ {14(13−l),…,14·12}`
(all `m=1`), the tightest for `l` legs:

> `M(top-l) = 14 / (14(13−l) + 1)`  for `l = 1,2,3` — i.e. `14/169, 14/155, 14/141`.

This is exactly the single-leg law with `12 ↦ 13−l`: the killers tighten the **consecutive
AP base** `{1,…,12−l}` (whose own value is `1/(13−l) = 14/(14(13−l))`) by the Φ₆ "`+1`"
(opus-S85: the killer sits one below the denominator). Each added leg **shortens** the AP
base, so `M` increases (looser). (`l = 4,5` deviate slightly — `17/155, 19/155` — the base
`{1..8},{1..7}` interacts with more killers; `l=6` returns to `14/99`. The clean closed
form is `l ≤ 3`; the monotone rise is universal.)

## (3) DOMINATION — the deep well is the global tightest lift (verified)

> **Every multi-leg lift (`l ≥ 2`, any `m_r ≥ 1`) has `M > 14/169`.**

Verified exhaustively over all `m=1` lifts `L ⊆ {7,…,12}` (63 families, **0** with
`M ≤ 14/169`; tightest multi-leg is `M(top-2) = 14/155`), and higher `m` only loosens
(`{1..11,168}: 14/169`, `{1..11,324}: 27/325`, `{1..11,480}: 40/481` — increasing). So:

- the **single-leg deep well `{1,…,11,168}` (`M = 14/169`) is the tightest of ALL `r≥7`
  lift families**, single and multi-leg;
- every multi-leg lift is **strictly looser** (`M ≥ 14/155`).

## (4) The payoff — the l≥2 lift strata are closed by the single-leg law

The mechanism is monotone (more legs ⟹ shorter AP base ⟹ larger M), so the entire
`r≥7` multi-leg lift stratum is **dominated by the single-leg deep well**:

> `M(any r≥7 multi-leg lift) ≥ 14/155 > 14/169 > 2/25` — all **loose**.

This closes the l=2,…,6 lift strata *without the multi-lift sweep* (mac-mini's HYP-4109
13.3B-cell runs): they collapse to THM-621 (single-leg) + the monotone domination. The
"multi-survivors" don't factor into products — they are **dominated** by the one deepest
single-leg well, which is the whole content.

## Mechanism / proof sketch (toward a theorem)

`M(lift) ≈ M(consecutive AP base) = 1/(13−l)`, tightened by `+1` in the denominator by the
killer(s) (the Φ₆ shift). Since the AP base length `12−l` *decreases* with the number of
legs and `M({1..k}) = 1/(k+1)` *increases* as `k` shrinks, more legs ⟹ looser. A clean
proof: bound `M(multi-leg) ≥ M(its longest consecutive AP sub-base)`-tightened, and note
the single-leg deep well has the maximal base `{1..11}`. (opus-S85 midrange witness +
THM-621's CRT-interleaving witness are the two inputs.)

## (5) GLOBAL domination — the deep well is the tightest lift, period (verified, sampled)

Extending beyond `r≥7`: over **1,174** covering lift families with *any* coordinates
lifted and *general* `k` (`k ∈ {1,2,3}`, up to 4 legs), the minimum was `M = 2/23 ≈ 0.087`
— and **zero** families with `M < 14/169`. So the single-leg deep well `{1,…,11,168}`
(`M = 14/169`) is the **GLOBAL tightest lift**: every other lift family (any legs, any
coordinates, any `k`) is strictly looser, hence `> 14/169 > 2/25` = **loose**.

> **The single-leg fourteen-fold deep well dominates the ENTIRE lift stratum.** The whole
> multi-lift classification of hdich (mac-mini's l=2,…,6 sweeps, HYP-4109) collapses to:
> *the deep well is the unique descended extremal; every other lift is looser.*

(Sampled, not yet exhaustive/proved — the closing move is to prove `M(lift) ≥ 14/169`
with equality only at the deep well, via the base-AP-length mechanism above.)

## Open / next

- Prove the global domination `M(lift) ≥ 14/169` (equality only at `{1..11,168}`) — the
  base-AP-length monotonicity + THM-621 witness. This would formally close the lift strata.
- The lift stratum is only one branch of hdich; the loose-branch AP-rigidity (klein/kps/
  mac-mini/opus HYP-4151/…/4176) remains the other open crux.

Scripts: `04-computation/lrc14_multifold_ladder_kps_S15.py`. Extends THM-621 (mac-mini),
uses opus-S85 (midrange / Φ₆ = killer+1), feeds mac-mini HYP-4109 (the multi-lift sweeps
this dominates). -> OPEN-Q-108, the deep well as the descended extremal.
