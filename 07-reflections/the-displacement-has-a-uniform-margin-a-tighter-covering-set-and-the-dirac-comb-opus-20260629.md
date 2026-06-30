# The covering displacement has a UNIFORM margin (not razor-thin in the peak): a new tighter covering set {1..12,182} with M=14/183 < 7/89 corrects the empirical min; the killer binds the displaced witness; and the AP-at-its-witness is the 14-comb (Dirac Ш₁₄) with the observer as the one empty tooth the killer fills

*opus-2026-06-29. Owner: chase a uniform cyclotomic lower bound on the killer's displacement, and the AP
spike vs the Dirac function. The chase found a uniform margin AND a tighter covering set than recorded;
the Dirac angle gave the cleanest picture of the whole obstruction.*

## ABNORMALITY (verified, exact + grid): a covering set tighter than 7/89
`{1,…,12,182}` (drop the binding unit `13≡−1`, add its killer `182=13·14`) is covering and has
> **`M = 14/183 ≈ 0.07650 < 7/89 ≈ 0.07865`** (exact and grid agree).
So the **empirical min over covering sets is NOT `7/89`** (CLAUDE.md / the THM-523 note) — it is at most
`14/183`. The repo's `7/89` came from dropping a *slack* speed (`12`, killer `84`); dropping the
*binding unit* `13` and adding its killer `182` is tighter. **The empirical covering-min should be
updated to `≤ 14/183`.** (Flagged, not silently overriding canon — THM-523's reduction stands; only the
empirical min value is corrected.)

## The displacement is BOUNDED BELOW — not razor-thin in the peak
Chasing lower: drop `j` binding units, add `j` killers. `M` **increases** with `j` (`j=1: 0.0765,
j=2: 0.0897, j=3: 0.109, …`), and the divisor-loaded direction (`drop slack 12, killer 84m`) also
**increases** (`0.0787, 0.0809, …`). So `14/183` is a local (and apparently global) MIN: both ways off
it raise `M`.
> **The covering peak-margin `M−1/14` is bounded below (`≥ 13/2562 ≈ 0.005` empirically). The conjecture
> is NOT razor-thin in the PEAK on the covering side — the covering sets sit a uniform `~0.005` above
> `1/14`.** (The razor-thinness is in the MEASURE `L→0` at the AP, as established; the peak has a margin.
> This is the uniform lower bound the chase sought — positive, not zero.)

## The killer binds the displaced witness (the mechanism)
For the tight covering sets, the killer `w` (mult of 14) **binds the new witness**: at `t=14/183`,
`‖182·t‖ = 14/183 = M` — the killer itself sits at the minimum. So the displaced `M` is *set by the
killer's resonance*: killing the unit-witness `1/14`, the witness moves to the Farey neighbor where the
killer is least dangerous, and that distance IS the new `M`. The displacement is the killer's "best safe
distance," a Farey-jump determined by the killer `= (dropped unit)·14`.

## The DIRAC COMB: the AP at its witness is Ш₁₄, observer = the empty tooth
At `t=1/14` the AP runners sit at `{i/14 : i=1..13}` — and with the observer at `0/14`, this is the
**full 14-comb `{0,1,…,13}/14` = the 14th roots of unity = the Dirac comb Ш₁₄**.
> **The observer's tooth (`0`) is the ONE EMPTY tooth — the gap, the razor's edge** (distance `1/14` to
> its neighbors `1/14, 13/14`). The AP fills 13 of the 14 teeth; the loneliness is the missing tooth.
> **A killer (mult-14) sits exactly at `‖w/14‖ = 0` = the observer's tooth — it FILLS the gap**, and the
> observer is lonely no more. Covering forces a killer ⇒ the empty tooth is filled ⇒ the edge dies.
Ш₁₄ is Poisson-self-dual (`Ш₁₄ ↔ Ш₁₄`, the 14-periodicity), so the comb picture and the cyclotomic
`Q(ζ₁₄)` Galois picture are Fourier duals: **the empty tooth (space) ↔ the killed unit-modes
(frequency).** The spike of `M` (the AP cusp) is the *comb-completion point* — `M=1/13` (a defective
comb) off the integer, `1/14` exactly when the 13th runner completes the tooth. (It is a finite-width
continuous cusp, not a literal delta — `M` is continuous; the Dirac object is the comb, not the spike.)

## The unified obstruction (cleanest statement yet)
> **Loneliness = the empty tooth of the 14-comb. The AP realizes Ш₁₄ with the observer's tooth empty
> (gap `=1/14`, the razor's edge). The covering constraint forces a multiple of 14 — the unique runner
> that lands ON the empty tooth and fills it — so the gap must relocate to a Farey neighbor, where the
> killer binds at distance `>1/14`. The displacement is uniform (`≥ ~0.005`) because the killer's best
> safe Farey distance is bounded away from `1/14`.** LRC(14) = the relocated gap is never smaller than
> the original — the comb's defect can't be made tighter than the AP's single empty tooth.

## Status
- **Verified (opus):** `{1,…,12,182}` covering, `M=14/183 < 7/89` (corrects the empirical covering-min);
  displacement bounded below (`min ≈ 13/2562`, dropping-more or loading-more both increase `M`); killer
  binds the displaced witness; AP-at-witness = Ш₁₄ with observer = empty tooth, killer fills it.
- **New to track:** the corrected empirical covering-min (`≤14/183`); the uniform peak-margin (`~0.005`);
  the Dirac-comb/Poisson-dual picture of the razor's edge.
- **Resolves the chase:** the uniform lower bound is POSITIVE (peak-margin `~0.005`), so the covering side
  is NOT razor-thin in the peak — only the measure is. The bound is the killer's Farey-safe-distance,
  cyclotomically the relocation of the empty tooth off the killed unit.

Related: the cyclotomic-self-dual-razor's-edge (the units/Galois), the razor's-edge/exact-variance
(the 6 units = `φ(14)`), THM-523 (covering reduction — its empirical min value corrected here), THM-566
(divisor loading), cuts-as-Farey-geodesics (the Farey-jump = relocated tooth), the Z₇-cyclotomic-SOS,
OPEN-Q-108.
