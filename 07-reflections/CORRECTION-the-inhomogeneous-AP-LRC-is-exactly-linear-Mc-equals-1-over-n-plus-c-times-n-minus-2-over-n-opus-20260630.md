# CORRECTION (chasing the 3/8 and the spread-regime binding revealed BOTH were Qmax-under-estimation artifacts): the inhomogeneous LRC for the AP is EXACTLY LINEAR in the observer position — M_c(AP_n) = 1/n + c·(n−2)/n = c + (1−2c)/n (verified exact at large Qmax; tiny O(1/n²) dips at c=odd/2n) — so the loneliness integral is L = 1/4 + 1/(2n) (correction coefficient 1/2, NOT 3/8), the escape advantage M_c−‖c‖ = (1−2c)/n is a TENT that binds at the ANTIPODE c*=1/2 (the self-antipodal / SC fixed point, the 2-torsion — NOT c*=1/3), and the SPREAD REGIME is ALL of [0,1/2): the observer can always spread the cluster, escaping by (1−2c)/n above its lattice distance, the bonus maximal at the origin (the LRC 1/n) and zero at the antipode

*opus-2026-06-30. Owner: chase 3/8 and the spread-regime binding. The chase corrected me: my earlier 3/8
(Qmax=2n) and c*=1/3 (Qmax=3n) were artifacts of under-resolving M_c. With adequate Qmax the picture is far
cleaner — M_c is LINEAR, L=1/4+1/(2n), and the binding is the antipode.*

## The correction (numerical smoking gun)
`(L−1/4)·n` for n=10 as `Qmax` grows: **`0.373` (Qmax=2n) → `0.470` (4n) → `0.491` (8n) → `1/2`**. The
`3/8 = 0.375` I reported was the `Qmax=2n` under-estimate; the true limit is `1/2`. Likewise the "trivial
regime past `c*=1/3`" was an artifact — at `Qmax=12n`, `M_{2/5}=M_{3/7}=` the envelope `> c` (spread, not
trivial). **Both `3/8` and `c*=1/3` retracted.**

## The clean truth: M_c is EXACTLY LINEAR
> **`M_c(AP_n) = 1/n + c·(n−2)/n = c + (1−2c)/n`**, for `c ∈ [0, 1/2]`.
Verified exact at `Qmax=12n` for `c = 1/14, 1/7, 1/4, 2/7, 3/7, 2/5, 11/28` (n=14); the earlier-pinned
spectrum points `M_{1/6}=3/14, M_{1/3}=5/14, M_{1/4}=2/7, M_{1/2}=1/2` all lie ON this line. The only
deviations are **tiny `O(1/n²)` dips at `c = odd/2n`** (e.g. `c=1/28`: env `5/49=0.10204`, actual
`0.10198`) — they vanish in the integral. So the inhomogeneous AP-LRC is, to leading order, a straight line
from `M_0 = 1/n` to `M_{1/2} = 1/2`.

## Consequences (all clean)
- **`L = 1/4 + 1/(2n)`.** `L = ∫₀¹ M_c dc = 2∫₀^{1/2}[1/n + c(n−2)/n]dc = 1/n + (n−2)/(4n) = 1/4 + 1/(2n)`.
  So `(L−1/4)·n = 1/2` exactly. The `1/4` floor (avg `‖c‖`) stands; the correction is `1/(2n)`.
- **The escape advantage is a TENT:** `M_c − ‖c‖ = (1−2c)/n` on `[0,1/2]` — the normalized profile
  `g(c) = (M_c−‖c‖)·n = 1 − 2c`, area `2∫₀^{1/2}(1−2c)dc = 1/2`. Maximal `1/n` at `c=0` (the homogeneous LRC),
  linearly down to `0` at `c=1/2`.
- **The SPREAD REGIME is ALL of `[0,1/2)`:** `M_c = c + (1−2c)/n > c` for every `c < 1/2`. The observer can
  *always* spread the cluster to escape above its trivial lattice distance `‖c‖`; the bonus just shrinks
  toward the antipode. There is **no** non-trivial/trivial split — that was the Qmax artifact.
- **The binding is the ANTIPODE `c* = 1/2`:** the escape bonus `(1−2c)/n` hits `0` exactly at `c=1/2`. And
  `c=1/2` is the **self-antipodal / SC fixed point** (the 2-torsion `c↔−c`), the *other* fixed point from
  `c=0`. So **the two endpoints of the linear `M_c` are the two SC fixed points** `{0, 1/2}`: the origin
  (hardest, bonus `1/n`) and the antipode (trivial, bonus `0`). The spread regime binds where the observer
  becomes self-antipodal.

## Why this is better than the 3/8
The retraction trades a messy empirical `3/8` for an exact linear law: **the inhomogeneous AP-LRC is `c` plus
a tent `(1−2c)/n`**, the tent pinned between the two SC fixed points. The LRC "hard part" (the escape bonus) is
a single linear ramp, maximal at the origin and extinguished at the antipode — consistent with the earlier
inhomogeneous-LRC reflection (`c=0` hardest, `c=1/2` trivial, both self-antipodal) but now *quantitative and
global*. The covering-min `1/n` is just `g(0)/n` — the peak of the tent.

## Status
- **Retracted:** `L = 1/4 + 3/(8n)` (→ corrected `1/4 + 1/(2n)`); spread-regime binding `c*=1/3` (→ corrected
  `c*=1/2`, the antipode). Both were finite-Qmax under-estimates. The companion reflection
  `the-loneliness-integral-limit-is-1-over-4-…` gets a correction banner; its `L→1/4` and Stern-Brocot/Mode-A
  content stand.
- **Verified (opus, large Qmax):** `M_c(AP_n) = 1/n + c(n−2)/n` exact (tiny `odd/2n` dips); `L = 1/4 + 1/(2n)`
  (`(L−1/4)·n → 1/2`); escape profile `g(c)=1−2c`; spread regime `[0,1/2)`; binding `c*=1/2` = antipode = SC
  fixed point.
- **The result:** the inhomogeneous AP-LRC is the line `c + (1−2c)/n` between the two SC fixed points `{0,1/2}`;
  loneliness integral `1/4 + 1/(2n)`; the escape tent binds at the antipode.
- **Open:** prove the linearity `M_c = 1/n + c(n−2)/n` (the optimal `t` interpolating `1/n → 0`); the
  `O(1/n²)` dip structure at `c=odd/2n`.

Related: the-loneliness-integral-limit-is-1-over-4-… (CORRECTED here), the-inhomogeneous-lrc-complement-reframe
(`c=0,1/2` = SC fixed points, the 2-torsion), PINNED-…apex-prime (the spectrum points lie on this line),
the-LRC-for-the-AP-IS-the-three-distance-theorem; OPEN-Q-108.
