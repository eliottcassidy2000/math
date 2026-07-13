# Support for opus's L² large-sieve route: the worst |core|=1 body is the multi-killer minimizer, and the margin is thin

*kind-pasteur-2026-07-11-S127 cont.69. Owner: "support the |core|=1 smooth-body discrepancy route." opus-S267
reduced LRC(14)-covering to a tight L² large-sieve energy bound `core·Σε² < (6/7)² = 0.735`, reporting the
energy `≤ 0.328`. Computing the true worst case pins the binding body, corrects that number to `0.60`, shows
the margin is **thin** (not comfortable), and confirms the crude Bessel bound genuinely fails — so the tight
large-sieve is necessary, and its target is now explicit.*

---

## opus's route (S267), in one line

LRC(14)-covering ⟸ `Σ|ε_v| < 6/7` (⟹ coreCover < 1). By Cauchy–Schwarz `Σ|ε_v| ≤ √(core·Σ_v ε_v²)`, so it
suffices to bound the **L² large-sieve energy** `core·Σε² < (6/7)² = 0.735`. opus has a rigorous Bessel bound
`core·Σε² ≤ (6/49)·core/|G'|` that is loose by a stated `3.1×`; a tight large-sieve estimate closes it. opus
reported the actual energy `≤ 0.328` (a comfortable `0.41` below the ceiling).

## The correction: the worst energy is 0.60, at the multi-killer minimizer

Searching `|core|=1` covering bodies, the maximum energy is **not** 0.328. It is:

> **`core·Σε² = 0.60`, at the body `{2,…,11, 13, 84}` — i.e. the family `{1,…,11, 13, 84}`, the multi-killer
> covering-min minimizer (`M = 7/89`).** There `coreCover = 0.92` (exactly mac-mini-S74's runner-1 figure, so
> that figure *is* `coreCover = |D₁∩G'|/|G'|`), `ε₁ = 0.92 − 1/7 = 0.777`, `|G'| = 0.0666`.

So the true margin to opus's ceiling is **thin, not comfortable**:

| quantity | worst-case value | ceiling | margin |
|---|---|---|---|
| `core·Σε²` (energy) | **0.60** | 0.735 | 0.13 |
| `Σ|ε_v| = √(energy)` | **0.777** | 6/7 = 0.857 | **0.08** |

opus's `0.328` was a non-worst sample; the tight large-sieve must clear `0.60 < 0.735` — with only `0.13` of
room — not `0.328` with `0.41`.

## The crude Bessel genuinely fails — the tight bound is necessary

At the worst body the crude Bessel bound is

> `(6/49)/|G'| = (6/49)/0.0666 = 1.84 ≫ 0.735` — **it fails**, and is `1.84/0.60 = 3.05×` loose (matching
> opus's stated `3.1×`).

So the tight large-sieve must gain a factor `1.84/0.735 ≈ 2.5×` over the crude Bessel to clear the ceiling at
this body. The route is viable — the *actual* energy `0.60 < 0.735`, so LRC holds with `Σ|ε| = 0.777 < 6/7` —
but there is no slack to spare, and the crude Bessel alone is genuinely insufficient.

## The key structural input: worst-discrepancy ≠ worst-loneliness

The most useful fact for opus's targeting:

> **The worst-discrepancy body is the multi-killer minimizer `{1..11,13,84}` (`coreCover = 0.92`, `M = 7/89`),
> NOT the deep well `{1..12,182}` (`coreCover = 0.72`, `M = 14/183`).**

The two extremal problems have *different* extremizers. The deep well minimizes *loneliness* (`M = 14/183`) but
is comfortably clear in *discrepancy* (`ε₁² = 0.33`); the multi-killer body is loose in loneliness (`7/89`) but
**maximizes the discrepancy energy** (`0.60`). So opus's tight large-sieve should be calibrated against
`{1..11,13,84}`, whose non-core body `{2,…,11,13,84}` has the small good set `|G'| = 0.067` and the
concentrated single arc — not against the deep well. Intuitively: the shorter core (`{2,…,11}` plus the split
`13, 84`) leaves a smaller, more arc-concentrated good set, pushing `coreCover` toward 1.

## Net — the target, sharpened

opus's L² route is sound and the reduction is rigorous up to the Bessel constant, but the true binding case is
`core·Σε² = 0.60` (margin `0.13`) at the multi-killer minimizer `{1..11,13,84}`, and the crude Bessel (`1.84`)
fails there by `2.5×`. The tight large-sieve estimate opus needs must (a) target that body's good set (`|G'| =
0.067`, single concentrated arc), and (b) gain `2.5×` over Bessel while the actual energy sits only `0.13`
below the ceiling. This is the concrete, quantified target for the last open analytic step — the worst case,
its margin, and the precision required, now pinned.

*Files: lrc14_L2energy_worstcase_kps_S127.py (+.out). Supports opus-S267 (L² large-sieve energy), S262
(near-orthogonality Cov ≤ 1/(3vv')); uses kps cont.62/63 (|core|=1 coreCover), cont.58 (multi-killer
minimizer). HYP-6246.*
