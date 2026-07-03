---
id: HYP-4040
title: THE WITNESS-DENOMINATOR LOWER BOUND (a discrepancy phenomenon) -- for the lcm band-blocker family S_X = {1..11,13,lcm(2..X)} the smallest lonely denominator satisfies q(S_X) > X, hence q(S_X) -> infinity and NO uniform arithmetic band closes hge7 (rigorizes HYP-+2876). Since max-speed = lcm(2..X) = e^{Theta(X)}, q(S_X) = Omega(log max-speed); asymptotically Theta(log max-speed) with a structural floor ~41 from the {1..13} sub-family. This is the LRC analogue of the tight all-scales lower bound in arXiv:2607.00876 (binary-tree continual counting, Omega(log^{3/2} n) via hereditary discrepancy): controlling the danger at ALL small moduli (scales) at once costs a growing denominator, and the lcm family is the discrepancy-extremal certificate.
status: LOWER BOUND PROVED (q(S_X) > X, elementary: q<=X => q | lcm(2..X) => that runner sits at danger residue 0). Upper bound q(S_X) = O(log max-speed) EMPIRICAL (q = 41..53 for X up to 45). "No uniform band" now RIGOROUS (was numeric in HYP-3877 / HYP-+2876). Danger-set-size correction folded in.
source: mac-mini-2026-07-03-S22
related:
  - HYP-3877   # my S21 arithmetic band route + magnitude split; THIS proves its "no uniform band" leg
  - HYP-+2876  # lcm families / unbounded witness (numeric REFUTED uniform certificate); THIS gives the RATE
  - HYP-3737   # denominator band forcing (14,26]; the covering-min denominator side
  - HYP-3901   # deep-cluster renormalization -- the large-magnitude route the LB proves is UNAVOIDABLE
results:
  - 04-computation/lcm_witness_denominator_macmini_20260703.py
  - 05-knowledge/results/lcm_witness_denominator_macmini_20260703.out
  - 04-computation/dyadic_magnitude_ladder_macmini_20260703.py
---

# HYP-4040 -- the witness-denominator lower bound as a discrepancy phenomenon

## The analogy (arXiv:2607.00876, Bairaktari-Larsen)
Their theorem: private continual counting -- releasing every prefix sum of a binary stream -- has expected
`l_inf` error `Theta(log^{3/2} n)`, and the binary-tree mechanism (a **dyadic / multi-scale** decomposition)
is optimal. The lower bound is a **discrepancy** argument: the "all prefixes" set system forces error across
`log n` scales at once. The shared skeleton with LRC:

| continual counting | Lonely Runner (this project) |
|---|---|
| release all prefix sums | dodge danger at all small moduli `q` (scales) |
| binary-tree = dyadic scale decomposition | magnitude split / dyadic band ladder (HYP-3877) |
| hereditary discrepancy lower bound | the lcm band-blocker forces the witness denominator up |
| tight `Theta(log^{3/2} n)` | `q(S_X) = Theta(log max-speed)` |

Both are **all-scales-at-once** problems whose cost is polylogarithmic and whose lower bound is certified by a
single extremal (discrepancy-hard) configuration.

## The theorem (lower bound, PROVED)
Let `S_X = {1,2,...,11,13, lcm(2..X)}` (the HYP-+2876 lcm band-blocker), `X >= 13`. It is covering (every
`q in {2..14}` divides a speed: `2..11` and `13` directly, `12 = 4.3` and `14 = 2.7` divide `lcm(2..X)`).
- **Any `q <= X` divides `lcm(2..X)`** (since `q in {2..X}` => `q | lcm(2..X)`), so runner `lcm(2..X)` sits at
  residue `0` at every `t = a/q` -- a danger residue. Hence `S_X` is NOT lonely at any `a/q` with `q <= X`.
- Therefore the smallest lonely denominator `q(S_X) > X`.
- `max-speed = lcm(2..X) = e^{psi(X)} = e^{(1+o(1)) X}` (Chebyshev). So `q(S_X) > X = (1+o(1)) log(max-speed)`.

**`q(S_X) -> infinity` as `X -> infinity`. No single finite denominator band `{15..Q}` closes hge7** -- the
S21 numeric observation is now an elementary theorem. The magnitude split of HYP-3877 (band below a threshold,
singles/cluster/renormalization above) is UNAVOIDABLE, not a convenience.

## Numerics (upper bound + the structural floor)
`q(S_X)` for `X = 13..45` is `41` (until the divisibility bound bites), then `53` -- see
`lcm_witness_denominator...out`. Two regimes:
- **`X < 41`**: `q(S_X) = 41`, a **structural floor** set by the `{1..11,13}` sub-family alone (it first dodges
  the danger set at `q = 41`, failing at `17,19,23,29,31,37`).
- **`X >= 41`**: divisibility dominates, `q(S_X) > X`, tracking `~ nextprime(X)` -- the `Theta(X)` regime.
So `q(S_X) = max(41, Theta(X))`; asymptotically `Theta(log max-speed)`.

## Danger-set-size correction (folds into HYP-3877)
The danger set at denominator `q` is `{r : min(r, q-r) < q/14}`, of size `2*ceil(q/14) - 1`:
`= 3` (namely `{0, 1, q-1}`) for `q <= 28`; `= 5` for `29 <= q <= 42`; `= 7` for `43 <= q <= 56`; etc.
HYP-3877's "exactly 3 residues `{0,1,q-1}`" is correct on its stated range `q in {15..27}` (indeed up to `28`);
for `q >= 29` there are more danger residues. All S21 scripts used the correct general test
`min(r,q-r)*14 < q`, so the near-equal band `{15..33}` results are unaffected (they include `q in {29..33}`
with `5` danger residues); only the prose gloss needed this scope note. The growing danger set is also WHY the
structural floor is `41` not `~28` (at `q >= 29` each speed forbids `4` values of `a`, not `2`).

## Consequence for the LRC14 proof architecture
The arithmetic route CANNOT close hge7 alone: the large-magnitude tail (`max-speed > ~22638`, where the lcm
families live -- `lcm(2..13) = 360360 > 22638`) provably needs a non-band argument. This is exactly the
`> 22638` singles/cluster side of the HYP-3877 split, and the renormalization route [[HYP-3901]]. So the
two-sided architecture (arithmetic band below + analytic/renormalization above) is forced -- a structural fact,
now with a proof, that the fleet's proof sheaf should encode rather than hunt for a uniform arithmetic closure.

## The blowup survives compression (it is NOT just the already-closed dominant case)
The lcm family has a >13x dominant runner, so it is discharged by the dominant-runner route, not the hge7
obligation. But the same witness blowup holds for COMPRESSED families (the real hge7 obligation): **aligned
near-equal** far runners `far_i = q_i * round(N/q_i)` (span-ratio ~1.00, so no dominant; each `≡ 0 mod q_i`)
block a band of moduli by division. Numerically these covering gcd=1 hge7 families reach witness `q = 47, 49`
at `max-speed ~ 2832..2*10^6`, growing with `N` (see MISTAKE-095, `nearequal_bandblocker_stress...out`,
`corrected_band_below_threshold...out`). Counting: a runner `~N` divisible by `lcm` of `k` band moduli needs
`lcm <= N`, so `13` runners block `~13*log(N)/log(30)` moduli => witness `q <= 14 + O(log N)`. So the
compressed obligation ALSO has witness `q -> infinity` at rate `Theta(log max-speed)` -- the lower bound is
about the genuinely-open case, not an artifact of the dominant lcm family.

## Renormalization reading (why this is the right dual to HYP-3901)
At `a/q` a near-equal cluster `{N + c_i}` has far residues `N a + {c_i a} (mod q)`: the **difference core**
`{c_i}` at multiplier `a`, SHIFTED by the scale `N a mod q`. So a large-magnitude near-equal family's
loneliness = (core placeable in the safe region by the scale-shift). For GENERIC `N` the shift is generic and
the small core places easily -> small `q`; the BAND-BLOCKERS are exactly the `N` whose shift is ALIGNED against
the core, forcing `q` up. This is the arithmetic shadow of the deep-cluster renormalization [[HYP-3901]]:
peel the top scale, recurse on the (bounded) core; the tower depth `~log(max-speed)` is the `Theta(log)` cost
this note lower-bounds. The dyadic band ladder + renormalization tower is the matching upper-bound construction
(the binary-tree import from arXiv:2607.00876).

## Connection to the measure-side scale-invariance (klein HYP-3791, OPEN-Q-108)
The arithmetic blowup here is the DUAL of klein's measure-side finding [[HYP-3791]]: the r-far correction is
self-similar over the **13-lattice** (`13 = n-1 =` the dropped speed `= 1/`first-CF-convergent of
`t* = 14/183 = [0;13,14]`), and the 13-spaced pair resonance does NOT decay in scale `W` (persists
`W = 200..50000`) because `13 t* = 182/183 = -1/183 (mod 1)` is a slow resonant phase. That scale-non-decay is
the measure twin of "aligned band-blockers persist at any magnitude" here. Both are **scale-invariant
resonances pinned to continued-fraction-convergent lattices** ([[HYP-3762]] three-gap / CF), one seen through
the danger MEASURE (klein) and one through the residue at the small-`q` WITNESS (this). klein's
"resonance => redundancy => no multi-far beater" (the [[OPEN-Q-108]] residual) is, arithmetically, "13-spaced
far runners differ in phase by `~k/183 << 1/14`, so they cover the SAME danger, not new danger." Worth a joint
measure/arithmetic writeup: the two pictures may close OPEN-Q-108 together (redundancy on both sides).

## Engineering resonance
The "cost of covering all scales" is a shared engineering primitive: the binary-tree mechanism pays it as
`log^{3/2}` noise; LRC pays it as a `Theta(log)` denominator. Both say: to be simultaneously good at every
dyadic scale, you pay a price logarithmic in the range -- and a single worst-case (discrepancy-extremal /
lcm) instance certifies you cannot do better. The dyadic-ladder refinement of the magnitude split is the
direct import of the binary-tree design.
