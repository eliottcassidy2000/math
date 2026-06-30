# The covering-min C(n) across all n: ONE transition at n=7 (drop-2 anomaly 2/(2n−1) for n≤6 → construction n/Φ₆(n) for n≥7, verified exact n=4..20); the margin C(n)−1/n = (n−1)/(n·Φ₆(n)) ~ 1/n² → 0 (the razor-thin covering margin, vanishing; n·C(n)→1⁻); the transition sits exactly at the first PG(2,n−1) FAILURE (n−1=6, Bruck-Ryser) and the apex prime 7; Φ₆(n) is built from Eisenstein-split primes (≡1 mod 6, or 3), with the apex abnormalities Φ₆(3)=7 and Φ₆(19)=7³

*opus-2026-06-30. Owner: push the optimality proof by understanding ALL transitions/abnormalities/structures
with nuance in n — see the bigger picture. Pulled klein HYP-3705/3706 + mac-mini HYP-3701; computed the full
C(n) landscape. The covering-min is one clean law with one transition and a vanishing margin.*

## The C(n) law (computed exact, n=4..20)
`C(n)` = covering-min of LRC(n) = min over covering (n−1)-sets of `M(S)`:
| regime | C(n) | family | n |
|---|---|---|---|
| small | `2/(2n−1)` | **drop-2** `{1,3,4,…} + tuned large` | n ≤ 6 |
| generic | `n/Φ₆(n)` | **construction** `{1,…,n−2,(n−1)n}` | n ≥ 7 |
> Verified: `C(4..6) = 2/7, 2/9, 2/11` (drop-2); `C(7..20) = 7/43, 8/57, 9/73, 10/91, 11/111, 12/133,
> 13/157, 14/183, 15/211, 16/241, 17/273, 18/307, 19/343, 20/381` — all `= n/Φ₆(n)` exactly (the construction,
> ζ₆-witness `M_exact` confirmed). **Exactly ONE transition, at n=7.**

## The margin: razor-thin and vanishing (`~1/n²`)
For `n≥7`, `C(n) = n/Φ₆(n)`, so
> `C(n) − 1/n = n/(n²−n+1) − 1/n = (n−1)/(n·Φ₆(n)) ~ 1/n² → 0`,
> equivalently `n·C(n) = n²/(n²−n+1) → 1⁻`.
**The covering-min approaches `1/n` from ABOVE, with a margin shrinking like `1/n²`.** At n=14 the margin is
`13/(14·183) = 0.00507`. The conjecture `C(n) ≥ 1/n` holds with this razor-thin, vanishing margin — the
"razor-thin" of LRC is precisely this `1/n²` covering margin (distinct from the apex measure-vanishing,
mac-mini HYP-3700). This is *why* the problem is hard at every n: the worst covering set sits just `~1/n²`
above the floor.

## The drop-2 anomaly is a FINITE-SIZE effect (washes out)
`2/(2n−1) < n/Φ₆(n)` for ALL `n>2` (since `n(2n−1)−2Φ₆ = n−2 > 0`), so drop-2 *would* always win **if it
scaled**. It doesn't: it achieves `2/(2n−1)` only while the interval `{3,…,n−1}` is short (n≤6); at n≥7 the
interval is too long and drop-2 degrades (n=7: `0.182 > 7/43`; n=14: `2/17, 9/83 ≫ 14/183`). So:
> **drop-2 is a small-n / finite-size anomaly; the construction (ζ₆-rotated AP, the hexagonal packing) is
> the scale-robust generic family.** The transition is the combinatorial threshold where the finite-size
> trick stops beating the asymptotic packing.

## The structural coincidences at n=7 (honest nuance)
The transition `n=7` coincides with three things — flagged, not claimed causal:
- **First PG(2,n−1) FAILURE:** `PG(2,m)` exists iff `m` is a prime power; `m=n−1=6` is the first
  non-prime-power, and `PG(2,6)` fails by **Bruck–Ryser**. For `n≤6` (`n−1=3,4,5` prime powers) a projective
  plane exists — the regime where the drop-2 *finite design* still wins; at the first PG failure the
  construction takes over. (Suggestive that the drop-2 anomaly lives where small projective planes exist;
  causality unproven.)
- **The apex prime 7.** Caveat: this transition is in `C(n)` for LRC(**7**); the apex-7 of the project is
  `14=2·7`. The numerical coincidence (n=7) may be shallow — the *deep* apex-7 is LRC(14)'s genus, a
  different role. Don't conflate.
- For the OPTIMALITY PROOF (klein's design bridge), the clean cases are PG-existence n's (`n−1` prime power):
  `n = 4,5,6,8,9,10,12,14,17,18,20`. **n=14 (n−1=13 prime, PG(2,13) exists) is a clean case** — klein's
  Steiner↔Kershner bridge applies directly. PG-failure n's (`7,11,13,15,16,19`) keep the construction value
  `n/Φ₆` but need a non-design optimality argument.

## Φ₆(n) abnormalities (the Eisenstein structure)
`Φ₆(n)=n²−n+1` is always odd; `3 | Φ₆(n) ⇔ n ≡ 2 (mod 3)`; every other prime factor is `≡ 1 (mod 6)` — the
**Eisenstein-split primes** (klein's Q(√−3) column). Factorizations: `13, 3·7, 31, 43, 3·19, 73, 7·13,
3·37, 7·19, 157, 3·61, 211, 241, 3·7·13, 307, 7³, 3·127`. Apex-7 abnormalities:
> **`Φ₆(3)=7`** (the apex prime itself — the "two Heegner fields meet at 7=Φ₆(3)" point) and **`Φ₆(19)=7³=343`**
> (the apex CUBED — echoing 14a's discriminant `−2⁶·7³`). The apex 7 recurs as a covering-denominator at
> `n=3,10,12,17,19` (every `Φ₆` divisible by 7).

## The bigger picture for the optimality proof
- **C(n) is one law with one transition.** The covering-min is `n/Φ₆(n)` for all `n≥7` (the hexagonal/ζ₆
  construction), with a vanishing `~1/n²` margin above `1/n`. The drop-2 anomaly is a finite-size artifact
  below the transition.
- **The proof target is uniform:** show the ζ₆-rotated AP (step-`n`) is the densest covering image for all
  `n≥7` (the continuous↔discrete hexagonal optimality, klein's Kershner bridge). The PG-existence cases
  (incl. n=14) get it from Steiner-design optimality; the PG-failure cases need the continuous (Kershner)
  side directly.
- **The transition at n=7 localizes the "where it gets hard":** below it, finite designs (drop-2/PG) win;
  at and above it, the asymptotic hexagonal packing is extremal. LRC(14) is firmly in the generic regime,
  with PG(2,13) existing — the best-case scenario for the bridge.

## Status
- **Computed/verified (opus):** `C(n)=2/(2n−1)` (n≤6) / `n/Φ₆(n)` (n≥7), exact n=4..20; one transition at
  n=7; margin `(n−1)/(nΦ₆)~1/n²→0`, `n·C(n)→1⁻`; drop-2 is finite-size; Φ₆ Eisenstein-split factors;
  `Φ₆(3)=7`, `Φ₆(19)=7³`; PG-existence map (n=14 clean).
- **Nuance (honest):** the n=7 transition coincides with the first PG failure (Bruck–Ryser) and the apex
  prime — flagged as suggestive, not claimed causal; the C(7) transition is LRC(7)'s, distinct from LRC(14)'s
  apex-7 genus.
- **Open (the optimality proof, klein's):** the ζ₆-AP is the densest covering image for n≥7 (hexagonal
  Kershner optimality); clean via Steiner when PG exists (n=14), via Kershner when it doesn't.

Related: the-covering-min-witness-is-kleins-zeta6 (the n=14 witness), klein HYP-3705/3706 (design/hexagonal),
mac-mini HYP-3701 (n-dependence — this maps it fully) + HYP-3700 (apex isolation), the-master-two-Heegner-
columns (Q(√−3) covering column), covering-min-Eisenstein-Φ₆; OPEN-Q-108.
