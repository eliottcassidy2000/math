# THM-619: The alignment-band confinement criterion — the geometry of hcomp's loose case

**Status:** PROVED (the criterion, parts i–iii) + VERIFIED-EMPTY on four loose bases (part iv)
**Author:** mac-mini-2026-07-04-S48 (HYP-4094)
**Verification:** `04-computation/lrc_alignment_band_confinement_macmini_S48.py` (+ `.out`).
**Context:** the sole open LRC(14) leaf is hcomp, re-scoped by klein-S131 to PRIMITIVE compressed covering ⟹ M ≥ 1/13; the tight-AP base case is done (CRT free-rider, S47); THIS theorem is the geometry of the remaining loose-base case (M(B) ≥ 2/25).

## The criterion

Let `V = B ∪ {w}` be a primitive compressed covering 13-set, base `B` (12 runners), killer `w ≤ 13·max(B)`. Then `M(V) < 1/13` REQUIRES all of:

**(i) [tooth-disjointness ⟹ one-tooth containment]** The killer's radius-1/13 comb has disjoint teeth (width `2/(13w)`, gaps `11/(13w)`); it fully covers an interval `J` iff `J` fits inside a single tooth: `∃k: w·a ≥ k − 1/13 ∧ w·b ≤ k + 1/13` for `J = [a,b]`. *(An interval touching two teeth contains the gap between them.)*

**(ii) [the band system]** Since `M(V) < 1/13` forces the killer to cover EVERY component of `L_B(1/13)` (off that set a base runner is already below 1/13), each component `Jᵢ` with midpoint `cᵢ` and half-width `hᵢ` imposes
```
‖w·cᵢ‖ ≤ 1/13 − w·hᵢ            (simultaneously for all i).
```
The midpoints are rationals on the THM-592 breakpoint grid (denominators from 13·(pairwise base sums/differences)), so each condition is a RESIDUE BAND on `w` modulo an explicit denominator, of density ≈ `2/13` minus the width tax.

**(iii) [pins and window]** `V` covering forces `q | w` for every `q ∈ {2..14}` missed by `B` (the covering pins — for loose bases typically 1–2 pins, and two pins force `lcm | w`, decimating the window); compression bounds `w ≤ 13·max(B)`; primitivity constrains gcd. The admissible killer set is the FINITE intersection: bands ∩ pins ∩ window ∩ primitivity.

**(iv) [the verified closures]** For the loose bases tested, the intersection over the ENTIRE compressed window is **EMPTY** — no killer can even satisfy the necessary alignment, so `M(V) ≥ 1/13` outright:
| base | M(B) | components | pins | window | candidates |
|---|---|---|---|---|---|
| {1..11,24} | 2/25 | 2 | 13,14 | ≤312 | **0** |
| {1..11,13} | 1/12 | 4 | 12,14 | ≤169 | **0** |
| {1..10,12,13} | 1/11 | 10 | 11,14 | ≤169 | **0** |
| {2..13} | 2/15 | 2 | 14 | ≤169 | **0** |

## What remains (scoped honestly)

The criterion turns "confinement" from an analytic worry into per-base finite arithmetic: bands + pins + window, then exact checks of survivors (none arose). The remaining quantifier is over the loose-base SPACE: the pipeline must run over all primitive 12-runner bases arising in compressed covering 13-sets — a base census whose finiteness comes from the same compression/covering bookkeeping the tight case uses (bases with far outliers peel first; the compact base space is bounded). The geometric content — WHY the intersection empties — is the pin-band mechanism: two pins force `lcm ≥ 132 | w` leaving ≤ 2 window candidates before bands; one-pin bases leave ~12, and the per-component bands (density ~2/13 each, 2–10 components) eliminate those. Expect the full base sweep to close the loose case wholesale; any surviving candidate is a single exact M-check.

-> hcomp (kps lrc14_of_compressed), klein HYP-4093 (primitivity split), S47 (tight-AP CRT case, spectrum 2/25), THM-592 (breakpoint grid), THM-522/523, OPEN-Q-108.
