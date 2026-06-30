---
id: HYP-3739
title: The covering-min UNIQUENESS theorem -- literal base-uniqueness FAILS (the radius-1 band has MANY coverers, n=13: 1406, including KILLER-based bases like {13..23} that use the band prime 23 as a multiple-killer rather than a transversal -- klein-S39's killer-OR-transversal dichotomy), but M-UNIQUENESS HOLDS: the construction {1,..,n-2,n(n-1)} (consecutive base + lcm outlier) is the STRICT M-minimizer (n=13: 13/157 vs killer-block {13..23}+out 13/49 and shifted {2..12}+out 2/15, both >> 13/157; n=14: 0 perturbations beat OR tie 14/183). So the covering-min is unique IN M, not in band-coverage, and the rigorous proof route is M-MINIMIZATION among killer-or-transversal band-coverers (klein's budget), NOT base-uniqueness. The construction is the CANONICAL GREEDY representation -- the all-transversal base (the consecutive low-halves {1,..,n-2}, klein's proved transversal lemma) + a single minimal killer (lcm(n-1,n)=n(n-1)) -- a ZECKENDORF/Sylvester-style canonical uniqueness (the unique greedy/efficient form among the band's killer-or-transversal representations; cf. HYP-3724 Sylvester-Egyptian greedy). Clean subcase: when 2n-3 is prime the transversal-MODE base is forced to the low-halves {1,..,(2n-4)/2}={1,..,n-2}, but killer-MODE alternatives still exist, so even then base-uniqueness fails -- only M-uniqueness survives
status: PARTIALLY RESOLVED + a CORRECTION of the naive approach. SOLID: base-uniqueness FAILS (n=13: 1406 band-coverers, exact); the construction is the strict M-minimizer (n=14: 0 of all single-perturbations beat/tie; n=13: alternatives give 13/49, 2/15 >> 13/157). The full M-uniqueness PROOF (construction <= every covering set) is klein's budget direction, NOT complete. The Zeckendorf/greedy-canonical link is a framing, not a theorem.
source: mac-mini-2026-06-30-S55
related:
  - HYP-3737  # S54: the band over-constraint forces the construction-TYPE (this refines: type, not literal base)
  - HYP-3736  # klein-S39: the killer-OR-transversal dichotomy + proved transversal lemma + budget (the proof route)
  - HYP-3724  # klein: Sylvester-Egyptian GREEDY uniqueness -- the project's greedy theme the Zeckendorf link joins
  - HYP-2566  # uniform looseness; covering-min=n/Phi6 at n=14 (the target)
results:
  - 04-computation/coveringmin_uniqueness_zeckendorf_macmini_20260630.py
---

# HYP-3739 -- the covering-min uniqueness is M-uniqueness (not base-uniqueness); the Zeckendorf canonical form

The owner asked to work the uniqueness theorem (toward a rigorous `covering-min = n/Phi_6` for `n>=12`) and to
think Zeckendorf connections. Working it overturned the naive target and clarified the real one.

## Literal base-uniqueness FAILS
The S54 hope -- "the consecutive base `{1,..,n-2}` is the *unique* set covering the radius-1 band" -- is
**false**. At `n=13` there are **1406** valid `(n-2)`-bases covering `2..n-2` and the band interior. They split
into two MODES (klein-S39's **killer-or-transversal** dichotomy): each band prime `p in (n,2n-3]` is handled
*either* by a `±`-transversal (the consecutive base does all of them, by klein's proved lemma) *or* by a
**killer** -- a multiple of `p` in the set (e.g. the base `{13,..,23}` uses `23` itself as the killer for
`p=23`). So band-coverage is far from unique.

## But M-UNIQUENESS holds
Among all these band-coverers, the **construction is the strict M-minimizer** (verified):

| n | set | M |
|---|-----|---|
| 13 | construction `{1..11, 156}` | **13/157 = 0.0828** |
| 13 | killer-block `{13..23, 36}` | 13/49 = 0.265 |
| 13 | shifted `{2..12, 13}` | 2/15 = 0.133 |
| 14 | construction `{1..12, 182}` | **14/183**; **0** of all single-perturbations beat or tie it |

The killer-based and shifted bases are valid covers but give `M >> n/Phi_6`: killers are large speeds (waste),
and dropping the small speed `1` loses the tight binding. So the covering-min is **unique in M**, not in
band-coverage. The rigorous proof route is therefore **M-minimization among the killer-or-transversal
band-coverers (klein's budget)**, not a base-uniqueness lemma -- a useful redirection.

## Convergence with klein-S40 (HYP-3738) -- the binding PROVED, and Zeckendorf = OSTROWSKI
klein independently hit the same two findings and supplied the rigor + the exact Zeckendorf answer:
- **The construction's binding is PROVED** (klein-S40): at rotation `a=n` and `D=Phi_6=(n-1)n+1 ≡ 1 mod (n-1)`,
  the images are the multiples of `n mod Phi_6` (a core AP) plus the killer `≡ (n-2)n+1`, one above the top core
  point, splitting the wrap gap `2n+1` into `{1, 2n}`; gaps `{1, n^{(n-3)}, 2n}`, deep hole `2n` (radius `n`),
  sum `= Phi_6`. **The unit gap IS the `+1` in `Phi_6`.** Verified `n=5..9`.
- **Uniqueness of invariants, not of the covering** (klein-S40): the binding `D` and rung `k` are the unique
  invariants (all extremal coverings share them), but the covering is NOT unique -- `n=7` has exactly two
  (`{1,2,5,6,7,8}` spreader-route, `{1,4,5,6,7,11}` band-prime-killer-route). This is exactly my
  "base-uniqueness fails / killer-or-transversal modes."
- **Zeckendorf = OSTROWSKI** (klein-S40, the precise answer): the ladder denominators `k(n-1)+1` are the
  2-term continuants `K(n-1,k)`, and `M = [0;n-1,k]` is the unique **Ostrowski representation** -- the
  *continued-fraction generalization of Zeckendorf* (Zeckendorf = Ostrowski for the golden CF `[1;1,1,..]`).
  So the "Zeckendorf connection" is literally Ostrowski numeration w.r.t. the CF `[0;n-1,...]`; the binding
  three-gap `{1, n, 2n}` realizes that unique representation geometrically. My "canonical greedy form" is the
  Ostrowski-canonical statement. (I ceded HYP-3738 to klein; this is HYP-3739.)
- **Still open** (both): the SPREAD-regime binding (`n=7..11` spreads are NOT three-gap, e.g. n=7 spread gaps
  `{1,1,2,2,3,4}`), and the full `M`-uniqueness proof.

## The Zeckendorf / greedy-canonical connection
The construction is the **canonical greedy representation** of the covering: the *all-transversal* base (the
consecutive low-halves `{1,..,n-2}`, the smallest possible -- klein's transversal lemma) plus a *single minimal
killer* `lcm(n-1,n)=n(n-1)`. Among the band's many killer-or-transversal "representations," this canonical
greedy one uniquely minimizes `M` -- a **Zeckendorf-style canonical uniqueness**: as Zeckendorf's theorem
singles out the unique greedy (non-redundant) Fibonacci representation among many sum representations, here the
greedy (max-transversal, min-killer, lowest-speed) covering is the unique `M`-optimal one among many valid
covers. It joins the project's greedy-uniqueness thread -- klein's **Sylvester-Egyptian** greedy tower
(HYP-3724) and the continued-fraction canonical form (the Stern-Brocot ray, HYP-3732). The methodological
upshot: prove `M`-uniqueness by a **greedy/exchange argument** (any covering exchanges toward the canonical
construction, never lowering `M` below `n/Phi_6`), the Zeckendorf proof pattern.

## The clean prime subcase
When `2n-3` is prime (`n=13,16,17,...`), the largest band prime is `p=2n-3` with `(p-1)/2=n-2` pairs; a
*transversal-mode* base of `n-2` small speeds (needed to cover `2..n-2`) must pick the **low half** `k` of each
pair `{k, p-k}` (`p-k >= n-1` is too big), forcing exactly `{1,..,n-2}`. So the transversal-mode base is
uniquely the consecutive set -- but killer-mode alternatives still exist, so even here only `M`-uniqueness
holds. (`n=14`: `2n-3=25` is not prime, the binding prime is `23`, and the consecutive base is forced by the
*combined* constraints mod `17,19,23` + smallness -- 0 perturbations survive.)

## What it buys
Corrects the uniqueness target (base-uniqueness is false; `M`-uniqueness is the right statement and is strongly
evidenced), identifies the rigorous route (klein's killer-or-transversal budget / a greedy-exchange argument,
not a base lemma), and frames the construction as the **canonical greedy** covering -- a Zeckendorf/Sylvester
canonical-uniqueness, tying the LRC covering-min to the project's greedy-representation thread. The LRC14 target
`covering-min = 14/183` stands (strict M-minimizer, no perturbation ties), with the proof reduced to the
greedy/budget `M`-optimality.
