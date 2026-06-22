---
id: THM-560
title: LRC(14) structured-tiler characterization -- the difference-closed exact tilers are exactly the dilated intervals d*{1..13}, all tight
status: PROVED (elementary)
author: kind-pasteur-2026-06-22-S31n
depends_on:
  - HYP-+2888   # mac-mini: exact-coverage extremal is scaling-invariant (this resolves its structured case)
  - OPEN-Q-108  # tight-locus finiteness (this isolates the residual SPORADIC case)
verified: computational (lrc_tiling_rigidity_kps.py, lrc_tight_vs_counterexample_kps.py) + proof below
---

# THM-560: the difference-closed exact tilers of LRC(14) are exactly d*{1..13}

## Setup
For a set `S` of 13 distinct positive integers (the runner speeds of LRC(14), relative frame),
define the **lonely constant** `M(S) = max_{t in [0,1)} min_{s in S} ||s t||` and the **safe set**
`Safe(S) = { t in [0,1) : ||s t|| >= 1/14 for all s in S }`. LRC(14) for `S` <=> `Safe(S)` nonempty.
Call `S` an **exact tiler** when `meas(Safe(S)) = 0` (the danger arcs `U_s={||st||<1/14}` cover
`[0,1)` up to measure zero). Call `S` **difference-closed** when `|s - s'| in S` for every
`s != s'` in `S`.

mac-mini (HYP-+2888) reframed the LRC(14) wide-bound crux as: *which `S` are exact tilers?* and
conjectured "only `d*{1..13}`." THM-560 resolves this **for the difference-closed sets** (and shows
the conjecture is INCOMPLETE -- a sporadic tiler exists).

## Statement
**(a)** A 13-element set `S` of positive integers is difference-closed **iff** `S = d*{1,...,13}`
for some integer `d >= 1` (a dilated interval).

**(b)** Every dilated interval `d*{1,...,13}` is an exact tiler, and is **tight**:
`M(d*{1..13}) = 1/14`, with `Safe = { j/(14 d) : gcd(j,14)=1 }` -- exactly 6 points, measure 0.

**(c)** Hence among difference-closed 13-sets, *exactly* the AP-dilates tile, all at `M = 1/14`.
The remaining exact tilers are **non-difference-closed (sporadic)**; one exists -- the
Goddyn-Wong-type set `S_GW = {1,...,11, 13, 24}` (verified `M = 1/14` exactly, witness `t = 5/14`,
`meas(Safe)=0`). Their classification is the residual hard core (OPEN-Q-108).

## Proof
**(a, <=)** `d*{1..13}` has nonzero differences `d*{1,...,12} ⊆ d*{1,...,13}`. Difference-closed.

**(a, =>)** Let `S` be difference-closed, `d = min(S)`. For `s in S` with `s > d`, the difference
`s - d > 0` lies in `S`. Iterating `s -> s-d -> s-2d -> ...` stays in `S` and `>= d` (minimality),
so the chain terminates at `d`; thus `s = (k+1)d` is a multiple of `d`. So `S = d*S'` with `S'`
difference-closed and `min(S') = 1`. For `s' in S'`, `s'-1 in S'` (difference with `1`), so
`s', s'-1, ..., 1` all lie in `S'`; hence `S' ⊇ {1,...,max(S')}`. With `|S'| = 13`, `max(S') = 13`
and `S' = {1,...,13}`. Therefore `S = d*{1,...,13}`. ∎

**(b)** Fix `S = d*{1..13}` and `t in Safe(S)`. Put `u = (d t) mod 1`. Consider the 14 points
`P = { k u mod 1 : k = 0,1,...,13 }`. Any two have circle distance
`|| (i-j) u || = || |i-j| u ||` with `|i-j| in {1,...,13}`. Now `||m u|| = ||m d t|| = ||(md) t||`
and `md in d*{1..13} = S` for `m in {1..13}`, so by `t in Safe(S)` every pairwise distance is
`>= 1/14`. Fourteen points on a circumference-1 circle, pairwise `>= 1/14`: the 14 cyclic gaps sum
to `1` and each (a distance between adjacent points) is `>= 1/14`, forcing **all gaps `= 1/14`**.
So `P` is the equally-spaced set `{0, 1/14, ..., 13/14}`, which forces `u = j/14` with
`gcd(j,14)=1` (else `{ku}` misses residues). Hence `d t ≡ j/14`, i.e. `t in {j/(14d) : ...}` --
finitely many points, so `meas(Safe) = 0`. At such `t`, `min_s ||st|| = min_k ||kj/14|| = 1/14`
(attained at the `k` with `kj ≡ ±1 mod 14`), and no `t` gives `min > 1/14` (that would be an open
safe set, contradicting `meas(Safe)=0`); therefore `M = 1/14`. ∎

**(c)** Immediate from (a),(b); `S_GW` checked computationally (not difference-closed: e.g. `24-1=23
∉ S_GW`; tight with `M = 1/14`). ∎

## The sporadic mechanism (part c, structural -- guides OPEN-Q-108)
At a tight witness `t`, the 13 speeds map to residues `||s t|| -> k/14`. Two regimes
(`lrc_sporadic_mechanism_kps.py`):
- **AP `{1..13}` at `t=3/14`**: the speed->residue map is a BIJECTION onto `{1,...,13}/14` --
  equal spacing, no gap, no collision. (The difference-closed case.)
- **GW `{1..11,13,24}` at `t=5/14`**: residue `4/14` is MISSING (a gap), residue `8/14` is DOUBLED
  (a collision -- speeds `10` and `24`, since `24 ≡ 10 mod 14`). So GW is exactly the AP with
  `12 -> 24`: because `24 ≡ 10 (mod 14)` the replaced speed lands on `10`'s residue at the witness,
  freeing `12`'s old residue `4`. A BALANCED gap+collision (one each) that still isolates `0` at
  exactly `1/14`.

So the sporadic tilers trade equal-spacing for a balanced gap+collision pattern, realized as an
AP-perturbation `s -> s'` with `s' ≡ s'' (mod 14)` for another speed `s''`, chosen so the witness
isolation of `0` survives. Crucially this is a **Diophantine (integer) phenomenon, not mod-14**: the
integer `24 != 10` supplies coverage at non-witness `t` that the residue `10` alone would not, which
is why GW tiles even though its residues mod 14 are NOT a complete system (missing `12`, doubling
`10`). **OPEN-Q-108 = count the balanced gap+collision perturbations that preserve `M = 1/14` over
ALL `t`** (tightness, not just witness isolation). This is a sharper, finite-combinatorial target
than the full measure-rigidity.

**GW is ISOLATED among single-replacements (supports finiteness).** Searching all AP single-speed
replacements `{1..13}\{rem} ∪ {v}`, `v<=60` (`lrc_sporadic_verify_kps.py`, with the RIGOROUS critical
points `t=k/(s_i +- s_j)`): the ONLY tight one is GW (`rem=12, v=24`). The tempting family
`{1..11,13,12j}` is LOOSE for `j>=3`: `M({1..11,13,36})=3/41`, `M(..,48)=4/53`, `M(..,60)=1/13`, all
`> 1/14` (positive safe measure, NOT tilers). [An earlier coarse search wrongly flagged these as tight
-- its candidate set missed the true argmax, underestimating `M`; the pairwise critical points fix it.
Disciplinary note: verify apparent extremal families with the exact critical set, not a breakpoint
sample.] So the witness-isolation condition (`min=1/14` at one `t`) is FAR from sufficient; global
tightness over all `t` is what isolates GW. This is real evidence the sporadic locus is small/finite.

## Finiteness of the sporadic locus (kps-S31o -- census + equidistribution mechanism)
**Single-swap census (RIGOROUS):** among all `{AP\{rem}} U {v}`, `v<=50`, with exact critical points
`t=k/(s_i+-s_j)`, the ONLY non-AP tight set is GW (`rem=12,v=24`)
(`lrc_sporadic_finiteness_kps.py`). So GW is the unique single-swap sporadic.

**Why (mechanism):** `M({1..11,13}) = 1/12 > 1/14` (the 12-subset is lonely at `1/12`). Adding `v`
lands `M` at exactly `1/14` only for tuned `v`: `v=24 = 2*12` kills the subset's `t=1/12` witness
(`24/12=2 in Z`) and lands precisely at `1/14`; nearby `v` keep `M=1/12` (don't kill the witness) and
far `v` give intermediate `M in {3/41,4/53,1/13,2/25}` -- ALL `> 1/14`, none equal.

**Finiteness heuristic (with identified gap):** an exact tiler (tight, meas-0 safe) must have BOUNDED
speeds. If `v=max(S)` is large, `U_v={||vt||<1/14}` (measure `1/7`) EQUIDISTRIBUTES (Weyl), so
`meas(U_v ∩ L) -> meas(L)/7` for the 1/14-lonely set `L` of `S\{v}`; thus `U_v` covers only `~1/7`
of `L`, leaving `~6/7 meas(L) > 0` safe -- NOT a tiling. So tilings need `max(S)` bounded => finitely
many (up to dilation). GAP: this needs `meas(L)` bounded below, which fails if `S\{v}` is itself
recursively near-tight (its own 1/14-lonely set tiny) -- the recursive-tightness gap is exactly why
OPEN-Q-108 stays open. **mod-27 = 2n-1 = 3^3 (HYP-2138):** GW swaps `12->24` WITHIN the gcd-3 shell
(both `≡ 0 mod 3`); sporadics exist iff `2n-1` is composite, finitizing the swaps to the non-unit
shells of `27`.

## Honest positioning (the search reframe, THM-523)
This whole tiling/tight-locus line is the **trivial boundary layer** of LRC(14): by THM-523 (q-witness
lemma) any speed set omitting all multiples of some `q<=14` is safe at `t=1/q`, so the tight sets here
(AP, GW -- both no multiple of 14) are TRIVIALLY lonely at `t=1/14`. The GENUINE crux is the
**covering sets** (containing a multiple of every `q in {2..14}`), where the easy witnesses all fail;
that is mac-mini's covering-system / prime-basis / single-atom over-determination (HYP-+2878). THM-560
fully resolves the structured boundary and bounds the sporadic boundary; it does not touch the covering
crux. The value is a clean, complete account of the boundary geometry + isolating where the real
difficulty is NOT.

## Significance
- **Resolves the structured half of mac-mini's crux**: the difference-closed (equivalently
  "interval-like") exact tilers are *completely* pinned down -- exactly the AP-dilates, all tight,
  by a one-paragraph pairwise-distance/equal-spacing argument. No additive energy, no convexity, no
  measure-LP -- just difference-closure forcing equal spacing.
- **Corrects HYP-+2888**: "only `d*{1..13}` tile" is FALSE -- `S_GW` tiles too. The tight locus is
  at least `{ AP-dilates } ∪ { GW-type }`, matching the known LRC tight locus `{AP, Goddyn-Wong}`.
- **Isolates the genuine hard core**: OPEN-Q-108 = "classify the SPORADIC (non-difference-closed)
  exact tilers." The structured case is no longer the obstruction; the sporadic finiteness is.
- Connects to the project frame: difference-closure = the speed set being its own difference set =
  the AP being the unique *self-similar* (Sidon-opposite) speed structure; the equal-spacing rigidity
  is the 1-D shadow of the regular/transitive extremality that runs through the whole tournament story.

-> OPEN-Q-108, HYP-+2888 (mac-mini), HYP-2885, `lrc-coverage-transcends-the-h-level-...md`.
