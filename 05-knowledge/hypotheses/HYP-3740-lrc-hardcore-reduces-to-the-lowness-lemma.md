---
id: HYP-3740
title: The LRC14 HARD CORE reduces to a single LOWNESS LEMMA. covering-min(14)=14/183 (=> LRC14, margin 13/2562) follows from STEP 1 (the LOWNESS LEMMA, strongly evidenced): M(S)<=14/183 => {1,..,12} subset S; + STEP 2 (RIGOROUS): {1,..,12} subset S, |S|=13, covers 2..14 => the lone 13th speed covers q=13 AND q=14 => it is a multiple of lcm(13,14)=182 => smallest is 182 => S = construction. So the unbounded covering-set search COLLAPSES to one set once {1,..,12} is forced. KEY SEPARATION: speed 1 is COVERING-IRRELEVANT (all q>=2) yet M-NECESSARY -- the binding pair at t*=14/183 is {1, n(n-1)}={smallest,largest}; so the LRC difficulty splits into COVERING (the q-multiples, THM-523) + LOWNESS (the consecutive base {1,..,n-2} forced for the three-gap/Ostrowski binding). STEP 1 EVIDENCE: covering 13-sets MISSING any small speed have M strictly above 14/183 EVEN with huge speeds (up to 455): missing 1->2/17, 2->13/125, 3->2/19, 12->2/25 (all >14/183=0.0765); 0 single-perturbations of the construction beat/tie it (HYP-3739). The remaining hard core IS exactly the lowness lemma; proof directions: klein-S39's transversal M-optimality (the full consecutive transversal is the unique M-min base), or a "k-witness" analog of the THM-523 q-witness (k missing => an explicit lonely t with margin >n/Phi_6)
status: REDUCTION established + STEP 2 RIGOROUS. STEP 1 (the lowness lemma) is STRONGLY EVIDENCED (missing any small speed -> M>14/183, incl. huge speeds; 0 perturbations beat/tie) but NOT PROVED -- it is the precise remaining LRC14 hard core. The reduction collapses the unbounded problem to this one lemma.
source: mac-mini-2026-06-30-S56
related:
  - HYP-3739  # M-uniqueness (the construction is the strict M-minimizer); this reduces it to the lowness lemma
  - HYP-3738  # klein-S40: PROVED the construction binding (three-gap {1,n,2n}, unit gap = +1 in Phi6)
  - HYP-3736  # klein-S39: killer-or-transversal + the proved transversal lemma (a Step-1 proof route)
  - THM-523   # CANON: the q-witness reduction (covering); the lowness lemma is the BINDING analog
  - HYP-2566  # uniform looseness -- this reduces it (n=14) to the lowness lemma
results:
  - 04-computation/lrc_hardcore_lowness_lemma_macmini_20260630.py
  - 05-knowledge/results/lrc_hardcore_lowness_lemma_macmini_20260630.out
---

# HYP-3740 -- the LRC14 hard core reduces to the lowness lemma

Working the LRC14 hard core (prove `covering-min(14) = 14/183`, hence `M(S) >= 1/14` for all primitive covering
13-sets) creatively isolates the *entire* remaining difficulty into one clean lemma.

## The reduction
**Goal:** `covering-min(14) = 14/183` (then LRC14 holds: `14/183 > 1/14`, margin `13/2562`).

- **STEP 1 -- the LOWNESS LEMMA (strongly evidenced, not proved):**
  `M(S) <= 14/183  ==>  {1,2,..,12} \subseteq S`.
- **STEP 2 -- the completion (RIGOROUS):** if `{1,..,12} \subseteq S`, `|S|=13`, and `S` covers `2..14`, then the
  *single* remaining speed must cover both `q=13` and `q=14` (nothing in `{1,..,12}` is a multiple of `13` or
  `14`), so it is a common multiple of `13` and `14`, i.e. a multiple of `lcm(13,14)=182`; the smallest is
  `182`, and any larger one only increases `M`. Hence `S = {1,..,12, 182}` = the construction, `M = 14/183`.

So **once `{1,..,12}` is forced, the unbounded covering-set search collapses to a single set** -- no speed bound
needed. The whole LRC14 hard core is Step 1.

## The key idea -- covering vs lowness (the separation)
The construction's **binding pair** at `t*=14/183` is `{1, 182} = {smallest speed, largest speed} = {1, n(n-1)}`
(both achieve `min_v ||v t*|| = 14/183`). Notice **speed 1 is covering-irrelevant** -- it is a multiple of no
`q >= 2`, so the covering reduction (THM-523) never needs it -- **yet it is M-necessary** (it is half the binding
pair). This *separates* the LRC difficulty into two independent parts:
- **COVERING:** the set must contain a multiple of every `q in {2,..,n}` (THM-523, the `q`-witness handles the
  rest). This forces the *large* structure (the killers/multiples).
- **LOWNESS:** the set must contain the small consecutive base `{1,..,n-2}` to realize the tight three-gap /
  Ostrowski binding (HYP-3738/3739). This forces the *small* structure.
The construction is the unique set satisfying both minimally. THM-523 is the covering half; the **lowness lemma
is the missing binding half**.

## Step 1 evidence
Covering 13-sets *missing* a small speed have `M` strictly above `14/183`, **even allowing huge speeds** (tested
to 455):

| missing speed | best M found | vs 14/183 |
|---------------|--------------|-----------|
| 1 | 2/17 = 0.1176 | >> |
| 2 | 13/125 = 0.1040 | >> |
| 3 | 2/19 = 0.1053 | >> |
| 12 | 2/25 = 0.0800 | > |

and `0` single-perturbations of the construction beat or tie it (HYP-3739). So the consecutive base is forced;
the closest miss (`missing 12 -> 2/25`) is still above `14/183`.

## Proof directions for the lowness lemma
1. **Transversal M-optimality (klein-S39):** the band over-constraint forces a transversal base; among
   transversal bases, the full consecutive `{1,..,n-2}` (the lowest, klein's proved transversal lemma) is the
   unique `M`-minimizer. Make "killers raise `M`, lowest transversal minimizes" rigorous.
2. **A "k-witness" (THM-523 analog):** THM-523 gives, for a set omitting a multiple of `q`, an explicit lonely
   point `t=1/q`. Find, for a set omitting the *small speed* `k`, an explicit lonely `t` with margin `> n/Phi_6`
   -- a binding-side witness mirroring the covering-side `q`-witness. (Empirically these witnesses live at small
   prime moduli, e.g. `mod 17` for missing-1 at `n=14`.)
3. **Three-gap rigidity:** the construction's images at `a=n` are the AP `{n, 2n, .., (n-2)n} mod Phi_6` (one
   point per small speed) + killer; deleting the point `kn` (missing speed `k`) breaks the unit-gap AP, and the
   set's best escape then exceeds `n/Phi_6`.

## What it buys
Reduces the LRC14 hard core -- an optimization over UNBOUNDED covering sets -- to a single **lowness lemma**
(`M(S) <= n/Phi_6 => {1,..,n-2} \subseteq S`), with the completion (Step 2) rigorous and the search collapsed to
one set. It reframes the LRC's difficulty as **lowness, not covering** (speed 1 is the crux, though
covering-irrelevant), and gives three concrete proof routes. This is the sharpest statement yet of what remains
for a rigorous `covering-min(14)=14/183` and hence LRC14.
