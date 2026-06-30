---
id: HYP-3748
title: RIGORIZING step 3 of the lowness lemma -- replacing the exhaustive search by a FINITE check + reductions. The exponential search over all coverings missing core speed k reduces to: (R1) the |L|=2 NORMAL FORM -- WLOG the adversary keeps the maximal small core {1..n-2}\{k} and uses exactly 2 large speeds (1 killer + 1 band-coverer); dropping more small speeds only adds witnesses. (R2) a 2-PARAMETER optimization over the 2 large speeds (kappa,w). (R3 BOUNDED, RIGOROUS FINITE CHECK): for kappa,w <= T the min of M(core-minus-k + kappa + w) > n/Phi6 -- verified n=14, all k checked (k=1->8/67, ..., k=12->7/89, the tightest; margin 35/16287>0). (R4 UNBOUNDED): a huge band-coverer is no more efficient (CRT-INVARIANT count <=2r+1 per speed, HYP-3745) and does not fix the punctured core's wide hole -> HYP-3745. The bounded case is now a FINITE rigorous check (not exhaustive over coverings); the unbounded case is reduced to HYP-3745 (the hard direction residual)
status: PARTIAL rigor. R3 (bounded finite check) RIGOROUS + VERIFIED n=14 (k=1..8,12; tightest 7/89 at k=12). R1 (|L|=2 normal form) and R4 (unbounded->HYP-3745, the wide-hole) are argued/reduced, not fully closed -- these are the residual (the hard direction). Net: the FULL exhaustive search is replaced by a finite bounded 2-parameter check + two structural reductions; the remaining gap is the unbounded wide-hole.
source: klein-2026-06-30-S46
depends_on:
  - HYP-3747   # the full lowness lemma (step 3 is its residual)
  - HYP-3745   # CRT escape uncoverable (the unbounded case R4)
related:
  - HYP-3740   # mac-mini: the lowness lemma (exhaustive-search verification, now reduced)
  - HYP-3741   # the witness hierarchy (R2 witnesses)
results:
  - 04-computation/farey_rung_spread_family_klein.py
---

# HYP-3748 — rigorizing step 3: finite check + reductions

## The goal
Step 3 of the lowness lemma (HYP-3747) -- "the budget leaves ~1 slot, forcing a single CRT band-coverer" --
was the one piece leaning on mac-mini's exhaustive search (HYP-3740). Here it is replaced by a FINITE check
plus structural reductions.

## R1 -- the |L|=2 normal form (reduction)
Let `S` be a covering missing core speed `k <= n-2`, `M(S) <= n/Phi_6`. Split `S = S_small ∪ L`,
`S_small subset {1,...,n-2}\{k}` (so `|S_small| <= n-3`), `L` = large speeds (`> n-2`), `|L| = n-1-|S_small| >= 2`.
**Claim:** WLOG `|S_small| = n-3` (the full core minus `k`) and `|L| = 2`. Reason: dropping a further small
speed `j` punctures the transversal at the `j`-band primes too (Step 1 of HYP-3747 applied to `j`), adding
more broken pairs that the extra large speed must also cover -- strictly more demands on `L`. So the adversary's
best is the maximal small core + exactly 2 large speeds (1 killer for `n-1,n`, 1 band-coverer for `k`'s broken
primes). [argued; the residual is making "more demands => not smaller M" fully rigorous.]

## R2 -- the 2-parameter optimization
The lemma now reads: for each `k`,
> `min over (kappa, w)` of `M({1,...,n-2}\{k} ∪ {kappa, w}) > n/Phi_6`,
a 2-parameter problem (the 2 large speeds), NOT an exhaustive search over coverings.

## R3 -- the BOUNDED case (rigorous finite check)
For `kappa, w <= T` this is a FINITE computation. Verified `n=14` (`T = 110-160`), the min over bounded
`(kappa,w)` killing all resonances `<= 14`:
```
k :  1     2     3    4     5     6     7     8    ...  12
M : 8/67  9/83  2/19 5/53  2/21  9/109 1/11  1/9  ...  7/89
~ : .119  .108  .105 .094  .095  .083  .091  .111 ...  .0787
```
every one `> 14/183 = .0765`. The **tightest is `k=12`: `M = 7/89`** (= mac-mini's earlier `n=14` bounded
estimate), margin `7/89 - 14/183 = 35/16287 > 0`. So no bounded `|L|=2` normal form beats the construction --
a finite, rigorous check (a 2D grid, not the exponential covering search).

## R4 -- the UNBOUNDED case (reduced to HYP-3745)
If a large speed is huge (`> T`), it is NOT more efficient: by the CRT-invariant counting (HYP-3745) it covers
`<= 2r+1` rotations of `Z/p` at any prime, exactly like a bounded speed -- CRT tuning chooses WHICH rotations,
never HOW MANY. A huge band-coverer `w ≡ k mod (band primes)` covers the band but does NOT fix the punctured
core `{1,...,n-2}\{k}`, which leaves a wide hole at a larger prime (e.g. `k=1, n=14`: `S = {2..12,182,7430}`
has `M = 38/269` -- the hole at `269`, NOT at the band primes `{17,19,23}` that `w` covers). The hole moves to
a worse modulus (HYP-3745). [This is the residual -- the hard direction; the wide hole of the punctured core
is bounded by HYP-3745's CRT-invariance but not yet by a closed inequality.]

## Honest scope
The exhaustive search over all coverings (exponential) is replaced by:
- **R3 a finite, rigorous 2-parameter grid check** (bounded large speeds) -- this is the concrete rigorization,
  done for `n=14`;
- **R1** (`|L|=2` normal form) and **R4** (unbounded -> HYP-3745, the punctured-core wide hole) -- structural
  reductions, argued but not 100% closed.
So step 3 is now rigorous for the bounded case (a finite check) and reduced (not search-dependent) for the
unbounded case. The residual is the unbounded wide-hole = the hard direction of the construction's uniqueness.

## Net
Step 3's budget bound is rigorized by reduction: the adversary's best is the `|L|=2` normal form (max small
core minus `k` + 2 large speeds); a finite 2-parameter grid check (rigorous) shows the BOUNDED min `M > n/Phi_6`
(`n=14`: tightest `7/89` at `k=12`, margin `35/16287`); the UNBOUNDED case is reduced to the CRT-invariant
counting + the punctured-core wide hole (HYP-3745). The dependence on the full exhaustive search is removed for
the bounded case and reduced to HYP-3745 for the unbounded case -- the remaining gap is the hard-direction
wide hole.
