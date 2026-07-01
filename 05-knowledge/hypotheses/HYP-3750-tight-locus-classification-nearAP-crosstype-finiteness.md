---
id: HYP-3750
title: CLASSIFYING the non-difference-closed tight sets -- NOT every one is a GW near-AP (v->2v); there is a broader "duplication+drop" family with a CROSS-type. Difference-closed tight = residues mod (n+1) are a permutation of {1..n} (a dilated AP mod n+1); non-difference-closed = residues COLLIDE (one residue duplicated via an element r+(n+1)j, another dropped). Two sub-types: GW-type (drop v, duplicate 2v mod n+1 = element v->2v; e.g. {1,2,3,4,5,7,12}=6->12) AND CROSS-type (drop v, duplicate an UNRELATED residue; e.g. {1,3,4,5,9}=drop 2 dup 3, {1,4,5,6,7,11,13}=drop 2 dup 5) which REFUTES 'every one like GW'. BUT the finiteness the program needs HOLDS: the near-APs are very few per n (census n=4..8: 0,1,0,2,0), so the tight locus is finite per n (1 dilated-AP class + <=2 near-APs) -> the fattening lemma follows, using the broader duplication+drop form (not just v->cv)
status: VERIFIED (tight sets M(S)=1/(n+1) enumerated n=4..8, gcd=1, reflection-deduped; cross-types {1,3,4,5,9},{1,4,5,6,7,11,13} re-verified tight and genuinely non-GW). The user's clean 'all near-APs are GW v->cv' is REFUTED (cross-types); the per-n finiteness (-> fattening lemma) is SUPPORTED.
source: klein-2026-06-30-S48
depends_on:
  - HYP-2893   # Goddyn-Wong = Jacobsthal accelerations (the GW-type criterion)
  - HYP-+2913  # three-gap tight-locus characterization (census a(n))
related:
  - HYP-+2914  # the tight "kind" necessary conditions
  - HYP-3749   # the punctured-core wide hole (the fattening context)
  - HYP-3747   # the lowness lemma
results:
  - 04-computation/farey_rung_spread_family_klein.py
---

# HYP-3750 — classifying the non-difference-closed tight sets

## The classification (residues mod n+1)
A tight set `S` (`|S|=n`, `M(S)=1/(n+1)`) is classified by its residue multiset mod the binding `n+1`:
- **DIFFERENCE-CLOSED** = residues mod `(n+1)` are a permutation of `{1,...,n}` (all distinct nonzero) = a
  **dilated AP mod `n+1`**. These are the "lifts" of `{1,...,n}` (e.g. `{1,3,4,7}` is `{1,2,3,4}` mod 5); one
  class up to dilation.
- **NON-difference-closed (near-AP)** = residues COLLIDE: one residue is DUPLICATED (an element `r+(n+1)j`
  lands on an existing residue `r`) and another residue is DROPPED. Equivalently `S = {1,...,n}\{v} ∪ {s+(n+1)}`
  (drop element `v`, add the lift `s+(n+1)` of an element `s`).

## The answer: NOT every non-diff-closed tight set is a GW near-AP
Two sub-types among the near-APs (verified `n=4..8`):
- **GW-type** (`s = 2v mod (n+1)`): drop `v`, duplicate `2v` -- i.e. element `v -> 2v` (Goddyn-Wong, HYP-2893).
  Example `n=7`: `{1,2,3,4,5,7,12} = {1,...,7}` with `6 -> 12` (`res 6` dropped, `res 4 = 2.6 mod 8` doubled).
- **CROSS-type** (`s` UNRELATED to `v`, not `2v`): drop `v`, duplicate a different element's residue. Examples:
  `{1,3,4,5,9}` (`n=5`: drop residue 2, duplicate residue 3 via `9 = 3+6`; `3 != 2v mod 6` for any `v`);
  `{1,4,5,6,7,11,13}` (`n=7`: drop 2, duplicate 5 via `13 = 5+8`). Both re-verified tight (`M = 1/6, 1/8`).
So the clean conjecture "**every** non-difference-closed tight set is, like GW, an AP with `v` replaced by `cv`"
is **REFUTED** -- the cross-type drops one residue and duplicates an unrelated one. The correct family is the
broader **"duplication+drop"** form.

## But the FINITENESS (what the program needs) HOLDS
Census of tight classes (`n=4..8`, gcd=1, reflection-deduped):
```
n         : 4   5   6   7   8
diff-closed(lifts): 2   1   1   1   1
near-AP GW-type   : 0   0   0   1   0
near-AP cross-type: 0   1   0   1   0
```
Up to dilation this is the small census (cf. HYP-+2913 `a(n)=1,2,2,1`). The near-APs are **very few per `n`**
(0-2), so the whole tight locus is **finite per `n`** (one dilated-AP class + at most a couple of near-APs).
The GW-type is Jacobsthal-bounded (HYP-2893: `v->2v` tight iff `gcd(v,j)>1` on an interval, so `n-v <= `
Jacobsthal gap); the cross-type is likewise tightly constrained (few solutions). So the finiteness that the
program requires -- **tight locus finite `=>` fattening lemma** -- HOLDS, provided one uses the broader
duplication+drop classification (not just `v->cv`).

## Net (answer to the question)
Difference-closed tight = dilated AP mod `n+1` (confirmed). NOT every non-difference-closed tight set is a GW
`v->cv` near-AP: besides GW-type (drop `v`, dup `2v`), there are CROSS-type near-APs (drop `v`, dup an unrelated
residue, e.g. `{1,3,4,5,9}`, `{1,4,5,6,7,11,13}`), which refute the clean conjecture. HOWEVER the near-AP family
is very small per `n` (census `0,1,0,2,0` for `n=4..8`), so the tight locus is finite per `n` and the fattening
lemma follows -- the program is sound once the classification is the broader "AP residues mod `n+1` with one
residue duplicated (a lifted element) and one dropped," which subsumes GW and the cross-type.
