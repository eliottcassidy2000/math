---
id: HYP-2040
status: SUPPORTED (verified exhaustively to speeds<=100); reduction is a THEOREM modulo the gap
source: monad-researcher-2026-06-01-S1
related:
  - HYP-2004
  - HYP-2039
  - THM-391
  - THM-369
---

# HYP-2040: LRC(n=4) reduces to a single measure gap + one boundary witness; {1,2,3} is the UNIQUE tight triple

## The methodological correction (important)

In the covering/Fourier framework, the safe measure
`|SAFE| = ∫_0^1 ∏_i 1[||s_i t|| >= 1/n] dt` is a product of nonnegative
indicators, so **`|SAFE| >= 0` is trivially true for every speed system and
every n.** Therefore a *measure lower bound* can NEVER prove LRC: the extremal
arithmetic-progression speed sets have `|SAFE|` **exactly 0**. HYP-2004 item (B)
("resonance correction `>= -(1-2/n)^{n-1}` for all systems `= LRC(n)`") is the
statement `|SAFE| >= 0`, which is trivial and is **not** equivalent to LRC.

The correct equivalence is:

```text
LRC(n)  <=>  every measure-ZERO ("tight") primitive speed system has a CLOSED
             boundary witness t* with ||s_i t*|| >= 1/n for all i.
```

(If `|SAFE|>0` the closed safe set is automatically nonempty; only the
measure-zero cases need the boundary argument.) This pinpoints HYP-2039's
"set-vs-measure gap (measure 0, set nonempty)" as the *entire* content of LRC,
not a residual.

## The n=4 result

At `n=4` (threshold 1/4, 3 runners), with the mod-4 character decomposition
`|SAFE| = 1/8 + R2 + R3` (R2 closed form = THM-391):

1. **Unique tight set.** Across ALL primitive triples with speeds `<= 100`
   (135,739 triples), the *only* measure-zero triple is the AP `{1,2,3}`.
   (Contrast n=5,6, where the AP is NOT the unique tight set — `{1,3,4,7}`,
   `{1,3,4,5,9}` are also tight. n=4 uniqueness is special.)

2. **Measure gap (CONJECTURE).** Every primitive triple `!= {1,2,3}` has
   `|SAFE| >= 1/28`, with equality at `{1,6,7}`. Verified exact to speeds<=24
   and by fast scan to speeds<=100; the minimum positive value is a stable
   `1/28` from speed-bound 20 through 100. The smallest-measure family is
   `(1, 4k+2, 4k+3)`: `|SAFE| = 1/28, 1/22, 1/20, 1/19, ...` (k=1,2,3,4).
   **S552 update:** THM-392 proves the exact formula for the whole adjacent
   family `(1,q,q+1)`, including that `(1,6,7)` is the unique positive minimum
   in that family.
   **S553 update:** THM-393 proves the exact formula for every primitive
   additive-return triple `(a,b,a+b)`.  This covers the first non-adjacent
   competitors `(2,3,5)`, `(2,5,7)`, `(3,10,13)`, `(3,5,8)`, and proves the
   same gap inside that family: only `(1,2,3)` is zero and only `(1,6,7)`
   attains the positive minimum `1/28`.
   **S553 second update:** THM-394 proves the exact formula for the next
   visible corridor `(1,3,q)`.  After excluding duplicates, the AP, and the
   adjacent/additive row `(1,3,4)`, this family bottoms at `(1,3,9)` with
   measure `1/18`, safely above `1/28`.

3. **Boundary witness at the AP.** For `{1,2,3}`, `t* = 1/4` gives
   `(||1/4||, ||2/4||, ||3/4||) = (1/4, 1/2, 1/4) >= 1/4`. Lonely. ∎ (for AP)

4. **Parity lemma (PROVED).** If all three speeds are odd then `t=1/4` is a
   witness (`||odd/4|| = 1/4`); hence every tight triple must contain an even
   speed. Verified: 0 all-odd tight triples.

### Reduction theorem (rigorous, modulo the gap)
LRC(n=4) follows from the **measure-gap conjecture (item 2)**: every non-AP
primitive triple then has `|SAFE| >= 1/28 > 0`, so its closed safe set is
nonempty; and the unique AP `{1,2,3}` has the explicit witness `t=1/4`. So
n=4 is reduced to proving the single inequality "`|SAFE| < 1/28 => triple is
{1,2,3}`" — a sharper, finite-flavored target than bounding the 3-term
resonance sum directly.

## Why the gap, not the raw measure, is the right object
- `R2` is known in closed form (THM-391), but `R3` (the genuine 3-term
  resonance, a triple `chi4`-character sum over the rank-2 lattice
  `{k : k.(a,b,c)=0}`) does not collapse to one character.
- Measure lower bounds via Bonferroni are too weak (e.g.
  `|S_a∩S_b∩S_c| >= |S_a∩S_b|+|S_c|-1 = 1/6+1/2-1 < 0`), confirming the
  difficulty is concentrated exactly at the measure-zero AP.

## OPEN / next
- (A) Prove the measure gap `|SAFE| >= 1/28` off `{1,2,3}` (closes LRC n=4 by
  this methodology). Likely route: bound `P_3 = |B_a∩B_b∩B_c|` away from
  `1/4 + R2` except at the AP, using `|SAFE| = 1/4 + R2 - P3`.
- (B) Prove `{1,2,3}` is the unique tight triple for ALL speeds (no bound).
- (C) Is there an analogous unique-tight + measure-gap statement at each n
  where the AP is the unique tight set? (n=4 yes; n=5,6 NO — extra tight sets.)
  What distinguishes n=4?

## Files
- `04-computation/lrc_n4_mod4_character_monad.py` (+.out): character identity,
  THM-391 verification (0 fails), full decomposition, AP tightness.
- `04-computation/lrc_n4_tight_characterization_monad.py` (+.out): unique tight
  set & measure floor to speeds<=100; parity lemma.
- `04-computation/lrc_n4_measure_gap_monad.py` (+.out): exact smallest measures,
  `1/28` gap confirmed.
- `04-computation/lrc_n4_adjacent_family_s552.py` (+.out): exact adjacent-family
  formula, proving the main near-tight family bottoms at `(1,6,7)`.
- `04-computation/lrc_n4_additive_return_s553.py` (+.out): exact
  additive-return formula for primitive triples `(a,b,a+b)`, proving the
  next obstruction family also bottoms at `(1,6,7)`.
- `04-computation/lrc_n4_13q_family_s553.py` (+.out): exact `(1,3,q)` formula,
  proving the next post-additive corridor bottoms at `(1,3,9)` with `1/18`.
- `01-canon/theorems/THM-392-lrc-n4-adjacent-family-measure-formula.md`.
- `01-canon/theorems/THM-393-lrc-n4-additive-return-measure-formula.md`.
- `01-canon/theorems/THM-394-lrc-n4-13q-family-measure-formula.md`.
- Reflection: `07-reflections/lrc-the-measure-is-trivially-nonneg-the-tight-set-is-everything-s1.md`.
