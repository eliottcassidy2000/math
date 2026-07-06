# Small statements on the gap members: the wall is the mediant, and defects fill the boundary

*kind-pasteur-2026-07-06-S34 — strengthening understanding of the open residual
through small, precise observations about the two nonempty cases (n=7, n=8), where
the gap members can actually be examined.*

Only n=7 and n=8 have nonempty gaps, so they are the only families we can dissect.
Here is what they show — a handful of small statements, each checkable.

## 1. The denominator wall IS the mediant's denominator

opus's formalised wall (`LRCFareyGap.lean`) says a gap value `p/q` has `q ≥ 3k+2`.
Observe: `3k + 2 = (k+1) + (2k+1)` is exactly the denominator of the **mediant**
`3/(3k+2)` of the gap `(1/(k+1), 2/(2k+1))` (numerators `1+2`, denominators
`(k+1)+(2k+1)`). So the wall is not an arbitrary bound — it *is* the statement
"no gap value is simpler than the mediant." The mediant is the unique
smallest-denominator rational in the gap, and the wall is its floor.

## 2. The wall is tight at n=8, slack at n=7

- **n=8 (k=7):** the *only* realised gap value is `3/23`, and `23 = 3·7+2` — the
  mediant, **at the wall**. The wall is achieved.
- **n=7 (k=6):** the mediant `3/20` (`20 = 3·6+2`, the wall) is **not** realised;
  the only realised value is `5/33` — a *deeper* Stern–Brocot descendant (the
  mediant of `3/20` and `2/13`), at denominator `33 > 20`.

So even among nonempty cases the mediant is not always the realised value: n=8
hits it, n=7 skips it for a level-2 descendant. The wall is tight exactly when the
mediant is realised.

## 3. The mediant realiser = base AP in the bulk + defects on the boundary

The n=8 realiser `{1,2,3,4,5,7,18}` at its witness `t = 4/23` splits cleanly:

- the **base AP** `{1,2,3,4,5}` maps to residues `{4,8,12,16,20}` — an arithmetic
  progression of gap `4 = a` (the witness numerator), evenly spanning the safe
  band `[3, 20]`;
- the two **defects** `{7, 18}` map to the *boundary* residues `{5, 3}` — the ones
  nearest the forbidden `[0, 3)` band — and `18 ↦ 3` is a **binding** runner
  (dist 3 = the mediant numerator).

So the near-tight structure is literally **bulk + boundary**: the base AP tiles
the interior like a coarse root-of-unity grid, and the defects are *tuned to the
boundary residues*, each contributing the minimal distance. The `M`-value is set
by how close a defect must sit to `0` — the mediant numerator `3` here.

## 4. Near-tight ⟹ small g (three-gap), and the numerator is the min-dist

The mediant realiser's residues `{3,4,5,8,12,16,20}` (with `0`) have gap multiset
`{3,1,1,3,4,4,4,3}` — **3 distinct gap lengths** (`g = 3`). This is the S29/S32
"spectral multiplicity" at work: a gap member is near-tight, so `g` is small
(here 3), consistent with mac-mini's three-gap route (HYP-4412). And in every
case the **`M`-numerator equals the minimum residue-distance**, carried by the
runner(s) pinned to the boundary residue `p` (or `q−p`).

## 5. The hidden connection: the density floor is a boundary-packing cost

Stringing these together: a gap member is a *base AP (bulk) plus defects (boundary)*,
and its `M`-rise above `1/(k+1)` is the **cost of placing the defects at the
boundary** without pushing any runner into the forbidden band. The window
`w(k) = 1/((k+1)(2k+1))` is the budget. At small `k` (n=7,8) the budget is large
enough to seat a defect at boundary distance `3/(3k+2)` (the mediant) — the
defects fit. As `k` grows the budget `~1/(2k²)` shrinks faster than the cheapest
boundary seat, and by `k=12` no defect placement stays under budget: the boundary
is too crowded. **The density floor is the statement that boundary-packing of
defects costs more than the window at `k=12`** — the same content as the
Selberg/Cohn–Elkies estimate (mac-mini HYP-4512/4532), read as a discrete
seating problem: the base AP occupies the roots of unity, and the defects have
nowhere cheap left to sit.

## Ledger

All five are checkable and small; none is a proof of the floor. Their value is
the *picture*: the wall = the mediant; realisation = boundary-defect seating; the
floor = the seating cost exceeding the window. It reframes mac-mini's Cohn–Elkies
LP (HYP-4532) as a concrete combinatorial packing — the base AP's roots of unity
plus a boundary that runs out of room — and pins the one number the proof turns
on: the cheapest boundary seat `3/(3k+2)` vs the window `1/((k+1)(2k+1))`, i.e.
`3(2k+1) vs 1` after clearing denominators — the mediant sits at the wall, and the
question is whether *any* seating beats it, which at k=12 it does not.

## Pointers

- `lrc_sternbrocot_realized_kps_S34.out` (realised gap values at n=7,8; mediant
  vs deeper descendant), the witness/residue dissection above.
- opus HYP-4456 (`LRCFareyGap`, the wall = `3k+2`), HYP-4466; mac-mini HYP-4502
  (metric half), HYP-4512 (Selberg `N~2k²`), HYP-4532 (Cohn–Elkies/Viazovska);
  kps S25 (Stern–Brocot), S29 (spectral multiplicity = g), S32/S33 (structure
  census, corrected).
