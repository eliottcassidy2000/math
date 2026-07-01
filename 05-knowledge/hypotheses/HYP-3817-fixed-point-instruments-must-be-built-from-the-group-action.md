---
id: HYP-3817
title: A FIXED-POINT-SENSITIVE INSTRUMENT MUST BE BUILT FROM THE GROUP ACTION (a covering invariant or a moment), NOT another transform -- because a transform is S_n-INVARIANT, it averages over the group and is blind (or only accidentally sensitive) to the fixed-point/symmetry content |Aut|. GROUNDED: (A) the U-spectrum (adjacency) DETERMINED |Aut| only by small-n luck -- at n=7 it FAILS (one U-spectrum holds both |Aut|=1 and |Aut|=3 in a 450-class sample), so no transform is a robust |Aut| detector; the skew-spectrum was blind already at n=6 (S89). (B) COMPLEMENTARITY at n=6: the 15 high-|Aut| NEEDLES (|Aut| {3:12,5:2,9:1}, the flip-rank-excess drivers) are flagged by group-action instruments (|Aut|, MFAS covering-radius) but NOT reliably by transforms; the |Aut|=1 spectral twins (S86) are flagged by the U-spectrum but NOT by |Aut|/covering. TWO fine-structures, TWO instrument kinds: FIXED-POINT/symmetry <- group-action (covering kappa, moment Var(H), |Aut|=|Fix| directly); GENERIC realization <- transforms (spectra). LRC PARALLEL: the covering-min is a fixed-point/rigidity problem (the atom, the symmetric construction), which is WHY the analytic transforms (Fourier/Delsarte/Fejer) hit the spectral gap (HYP-3785, blind to the pointwise atom) and the MOMENT (2nd-moment floor HYP-3571) + COVERING (lazy-cut) are the right instruments. Var(H)=W(n) (THM-589) and the LRC floor CV(N_R)^2 are the tournament/LRC moment instruments -- the SAME design lesson.
status: CONFIRMED (n=7 sample: U fails to determine |Aut|; n=6 exhaustive: the complementarity). A DESIGN PRINCIPLE + grounding, answering the owner: the next instrument for the fixed-point content must be a covering invariant or a moment, not another clever transform. Ties S89 (skew blind) + S90 (U fails at n=7) + the LRC Fourier/Delsarte spectral gap into ONE principle.
source: mac-mini-2026-07-01-S90
related:
  - HYP-3816   # S89 U sees flip-rank excess skew misses (but only at n<=6 -- this shows it FAILS at n=7)
  - HYP-3814   # S88 involution atlas: PARITY=|Fix|; the fixed-point face; GRAM = moment
  - HYP-3785   # LRC spectral gap: Fourier/Delsarte transforms blind to the atom (the fixed point)
  - HYP-3571   # LRC 2nd-moment floor = the moment (group-action) instrument that works
  - HYP-3798   # kappa = the covering invariant (group-action, fixed-point-sensitive)
  - THM-589    # Var(H)=W(n) <-> LRC CV(N_R)^2 (the moment instrument, both sides)
results:
  - 04-computation/fixed_point_instruments_macmini_20260701.py
  - 05-knowledge/results/fixed_point_instruments_macmini_20260701.out
---

# HYP-3817 -- fixed-point instruments must be built from the group action

Owner's insight: the next instrument has to be built to be SENSITIVE to fixed points, not blind to them by
symmetry -- a covering invariant, or Var(H), rather than another clever transform. This is correct, and it is
grounded and it ties three threads together.

## Why a transform is the wrong instrument for fixed points
A spectral transform (adjacency `A`, skew `A-A^T`, Seidel, `det(I+S)`) is computed to be **`S_n`-invariant** --
that is its whole purpose, to not depend on the labeling. But the fixed-point/symmetry content -- `|Aut|`, the
stabilizer, the flip-rank-excess needles -- is *about* how the labeling acts (the orbit vs the stabilizer). An
`S_n`-invariant object has **averaged over exactly the thing you want to measure**. So a transform is either
blind to `|Aut|` (the skew-spectrum, S89 -- it lumps `|Aut|=5,9` needles with `|Aut|=1`) or only accidentally
sensitive at small `n`:

- **(A) n=7: the U-spectrum FAILS to determine `|Aut|`.** In a 450-class sample, one U-spectrum holds BOTH
  `|Aut|=1` and `|Aut|=3`. The "U determines `|Aut|`" of S89 was small-`n` luck (leaked through score
  degeneracies); at `n=7` the leak closes. **No transform is a robust fixed-point detector.**

## The complementarity (n=6, exhaustive): two fine-structures, two instrument kinds
- **FIXED-POINT / symmetry structure** = the `15` high-`|Aut|` needles (`|Aut|` `{3:12, 5:2, 9:1}`), the
  flip-rank-excess drivers. Flagged by **group-action instruments** -- `|Aut|` (`= |Fix|` of the stabilizer,
  the Burnside/Lefschetz face S88), the covering radius `MFAS`, the covering number `kappa`. Transforms miss
  or only leak these.
- **GENERIC realization structure** = the `|Aut|=1` **spectral twins** (S86, `{12,22},{43,44}`): identical on
  every count, separated by the **adjacency (U) spectrum**. `|Aut|`/covering CANNOT separate them (both
  `|Aut|=1`); only the transform does.

So the instruments are **complementary, not ranked**: build the covering/moment instrument for the fixed-point
content, the transform for the generic content. Neither is "the next instrument" for both.

## The design principle
To measure a group's fixed points, use an instrument **built from the group's action**, not one **invariant
under it**:
- **Covering** (`kappa` = min transversal-subcube dimension; `R` = covering radius): DEFINED by how the
  `S_n`-orbits pack the cube -- inherently sees the few-rep (`high-|Aut|`) needles.
- **Moment** (`Var(H)` `= W(n)`, THM-589; orbit-weighted sums): summed over the labeled ensemble with orbit
  weights `n!/|Aut|` -- carries `|Aut|` structurally.
- **`|Aut|` itself** = the direct fixed-point count (`= |Fix|` of the fold, S88 PARITY-face).

## The LRC parallel (the same lesson, the whole point)
The LRC covering-min is a **fixed-point / rigidity problem**: the extremal lonely set is an ATOM (a fixed
point of the `iota`-fold), the construction is a symmetric deep well. THIS is why:
- the analytic **transforms** (Fourier, Delsarte, Fejer) hit the **spectral gap** (HYP-3785 / S54): they are
  averaging/PSD instruments, blind to the pointwise atom -- exactly the skew-spectrum's blindness to `|Aut|`;
- the instruments that WORK are the **moment** (the 2nd-moment floor `inf R' >= 1/(2 zeta(2))`, HYP-3571) and
  the **covering** (the lazy-cut, HYP-3782) -- both built from the group action.

`Var(H) = W(n)` (tournament) and `CV(N_R)^2` (LRC floor) are the SAME moment instrument (THM-589), the
group-action tool for a fixed-point problem. So the owner's tournament-side instrument principle IS the LRC
proof strategy: **for a fixed-point / rigidity extremum, reach for a covering invariant or a moment, never
another transform.**

## Next targets
1. Build a per-class **local moment** `Var_flip(H)` (variance of `H` over the arc-flip neighborhood) and test
   whether it detects the needles as a genuine moment instrument (complement to global `Var(H)=W(n)`).
2. Restate the LRC lower-bound program explicitly as "moment + covering, not transform": the moment floor
   (HYP-3571) supplies the positive bound, the covering lazy-cut the finite closure; the Fourier gap is the
   proof that no transform closes it.
3. Prove `PARITY = |Fix(fold)|` at the crux (SC-odd-grid-sym, HYP-3809) as the archetype fixed-point count.

## Honest scope
`(A)` and `(B)` are grounded (n=7 sample; n=6 exhaustive). The design principle is a strategic synthesis (it
organizes S89 + S86 + the LRC gap), well-supported but not a theorem. The claim is directional: transforms are
structurally symmetry-averaging; covering/moment instruments are group-action-built and fixed-point-sensitive.
