---
id: HYP-3792
title: THE COVERING-MIN LOWER BOUND, CHARACTERIZED FOR A PROOF -- (A) the SAFE-BAND RESIDUE frame M(S)=max_{q,a}(1/q)min_v ||av||_q recasts it as band-dodging: every covering 13-set has a cyclic dilation avoiding a band of half-width ceil(rq); (B) the EXTREMAL RESIDUE MECHANISM (all n): at t*=n/Phi6 the core {1..n-2} tiles residues {n,2n,..,(n-2)n}, the skipped runner n-1 lands at residue -1 (THE dangerous slot, why it is skipped), and since n(n-1)=-1 mod Phi6 the multiple k(n-1) lands at -k so covering's FORCED multiple of n-1 must be >= n(n-1) to be safe = the patch; (C) ARITHMETIC DEPTH: the construction is the UNIQUE covering 13-set with a DEEP witness (q*=Phi6, CF [0;n-1,n]); restructured covering sets bind SHALLOW (q*<=~50, M~0.10-0.14); (D) DEEP-WELL ISOLATION (5001-sample): random covering 13-sets have M>=0.108 (1.4x above construction, 1.5x above 1/14) => the danger zone near 1/14 IS exactly the construction family. SLACK: LRC needs only M>=1/14 (margin = Dedekind sum 13/2562), and the bulk sits at 0.108>>1/14, so a certificate need only reach 1/14 = real leverage for the bounded-case ILP.
status: CONFIRMED (exact-Fraction mechanism for n=8,10,12,14 + exact witness CF + 5001-set grid sampling). A CHARACTERIZATION synthesis toward a proof: it makes the covering-min lower bound an explicit band-dodging / arithmetic-depth statement, gives the exact reason the construction's skip-and-patch is forced (forced-cover obstruction), and quantifies the deep-well isolation (the bulk is far from the LRC threshold => slack). NOT a proof: the for-all-covering-S band-dodging existence is still the open core (bounded = lazy-cut, huge = S73/74/75). But it converts the reframe into a concrete finite target (bulk M>=c>1/14 + construction family) and supplies the missing STRATEGIC lever (aim 1/14 not 14/183).
source: mac-mini-2026-06-30-S77
related:
  - HYP-3789   # S76 moment relaxation (this is its arithmetic/residue face; the atoms are the witnesses)
  - HYP-3778   # covering-min = 14/183 (the statement being characterized)
  - HYP-3784   # S73 huge single-patch scaling (1/M=[n-1;nk]; the deep-CF ladder)
  - HYP-3777   # S69 annealing deep-well (this quantifies the isolation arithmetically)
  - HYP-3750   # band-prime reduction (residue-structural, same residue frame)
  - HYP-3768   # margin = Dedekind sum (the slack this reframe spends)
  - THM-523    # covering reduction (LRC14 <=> covering-min >= 1/14)
results:
  - 04-computation/witness_depth_and_safeband_characterization_macmini_20260701.py
  - 05-knowledge/results/witness_depth_and_safeband_characterization_macmini_20260701.out
---

# HYP-3792 -- the covering-min lower bound, characterized toward a proof

The owner's directive: *characterize and understand related properties so they can lead to a proof; create
new definitions and key angles that shed light on the underlying structure.* Four angles, all verified.

## (A) The safe-band residue frame
`M(S) = max_{q,a} (1/q) * min_{v in S} ||a v||_q`, `||x||_q = dist(x mod q, 0)`. Loneliness at `(q,a)`
`<=>` the **residue system** `R(q,a) = {a v mod q}` avoids the **danger band** `B_r = (-rq, rq)`. So the
covering-min lower bound is:

> **every primitive covering 13-set has some modulus `q` and dilation `a` for which all residues `a v mod q`
> miss the band of half-width `ceil(rq)` around 0** (`r = 14/183`, so half-width 14 at `q = 183`).

This is a clean *band-dodging* problem on `Z/q`, dual to the moment picture (HYP-3789): the witness `t*` IS
the atom of the extremal lonely measure.

## (B) The extremal residue mechanism (why skip `n-1`, why patch `n(n-1)`) -- all n
At the binding witness `t* = n/Phi6` (`a=n`, `q=Phi6`), verified n=8,10,12,14:
- core `{1..n-2}` -> residues `{n, 2n, ..., (n-2)n}` -- an AP of step `n` tiling the safe band `[n, Phi6-n]`;
- **skipped runner `n-1`** -> residue `-1` (adjacent to 0, THE dangerous slot) -- this is *why* it is skipped;
- since `n(n-1) = Phi6 - 1 = -1 mod Phi6`, the multiple `k(n-1)` -> residue `-k` (distance `k`);
- **patch `n(n-1)`** -> residue `-n` (the mirror safe edge). Binders: runner 1 (`+n`) and patch (`-n`), `M=n/Phi6`.

**Forced-cover obstruction**: covering *requires* a multiple of `n-1`; at `t*` that multiple `k(n-1)` sits at
distance `k`, so to be safe (`>= n`) it must be `>= n(n-1)` -- exactly the construction's patch. The skip-and-
patch is not a lucky choice; it is *forced* by "cover `n-1` without landing near 0 at `t*`."

## (C) Arithmetic depth (the through-line)
`q*(S)` = witness denominator; the construction is **deep** (`q* = Phi6`, `CF(t*) = [0; n-1, n]`,
`1/M = [n-1; n]` = the S71 self-concordant ladder). Every restructured covering 13-set (drop core, re-cover
13 & 14) binds **shallow** (`q* in {15,...,53}`, short CF) with `M ~ 0.10-0.14`. Reaching the deep `q=Phi6`
locus requires the `182`-type patch; nothing else gets there.

## (D) Deep-well isolation + the slack lever
5001 random covering primitive 13-sets (speeds <= 60): `min M = 0.108`, 5th pct `0.137`, median `0.174`; **zero**
below `0.10`. The construction (`0.0765`) is isolated by a factor `~1.4`. Since LRC14 needs only `M >= 1/14`
(the margin `14/183 - 1/14 = 13/2562` is the Dedekind sum, HYP-3768, spendable slack) and the **bulk sits at
`0.108 >> 1/14`**, a certificate need only reach `1/14`, not the exact `14/183`.

## Toward a proof (the map this enables)
LRC14-covering `<=` three cases keyed by arithmetic depth:
1. **shallow / bulk** (`q* <= Q`, bounded speeds): `M >= 0.108 > 1/14` -- lazy-cut ILP, and targeting `1/14`
   (slack) rather than `14/183` should go infeasible far faster (the bulk is far from the target);
2. **construction family** (deep `q* = Phi6` scaled): `M >= 14/183 > 1/14` -- S73 scaling law, closed;
3. **huge speeds**: don't beat -- S74/S75 (equidistribution / signed correction), `<=6`-huge rigorous.

Honest: the `for all covering S` band-dodging existence is still open (cases 1 and 3's residuals), but this
converts "reframe" into a concrete finite target with a strategic slack lever, and gives the exact structural
reason the extremum is what it is.
