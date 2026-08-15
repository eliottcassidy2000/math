# AMM 12592 cross-frontier: an exact fractional vertex and the finite golden profile cell

Status: **PROVED structural lemma + FINITE-EXACT hostile computation.**  This
does not prove an integer point at `R=512`, an all-epoch construction, a new
deadline bound, or the general golden lower bound.

## Current truth surface at `origin/main e6ec5465c`

The newest AMM chain has separated three questions which older notes often
blended.

1. [THM-3340](../01-canon/theorems/THM-3340-single-donor-cyclic-rotation-proves-all-pointwise-AMM-floors.md)
   proves the pointwise optimum `T_opt(n)=n+1` for every `n`, but the
   extractor depends on the target horizon.
2. [THM-3342](../01-canon/theorems/THM-3342-sublinear-deadline-excess-is-impossible-for-fair-critical-run-extractors.md)
   proves that no one fixed extractor has `T(n)=n+o(n)`.  This is
   nonattainment at slope one, not a positive uniform gap above one.
3. [THM-3343](../01-canon/theorems/THM-3343-shifted-donor-rotation-bisects-exactly-the-dyadic-annuli.md)
   and [THM-3344](../01-canon/theorems/THM-3344-orientation-splitting-saves-exactly-one-dyadic-donor-bit.md)
   give one exact uniform donor construction, with every nonpower at its
   pointwise floor and sparse dyadic donors at `2n-1`.

The closest proved mechanism to the kps-S169 real LP is instead
[THM-3329](../01-canon/theorems/THM-3329-r128-triple-closure-attractor-rule-and-golden-attainment-to-1023.md):
at dyadic epochs the substitution to halved transportation variables makes
parity automatic, while the evaluation-at-one ballot cone remains a genuine
global sidecar.  [THM-3332](../01-canon/theorems/THM-3332-s-cone-one-row-certificate-32768-closures-and-the-256-floor-point.md)
supplies the nearest positive control: the exact-floor `R=256` polytope has
an integer point.  The inherited hostile is that plain clamp dies even on
known-feasible controls; policy death is not polytope infeasibility.

Thus the exact question complementary to the sparse real LP is not “does
parity round?”  At dyadic `R`, parity has already disappeared.  It is:

> Is the integer lattice section of the active sparse polytope forced by
> real feasibility together with the ballot cone?

The computation below does **not** answer that existence question.  It gives
a sharp negative only to the stronger shortcut “the whole real polytope is
integral / every vertex rounds automatically.”

## Exact hostile: the sparse polytope is nonintegral already at `R=8`

Put `y=j/2` in the sparse junk-flow system of kps-S169.  At a dyadic epoch
the initial load and every feed coefficient are even, so all coefficients,
bounds and right-hand sides in `y` are integers.  For `R=8` the exact golden
profile is

```text
(4,5,5,6,7,7,8,8).
```

The resulting polytope has `42` variables and `106` displayed causal and
terminal inequalities.  The companion script verifies a feasible point
with `20` nonintegral coordinates, all of denominator `103`.  It has `45`
active constraints.  A displayed subset of `42` active integer rows has

```text
det(B)=1648=16*103 != 0.
```

Therefore those rows meet in exactly one point, so the feasible point is a
vertex.  It is fractional.  Consequently:

* the sparse real polytope is not integral, even at `R=8`;
* the full displayed constraint matrix, including coordinate-bound rows,
  is not totally unimodular (it has a square minor of determinant `1648`);
* an exact active-face solve need not return an integer point; and
* this failure coexists with the already known integer witness at `R=8`.

The last item is the important hostile control.  Nonintegrality of the
polytope does **not** mean integer infeasibility.  It means that generic
vertex selection, rational reconstruction, and total-unimodularity
arguments cannot supply the missing lift.  A successful lift must exploit a
typed entry section, a circuit-cancellation operation, or the ballot sidecar
rather than arbitrary LP vertices.

The prime `103` comes from this active lattice basis.  It has no visible
relation to the golden constant.  This sharply separates the finite
integrality obstruction from the asymptotic entropy constant.

## Exact role of `C=1+2 log(phi)/log(5)`

Write

```text
gamma_* = C-1 = log_5(phi^2).
```

For a fixed epoch `R`, every coefficient and bound of the sparse polytope
depends on `gamma_*` only through the finite degree word

```text
d_n=floor(gamma_* n),                 R<=n<2R.
```

This gives an elementary profile-cell lemma.  For any prescribed word
`(d_n)`, the set of slopes producing it is exactly

```text
[ max_n d_n/n,  min_n (d_n+1)/n ).                    (1)
```

The endpoints are rational.  Hence no finite epoch can detect whether its
slope is logarithmic, algebraic, or rational; it detects only the cell (1).

At `R=512`, exact Fibonacci--Lucas comparisons with
`phi^(2m)=(L_(2m)+sqrt(5)F_(2m))/2` give

```text
534/893 < gamma_* < 119/199,

profile cell = [534/893,119/199),
width        = 1/177707.
```

The binders are `n=893` and `n=597`.  Moreover

```text
119*893-534*199=1,
```

so the endpoints are Farey neighbours.  The standard determinant-one
lemma says that every reduced rational strictly between `a/b<c/d` has
denominator at least `b+d`; its proof writes

```text
bx-ay >= 1,       cy-dx >= 1
```

for `a/b<x/y<c/d`, and then uses
`y=(bc-ad)y=b(cy-dx)+d(bx-ay)>=b+d`.  Thus the smallest
denominator interior replacement is the mediant

```text
(534+119)/(893+199)=653/1092.
```

The script independently exhausts every smaller denominator and verifies

```text
floor((653/1092)n)=floor(gamma_* n),       512<=n<1024.
```

Therefore the `R=512` sparse LP built with the golden logarithm is literally
the same integer matrix, right-hand side and bounds as the LP built with the
rational slope `653/1092`.

This pins the golden constant's role precisely:

* **finite epoch:** it only selects a rational Farey cell / mechanical word;
* **all scales:** the nested cells shrink to the slope because
  `d_n/n <= gamma_* < d_n/n+1/n`;
* **continuum threshold:** THM-3009's entropy tangency
  `delta^2=1-delta` gives `delta=1/phi` and then
  `gamma_*=2 log(phi)/log(5)`.

So the golden constant is not merely an analogy, but neither is it a
coefficient-field obstruction inside the finite LP.  It is the exact
all-scale entropy/tangency rate.  The finite lifting problem is an integer
entry/lattice-section problem whose first hostile vertex denominator is
instead `103`.

## Sharp next target

The fractional vertex removes “prove the whole matrix is TU/integral” from
the live list.  THM-3371 now replaces the former **distinguished-face** target
by a more faithful **entry-section compiler**:

> find an integer causal prefix whose first feed-free state lies in the
> linear THM-3332 cone; the theorem then constructs the tail.

At `R=8`, the exact integer point has `25` active constraints of rank `22`.
Only `8` hostile-vertex active rows, of rank `7`, remain active there; the
other `37` become strict.  Its state entering the first feed-free row is
`(-2,0,0,0,0)`, which satisfies `F1--F3` and captures in one clamp.  Thus the
positive certificate leaves the hostile vertex's minimal face; no circuit
inside that face can reach it, and no theorem says the entry section itself
is a face.

At `R=512`, the same reframe is computationally material.  The last feed row
is `129`, the first feed-free row is `130`, and the degree is `383`.  Searching
for an integer prefix landing in the divided `F1--F3` cone uses `44,750`
variables and `89,146` causal inequalities, versus `234,117` variables in the
full sparse system.  The rational profile selector `653/1092` removes all
transcendental arithmetic from this finite sufficient search.  Feasibility of
that prefix cone would close the exact-floor `R=512` case; infeasibility would
exclude only first-feed-free cone entry, not every witness.

THM-3373 supplies the first exact operation law between the two `R=8` points.
After adjoining all inequality slacks, their denominator-cleared displacement
has conformal causal-row locality width exactly five: integer Farkas
certificates rule out widths at most four even over the reals, while twelve
integer width-five kernel moves give a monotone path in the `103`-dilated
polytope.  This is more rigid than “find a circuit” and less than a Graver
compiler.  The next lawful transfer test is to translate the twelve labelled
atom shapes through the two `delta=0,1` Pascal steps at `R=512`, retaining cap
slacks and the entry margin; the naked kernel vectors are insufficient.

THM-3374 supplies the hostile transfer test.  The exact `R=512,D0=0` rule-A
trajectory reaches its row-107 death with a 277-bit fatal constant.  The total
coefficient capacity of arbitrary admissible changes in the five immediately
preceding rows is only 42 bits.  More generally, the capacity of the preceding
57 rows is still too small, while the coarse bound first becomes inconclusive
at 58 rows.  Therefore any successful prefix must already differ from rule A
by row 49.  Width five remains an operation scale, but it is not an end-of-life
repair scale: the live experiment must translate and compose such moves early,
while slack is still being accumulated.

## Reproduction

```bash
python 04-computation/amm12592_crossfrontier_20260812_exact_vertex_and_golden_cell.py
python -O 04-computation/amm12592_crossfrontier_20260812_exact_vertex_and_golden_cell.py
```

Expected output:
`05-knowledge/results/amm12592_crossfrontier_20260812_exact_vertex_and_golden_cell.out`.
The script uses only the Python standard library and exact integer/Fraction
arithmetic.  No solver output is trusted or imported.
