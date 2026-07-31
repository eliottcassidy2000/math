# The integrated-row service law independently audits THM-2687

**Status: FINITE-EXACT ALTERNATE PROOF + INDEPENDENT AUDIT.**  This reflection
does not create a new theorem ID.  THM-2687's canonical proof uses a finer
future-coordinate envelope; the argument below is a simpler integrated-row
projection and is valid only in the necessary no-go direction.

## Statement

Use the physical THM-2640 carrier with

```text
p=13,       R=13^6,       S=13^5,       tau_delta=7 delta/R.
```

For each target label `delta in F_13`, allow an arbitrary THM-2640
configuration

```text
(rail j, delayed sector sigma, deep half-edge epsilon,
 predecessor carry c, future half-digit (h,kappa)),
```

and allow the rail, sector, edge, and seven-clock row to depend on `delta`.
Require only that the pulled-back chart is a retained THM-2640 unit chart at
`x+tau_delta`.  Then no positive-measure common intersection of all thirteen
labels exists.  More precisely, the unrestricted union-labelled
slope-seven chart family has maximum positive-overlap multiplicity exactly
`12`, hence positive-overlap nerve dimension exactly `11`.

The upper bound permits all `162` rails, both delayed sectors, both deep
half-edges, all base cells, and an unrelated configuration at every label.
The lower bound is the explicit positive twelve-fold physical component in
THM-2672, Section 3.

## 1. Pointwise invariants of the slope-seven lift

Suppose a common open component existed and choose a point away from its
finite endpoint set.  Put

```text
y={R x},
d=floor(26y)=2h+kappa,       0<=h<=12, kappa in {0,1},
c_0=floor(Rx) mod 13.
```

For every label,

```text
R(x+tau_delta)=Rx+7delta.
```

Consequently `y`, `d`, `h`, and `kappa` are identical at all thirteen
labels, while the predecessor carry is

```text
c_delta=c_0+7delta mod 13.                                  (1)
```

Since `7` is a unit modulo `13`, `(c_delta)_delta` is all of `F_13`.
Thus arbitrary configuration switching can help only if, at this one fixed
physical half-digit `d`, the union of every THM-2640 configuration services
all thirteen carries.

## 2. Exact full-configuration service law

For fixed `(h,kappa,c)`, define `c in M_(h,kappa)^raw` when at least one of
the `162*2*2` seven-clock raw rows

```text
(j,sigma,epsilon,c,kappa,h)
```

is nonzero.  Define `M_(h,kappa)^unit` analogously using the THM-2640
cyclotomic unit test.  The exhaustive carrier calculation gives equality of
these two sets and the following exact law, written in terms of
`d=2h+kappa`:

```text
d=0:       F_13 \ {0}
d=1..11:   F_13 \ {0,6}
d=12..13:  F_13 \ {6}
d=14..24:  F_13 \ {6,12}
d=25:      F_13 \ {12}.                                  (2)
```

Hence the all-configuration service-size histogram over the twenty-six
halfdigits is

```text
11^22 12^4.                                                (3)
```

In particular, every `M_(h,kappa)` is proper.  The missing carries are
always among the three seam values `{0,6,12}`.  Equality of raw and unit
service means the obstruction precedes cyclotomic cancellation: at the
missing seam carry there is no raw THM-2640 row in the full configuration
union.  Each raw clock entry is the exact integral of a nonnegative physical
Boolean product.  Therefore a positive open chart component would make at
least one raw clock entry positive.  This gives the necessary implication
used below; no converse from an integrated raw row to a fixed `y`-component
is asserted.

## 3. Exact audit universe and an independent unit classifier

The referee rebuilds the four canonical THM-2640 shards and their global
primitive content `26`.  Its universe has exactly

```text
162 rails * 2 sectors * 2 edges * 13 carries
          * 2 kappa values * 13 h values = 219,024         (4)
```

seven-clock rows.  It does not reuse the canonical multiplication
determinant as its only unit test.  If a row is `(A_0,...,A_6)` and its
private root is `r`, it forms

```text
Y(z)=sum_(ell=0)^6 (A_ell/26) r^(-1) z^ell
```

over `F_13`, reduces `Y` by
`Phi_7(z)=1+z+...+z^6`, and computes the Euclidean gcd with `Phi_7`.
All `219,024` classifications agree with the canonical determinant flags.
The exact row counts are

```text
raw nonzero seven-clock rows = 92,696,
unit seven-clock rows        = 80,955.                     (5)
```

Their configuration unions nevertheless have the same missing table (2).

## 4. Exclusion and exactness

At the common point, (1) visits every carry.  For its fixed half-digit `d`,
(2) supplies a carry `c_*` with no raw row in any configuration.  The unique
label satisfying `c_delta=c_*` therefore has no THM-2640 chart at that point,
contradicting thirteen-fold intersection.  This proves multiplicity at most
`12` without invoking pairwise overlaps, gain balance, or a fixed-configuration
twelve-subface.

THM-2672 gives an honest positive open component for twelve slope-seven
labels, with all physical factors and one common delayed clock retained.
Therefore the maximum is exactly `12` and the unrestricted union-labelled
positive-overlap nerve has dimension exactly `11`.

## 5. Rainbow and holotopy boundary

The proposed THM-2648/2656 two-rainbow bridge fails before physical
selection.  At fixed `(h,kappa)` there are `648` possible
`(j,sigma,epsilon)` configurations, but even their total service union is
the proper set (2).  Exhausting all two-configuration pairs gives

```text
full thirteen-carry pairs = 0,
best pair service histogram over halfdigits = 11^22 12^4. (6)
```

Thus there is no source relation on which a same-base rainbow selector could
act: all configurations share at least one seam puncture.  The first failed
predicate is raw service, not rail overlap, present-factor overlap, delayed
word compatibility, integer gain balance, or connected-component Helly.
This is precisely when the component-refined THM-2658 machinery should *not*
be invoked: a lower-dimensional necessary projection is already empty.

## 6. Scope

This result closes only the full union of the promoted THM-2640 charts under
the canonical slope-seven lifts.  It does not exclude a different physical
translation, a broader carrier, a nonunit/signed coefficient construction,
or the `165` live LRC rows.  It supplies no semantic owner, endpoint current,
component transition, physical sphere, row exclusion, or LRC(14) conclusion.
It also does not assert that every twelve-label subset intersects; it asserts
only the sharp maximum, with one lower witness supplied by THM-2672.

## 7. Reproduction

From the repository root run

```bash
python 04-computation/lrc14_slope7_global_configuration_switching_integrated_referee.py
python -O 04-computation/lrc14_slope7_global_configuration_switching_integrated_referee.py
```

Both modes must byte-match
`05-knowledge/results/lrc14_slope7_global_configuration_switching_integrated_referee.out`.  The completed normal and
optimized full-carrier rebuilds match each other and the stored transcript
after LF normalization.  The LF-normalized SHA-256 hashes are

```text
script  fe92e886874ed10cebf5d877d4e45d7a4c72dbc60c3050b2d5cc775aa047ce2e
output  e16bb1f643e6c3ec85b494de9b24fec57db3548cc247e6433ea95c5721f0512f
```

QED.
