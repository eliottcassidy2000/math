---
id: THM-2624
title: "Two-clock root tomography and disjoint-carrier holotopy boundary"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.
  For each of THM-2614's 84 base cells, form the globally primitive integer
  matrix whose active target-section rows q record the twelve nonzero
  deep-probe-root weights.  Every one of the 144 mixed C13 characters is
  nonzero, but no one-clock matrix has full root rank: the exact rank
  histogram is 5^6,6^12,7^4,8^22,9^6,10^2,11^32.  The 32 rank-eleven
  kernels have two exact projective coordinate-pattern classes.  For every
  fixed source shift s, some pair of owner clocks stacks to rank twelve;
  60 of 84 adjacent cyclic pairs already do so, and consecutive windows have
  sharp width histogram 2^60,3^22,4^2.  Thus two charts give a rational
  signed left inverse.  They do not give physical descent: the clock strata
  are disjoint, every available row is strictly positive on all twelve roots,
  and hence no nonnegative left inverse exists.  A labelled guard wall cospan
  refines each chart but supplies neither a common clock overlap nor a
  positive root transport.  No ancestry connector, row exclusion, or proof
  of LRC(14) follows.
source: wild-holotopy-mining-2026-07-28-two-clock-tomography
depends_on:
  - THM-2614-punctured-target-root-cosupport-and-principal-deck-no-go
related:
  - THM-2356-finite-field-chirp-gram-tomography-and-bockstein-pairing
  - THM-2616-cross-time-target-future-diagonal-and-principal-action-no-go
  - THM-2622-affine-torsor-holonomy-fixed-section-spectrum-and-v4-c13-dictionary
script: 04-computation/lrc14_two_clock_root_tomography_thm2624.py
output: 05-knowledge/results/lrc14_two_clock_root_tomography_thm2624.out
script_sha256: dca5618da2bb7165e9cb84cd8af0583b4ba4c27272c7e6e1035c7c10d0d85a0b
output_sha256: 64dc998822421647d7d53762eb3568a314a97fc9de869d480ebb390fff07cc76
hash_basis: LF-normalized bytes
---

# THM-2624 -- two clocks see every root direction, but do not glue a root

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.**

THM-2614 finds maximal support: in each retained base cell the same-event
target/deep-probe relation is `Q_(s,ell) x F_13^*`.  Support, however, forgets
the weights.  Restoring those weights reveals a sharp intermediate object.
No single clock chart determines all twelve root coordinates, even though
every mixed Fourier character survives.  Two clock charts do determine them
linearly.  The inverse is necessarily signed and the two charts live on
disjoint physical strata, so this is tomography rather than transport.

That distinction is the finite holotopy boundary needed after THM-2614 and
THM-2616: local separation of states is not descent of states.

## 1. The canonical target/root weight matrices

Use THM-2614's exact notation.  Its globally primitive content is

```text
G=4,244,240.                                                (1)
```

For a base cell `c=(s,ell_4)`, let `E_c` be its one or two retained rails and
let `Q_c` be its active unit target-section set.  Define the integer matrix

```text
W_c(q,r)
 =G^(-1) sum_(e in E_c, e unit at q) sum_(ell_5=0)^6
             A_(e,q;ell_5,r),

q in Q_c,                    r in F_13^*={1,...,12}.         (2)
```

Here `A` is the exact nonnegative common-`x` numerator rebuilt in THM-2614.
Division by `G` is entrywise integral.  The aggregate unit predicate still
belongs to the whole `(ell_5,r)` slice; equation (2) does not assign it to a
hidden Perron sheet.

THM-2614's product support has a useful weighted consequence:

```text
W_c(q,r)>0 for every q in Q_c and every r in F_13^*.        (3)
```

There are exactly `1,080` such strictly positive rows over all `84` cells.

## 2. Every mixed character survives, but every one-clock matrix is singular

For `alpha,beta in F_13^*`, form

```text
F_c(alpha,beta)
 =sum_(q in Q_c) sum_(r=1)^12 W_c(q,r) zeta_13^(alpha q+beta r).
                                                                    (4)
```

Exact reduction modulo
`Phi_13(X)=1+X+...+X^12` gives

```text
F_c(alpha,beta) !=0
for all 84*12*12=12,096 choices.                            (5)
```

Nevertheless the exact rational ranks are

```text
rank          5   6   7   8   9  10  11
cells         6  12   4  22   6   2  32.                  (6)
```

In particular no `W_c` has column rank twelve.  Thus even simultaneous
nonvanishing of every mixed cyclotomic character does not imply linear root
reconstruction.

The `32` rank-eleven kernels split into two exact coordinate-pattern classes,
each on `16` cells.  Under the following **display compression only** --
primitive projectivization, reversal, and cyclic rotation of the printed
ordered column list `(1,...,12)` -- representatives are

```text
v_small=(0,0,0,0,0,0,0,0,0,1,0,-1),

v_large=(0,0,0,0,0,84919023,-829466001,7730416511,
         829466001,-15630671068,0,7815335534).             (7)
```

The compression in (7) is not asserted to be a physical `C_13` deck
symmetry; it only records the exact finite atlas compactly.  The small class
has primitive height one and the large class height `15,630,671,068`.

## 3. Two-clock signed tomography is exact and sharp

Fix `s` and stack two clock matrices vertically:

```text
W_(s;ell,ell') = W_(s,ell) direct-sum W_(s,ell').          (8)
```

For every `s in F_13^*`, at least one pair `{ell,ell'}` has exact rational
rank twelve.  Since no one-clock matrix has rank twelve, the minimum number
of clock charts is exactly two for every fixed `s`.

Equivalently, with

```text
K_(s,ell)=ker_Q W_(s,ell) subset Q^12,                     (9)
```

a winning pair is exactly a transverse pair

```text
K_(s,ell) intersect K_(s,ell')={0}.                        (10)
```

The `12` labelled two-clock graphs on vertex set `F_7` have only three exact
types.  For

```text
s in {1,3,4,5,8,9,10,12},                                 (11)
```

there are `18` winning pairs.  For each of the two classes

```text
{2,7},                 {6,11},                             (12)
```

there are `6` winning pairs, with the exact edge lists printed by the
companion.  Over all `12*binom(7,2)=252` fixed-`s` clock pairs, the exact
rank histogram is

```text
rank                         8   9  10  11  12
pairs                        6  14  26  38 168.            (13)
```

In the cyclic adjacent subatlas,

```text
rank of W_(s;ell,ell+1)   9  10  11  12
oriented pairs             8   6  10  60.                 (14)
```

Starting from a labelled clock and adding consecutive cyclic clocks, the
first full-rank width has histogram

```text
width                       2  3  4
starting cells              60 22  2,                     (15)
```

and width four occurs exactly at `(s,ell)=(2,5),(11,0)`.

Full column rank in (8) gives a rational left inverse

```text
L_(s;ell,ell') W_(s;ell,ell')=I_12.                        (16)
```

This is the precise positive result: if one independently supplies a single
root vector that is common to both charts, its two response tables determine
it exactly.

## 4. The inverse is necessarily signed

Equation (3) makes the first obstruction immediate.  A nonzero nonnegative
linear combination of available rows is strictly positive in all twelve root
coordinates.  It cannot equal a coordinate row of `I_12`, which has eleven
zeros.  Therefore

```text
no stack of the available W_c has a nonnegative left inverse.              (17)
```

In particular every left inverse in (16) uses cancellation.  It is not a
stochastic kernel, a Markov transport, a Boolean owner selector, or a positive
physical correspondence.  This remains true even if all seven owner clocks
are stacked.

This is stronger than merely saying that the inverse is noncanonical: its
sign is an invariant obstruction to interpreting the reconstruction as a
positive root transition on the present carrier.

## 5. Why disjoint clock charts do not form descent data

The clock label `ell_4` in THM-2614 selects disjoint owner-clock strata of the
physical base.  The direct sum in (8) therefore has no common physical overlap
on which to compare two root origins.  It proves injectivity of the formal map

```text
x |-> (W_(s,ell)x, W_(s,ell')x)                            (18)
```

for one vector `x` *assumed beforehand* to be common.  It does not prove that
the two physical strata carry the same `x`, nor provide a transition map that
could make them so.

In Cech language, the detection presheaf is separated after two charts, but
the physical atlas has no supplied overlap isomorphism and hence no descent
datum.  Empty intersections make compatibility vacuous; they do not identify
fibres.  In THM-2622's affine-holonomy language, (16) recovers coordinates
after a frame has been imposed, while the affine translation/cocycle between
frames remains absent.

This is the promised holotopy distinction:

```text
two-chart coefficient tomography:  PROVED, exact, signed;
common-carrier root transport:      NOT SUPPLIED.          (19)
```

This is the linear root analogue of THM-2356 and MISTAKE-261.  There, planar
chirp data reconstruct a refined graph signal but independent graph phases
prevent descent to the coarse physical target.  Here, two clocks reconstruct
a root vector only after one declares it common across the clock label; the
missing clock-overlap map is the corresponding phase/reference debt.  In both
cases tomography is valid on the refined carrier and invalid after forgetting
the label that made the reconstruction meaningful.

## 6. Comparison with the guard wall cospan

THM-2616's delayed guard-safe word and its guard-danger complement are
disjoint Boolean sectors with an exact guard-omitted common parent before the
unit test.  This is a genuine wall cospan *inside each fixed clock chart*.
It does not change the conclusion above.

First, the wall bit refines a fixed `(s,ell_4)` stratum; it creates no overlap
between different `ell_4` strata in (8).  Second, every nonzero fixed
`(rail,q,ell_5,h=q)` row in either labelled sector still has full support on
all `r=1,...,12`.  Hence the nonnegative-left-inverse obstruction (17)
survives sectorwise.  Third, raw Boolean tensors add across the cospan, but
the predicate "the seven-clock Bockstein is a unit" is nonlinear.  A union of
unit-support atlases is not the unit atlas of the summed parent tensor.

Thus a common wall-edge name is not yet affine descent data.  A successful
sidecar must supply at least one of:

1. a positive root-atomic row or selector, rather than another full-support
   positive row;
2. a lawful chronological overlap map identifying the same root state on two
   clock strata; or
3. an affine transition cocycle whose translations are physically typed and
   satisfy the cycle law.

The guard wall is useful retained structure, but carrier absence remains
absolute for the present two-clock stack.

## 7. Exact evidence and scope

Run

```text
python 04-computation/lrc14_two_clock_root_tomography_thm2624.py
python -O 04-computation/lrc14_two_clock_root_tomography_thm2624.py
```

The companion rebuilds THM-2614's complete exact bank, divides only by its one
global primitive content, checks all `12,096` mixed characters by exact
cyclotomic reduction, and uses `Fraction` row reduction for every rank,
kernel, two-clock edge, and consecutive window.  It checks every logical
claim with explicit optimized-mode guards.  Normal and optimized executions
must byte-match the stored transcript after LF normalization.

The theorem concerns the one canonical typed row and the collapsed integer
matrices (2).  It does not assign unitness to a hidden sheet, retain one old
root across clocks, build a principal `C_13` action, identify THM-2613's local
root with THM-2585's next target, preserve semantic owner/repair provenance,
exclude a scalar profile, or prove LRC(14).

QED (candidate; independent hostile audit pending).
