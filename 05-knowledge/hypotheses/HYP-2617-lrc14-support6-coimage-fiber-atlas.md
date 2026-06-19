---
id: HYP-2617
title: LRC(14) support-6 coimage fiber atlas
status: OPEN
source: codex-2026-06-19-S14
depends_on:
  - THM-538
  - HYP-2608
  - HYP-2614
  - HYP-2615
related:
  - HYP-2616
  - HYP-2612
  - HYP-2613
  - MISTAKE-078
  - OPEN-Q-108
---

# HYP-2617 - LRC(14) Support-6 Coimage Fiber Atlas

## Claim

The exact support-6 LRC(14) tail has a finite mod-7 coimage atlas sitting
between the large relation lattice and the reciprocal analytic tail.

For a six-speed support with speed residues

```text
a = (a_1,...,a_6),   a_i = e_i mod 7,
```

and nonzero Fourier coefficient residues

```text
r = (r_1,...,r_6) in (F_7^*)^6,
```

the mod-7 shadow of the relation hyperplane is

```text
a_1 r_1 + ... + a_6 r_6 = 0 mod 7.
```

The leading coimage fiber coefficient is therefore

```text
S_d(a) = sum_{r in (F_7^*)^6, a.r=0} C_d(r),
```

where `C_d` is the residue coefficient from HYP-2614.  The speed residues
`a_i` are allowed to be `0`.  This is essential: a speed divisible by `7` is
still Fourier-live when its coefficient residue is nonzero, but that coordinate
is invisible to the mod-7 relation.  This is the concrete coimage degeneration
behind the user's prompt.

Modulo scalar multiplication by `F_7^*` and coordinate permutation, there are
only `159` projective speed-residue coimage classes.  Any final support-6
reciprocal-tail theorem can target this finite address table rather than an
unstructured family of raw supports.

## Computation

Script:

- `04-computation/lrc14_support6_coimage_fiber_codex_s14.py`
- output: `05-knowledge/results/lrc14_support6_coimage_fiber_codex_s14.out`

The class inventory is:

```text
projective classes modulo F_7^* and S_6: 159
zero-speed-residue histogram: {0:80, 1:42, 2:22, 3:10, 4:4, 5:1}
```

The largest coimage fibers by ambient dimension are:

```text
d=6   max |S_d| = 2.1852492   class (1,1,2,2,3,3)
d=7   max |S_d| = 1.4048031   class (1,1,2,2,3,3)
d=8   max |S_d| = 0.82504308  class (1,1,2,2,3,3)
d=9   max |S_d| = 0.43959824  class (1,1,2,2,3,3)
d=10  max |S_d| = 0.27850323  class (1,1,1,1,4,4)
d=11  max |S_d| = 0.26426602  class (1,1,1,1,2,2)
d=12  max |S_d| = 0.22155439  class (1,1,2,2,4,4)
d=13  max |S_d| = 0.2106327   class (1,1,1,1,1,1)
```

Many classes are coimage-null or coimage-small:

```text
d     <1e-12  <0.001  <0.01  <0.1
6        142      142     142   142
7        113      113     113   121
8         80       80      80   118
9         43       43      77   117
10        34       34      43   140
11         3        3      35   148
12         3        3      40   149
13         3        3      25   150
```

The three classes null for every `d=6..13` are:

```text
(0,0,0,0,0,1)
(0,0,0,1,1,1)
(0,1,1,1,1,1)
```

## Named Support Readout

The S12 hard supports land in the coimage atlas as follows:

```text
AP / dissociated 211, d=7:
class (1,2,3,4,5,6), z=0,
S_d=-0.93653539, abs/signed=15.9731.

resonant 21 support, d=7:
class (0,1,2,3,4,5), z=1,
S_d=-0.11706692, abs/signed=125.46.

k=9 wide 68 support, d=8:
class (1,1,2,4,5,6), z=0,
S_d=-0.21741, abs/signed=38.0849.

k=10 wall 22 support, d=9:
class (0,1,1,1,2,4), z=1,
S_d approximately 0, abs/signed approximately 1.13e16.
```

The last row is the strongest new clue.  The `k=10` height-one wall has large
absolute fiber mass, but its leading coimage fiber is numerically zero in its
own ambient dimension `d=9`.  That means the wall's danger belongs in the
finite low-height ledger, not in the infinite analytic tail.  This is exactly
compatible with HYP-2616, which independently clears the one-large height-one
wall ledger.

## Interpretation

HYP-2614 showed that the exact support-6 term is a residue-addressed reciprocal
hyperplane sum.  HYP-2615 reframed the large absolute/small signed split as a
coimage phenomenon.  HYP-2617 makes that coimage finite:

```text
raw relation volume
-> projective speed-residue class
-> finite mod-7 coimage fiber
-> reciprocal tail by class
```

The quotient preserves the analytic predicate "support-6 correction after
low-height deletion" but discards witness-time geometry.  That is the right
loss for this proof obligation: the remaining problem is not to find a lonely
time, but to show the signed correction cannot cross the cap margin.

## Tournament Analysis

The tournament vertices are quotient stages, not runners, arcs, or residue
tuples:

- raw relation volume;
- speed-residue projectivization;
- coimage fiber sum;
- fixed-boundary residue;
- named wall nullity;
- low-height wall ledger;
- reciprocal tail estimate.

The proof-quotient tournament is transitive, with Hamiltonian path:

```text
named_wall_nullity
> coimage_fiber_sum
> low_height_wall_ledger
> fixed_boundary_residue
> speed_residue_projectivization
> reciprocal_tail_estimate
> raw_relation_volume.
```

Fingerprints:

```text
score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1}
directed_3_cycles = 0
```

Challenged assumption: the relevant tournament vertices are not runners.  They
are proof-obligation quotient stages.  The quotient destroys time-location
data but preserves exactly the support-6 signed-tail predicate.

## Proof Target

HYP-2617 does not prove LRC(14).  Combined with HYP-2616's height-one ledger,
it narrows the next theorem to:

```text
delete low-height wall ledger
then bound the non-null projective coimage classes
by a signed reciprocal-tail estimate.
```

The most promising route is a finite class-by-class cotangent/Dedekind or
summation-by-parts bound over `sum e_i n_i=0`, with the coimage-null and
coimage-small classes removed from the analytic tail and recorded in the finite
wall ledger.
