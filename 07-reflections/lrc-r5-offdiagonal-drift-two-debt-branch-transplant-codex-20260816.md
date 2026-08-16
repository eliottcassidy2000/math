# The typed branch transplant must puncture the diagonal, and that still does not create the K4 connection

**Status: FINITE-EXACT SUPPORT/REPRESENTATION SIDECAR, NOW INTEGRATED WITH
THE CONCURRENT COMMON-GAUGE POSITIVE; promotion recommended only as a scoped
stopping theorem.**  The exact companion is
`04-computation/lrc_r5_offdiagonal_drift_branch_transplant_probe_20260816.py`.
It reconstructs the independently audited THM-2594 source and U_full endpoint
bank.  It does not construct U_full ancestry, a physical current, an absolute
graph `H^1` class, a bispectrum, a scalar-row exclusion, or LRC(14).

## Verdict

The common torsor skeleton fixes the source and target drift coordinates up
to orientation:

```text
source:  s=u-q,
target:  d=b-a,
gauge:   u=a+c, q=b+c,
hence:   s=-d.                                         (1)
```

THM-2471 same-root disjointness makes the entire source fibre `s=0` vanish.
The Cartesian U_full endpoint bank instead has four nonzero `d=0` chamber
corners.  Therefore the full endpoint Cartesian pairing cannot itself be a
drift-preserving THM-2471 support relation.

Deleting `d=0` repairs this support mismatch exactly and preserves maximal
endpoint signal.  It does **not** repair the connection: the same augmented
projective system remains annihilator-only, and even the reduced three-channel
K4 potential has no map at all from the six source edge channels.

At the time of this probe the branch-transplant problem separated into two
independent debts:

```text
support debt:     remove/reweight the Cartesian same-sheet fibre;
connection debt: build the punctured chamber/guard incidence.   (2)
```

Concurrent work subsequently paid the first debt at the ancestry-service
level.  The exact companion
`lrc_r5_common_ancestry_guard_atom_root_drift_probe_20260816.py` evaluates the
same 39 guard atoms on the two linked THM-2471 nodes and imposes
`a-u=b-q` before marginalization.  It realizes precisely the 48 buckets
predicted here.  Pulling the actual endpoint pair function through that
address measure is nonzero at all thirteen owner characters, while its
post-marginal projective systems still have nullity zero.  Thus the remaining
debt is now temporal/Fubini realization of the nonlinear address pullback,
not discovery of the 48-bucket support relation.

## 1. The diagonal witness is forced by typing

The audit rebuilds the full common-base table and checks all

```text
13 roots x 7 cells x 3 nonempty deep windows = 273             (3)
```

same-root entries.  Every one is zero.  After summing only the lawful
indices, the root-difference support is

```text
ell=0:       0/13,
ell=1,...,6: 12/13 each,
total:       72=6*12.                                         (4)
```

This is the pointwise shadow of `A(y,u)F(y,u)=0` in THM-2471, not a numerical
cancellation after Fourier transform.

At the target, inversion of the four Walsh rows recovers the drift-zero
corner values

```text
LL  542919271526962761242
LR  564956562728457577520
RL  281582681509352363753
RR   79545704041698666254                              mod p, (5)
```

all nonzero.  Their sum is the independently audited same-sheet bridge

```text
324498447313453607031 mod p.                           (6)
```

No common translation can move this discrepancy to another drift: the
translation `c` cancels in both differences in (1).  Orientation reversal
changes `d` to `-d` and still fixes zero.  Thus deleting another drift is an
algebraic hostile, not a typed alternative.

## 2. The punctured carrier has maximal spectral signal

Set all four target corner values at `d=0` to zero and leave every other
bucket unchanged.  The resulting off-diagonal endpoint has

```text
corner support:          12/13 in each of four rows,
Walsh support:           12/13 in each of four rows,
drift Fourier support:   13/13 in each of four rows,
row rank:                 4.                             (7)
```

Hence same-sheet removal does not destroy endpoint weighted-tree spectral
closure.  It produces the clean common carrier

```text
source: 6 x F_13^*,       support 72,
target: 4 x F_13^*,       support 48.                    (8)
```

This is the strongest positive survivor.  The target still lives in endpoint
`B^1`; puncturing does not manufacture a graph cycle or absolute flux.

## 3. Support repair is not connection repair

For the six-dimensional source root-difference curve and four-dimensional
target Walsh curve, the augmented system retains all thirteen common drift
multipliers.  The exact certificate fields are

```text
(source rank, wedge rank, wedge nullity, wedge excess,
 augmented rank, augmented nullity,
 lambda-projection rank, nonannihilator entries).
```

Both the full and punctured targets give exactly

```text
(6,24,4,0,37,4,0,0).                                  (9)
```

Every one of the twelve nonzero torsor dilations gives (9).  Deleting any
single target drift also gives (9), while only deletion of zero satisfies the
typing in (1).  Thus the diagonal mismatch is real but is not the algebraic
cause of the connection obstruction.

## 4. Why six source channels are not automatically the six K4 edges

The source table has the exact simultaneous reflection

```text
B(ell,s)=B(-ell,-s)                                    (10)
```

on all `7*13` entries.  Its six nonzero `ell` channels can therefore be
drawn as a reversal-equipped six-edge carrier.  This matches several exact
combinatorial counts:

```text
|E(K4)|                    = 6,
|F_13^*/{+/-1}|            = 6,
|F_7^*/{+/-1}|             = 3 = number of K4 matchings. (11)
```

The count is useful, but a carrier is not a connection.  Remove the constant
Walsh channel from the punctured target, leaving the three reduced K4
potentials.  Allow an arbitrary fixed `3x6` map and one common drift operator.
The `39x31` augmented system has

```text
source rank 6, target rank 3,
wedge rank 18, wedge nullity 0,
augmented rank 31, augmented nullity 0.                (12)
```

There is no formal map—not even a source annihilator.  Since target edge
gradients are an invertible presentation of the same reduced potential
space, choosing a K4 incidence matrix, tournament section, or edge ordering
cannot evade (12).

This also clarifies the tournament language.  `F_13^*` contains both `d` and
`-d`, so the drift carrier is bidirected and loopless after puncturing, not a
tournament.  A tournament is a section choosing one orientation from each of
the six reversal pairs.  Such a section supplies gauge data but cannot alter
the zero-nullity certificate (12).

## 5. The successful support map and the remaining temporal debt

The concurrent common-gauge relation is not the Cartesian product with a
cosmetic filter.  On the unmarginalized coordinates

```text
(u,s,ell,theta; a,d,C,D),                              (13)
```

it satisfies both structural conditions isolated by this probe:

1. `s=-d` under one common root/sheet gauge, hence zero ancestry-service mass
   on `d=0`; and
2. a chamber-dependent Boolean coupling on `d!=0`, formed before the root and
   guard-sheet sums.

The realized common-gauge service has support `48/52`, rank four, and—after
retaining its common offset—full owner-frequency/Walsh/drift spectral support
at every nonzero owner character.  The separate endpoint-weight companion
then constructs the finite address pullback

```text
P_k=sum_(omega,nu,c) M_(omega,nu,c) E_(omega,nu) zeta^(-kc),       (14)
```

and finds all thirteen `P_k` nonzero.  This is exactly the nonlinear incidence
that equations (9) and (12) said could not factor through a fixed channel
map.  Indeed all of the new full and punctured `29`-unknown projective systems
have rank `29` and nullity zero.

What remains is not another support or spectral census.  The scalar endpoint
pair function `E_(omega,nu)` is still computed by a separate endpoint
integration and then pulled back on the finite address set.  A temporal/Fubini
theorem must lift the factors defining `E` onto the same linked-node fibre
before integration, or identify (14) as the pushforward of a lawful relation
current with chronology and the grouped coefficient.  Until then the positive
pullback is an auxiliary address-weighted contraction, not a physical current.

## Scope boundary

The off-diagonal projection in this script is a necessary typed hostile, not
by itself a U_full ancestry selector.  The separate common-gauge companion is
a lawful THM-2471 arrival-address refinement but still respects the theorem's
arrival/source temporal-atom warning.  The six-edge/K4 correspondence in (11)
is a representation atlas, not a canonical edge labelling or `H^1` class.  No
physical current, chronology, grouped coefficient, scalar exclusion, or
LRC(14) consequence follows.

## Reproduction

Run

```text
python -B 04-computation/lrc_r5_offdiagonal_drift_branch_transplant_probe_20260816.py
python -B -O 04-computation/lrc_r5_offdiagonal_drift_branch_transplant_probe_20260816.py
```

The pinned semantic digest is

```text
84a21447326ac6e86ff7743704653d1c280eee265cfdbe357735f5ad6d80915a.
```
