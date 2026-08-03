# LRC(14) `z1=216`: threshold-layer multicovers and support-minimal status circuits

**Status:** structural synthesis for
[THM-3308](../01-canon/theorems/THM-3308-threshold-chain-modular-multicovers-and-three-layer-status-circuit.md),
with projected-atlas scope only.  The `271` exact common-status eliminations in the second intrinsic
ruler prefix all admit deterministic modular-cover certificates which do not
use a solver-selected Farkas basis or its normalization.  Their exact
tail-support spectrum is

```text
260  one load layer,
 10  two load layers,
  1  three load layers.                                      (1)
```

Thus the eleven states previously left as a “weighted exact-Farkas core” are
not an unstructured LP residue.  They are exactly the eleven cases in which
one threshold is insufficient.  Ten are two-layer circuits and one is a
genuinely three-layer circuit.

Companions:

- [`lrc14_j7_k3_z216_threshold_layer_multicover_dual_anatomy_scout_20260803.py`](../04-computation/lrc14_j7_k3_z216_threshold_layer_multicover_dual_anatomy_scout_20260803.py);
- its [frozen exact output](../05-knowledge/results/lrc14_j7_k3_z216_threshold_layer_multicover_dual_anatomy_scout_20260803.out).

This recertifies an existing projected closure.  It changes no row ledger and
does not restore a physical speed, endpoint, owner, phase, current, arbitrary
`k<=1`, the rung, or LRC(14).

## Inheritance pass

**Closest proved mechanism.**  [THM-2928 -- critical seven-comb grid
tensorization and drift tariff](../01-canon/theorems/THM-2928-critical-seven-comb-grid-tensorization-and-drift-tariff.md)
proves that the maximum mass of one upward Boolean event with prescribed
one-marginals is its fractional-cover value.  More importantly here, its
all-threshold condition `(37tg+)` warns that all target-load tails must use
one common status table.  Optimizing a different table at every threshold is
unsound.

**Canonical hostile.**  [THM-3264 -- projected-k3 `z216` low-cost gcd-eight
terminal descent](../01-canon/theorems/THM-3264-projected-k3-z216-low-cost-gcd8-seventeen-row-terminal-descent.md)
contains a screen-stage hostile: row `64` has `57` exact status eliminations
but also eight status survivors.  A proposed status template must therefore
distinguish infeasibility from a merely complicated capacity function.

**Corrected near miss.**  MISTAKE-331 and MISTAKE-333 identify exact validity
of one returned Farkas vector with neither uniqueness of its basis nor
canonicity of its magnitude.  The new computation reconstructs only the
deterministic instance

```text
(divisors, q, marginals, capacities, load histogram)          (2)
```

and derives a lexicographically deterministic integral cover.  It neither
uses nor hashes the inherited dual vector.

**Least-used relevant sidecar.**  The one-table clause in THM-2928 was more
useful than another one-threshold cover formula.  The obstruction below is
precisely that every selected threshold is individually—or, in one case,
every selected pair is jointly—feasible, while all selected layers cannot be
realized by one table.

## The common-table object

Let `x_P>=0` be the mass of status pattern `P subseteq {0,1,2,3}`.  For one
deterministic status instance, the inherited relaxation prescribes

```text
sum_P x_P = q,
sum_(P contains i) x_P = m_i.                                (3)
```

The pattern has capacity `c(P)`.  If `H_t` denotes the target histogram mass
at load at least `t`, the same table must obey every relevant tail inequality

```text
sum_(P : c(P)>=t) x_P >= H_t.                                (4)
```

For a finite threshold set `T`, define its layer count

```text
f_T(P) = sum_(t in T) 1[c(P)>=t].                             (5)
```

Suppose nonnegative integral coordinate weights `w=(w_0,...,w_3)` give the
pointwise modular majorant

```text
f_T(P) <= sum_(i in P) w_i             for every P.           (6)
```

Summing `(6)` against `x_P` and using `(3)--(4)` forces

```text
sum_(t in T) H_t <= sum_i w_i m_i.                            (7)
```

Therefore a strict reverse inequality in `(7)` is an exact contradiction.
It is a Farkas certificate, but one presented as a **threshold-chain
multicover**: selected nested load layers on the left, a zero-constant modular
coordinate cover on the right.

For `k=|T|`, `0<=f_T<=k`, so every nonnegative integral `w_i` can be truncated
to `k` without losing `(6)`.  The finite search

```text
T subset of the nontrivial thresholds,
w in {0,...,|T|}^4                                           (8)
```

is therefore exhaustive for binary-layer integral modular covers.

## Complete `271`-state taxonomy

The earlier anatomy was

```text
248 coordinate-union identities,
  1 zero-marginal-reduced union,
 11 two-fans,
 11 weighted exact-Farkas states.                             (9)
```

The multicover compiler sharpens `(9)` to

| old branch | one layer | two layers | three layers |
|---|---:|---:|---:|
| coordinate union | 248 | 0 | 0 |
| zero-reduced union | 1 | 0 | 0 |
| two-fan | 11 | 0 | 0 |
| weighted core | 0 | 10 | 1 |

Thus all `271` states have an elementary deterministic proof of the form
`(6)--(7)`.  A two-fan is not a separate LP species: its best single-layer
cover chooses either the fan stem or its two leaves according to the given
marginal costs.  The only genuinely new phenomenon in the old weighted core
is simultaneous use of several load layers.

## The eleven support-normalized circuits

There are eight distinct deterministic instances among the eleven addresses.
Each row below gives the selected load layers `T`, primitive coordinate
weights `w`, and strict gap

```text
Delta = sum_(t in T) H_t - sum_i w_i m_i > 0.                (10)
```

| atlas row / divisor address | multiplicity | `q` | `T` | `w` | `Delta` | `f_T` submodular? |
|---|---:|---:|---|---|---:|---|
| `134 / (385,1386,1980,3080)` | 1 | 3080 | `(3,4)` | `(2,1,1,2)` | 66 | yes |
| `121 / (45,280,455,3640)` | 1 | 4095 | `(3,4)` | `(2,1,2,1)` | 194 | yes |
| `17 / (49,63,196,252)` | 1 | 441 | `(1,2)` | `(2,2,1,1)` | 18 | yes |
| `27 / (49,126,196,252 or 1764)` | 2 | 294 | `(3,5)` | `(2,1,1,0)` | 6 | **no** |
| `56 / (49,126,196,252)` | 1 | 294 | `(2,5)` | `(2,1,1,0)` | 12 | **no** |
| `138 / (18,49,196,882)` | 1 | 294 | `(1,4,5)` | `(1,3,1,1)` | 11 | yes |
| `138 / (49,63,196,252 or 1764)` | 2 | 441 | `(1,2)` | `(2,2,1,1)` | 18 | yes |
| `298 / (49,63,196,252 or 1764)` | 2 | 441 | `(1,2)` | `(2,2,1,1)` | 18 | yes |

Every displayed certificate has:

- binary coefficient one on each selected tail row;
- zero constant term;
- nonnegative integral coordinate weights;
- gcd one;
- coordinatewise minimal `w`; and
- an explicitly checked set of tight Boolean patterns.

This is a normalization derived from the Boolean incidence problem, not a
rescaling of whatever basis the inherited optimizer happened to return.

## Why this is a compatibility obstruction

The row-`27` circuit makes the mechanism visible.  Its marginal vector is

```text
m=(42,126,84,252).                                           (11)
```

At thresholds `3` and `5`, the required tails are `252` and `48`.  The
pointwise inequality is

```text
1[c(P)>=3] + 1[c(P)>=5]
    <= 2*1[0 in P] + 1[1 in P] + 1[2 in P].                 (12)
```

Hence a common table would force

```text
252+48 <= 2*42+126+84 = 294,                                 (13)
```

which fails by `6`.  Neither threshold alone is contradictory.  The gain is
not “a better union bound at threshold `3`” or “a better union bound at
threshold `5`”; it is the impossibility of using two separately favorable
joint laws while pretending they are the same law.

This connects directly to THM-2928's distinction between separate
one-threshold fractional-cover optima and the all-threshold common table.
The new object is the small dual circuit witnessing that distinction.

## Exact minimality, not absence by search failure

Minimum tail support was audited in the full common-table LP, not only inside
the integral modular template.

For each of the eight distinct multi-layer instances, the script takes every
threshold subset of cardinality strictly below the claimed minimum.  It then
constructs an exact nonnegative table satisfying `(3)` and all selected
instances of `(4)`.  The tables are found by lexicographic rational vertex
enumeration over the sixteen Boolean columns.

The census is

```text
59 addressed proper tail supports,
47 distinct proper tail supports after identifying repeated instances,
47/47 with exact feasible common tables.                     (14)
```

Primal feasibility of every proper subsystem rules out **any** Farkas
certificate supported there, including certificates with arbitrary rational
tail weights, signed affine coordinates, or a nonzero constant.  Consequently
the ten support-two claims and the one support-three claim are minimum support
statements in the full dual cone.

The exceptional instance is

```text
row 138, divisors (18,49,196,882), q=294,
m=(147,42,84,126),
histogram=((0,42),(1,36),(2,6),(3,36),(4,106),(5,56),(6,12)). (15)
```

Every singleton and every pair among its six nontrivial tail rows has an exact
feasible common table.  The three layers `(1,4,5)` obey

```text
H_1+H_4+H_5 = 252+174+68 = 494,
(1,3,1,1).m = 147+126+84+126 = 483,                          (16)
```

giving gap `11`.  This records the first failed implication of the bold
two-threshold reframe:

```text
nested threshold events
    does not imply
every common-table obstruction has a two-layer certificate. (17)
```

Three layers are genuinely necessary here.

## Submodularity is not the mechanism

The selected layer-count functions are submodular at eight of the eleven
addresses, but not at the two copies of row `27` or at row `56`.  Exact first
violations include

```text
row 27: A=2, B=8,
f(A)+f(B)-f(A union B)-f(A intersection B) = -1;

row 56: A=6, B=8,
f(A)+f(B)-f(A union B)-f(A intersection B) = -1.             (18)
```

None of the eleven selected functions is supermodular.  Yet all eleven have
valid modular majorants.  The source-to-target map is therefore

```text
nested capacity tails -> layer-count function -> modular majorant,          (19)
```

not

```text
nested capacity tails -> submodular rank function.                          (20)
```

The submodular/polymatroid idea was a useful hostile probe, but demanding
submodularity would discard three valid certificates.

## Hostile survivor

Replaying row `64` gives

```text
115 states = 50 crude + 57 status + 8 residual.              (21)
```

All `57` genuine status eliminations have one-layer modular covers.  For the
first residual divisor state

```text
(80,2156,4312,5390), q=43120,
m=(6468,6160,6160,6160),
c(empty)=0, c(P)=1 for P nonempty,
histogram=((0,28046),(1,15074)),                              (22)
```

the sparse table

```text
x_0=28046,
x_3=3234, x_5=3234, x_6=2446, x_8=5680, x_14=480,
all other x_P=0                                              (23)
```

has exactly the required total, marginals, and nonempty-event tail `15074`.
The multicover search correctly returns no contradiction.  This is a hostile
negative control against the claim that the new compiler simply kills every
four-bit packet.  THM-3264 later closes row `64` by its terminal machinery;
`(23)` concerns only the earlier common-status stage.

## Beach-search synthesis

Three older motifs meet here without being identified:

1. **Fractional covers.**  THM-2928 gives the exact one-threshold marginal
   transport optimum.  The new certificate is its common-table,
   multiple-threshold dual companion.
2. **Laminar cuts.**  The selected load events form a genuine nested chain.
   The laminarity lives in threshold space, not necessarily in the minimal
   true clutters on status coordinates.
3. **Small Boolean circuits.**  As in low-depth Boolean/Mobius cut work, one
   should form the full joint law before selecting a few features.  Here the
   surviving features are two or three tail layers, and the lost information
   is certified irrelevant only after the pointwise sixteen-pattern audit.

The reusable procedure is:

```text
capacity table c(P)
 -> nested tail indicators
 -> enumerate support-minimal layer sums
 -> compute a modular multicover
 -> prove every smaller tail subsystem feasible
 -> retain the circuit, tight patterns, and physical-loss sidecar.           (24)
```

This replaces a normalization-sensitive generic dual by a finite certificate
compiler with an intrinsic complexity statistic: **tail-circuit arity**.

## Next operations

1. Apply `(24)` to the next intrinsic `z216` ruler prefix and record capacity
   signatures before solving each status LP.  Reused signatures may compile
   whole status families at once.
2. Run the circuit compiler across older projected `k=3` status populations.
   The observed arity bound `<=3` is proved only for these `271` instances;
   it is not yet a universal four-bit theorem.
3. Classify which capacity tables force integral modular multicovers and which
   need fractional coordinate weights.  THM-2928's Fano-plane caveat warns
   that integrality cannot be assumed at larger arity.
4. Keep the physical-loss ledger explicit.  A compact status circuit is still
   only a proof inside the necessary quotient until endpoint, owner, phase,
   and current data lift it lawfully.

## Reproduction and hashes

```text
python3 04-computation/lrc14_j7_k3_z216_threshold_layer_multicover_dual_anatomy_scout_20260803.py
python3 -O 04-computation/lrc14_j7_k3_z216_threshold_layer_multicover_dual_anatomy_scout_20260803.py
```

Normal and optimized transcripts are byte-identical.

```text
script  e7462eabab773133688079141708d2377742e2d6f9a9480089666956f472020c
output  f83d7ca1a3e2557330ec7ed655dc756cef21d6abdb52326d3623563c2abc407d
semantic bd7ac011d960afe7c4566ac9c17e9df3d8c2cf710e6023a715cb6a6d1ce1408d
```
