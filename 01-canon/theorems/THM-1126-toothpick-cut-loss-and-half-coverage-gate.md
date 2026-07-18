---
id: THM-1126
title: Toothpick cut-loss duality and the half-coverage gate for the sharp four-comb target
status: PROVED — elementary one-dimensional topology; reduces failure of a component longer than one finest tooth to at least half-coverage, and gives an explicit spanning-tree overlap criterion
source: codex-2026-07-18-S73 r5 continuation
depends_on: [THM-1097, THM-1111]
related: [THM-1094, THM-1101, MISTAKE-164]
---

# THM-1126 — toothpick cut-loss duality and the half-coverage gate

For a positive integer `k`, write

```text
D_k={t in [0,1] : ||kt||<1/14}.
```

Each complete tooth of `D_k` has length `1/(7k)`.  Endpoints may be made
open or closed throughout: they do not affect component lengths or measure.

## 1. The cut-loss lemma

Let `G` be a finite union of intervals, let `C(G)` denote its number of
positive-length components, and put

```text
delta = |G intersect D_k|.
```

Then

```text
C(G minus D_k)-C(G) <= 7k delta.                       (1)
```

Indeed, a positive increase in component count occurs only when a tooth is
wholly contained in the interior of a component of `G`.  Such a split costs
the complete tooth length `1/(7k)`.  Distinct new splits use distinct teeth.
Trimming an endpoint, deleting a component, or joining the two boundary
representatives of the circular tooth does not increase the count.  Charging
each new component to its complete contained tooth proves (1).

This is the coupling discarded by the separate discrepancy/component bounds
in THM-1097: fragmentation is not free.  Every added cut carries an exact
amount of lost mass.

## 2. The half-coverage gate

Let `G` be a finite union of `c` positive-length intervals, of total measure
`ell`, and let

```text
E = G intersect (D_k1 union ... union D_km),
A = G minus E,
K = max_i ki,
tau = 1/(7K),
delta = |E|.
```

Then

```text
delta < (ell-c*tau)/2                                  (2)
```

implies that `A` has a component of length strictly greater than `tau`.
Equivalently, failure of the sharp component target forces

```text
delta >= (ell-c*tau)/2.                                (3)
```

To prove this, count separately inside each component of `G`, and let `b` be
the total number of components of `E` which meet neither endpoint of their
ambient `G` component.  Every such internal component contains a complete
tooth of at least one `D_ki`: take a point in the component and a tooth
containing it; if that tooth met the ambient boundary, connectedness would put
the danger component at the boundary too.  Consequently every internal component has
length at least

```text
1/(7ki) >= 1/(7K)=tau,
```

and hence `b tau <= delta`.  The complement `A` has at most `b+c`
positive-length components.  If all of them had length at most `tau`, then

```text
|A| <= (b+c)tau <= delta+c*tau,
ell=|A|+delta <= 2delta+c*tau,
```

which is (3).  Taking the contrapositive proves (2).

The single-interval form used below is the case `G=J`, `c=1`.  For the open
`r=5` stratum, `m=4` and `K=k4`.  Thus it is enough to prove that the union of
the first four killer combs occupies less than

```text
ell/2 - 1/(14k4)                                       (4)
```

inside one core-safe interval.  Alternatively, one may take `G=S(P)` itself:
if it has `c(P)` components and measure `mu(P)`, the sufficient whole-core
bound is

```text
|S(P) intersect union_i D_ki|
  < mu(P)/2 - c(P)/(14k4).                              (4')
```

The whole-core form can recover a long component that is not descended from a
preselected largest core interval.  Both forms are strictly more faithful
than bounding survivor mass and cut incidence independently.

## 3. A spanning-tree overlap certificate

Put `Ai=J intersect D_ki`.  For any spanning tree `T` on the comb indices,
the elementary tree-union inequality is

```text
|union_i Ai|
 <= sum_i |Ai| - sum_{ij in T}|Ai intersect Aj|.       (5)
```

Root `T` and add its vertices parent-before-child.  When `Ai` is added, its
overlap with the preceding union is at least its overlap with its parent;
summing the four incremental bounds proves (5).  Maximising the subtracted
term gives exactly the maximum-spanning-tree observable of THM-1111.

The sharp one-comb discrepancy from THM-1097 gives

```text
|Ai| <= ell/7 + 6/(49ki).                              (6)
```

Therefore, for four combs `k1<k2<k3<k4`, the sharp target follows whenever
some spanning tree `T` satisfies

```text
sum_{ij in T}|J intersect D_ki intersect D_kj|
  > ell/14 + (6/49)sum_i(1/ki) + 1/(14k4).             (7)
```

Indeed, (5)--(7) make the union measure strictly smaller than (4), and the
half-coverage gate supplies a component longer than `1/(7k4)`.  A fifth
killer `k5>k4` then cannot cover that component by THM-1097's sharp final-comb
argument.

Criterion (7) is sufficient, not asserted here to hold for every tuple.  Its
value is that it states the missing r=5 information exactly: the coarse
four-comb calculation misses the target by only `ell/14` plus endpoint terms,
and three labelled pair-overlap edges are enough to pay that deficit.

### Exact negative audit of the sufficient gates

THM-1123 subsequently evaluated all `495*binom(16,4)=900900` tuples formed
from the first sixteen legal killer speeds above each eight-core.  The sharp
component target itself has zero failures, but the sufficient gates here do
not explain it:

```text
largest-J half-coverage failures     441576 / 900900,
largest-J MST-criterion failures     823588 / 900900,
whole-core half-coverage failures    671455 / 900900.
```

Exact counterrows include

```text
half gate worst:
  P={1,2,3,4,5,8,9,10}, k=(131,132,140,141),
  gate slack=-539071/106195936,       7k4L=438/131;

MST gate worst:
  P={2,3,4,5,6,7,8,9}, k=(121,124,127,130),
  gate slack=-466976387/54621386820,  7k4L=813/242;

whole-core gate worst:
  P={1,2,3,4,6,10,11,12}, k=(158,161,162,164),
  gate slack=-29245585/929275578,     7k4L=55/18.
```

Thus THM-1126 is a genuine reduction and diagnostic, not a hidden uniform
closure.  The surviving sharp margin is carried by the nonuniform distribution
of individual gap lengths, which scalar union mass and comb-level pair-overlap
trees erase.

## 4. The gap-energy / truncated-covariogram gate

The exact negative audit shows what the scalar gates erase: the survivor
lengths are highly nonuniform.  There is a natural quadratic carrier for that
dispersion.  If the components of the final safe set `A` have lengths
`lambda_1,...,lambda_s`, then

```text
H(A)=sum_j lambda_j^2.                                    (8)
```

The elementary implication

```text
H(A)>tau |A|  ==>  max_j lambda_j>tau                    (9)
```

follows termwise: if every `lambda_j<=tau`, then
`lambda_j^2<=tau lambda_j`.  The arithmetic statement is kernel-checked as
`gapSquareEnergy_forces_long` in `LRCSharpCombArithmetic.lean`.

There is also a component-free functional form.  Consecutive components of
`A` are separated by internal components of the danger union, each of length
at least `tau`.  For translation on the ambient line put

```text
Gamma_tau(A)=integral_(0<=h<=tau) |A intersect (A-h)| dh. (10)
```

No translate of size `h<tau` can jump from one safe component to another, so

```text
|A intersect (A-h)|=sum_j max(lambda_j-h,0),

Gamma_tau(A)=sum_j phi_tau(lambda_j),

phi_tau(x) = x^2/2                    if x<=tau,
             tau*x-tau^2/2            if x>=tau.          (11)
```

Consequently

```text
Gamma_tau(A)>tau |A|/2  ==>  max_j lambda_j>tau.         (12)
```

Unlike mass and component count, `Gamma_tau` is a two-point functional of the
final safe indicator:

```text
|A intersect (A-h)|=integral 1_A(t)1_A(t+h) dt.
```

This is the useful functional form of an `H`-drift: a comb cut changes a
short-shift autocorrelation, not merely one scalar mass.  It is accessible to
Fourier/covariance methods while retaining the gap-length dispersion that the
comb-level MST loses.

On the hardest bottom-bank row

```text
P={1,2,4,5,7,9,11,12}, k=(158,160,162,164),
L=41/25920,
```

the mean-length test fails (`|A|-16tau<0`), whereas the exact covariogram
ratio is

```text
Gamma_tau(A)/(tau |A|/2)
  = 336218180863/317189087460
  = 1.05999... > 1.                                      (13)
```

Thus (12) detects precisely the nonuniformity carrying this row.  No uniform
lower bound for (10) is claimed here; that is a sharper next target.

## 5. Structural and tournament reading

The useful vertices are not automatically runners.  Three alternate carriers
appear naturally:

1. complete teeth, for the split-to-lost-mass charging in (1);
2. connected danger components inside `J`, for the half-coverage gate; and
3. the four killer combs, for the overlap tree in (5)--(7); and
4. a safe component paired with its translate, for the covariogram (10).

On the third carrier, the pairwise observable is exact intersection measure.
Orient `i -> j` when the chosen deterministic tie-break regards edge weight
from `i` as larger; equal weights are ties, and sorting by that gauge gives a
tie Hamiltonian path.  The maximum spanning tree preserves the precise datum
needed by (7), but destroys endpoint order and higher intersections.  The
tooth/component carrier preserves the sharp LRC predicate (existence of a
safe gap longer than one finest tooth), while a scalar union measure alone
forgets where its mass lies.  This challenges the default assumption that
runner or residue vertices are the faithful tournament object: here the
proof obligation lives on labelled overlap edges plus a one-dimensional
component sidecar.

## 6. What this proves and what remains

THM-1126 does **not** close uniform `r=5`.  It proves two reusable exact
lemmas and replaces the failing coarse quotient by the following narrower
task:

> control the labelled endpoint word strongly enough to prove the sharp gap
> directly, or obtain a uniform lower bound for the truncated covariogram
> (10).  The half-coverage and comb-level overlap gates are useful on large
> subbanks but are not uniform.

The reduction also explains why the near-equal toothpick rays can recover at
large scale even though the uncoupled sufficient inequality in THM-1101 gets
worse: overlapping teeth simultaneously reduce occupied union mass and the
number of effective cuts.
