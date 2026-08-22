# The pointed six is a bidirected-path arc module, not a global ceiling

**Status: DERIVED FINITE-EXACT REPRESENTATION SIDECAR WITH AN INDEPENDENTLY
AUDITED FIVE-COORDINATE PARENT.  THE POINTED RESPONSE TABLE IS RECONSTRUCTED
THROUGH THE INDEPENDENT IMPLEMENTATION; THIS CHEAP SIDECAR DOES NOT REPEAT THE
FULL ENDPOINT SWEEP.  COMMIT `1d4b7f1ee` HAS SINCE RECONSTRUCTED THE COMPLETE
DIAGONAL BUNDLE WITHOUT IMPORTING THE CANDIDATE AND ACCEPTED EVERY PARENT FACT
USED BELOW.  EVERYTHING HERE IS STATIC.  THERE IS NO CHRONOLOGY, LOCAL SYSTEM,
COCYCLE, PHYSICAL CURRENT, ROW EXCLUSION, OR LRC(14) CONCLUSION.**

The six pointed incidences have a canonical combinatorial model.  Put the
owner states in Gray-path order

```text
A=0, B=1, C=3, D=2.
```

The three roots occur in exactly the adjacent state pairs

```text
root 0  : {A,B},
root 6  : {B,C},
root 12 : {C,D}.
```

Sending `(S,u)` to the arc from `S` to the other state containing `u` gives

```text
(0,0)   -> A->B,     (1,0)  -> B->A,
(1,6)   -> B->C,     (3,6)  -> C->B,
(3,12)  -> C->D,     (2,12) -> D->C.               (1)
```

Thus the carrier is simultaneously:

- the six edge labels of the alternating incidence path

  ```text
  A--root0--B--root6--C--root12--D = P7;
  ```

- the six directed arcs of the bidirected path

  ```text
  A<->B<->C<->D;
  ```

- the square `Q2` in Gray order `0--1--3--2--0` after deleting the closure
  edge `2--0`, then retaining both directions on each surviving edge.

This is the exact size-four generalized-relation interpretation.  Among the
six unordered state pairs its profile is

```text
(both directions, one direction, missing) = (3,0,3).             (2)
```

It is not a tournament.  A tournament on four vertices also has six arcs,
but its profile is `(0,6,0)`.  Any cardinality-only bijection with the six
unordered edges of `K4` destroys arc endpoints, root pairing, tail state,
incidence, and both involutions below.

## Exact edge-response theorem

Let `k=F_p`, with

```text
p=755373809845391722745761,
```

and let `A=k^6` have basis `e_0,...,e_5` in the order (1).  The independently
reconstructed pointed response rows define

```text
rho : A -> Fun(F_13,k),
rho(e_i) = the relation row of pointed incidence i.               (3)
```

The exact rank of (3) is six.  Hence `rho` is injective: the relation rows do
not merely have the same dimension as an edge module; they are a labelled,
response-realized copy of that module.  Their canonical row-space hash is

```text
6e9083f15408f6d2d85fb3f2747ba0bd1f987e83ce4b836cb7298aaccc84e0c4.
```

For a state `S`, let `t_S` be the sum of arcs with tail `S`.  Then

```text
rho(t_A), rho(t_B), rho(t_C), rho(t_D)
```

agree entrywise with the four Boolean-parent relation rows reconstructed by
summing the pointed root-difference tensor.  Their rank is four and their
canonical row-space hash is

```text
1b4fef00a23dcb79dc52ace662bae2f858ce3da27b6ef19b902ae40f5a79e755.
```

This explains the statewise carrier counts.  In path order `A,B,C,D` the
outdegrees are `(1,2,2,1)`; in the repository's numeric state order
`0,1,2,3` they are `(1,2,1,2)`.  The rank-six channel is therefore the full
six-arc carrier, while the old rank-four Boolean parent is its tail-star
submodule.

## Two involutions, not one cosmetic reflection

There are two commuting typed involutions on `A`.

Arc reversal acts within a root pair:

```text
iota=(0 1)(2 3)(4 5).
```

Path reflection comes from `state -> state XOR 2` and `root -> 12-root`:

```text
j=(0 5)(1 4)(2 3).
```

Both split the full arc module as `3+3`, but they mean different things.  The
previously audited pointed-path reflection is `j`, not `iota`.  Under `j`, the
tail-star parent splits `2+2`; the three root-pair sums split `2+1`.

Arc reversal gives the useful representation decomposition

```text
A = A^+ direct_sum A^-,

A^+ = span(e0+e1, e2+e3, e4+e5),
A^- = span(e0-e1, e2-e3, e4-e5).                  (4)
```

The even module `A^+` is exactly the three-root quotient.  The odd module
`A^-` is the formal oriented-edge-current module of the underlying path.

## The period-three conditional ranks are now explained exactly

Write the diagonal transition first exposed at `b1baa781a` as

```text
K_(r0,r1)=diag(k_0(r0,r1),...,k_5(r0,r1)).         (5)
```

The independent hostile audit at `1d4b7f1ee` reconstructs the ordered
`(u,q)` contributions by exact OWNER-lawful source cell before expanding the
169 addresses.  It matches all four candidate tensor digests and proves that
all 169 matrices (5) exist uniquely, are diagonal, and satisfy

```text
sum_r1 K_(r0,r1)=I_6
```

for every fixed `r0`.  Its semantic digest is

```text
b58d98c283a2dc42111365a3d61af0948757feda3b1ff11ca65cd7d562b15a56.
```

The smaller source-profile audit here independently reconstructs all
coefficients for which child/parent proportionality holds and uses the now
independently accepted exact `954/954` source-to-output agreement.

For every interior nonmultiple of three,

```text
r0 in {1,2,4,5,7,8,10,11},
```

all source ratios are defined, the two coefficients on each reversed arc
pair are equal, and the thirteen coefficient rows span exactly

```text
A^+.                                                       (6)
```

Every one of these eight row spaces has canonical hash

```text
7016acfc7f118cac9d673713ebd4139b2d57800d72a2527007322f81a609a4c6.
```

At `r0=6`, every source ratio is again defined.  The coefficient rows span
exactly

```text
Q = A^+ direct_sum k*(e2-e3),                         (7)
```

whose canonical hash is

```text
a3ca24b41fb507941de6e274058bf4658918a13d9ca9902e6483fb52be444446.
```

At `r0=3` and `r0=9`, the only undefined source ratios lie respectively on
one of the two middle arcs.  Every coefficient on the outer root pairs is
defined and pair-equal, so the completed output rows from (5) lie in `Q`.
The independently audited conditional rank is four.  Because the six pointed
state-embedded response rows are independent, the map from a coefficient row
`(k_0,...,k_5)` to its `(state,relation)` response is injective.  Therefore the
coefficient-row rank is also four, and containment in the four-space `Q`
forces equality.

Consequently the exact interior theorem is

```text
span_r1 diag K_(r0,r1) =
  A^+                         if 3 does not divide r0,
  A^+ + k*(B->C - C->B)      if r0 in {3,6,9}.        (8)
```

Equation (8) gives a representation-theoretic explanation of the interior
conditional ranks

```text
3,3,4,3,3,4,3,3,4,3,3.
```

It is not a `C3` action.  Divisibility by three selects whether the middle
arc-odd line occurs; translation by three is still a 13-cycle and the rank
record remains nonconstant along that orbit.

The boundary digits `r0=0,12` also have conditional rank four, but they arise
from the previously audited boundary projection `6 -> 4`.  The present
source constraints do not identify their completed four-spaces with `Q`, so
(8) deliberately excludes them.

## Why the extra line is not an H1 flux

Call

```text
j_mid=e2-e3=B->C-C->B.
```

This is a useful formal current coordinate, not a physical current.  In the
ordinary oriented chain complex of the state path, its boundary in state
order `A,B,C,D` is

```text
partial(j_mid)=2(C-B) != 0.                           (9)
```

The boundary map on all three generators of `A^-` has rank three, so its
kernel is zero.  Equivalently, `P4` is a tree.

The original alternating-incidence complex is even more direct: its six P7
edges have a `7 x 6` boundary matrix of rank six, hence

```text
H1(P7;k)=0.                                           (10)
```

There is a subtle typing trap here.  If the six independent arcs of the
bidirected path are fed naively into a directed-multigraph boundary, each
forward-plus-backward pair looks like a length-two cycle.  Those are artifacts
of changing the chain complex.  The bijection (1) is a carrier-module
bijection, not a chain isomorphism from P7 to a multigraph with parallel
opposite edges.  Equations (9)--(10), not the artificial backtracks, govern
the lawful static topology.

## The exceptional tent is a section, not yet a divisor

The 60 failures of source-cell proportionality occur exactly on

```text
D = {(e_i,r0=h_i,r1): 0<=i<6,
     r1 in F_13 \ {0,6,12}},

h=(12,12,9,3,0,0).                                  (11)
```

The affine reflection law is

```text
h_(5-i)=12-h_i,
```

and the successive drops along the arc/P7-edge order are

```text
h_i-h_(i+1)=(0,3,6,3,0)=3*(0,1,2,1,0).              (12)
```

Thus (11) is exactly a reflection-equivariant edge section whose discrete
derivative is a tent on the five internal adjacencies.  Its line-graph
Laplacian is

```text
Delta h=(0,3,3,-3,-3,0).                             (13)
```

Arc reversal gives the further exact decomposition

```text
h_even=(12,12,6,6,0,0),
h_odd =(0,0,3,-3,0,0)=3*j_mid.                       (13a)
```

The [typed C4 cospan audit](lrc-r5-tent-location-c4-cospan-hostile-audit-codex-20260816.md)
shows both the value and the limit of (13a).  It lives in the
exception-**location** module.  Under a marked normalized odd-edge basis,
middle coefficient three has the unique formal `C4` completion
`3*(1,1,1,1)` and mod-13 seam `12=-1`; under literal directed-arc chain
realization `j_mid` maps to twice the middle edge, giving coefficient six and
seam `11`.  The formal cycle restricts to `(3,3,3)`, not the observed
`(0,3,0)`, so it completes only one selected coordinate.  No lawful closure
edge or location-to-response coefficient transport is present.  The exact
diagonal-profile hostile is stronger: on the same middle support the actual
odd response takes eight distinct values in each reflected `r0` row, and the
two rows have rank two, so the location coefficient three admits no
`r1`-blind scalar transport.

What (11) is not yet is a divisor of the diagonal bundle.  No valuation,
zero/pole order, or section whose zero locus is `D` has been supplied.  The
observed property is failure of source proportionality, not vanishing of
`K`.  In the purely combinatorial line-graph sense, the principal divisor of
`h` would be (13), which is not the 60-point support (11).

There is also an exact hostile to a simple P7-incidence explanation of the
endpoint-killed remainders.  Quotient each exceptional source stalk by its
one-digit parent line.  The six stalk dimensions are

```text
(2,1,1,1,1,2),                                       (14)
```

so the global quotient defect has dimension eight and reflection split
`4+4`.  Before quotienting, the 60 source rows span dimension 14 and the six
parent lines span dimension six.  Ordinary `C1(P7)` has dimension six and
`H1(P7)=0`; neither can be the eight-dimensional kernel mechanism.  The
independent five-coordinate audit now additionally verifies that every one of
the 60 nonproportional source remainders is annihilated by the common endpoint
operator.  Thus the eight-dimensional quotient defect is an exact subspace of
the induced endpoint kernel, not merely a source-side dimension count.

The strongest exact survivor is therefore a decorated edge sheaf: rank two
at the terminal arcs and rank one at the four inner arcs.  The endpoint map
annihilates its classes after quotienting the target by the corresponding
pointed response line.  What remains open is a basis-level factorization of
that eight-space through a decorated incidence operator; ordinary `P7`
incidence is already refuted.

## Correct static formalization of the diagonal bundle

Because every `K_(r0,r1)` is diagonal, write its six diagonal entries as

```text
k_e : F_13(r0) x F_13(r1) -> k,   e in E(P7).        (15)
```

The projective sum law is precisely

```text
sum_r1 k_e(r0,r1)=1                                  (16)
```

for every `(e,r0)`.  The best current language is therefore:

> six normalized scalar cylinder sections, or equivalently an
> endomorphism-valued partition of unity on the trivial rank-six arc bundle.

“Measure” is acceptable only in the algebraic finite-additive sense; the
field has no positivity order.  “Diagonal local system” is too strong:
35 matrices are zero, other matrices have ranks two or four, no stationary
`K_r1` exists, and there are no typed base edges or composition maps.
“H1/cocycle” is a type error at this stage: (15) is static function data, not
parallel transport.

The complete-profile companion at `7f88e61af`, now checked against the
independent reconstruction at `1d4b7f1ee`, sharpens this classification
without changing it.  Its six `13 x 13` matrices
`(k_e(r0,r1))` have ranks

```text
(11,11,12,12,11,11),
```

each has right and left eigenvalue-one nullity one, and their direct-sum rank
is `68/78`.  Arc-reversal equality holds at exactly `143/169` addresses.  The
26 failures all split the middle pair: 20 are the nonroot digits at
`r0=3,9`, and six are the endpoint residuals

```text
(0,11),(0,12),(6,5),(6,7),(12,0),(12,1).
```

Chamber reflection holds on all 1,014 scalar entries.  These complete-profile
facts support (8) through two independent construction routes, but row-sum
one over a finite field still supplies neither positivity nor composable
transport.

## Third-digit hostile: diagonal is not a cocycle

The incoming exact `r2` refinement at `840cec984` retains

```text
a=r0+13*r1+169*r2+2197*c.
```

It gives the cleanest falsification of a hidden local-system interpretation.
Every one of the 2,197 third-digit children remains on the same six pointed
arc lines, and every canonical map is diagonal there.  Full `Mat_6`
uniqueness holds only for the `130*13=1690` children above rank-six parent
addresses; over the 39 deficient parents only the action on live arc lines is
observable.  Accordingly,

```text
sum_r2 L_(r0,r1;r2) = I_6
```

on the 130 full-rank parents, while on every parent it is exactly the
projector onto the live arc lines.  This is a partition-of-unity law on a
varying observable support, not parallel transport.

Both plausible adjacent-pair laws fail nonvacuously.  The literal equation

```text
L_(r0,r1;r2) = K_(r1,r2) K_(r0,r1)
```

and its type-correct cumulative version hold at only `455/2197` triples.
Those are precisely the 13 children of each of the 35 zero parents.  Thus
there is no surviving adjacent-digit cocycle on a live carrier.

Nor does adjoining the previous digit repair the defect.  The 78 rows indexed
by `arc x r1` have parent rank 68, but for every fixed `r2` the parent/child
union has rank 78.  Hence no linear operator on the old 68-dimensional
observable rowspace produces the next digit.  The ten new directions are
address-amplitude/history data: they do not mix the six arcs and they do not
constitute chronology.  This sharply separates a stable six-line carrier
from a growing amplitude state and rules out interpreting “rank six” as a
global endpoint-response ceiling.

The same refinement finds 104 live source profiles with two source-cell
ratios, all endpoint-annihilated enough to preserve a diagonal response.  That
is consistent with the decorated-kernel diagnosis above and further hostile
evidence against ordinary `P7` incidence as its mechanism.

A future D5/H1 interpretation needs a lawful clock/address graph `Gamma`, a
typed transport on each clock edge, and a composition law.  Only after
restricting to invertible transports can one test cycle holonomy or a
`Z^1(Gamma,(k^x)^6)` class.  A clock edge closing the missing `D--A` edge is a
possible test object, not a current theorem.  If such an edge is lawfully
constructed, the state path closes to `C4`, whose first (co)homology is
one-dimensional; in a consistent orientation the candidate class coordinate
is the sum of the four edge values.  The present tensor supplies only three
of those four typed edges.

## Connection contract

| field | exact answer |
|---|---|
| source | independently reconstructed pointed root-difference response plus exact two-digit source profiles |
| target | six labelled pointed relation rows and diagonal coefficient sections |
| carrier map | `(S,u) -> S -> other state containing u` |
| preserved | tail owner state, absolute root, root-paired arc reversal, path reflection, relation response, `r0,r1` labels |
| destroyed by arc view | explicit root vertex as a chain-complex vertex; source cells; root difference; upper address digits; chronology |
| exact carrier | injective six-dimensional P7-edge/bidirected-P4 arc module |
| Boolean parent | four-dimensional tail-star submodule |
| root quotient | arc-reversal-even module `A^+`, dimension three |
| interior rank-four extension | `A^+` plus middle arc-odd line |
| exceptional source object | graph section (11), decorated quotient stalks `(2,1,1,1,1,2)` |
| third-digit hostile | all six arc lines preserved; adjacent cocycle dead-only; old `68`-space grows to `78` |
| lost claims | no boundary-space classification, local system, cocycle, physical current, row exclusion, or LRC(14) |

## Cheapest next decisive computations

1. Add postprocessing to the independent full-tensor audit: RREF the six
   diagonal coefficient vectors for each `r0` against `A^+` and
   `A^+ + k*j_mid`.  The matrices are already resident, so this directly
   confirms (8) and identifies the two boundary four-spaces without another
   endpoint sweep.
2. In that same resident tensor, retain a basis of the 60 absolute remainders
   and test every natural decorated-stalk incidence matrix with stalk ranks
   (14).  This decides whether the exact eight-dimensional kernel subspace has
   a sheaf-boundary factorization or is only endpoint-specific.
3. Decompose the exact ten-dimensional third-digit quotient
   `span(parent,child)/span(parent)` under arc reversal and path reflection.
   This is a small postprocessing calculation on the resident 78-row banks
   and decides whether the missing history state has a canonical
   representation, without rerunning the endpoint sweep.
4. If a lawful address clock is constructed, add one typed clock edge at a
   time and test invertibility plus composition before computing holonomy.
   Without those gates, no D5/H1 promotion is admissible.

## Reproduction

```text
python -B 04-computation/lrc_r5_p7_arc_carrier_tent_section_audit_20260816.py
python -B -O 04-computation/lrc_r5_p7_arc_carrier_tent_section_audit_20260816.py
```

Normal and optimized replays are byte-identical.  The pinned record and
semantic SHA-256 values are respectively

```text
c1e17b60d62fdafe6a2b99fd5df6d097db1d327590ab0995bda182466d609ec4
54f2fd0702cfee4b3af954ca3719b9966cb06e9227a94bdc8b14651c452d5feb.
```
