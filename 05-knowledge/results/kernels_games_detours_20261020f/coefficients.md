# Completing the deficiency-three kernel sum

**PROVED:** the structural shortcuts for the two largest primitive-kernel
budgets, and the explicit all-k formula below obtained from the inherited
finite-kernel theorem. **FINITE-EXACT:** the C4 owner census, the small
placement counts for C5, and the stated verification receipts. These are
ordinary finite mathematical computations, not a formal proof-assistant
verification. No claim of literature novelty is made.

The closest mechanism is the external-port extraction and matching-completion
theorem in [the preceding kernel packet](../rotation_owners_spectrum_20261020e/owners.md).
Its hostile is that two matchings with identical induced light gaps can have
different outside partners and different owners. Its corrected near miss is
that contracting clones can change deficiency. This computation therefore
retains every defect edge and uses the pair-free predicate before ownership;
it never identifies ports by their local gap type. The underused comparison
is the SCC-word recognition grammar discussed in that packet: a bounded
exception can be counted around unbounded bulk, provided its interface is
retained. Here the bulk consists of matching endpoints, not SCCs. Triangle
modules below are never treated as transitive chain weights.

The live objects are primitive defect components, exceptional gaps, exact
outside partners, clone pairs, and least-good-record ownership. The only new
quotient is deletion of inactive matching components, already proved to
preserve the owner predicate. This note finishes the finite arithmetic that
the preceding packet explicitly left open.

## 1. Explicit result and scope

Let N3(k) be the number of capped weighted, pair-free chain-core isomorphism
classes with feedback budget at most k and exactly `h=2k-2` heavy coordinates.
The cap is k+1. This is the deficiency-three layer of the fixed-feedback
counting grammar; it is **not** A000568 itself.

For every integer k>=2,

```text
N3(k) = (2k-3)!! / 6 *
        (16k^7 + 64k^6 + 8k^5 - 140k^4 - 22k^3
         + 135k^2 + 61k - 110).                              (1)
```

At k=1, the separate zero-heavy calculation gives N3(1)=5. Formula (1)
would give 2 there; its k>=2 restriction is essential.

The computed initial values are:

| k | N3(k) |
|---:|---:|
| 1 | 5 |
| 2 | 756 |
| 3 | 36,473 |
| 4 | 1,243,815 |
| 5 | 38,295,600 |
| 6 | 1,157,009,490 |
| 7 | 35,523,051,795 |
| 8 | 1,126,532,792,385 |
| 9 | 37,192,164,159,150 |
| 10 | 1,282,993,311,600,000 |

The values at k=2,3 agree with the earlier independent full capped-class
census. The later values are consequences of the finite-kernel theorem and
the computations here; they are not additional full-class censuses.

For an alternative form with nonnegative coefficients, set x=k-2. Then

```text
N3(k) = (2k-3)!! P(k-2),

6 P(x) = 16x^7 + 288x^6 + 2120x^5 + 8260x^4
         + 18378x^3 + 23395x^2 + 15953x + 4536.              (2)
```

The useful factored-basis form is

```text
P(x) = C2(x) + x C3(x)/3 + x(x-1) C4(x)/30
       + x(x-1)(x-2) C5(x)/630
       + x(x-1)(x-2)(x-3) C6(x)/22680.                    (3)
```

This basis displays exactly which primitive budget first contributes, so it
is more useful here than an attempted factorization of the expanded septic.
The inherited eighth-order recurrence is now numerically executable, since
(1) supplies all its initial values. A layer-by-layer structural recurrence
has advanced, but an efficient scalar recurrence for unrestricted A000568
has not been obtained.

## 2. The five finite coefficient polynomials

An accepted primitive symbolic kernel contributes y for each exceptional
gap of symbolic weight L. The symbols are `0,1,2,L,H`, with L representing
any light weight from 3 through k and H a heavy coordinate. Thus y=k-2.
The complete polynomials, in ascending coefficient order, are

| j | 1 | y | y^2 | y^3 |
|---:|---:|---:|---:|---:|
| 2 | 756 | 792 | 270 | 30 |
| 3 | 13,380 | 12,974 | 4,140 | 435 |
| 4 | 142,335 | 130,750 | 39,810 | 4,020 |
| 5 | 893,760 | 788,970 | 231,840 | 22,680 |
| 6 | 2,593,080 | 2,222,640 | 635,040 | 60,480 |

In particular,

```text
C5(y) = 210 (2y+7)(6y+19)(9y+32),
C6(y) = 7560 (2y+7)^3.                                    (4)
```

The preceding packet proved, independently of these numerical values,

```text
N3(k) = sum_(j=2)^min(6,k)
        Cj(k-2) binom(2k-3,2j-3) (2k-2j-1)!!,             (5)
```

with `(-1)!!=1`. Cancelling factorials in (5) gives (3); multiplying out gives
(2) and substituting x=k-2 gives (1). No interpolation or polynomial fit is
used at any stage.

## 3. Exhaustive small graph grammar and C4

For a primitive kernel write `q=j-e`, `delta=2e-t`, and let u be its number
of exceptional gaps. The inherited six-menu identity is

```text
2q + delta + u = 3,  q in {0,1},  t=2e-delta.
```

All defect endpoints have weight one, except that q=1 permits a single
weight-two endpoint of degree one. It also permits the all-unit endpoint
word. Every endpoint is activated by a nonmatching component, a nonunit
matching endpoint, or adjacency to an exceptional gap, followed by closure
under defect adjacency.

A connected defect component with v vertices and e edges has charge 2e-v.
For charge at most three its complete possibilities are an edge (charge 0),
P3 (charge 1), a four-vertex tree (charge 2), a five-vertex tree (charge 3),
or a triangle (charge 3). Connectedness proves v<=charge+2; enumerating the
edges on those bounded vertex sets then proves the list exhaustive.

The generator chooses the component containing the least unused endpoint,
then recurses. This gives each ordered labelled defect graph exactly once.
For all edge budgets through four, a separate enumeration of every subset
of the complete graph, filtered only by absence of isolated endpoints,
checks equality of the graph universes and absence of duplicates. Thus the
component grammar is checked independently of its own recursion.

At j=4 the universe consists of 930 defect graphs across the six menus and
337,425 primitive symbolic records before the clone and owner filters.
Every record is expanded to its exact weighted tournament. A homogeneous
pair rejects it. Otherwise the existing exhaustive good-order algorithm
permutes only intrinsic light gaps, computes exact weighted backward cost,
and accepts exactly the least numerical good record. H is represented by
100, above every kernel budget; L is represented by 3, using the inherited
proof that its value cannot affect ownership. No global isomorphism seen-set
or unbounded cache is needed.

The independently reproducible C4 split is:

| (q,delta,u) | 1 | y | y^2 | y^3 |
|---|---:|---:|---:|---:|
| (0,0,3) | 100,080 | 102,840 | 35,220 | 4,020 |
| (0,1,2) | 37,560 | 26,280 | 4,590 | 0 |
| (0,2,1) | 4,150 | 1,510 | 0 | 0 |
| (0,3,0) | 125 | 0 | 0 | 0 |
| (1,0,1) | 360 | 120 | 0 | 0 |
| (1,1,0) | 60 | 0 | 0 | 0 |

Repeating the entire procedure for j=2,3 recovers the independently computed
C2,C3 in the preceding packet. C4 is a finite-exact computation with the
same justified owner predicate, rather than a hand classification of every
small collision.

## 4. C5: only one possible local clone obstruction

The activation bound `j<=6-2q` forces q=0 at j=5. Thus all endpoints have
weight one, and `(delta,u)=(0,3),(1,2),(2,1),(3,0)`. Call an endpoint marked
if it lies immediately beside an exceptional gap. Every isolated matching
edge must have a marked endpoint; all nonmatching components are activated
automatically.

**Proposition.** Among activated j=5 records, the only possible defect edge
internal to an intrinsic light gap is one isolated matching edge joining
the two endpoints immediately beside one exceptional gap. There is at most
one such gap. All records are pair free and self-owned, except exactly
those in which this matched gap is empty.

**Exhaustion of the first assertion.**

* At delta=0 there are five matching edges and at most six marked endpoints.
  With five marked endpoints every edge crosses from a marked to an
  unmarked endpoint. With six, the exceptional gaps are three disjoint
  internal endpoint pairs; precisely one matching edge has both endpoints
  marked. It is internal to a light gap only when it joins one such pair.
* At delta=1 there are a P3 and three matching edges. At least three of the
  at most four marked endpoints belong to the matching edges. If exactly
  three do, every matching edge has an outside partner and at most one P3
  vertex is marked. If four do, no P3 vertex is marked and only one matching
  edge has both endpoints marked. Again only its own exceptional gap can
  make that edge internal.
* At delta=2 there are either one four-vertex tree and two matching edges,
  or two P3 components and one matching edge. There is one exceptional gap.
  In the former case its two neighbors must belong to different matching
  edges. In the latter it must touch the matching edge, so it cannot touch
  two P3 vertices; only the two ends of that matching edge can form an
  internal defect.
* At delta=3 there are no exceptional gaps. Primitivity excludes isolated
  matching edges, leaving a four-vertex tree and a P3. Every light gap is a
  singleton.

All remaining defects cross distinct heavy gaps. If there is no internal
defect, this fixed crossing cost is the entire budget five, and every light
gap is transitive in its displayed order. Any other order adds positive
cost. Hence the good order is unique. The exact homogeneous-pair criterion
also excludes every clone pair: a potentially adjacent pair in a light
cluster has distinct external defect partners; different light gaps are
separated by a heavy witness.

If the one matched exceptional gap is empty, its two neighbors are adjacent
homogeneous vertices, so the record is rejected. If it has weight v>0, the
three positions form an unported cyclic triangle with weights `(1,v,1)`;
all other crossing defects cost four. For v=1 its three cyclic good orders
produce the same numerical endpoint-gap record. For v>=2, the only good
cyclic order keeps v off the defect edge, again producing that record.
Noncyclic orders cost more than one inside the cluster. This triangle is
pair free, and its three vertices have identical external comparisons.
Consequently it introduces neither an owner duplicate nor a new clone.
This proves the proposition.

Write b=0 or 1 for the number of internally matched exceptional gaps.
An activated placement therefore contributes exactly

```text
(3+y)^(u-b) (2+y)^b.                                      (6)
```

The factor 3+y allows gap values 0,1,2,L; the factor 2+y excludes only 0.
The small exact placement census is:

| delta | u | b=0 | b=1 |
|---:|---:|---:|---:|
| 0 | 3 | 20,160 | 2,520 |
| 1 | 2 | 27,720 | 2,520 |
| 2 | 1 | 12,180 | 630 |
| 3 | 0 | 1,680 | 0 |

These integers enumerate component graphs and exceptional-gap positions,
without any permutation or ownership search. Substitution in (6) proves
the stated C5 polynomial and its factorization in (4). Selected records
from every `(delta,b)` sector, with every assignment of 0,1,2,L to their
exceptional gaps, are independently passed through the full owner filter;
all 1,352 tests agree. A separate complete C5 owner replay also passed:
15,435 defect graphs and all 1,988,280 primitive records produce exactly the
same polynomial, in 215.0 seconds. Thus the shortcut is checked exhaustively
without relying on its proof.

## 5. C6: maximal activation makes every order rigid

At j=6 equality holds in the activation bound. Thus q=0; every nonmatching
component is P3; no exceptional gap touches a P3 vertex; and the u=3-delta
exceptional gaps are disjoint internal pairs of endpoint positions. Each
of their 2u marked endpoints has a distinct matching partner outside the
marked set. The other 3delta endpoints form delta P3 components.

Every defect goes between different intrinsic heavy gaps. Their fixed cost
is already six, so all transitive light clusters have a unique good order.
Distinct outside partners prevent clones. Every assignment of 0,1,2,L to
the u exceptional gaps is consequently accepted, with weight `(3+y)^u`.

The exact placement counts follow directly:

| delta | u | count | derivation |
|---:|---:|---:|---|
| 0 | 3 | 60,480 | `binom(9,3) 6!` |
| 1 | 2 | 90,720 | `binom(9,2) binom(7,3) 3 4!` |
| 2 | 1 | 45,360 | `9 binom(8,6) 90 2!` |
| 3 | 0 | 7,560 | `9!/(3!^3 3!) 3^3` |

Here 90 partitions six chosen vertices into two unordered triples and
chooses the center of each P3: `10*3^2`. The factors choosing disjoint
internal endpoint pairs are elementary path-matching counts. All labels
are ordered endpoint positions, so no hidden graph-isomorphism divisor
is present.

It follows that

```text
C6(y) = 60480(3+y)^3 + 90720(3+y)^2
        + 45360(3+y) + 7560
      = 7560(2y+7)^3.
```

The script additionally constructs 2,787 explicit full-owner controls,
covering every delta sector, every center choice in the selected P3
placements, multiple exact outside-partner bijections, and every gap
symbol. No local-gap isomorphism replacement is used in these controls.

## 6. Reproduction, independent checks, and remaining limit

Run from the repository root:

```text
python 05-knowledge/results/kernels_games_detours_20261020f/coefficients_check.py
python -O 05-knowledge/results/kernels_games_detours_20261020f/coefficients_check.py
python 05-knowledge/results/kernels_games_detours_20261020f/coefficients_check.py audit-five
```

The normal receipt is [coefficients_check.out](coefficients_check.out).
The full C5 replay is [coefficients_five.out](coefficients_five.out).
No removable assertions are used. Normal and `-O` standard runs passed in
30.5 and 32.9 seconds respectively, with identical 1,855-byte output.

In addition to the graph-universe, old C2/C3, clone, and owner controls,
the script checks the exact rational polynomial algebra and equality of
(1) and (5) for k=2,...,40. It also checks the independent all-L matching
universe identity

```text
sum_(j=2)^min(6,k) [y^3]Cj binom(2k-3,2j-3)(2k-2j-1)!!
    = (2k-1)!! binom(2k+1,3).
```

On the right, select an arbitrary matching of 2k ordered endpoints and
three L gaps. Their positive weights make all endpoint positions rigid.
This checks the leading coefficient of every kernel polynomial without
using the owner algorithm. The k=1 discrepancy is explicitly required as
a hostile boundary control rather than silently extrapolated away.

The remaining unrestricted obstacle is structural: fixed deficiency admits
a bounded exceptional interface, but unrestricted counting requires all
deficiencies or another exact way to sum their interfaces. Completing C4,
C5,C6 removes the numerical gap for D3; it does not justify erasing ports
or claiming a scalar closure for the full tournament sequence.
