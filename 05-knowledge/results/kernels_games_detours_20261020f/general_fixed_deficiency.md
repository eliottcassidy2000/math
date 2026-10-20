# Every fixed feedback-deficiency layer has a finite kernel sum

**PROVED below from the inherited canonical-chain and exact weighted
feedback mechanisms:** a finite symbolic-kernel formula for every fixed
deficiency D in the stable positive-heavy range. Its double-factorial
normalization is a polynomial of degree exactly `floor(5D/2)`, with explicit
leading coefficient. This is a structural theorem, not a fit to D=0,1,2,3.
The two small-budget boundary strata are stated separately. No unrestricted
scalar recurrence for A000568 is claimed. **FINITE-EXACT:** the new-owner
comparisons and port-preserving completion controls in section 8.

Dependencies are the unique transitive-chain core and weighted feedback
formula in [the sparse theory](../defect_boundary_20261020/sparse_theory.md),
the exact endpoint/gap deficiency identity and least-record counting
principle in [the ownership packet](../deficiency_three_20261020/ownership.md),
and the matching-completion pair-free argument in
[the preceding kernel packet](../rotation_owners_spectrum_20261020e/owners.md).
The extension below replaces the old endpoint-word owner with a full-position
owner. It does not assume the old lexicographic insertion argument extends
to arbitrary numbers of exceptional gaps.

## 1. Statement

Fix an integer D>=0. Let `m=floor(D/2)` and `M=m+1`. For each k>=M, let N_D(k)
count capped weighted pair-free chain cores of weighted feedback at most k
with exactly `h=2k+1-D` heavy coordinates, where the cap is k+1. There are
finite polynomials `C_(D,j)(y)` with nonnegative integer coefficients such
that

```text
N_D(k) = sum_(j=ceil(D/2))^min(2D,k)
         C_(D,j)(k-M) binom(2k-D,2j-D)(2k-2j-1)!!.          (1)
```

Use `(-1)!!=1`. For D=0 include the empty-endpoint kernel j=0, consisting
of one heavy vertex. It has `C_(0,0)=1`, so `N_0(k)=(2k-1)!!` for k>=1;
the separate k=0 value is 1 as well.

Each coefficient polynomial can be computed from kernels with at most 4D
endpoints and 2D defect edges. In particular it is a finite object whose
size does not depend on k. Define

```text
B_D(k) = (2k-2m-1)!!.
```

There is a polynomial P_D(k) such that, for every k>=M,

```text
N_D(k) = B_D(k) P_D(k),
degree(P_D) = floor(5D/2),
leading_coefficient(P_D) = 2^floor(3D/2) / D!.              (2)
```

Thus each fixed-D layer has an effective finite polynomial-coefficient
recurrence. Equation (2) does not claim the finite polynomial coefficients
are easy to evaluate for large D.

For D=3 the numerical N3 formula in [coefficients.md](coefficients.md) is
consistent with (2): degree seven and leading coefficient 8/3. **The new
owner can redistribute primitive kernels between budgets j.** Therefore
the C2,...,C6 table computed using the old owner is not asserted to be the
table `C_(3,j)` for the new owner. Both constructions count the same N3.

## 2. Exact record identity and bounded active interface

Choose any good total order of a pair-free weighted core: its backward
weighted cost is at most k. Let r be its number of distinct backward edges,
t their endpoint count, and let u count exceptional endpoint gaps (empty
or light, instead of one heavy vertex). Consecutive vertices untouched by
backward edges would be a homogeneous pair, so each endpoint gap contains
at most one core vertex. Consequently

```text
q=k-r >=0,  delta=2r-t >=0,
h=t+1-u,
D=2q+delta+u.                                             (3)
```

A defect endpoint of weight w has an incident edge of cost at least w.
Compared with the unit cost of each of the r distinct edges this consumes
at least w-1 of the excess budget q. Hence

```text
w <= q+1 <= floor(D/2)+1 = M.                             (4)
```

Activate every nonmatching defect component, every matching component
having a nonunit endpoint, and every defect component touching an endpoint
immediately beside an exceptional gap. Keep every vertex and edge of an
activated component. The remaining components are unit matching edges.
Delete their endpoints and compress all-heavy intervals to one heavy slot.
This retains all external partners of the exceptional interface.

Let e and a be the numbers of retained edges and endpoints and set j=e+q.
For a connected nonmatching component of charge d=2e-v, connectedness and
d>=1 imply `v<=d+2<=3d` and `e<=2d`. Each decorated matching edge consumes
at least one unit of q. Each exceptional gap touches at most two other
matching edges. Therefore

```text
e <= 2delta+q+2u,
j <= 2delta+2q+2u = 2D-2q <= 2D,
a <= 3delta+2q+4u <= 4D.                                 (5)
```

The kernel heavy-slot count is

```text
b=a+1-u=2j+1-D.                                          (6)
```

Because the original good order has a heavy vertex, its compressed kernel
has one too, so b>=1 and j>=ceil(D/2). D=0 allows a=e=j=0. If s=k-j unit
matching edges were removed, reinsert their 2s ordered endpoints into the b
heavy slots. Inserting v endpoints expands one H into

```text
H, endpoint, H, endpoint, ..., endpoint, H.
```

The endpoint allocation is a weak composition of 2s into b slots. Pair the
inserted endpoints arbitrarily. Thus the exact completion count is

```text
binom(2s+b-1,b-1)(2s-1)!!
  = binom(2k-D,2j-D)(2k-2j-1)!!.                          (7)
```

This is a count of ordered numerical presentations. Sections 3-5 prove it
is also the correct count after rejecting clones and retaining one owner.

## 3. Pair-freeness and good orders under completion

The pair-free completion lemma in the preceding packet used only that
inserted endpoints have unit weight, have exactly one defect partner among
other inserted endpoints, and are separated by heavy vertices. It did not
use D=3. Its proof applies unchanged here: a homogeneous pair in the old
kernel survives using the appropriate boundary heavy vertex; an old
distinguishing witness survives; and any attempted new pair involving an
inserted endpoint would require its matching partner inside the interval
while a separating heavy vertex in that interval gives a contradiction.
Therefore a kernel is pair free iff every one of its completions is.

Every heavy vertex has all its incident arcs forward in every good order.
Their order and the intrinsic light gap occupied by every light vertex are
fixed. The inserted endpoints each occupy a singleton gap between two
heavy vertices, so they too are fixed in every good order. Good orders of
a completion correspond bijectively to good orders of the kernel. Every
order gains exactly s in backward cost from the s new unit matching edges;
all other new comparisons are forward. Its residual budget is j=k-s.

In particular the choice of the bulk matching cannot change which kernel
orders are good. It nevertheless remains in the completed edge data and
is counted, not discarded.

## 4. A different owner removes the arbitrary-D insertion problem

For a good order pi define its numerical record to be

```text
R(pi) = (q, complete ordered vertex-weight word,
         sorted backward-edge list on complete vertex positions). (8)
```

Take the lexicographically least such record as owner. This records every
vertex, including nonendpoint gap vertices; it is not the old tuple of
endpoint weights and endpoint gaps. The set of records is invariant under
weighted-core isomorphism, and a record specifies its complete ordered
weighted tournament. Consequently exactly one numerical owner record
corresponds to each weighted-core class, even when several orders attain
the same record.

**Insertion-invariance lemma.** For any fixed bulk completion, the relative
order of every pair of good kernel records under (8) is preserved.

**Proof.** The candidate orders can only permute vertices inside fixed
intrinsic heavy gaps. Hence each heavy occurrence has the same absolute
full position in every candidate order. Replacing the i-th heavy vertex by
its prescribed insertion string therefore gives one common increasing
embedding of old full positions, identical for every candidate order.

The q coordinate is unchanged: both budget and distinct backward-edge
count increase by s. If two q values differ, their comparison is preserved.
Otherwise the two full weight words have the same length, and the same
strings are inserted in the same positions. Their first differing entry,
if any, is preserved. If the weight words tie, the edge lists have equal
cardinality, because q is tied. The common position embedding preserves
the lexicographic order of old edge pairs, and the same new matching-edge
list is merged into both lists. Equal-cardinality sorted lists compare by
their first differing element, equivalently their least symmetric
difference; merging common elements cannot change this comparison. QED.

The q coordinate is a substantive guard. Without it, lexicographic list
order can reverse when one list is a proper prefix: `[1] < [1,3]`, whereas
adding common element 4 gives `[1,4] > [1,3,4]`. Keeping q prevents that
failure. A full adjacency bitword would be another possible repair, but
is not the owner convention used in this theorem.

Thus a record is an owner iff its primitive kernel is an owner at residual
budget j. No global isomorphism seen-set is required. Conversely every
completed owner has a unique extraction: the activation rule determines
exactly which endpoints were removed, their slots, and their matching.
This proves the quotient handling required for multiplying by (7).

## 5. Finite light alphabet and independent large-value choices

Use the alphabet

```text
gap values 0,1,...,M,L,H; endpoint values 1,...,M.
```

Here L stands for a light weight in `M+1,...,k`; H is a prohibitive heavy
symbol. A vertex with weight above M cannot be a defect endpoint in any
good order by (4). Every one of its incident arcs is consequently forward.
Its absolute position in every good order equals its indegree in the
unweighted core. Thus each L vertex is fixed at one full position in all
candidate orders.

Changing its weight to any other value above M leaves the good-order set
unchanged in both directions. Such a vertex remains excluded from every
backward edge, so no good cost depends on its value. The pair-free predicate
depends only on the underlying tournament. In owner comparisons its entry
lies in the same absolute position with the same value in both candidate
words, so it cannot decide their first difference. This proves invariance
under independent changes of all L values, including changes of their
relative sizes.

Every accepted symbolic kernel with ell L vertices therefore gives exactly
`(k-M)^ell` numerical owner records. Different assignments cannot be
identified by an automorphism: each L vertex has its own fixed position
in every good order. All endpoint values lie in the finite alphabet by
(4), and `ell<=u<=D`.

For the finite calculation one may represent L by M+1 and H by any larger
value exceeding 2D. Pair-free symbolic kernels at residual budget j can
be enumerated finitely from (5), and all good orders tested exactly. Define
`C_(D,j)(y)` as their owner count, with y per L vertex. It has nonnegative
integer coefficients and degree at most D. Combining this count with (7)
proves (1).

## 6. Polynomial degree and the unique top sector

For even D=2m, exact factorial cancellation gives

```text
binom(2k-D,2j-D)(2k-2j-1)!! / B_D(k)
  = 2^(j-m) (k-m)_(j-m) / (2j-2m)!.
```

For odd D=2m+1 it gives

```text
binom(2k-D,2j-D)(2k-2j-1)!! / B_D(k)
  = 2^(j-m-1) (k-m-1)_(j-m-1) / (2j-2m-1)!.
```

The falling factorial makes every j>k term vanish automatically in the
stable range. Since a kernel with parameters q,delta,u has degree at most u
and `j<=2delta+2q+2u`, its contribution to the normalized degree is at most

```text
floor(5D/2) - 4q - delta.
```

To attain the upper bound requires q=delta=0, u=D, j=2D, and all D
exceptional gaps of type L. The defect graph is a matching on 4D endpoints.
Maximal activation forces D disjoint internal adjacent endpoint pairs and
a bijection from their 2D marked endpoints to the other 2D endpoints.
There are exactly

```text
binom(3D,D)(2D)!
```

such records. They are pair free, and all 2D defects cross distinct heavy
gaps. The fixed crossing cost exhausts the budget 2D; every light cluster
is transitive, so the good order is unique and every record is its owner.
This includes the empty record at D=0. Thus the top coefficient is positive.
Dividing by the factorial in the appropriate normalization gives

```text
2^floor(3D/2) * binom(3D,D)(2D)! / (3D)!
  = 2^floor(3D/2) / D!.
```

No other sector attains this degree, so no cancellation can remove it.
This proves (2), including exact degree and leading coefficient.

As independent low-deficiency controls, the preceding counts give
`P_0=1`, `P_1=2k^2+4k-1`, and

```text
P_2=4k^5+12k^4-2k^3-(40/3)k^2+5k+10/3.
```

The last expression is a direct rewrite of the previously proved N2
formula in [the deficiency-two audit](../defect_transitions_20261020/strata_audit.md).
Together with the explicit P3 in [coefficients.md](coefficients.md), their
degrees are 0,2,5,7 and leading coefficients 1,2,4,8/3, exactly as (2)
predicts. These controls do not enter the general proof.

## 7. Boundaries and the unrestricted obstruction

The stable range k>=M is necessary for the alphabet specialization:
`k-M` counts actual available large-light values. For even D=2m, the
positive-heavy boundary k=m has only one heavy vertex and would give
`k-M=-1`; its value must be calculated separately. For odd D=2m+1, the
boundary k=m has no heavy vertex, no forced heavy-gap order, and no bulk
slot. It is also excluded. In D3 this is exactly the separately computed
N3(1)=5 boundary.

This is not a minor issue for unrestricted enumeration. Any fixed unweighted
pair-free core, with all weights one and k chosen at least its feedback
number and at least one, has h=0. Its deficiency is D=2k+1, and it lies at
`k=M-1`, outside the theorem. Thus arbitrary core complexity can remain in
the zero-heavy boundary. At a fixed budget k, the theorem covers the strata
with h>=2 and leaves h=0,1 as separate boundaries. Increasing k and D together
does not turn the fixed-D law into a bounded-interface formula for arbitrary
tournaments.

The positive conclusion is nevertheless stronger than a list of numerical
layers: every fixed-D stratum has a finite, computable exceptional grammar,
an exact bulk factor, a polynomial normalization of known degree, and a
known leading coefficient. The unresolved unrestricted task is to control
or sum the growing interfaces and the zero-heavy boundary, rather than to
repeat a separate all-k argument for each fixed deficiency.

For an explicit recurrence, set `d=floor(5D/2)` and `r=d+1`. For every k>=M,

```text
sum_(i=0)^r (-1)^(r-i) binom(r,i)
    [product_(v=i)^(r-1) (2k+2v-2m+1)] N_D(k+i) = 0.       (9)
```

This is the vanishing r-th forward difference of P_D, multiplied by
`B_D(k+r)`. The empty product is one. The finite kernel sum determines its
initial values, so (9) is an effective recurrence for a fixed layer even
when those initial values have not yet been numerically evaluated.

## 8. Exact new-owner controls

The independent ownership convention changes actual choices: among all
36,942 primitive D3 records at j=2,3, the old and new owner tests disagree
on 392 accept/reject decisions. Nevertheless their resulting C2,C3 tables
agree in this finite universe, and both reconstruct the independent
full-class counts N3(2)=756 and N3(3)=36,473. This finite agreement does not
assert all primitive coefficient tables are owner-independent.

The checker then retains accepted and rejected records from every available
`(q,delta,u,L-count)` sector. It enumerates all one- and two-edge bulk
completions of these records, including every matching of the inserted
endpoints, and checks the new owner decision directly. It independently
assigns each L vertex a value from 3,4,7 and also checks joint substitution
and insertion. The 7,105 completion/joint-operation controls and 313
substitution controls pass; 168 substitution cases have unequal
L values. The empty D0 kernel and all three of its two-edge matching
completions also pass. The failed proper-prefix simplification is a
required hostile control.

Reproduce from the repository root:

```text
python 05-knowledge/results/kernels_games_detours_20261020f/general_fixed_deficiency_check.py
python -O 05-knowledge/results/kernels_games_detours_20261020f/general_fixed_deficiency_check.py
```

See [the exact receipt](general_fixed_deficiency_check.out). The code uses
explicit exceptions, not removable assertions. Normal and `-O` runs passed
in 7.5 and 8.1 seconds, with byte-identical 468-byte output. These controls
exercise the new owner and its interfaces; the all-D conclusion follows
from the proof above, not from these bounded tests.
