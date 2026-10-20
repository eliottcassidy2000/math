# Summing deficiency-three owners through bounded external-port kernels

**PROVED below:** a finite port-kernel sum for the entire deficiency-three
stratum, and consequently an eighth-order polynomial-coefficient recurrence.
The normalized count is a polynomial of degree exactly seven. This is
derived from an exhaustive structural grammar, not fitted to small counts.
**FINITE-EXACT:** the first two kernel polynomials and explicit hostile
completion controls. The remaining three kernel polynomials are finite
uncomputed objects; a fully numerical all-k formula is not claimed.

## Inheritance, scopes, and the interface retained

Use the pair-free chain cores, exact clone-merger criterion, six D3 menus,
and least-good-record ownership theorem in
[`deficiency_three_20261020/ownership.md`](../deficiency_three_20261020/ownership.md).
The current hostile is two matchings with identical induced light gaps but
different external partners: local gap types alone cannot determine the
owner. The repaired near miss is that clone contraction can decrease
deficiency, already on two light vertices at k=1. We use the direct
pair-free owner test, so no merged presentation is counted twice.

The least-used incoming comparison is the finite recognition grammar in
[`prime_batch_triangle_20261020`](../prime_batch_triangle_20261020/README.md):
an exceptional component must retain its actual outside incidences and
ports. Its triangle modules are **not** transitive-chain weights. The
rational SCC strata in
[`scc_basins_20261020d/strata.md`](../scc_basins_20261020d/strata.md) show a
parallel strategy: a bounded exceptional alphabet can be summed around an
unbounded background. Here the background is a matching, not an SCC word,
and the resulting factor is a double factorial, not a geometric series.
The strong-centre bounds in that packet do not imply any bound used here.

Our live objects are canonical weighted chain cores, heavy-gap owners,
exceptional endpoint gaps, external matching ports, and matching insertion
slots. The map below erases unmarked matching components only for the
purpose of counting completions. It does not preserve the full tournament,
Hamiltonian paths, modules, or blue/black destinations. It preserves the
pair-free predicate and the comparison between all good owner records.
All erased endpoint positions and matching partners are counted explicitly
when the map is inverted.

Throughout D means feedback deficiency `2k+1-h`, not black half-edge
deficit, symmetry slack, or a strong-centre generating series.

## 1. The resulting finite sum and recurrence

Let N3(k) count capped weighted-core isomorphism classes with feedback at
most k and exactly `h=2k-2` heavy coordinates, where the cap is k+1.
The inherited small values are `N3(1)=5`, `N3(2)=756`, and `N3(3)=36473`.

There are five finite polynomials C2,...,C6 with nonnegative integer
coefficients and degree at most three such that, for every k>=2,

    N3(k) = sum_(j=2)^min(6,k) Cj(k-2)
               binom(2k-3,2j-3) (2k-2j-1)!!.                 (1)

Here `(-1)!!=1`. The polynomials count the explicitly defined finite
primitive port kernels in section 2, so (1) is an effective terminating
description independent of k, even before every coefficient is computed.

Using falling factorials `(x)_d=x(x-1)...(x-d+1)`, (1) becomes

    N3(k)/(2k-3)!! = P(k-2),

    P(x) = sum_(j=2)^6 [2^(j-2)/(2j-3)!] (x)_(j-2) Cj(x).   (2)

Terms j>k vanish automatically. The polynomial P has degree exactly seven
and leading coefficient `8/3`. Consequently, with a_k=N3(k),

    sum_(i=0)^8 (-1)^(8-i) binom(8,i)
         [product_(v=i)^7 (2k+2v-1)] a_(k+i) = 0,            (3)

for every k>=2, with the empty product equal to one.

Thus an exact recurrence is obtained for this whole fixed-deficiency layer.
It is not a recurrence for unrestricted A000568, and it does not yet supply
all initial values needed to execute (3). Computing C4,C5,C6, or enough
independent initial values once the degree theorem is known, is a finite
remaining arithmetic task. N3(1) is handled separately because it has no
heavy vertex and no bulk insertion slot.

## 2. Extracting the primitive external-port kernel

Use the inherited good endpoint record

    (q, delta, u, endpoint weights, gap weights, defect edges),
    q=k-r, delta=2r-t, 2q+delta+u=3.

The defect graph has no isolated endpoint. Activate the following entire
defect components:

1. Every nonmatching component.
2. Every matching edge with a nonunit endpoint weight.
3. Every component meeting an endpoint immediately beside an exceptional
   endpoint gap (an empty or light gap, rather than a heavy one).

Equivalently mark those endpoints and close under defect adjacency. All
unactivated components are unit matching edges, and neither endpoint is
beside an exceptional gap. This definition keeps all outside partners of
an exceptional cluster; none of its ports is replaced by a mere count.

Delete the unactivated endpoints and compress each remaining all-heavy
endpoint interval to one heavy gap. Retain the order of active endpoints,
every defect edge between them, endpoint weights, and every exceptional
gap with its weight. Call this ordered record the primitive kernel.
Its defining activation procedure activates every endpoint, so there are
no arbitrary extra marked matching edges.

Let e be the number of active defect edges, a the active endpoint count,
and put `j=e+q`. Deleting a matching edge removes two endpoints, so delta
is unchanged. The kernel therefore has exactly

    b=a+1-u=2j-2                                             (4)

heavy gap slots. If `s=k-j` unit matching edges were deleted, their
`2s` endpoints are inserted into these b slots. Each insertion of v
endpoints expands one heavy gap into

    H, endpoint, H, endpoint, ..., endpoint, H.

Thus an insertion allocation is a weak composition of 2s into b parts.
The inserted endpoints are then paired by an arbitrary matching. There
are precisely

    binom(2s+b-1,b-1) (2s-1)!!                              (5)

such completions when b>0. They reconstruct numerical ordered records
uniquely. Exceptional gaps force zero inserted endpoints at their own
positions; every other compressed gap is exactly one allowed slot.

### The kernel budget is bounded by six

For a connected nonmatching component with defect charge d=2e-v,
connectedness gives v<=d+2<=3d and e<=2d. Hence the total edges in all
nonmatching components are at most 2delta. At D3, q is zero or one;
at most q decorated matching edges occur. Each of the u exceptional gaps
touches at most two other matching edges. Consequently

    e <= 2delta + q + 2u = 6-3q,
    j=e+q <= 6-2q <=6.                                    (6)

All active endpoints together number at most `3delta+2q+4u<=12`.
The lower possibility j=1 has b=0, so has no nontrivial bulk completion
and contributes only at k=1. For k>=2 only j=2,...,6 occur.

This is a finite ordered-port grammar: at most twelve endpoints, at most
six defect edges, at most three exceptional gaps, and bounded endpoint
weights. No unbounded core bank is an input.

## 3. Pair-freeness is unchanged by bulk completion

Recall the exact homogeneous-pair criterion for ordered defect positions
i<j:

    (N_E(i) symmetric_difference N_E(j)) minus {i,j}
          = {i+1,...,j-1}.                                (7)

**Lemma.** A primitive kernel is pair free iff any one of its bulk
completions is pair free. The choice of bulk matching does not affect this.

**Proof.** Every added endpoint has one backward neighbor, its added
matching partner; all added endpoint gaps are heavy. No old vertex has
a defect incidence to an added endpoint. An old heavy gap is replaced
by an alternating string of new heavy vertices and matching endpoints.

An old homogeneous pair not meeting the replaced heavy gap stays
homogeneous: it cannot straddle that gap, since the old unmarked heavy
vertex would contradict (7). A homogeneous pair meeting the old heavy
vertex survives using the new heavy vertex nearest its old partner.
Conversely, any old witness distinguishing two old vertices remains;
the old heavy witness can be replaced by a new heavy in the same slot.

There is no new homogeneous pair involving an added matching endpoint.
For a pair with its own partner, (7) would require adjacency, but a heavy
vertex separates their positions. For a pair with a different vertex,
its matching partner would have to lie inside the interval in (7), and
then a separating heavy vertex also lies inside that interval. This
unmarked interior vertex violates (7). Two new heavy vertices likewise
cannot merge across an added endpoint. A new heavy paired with an old
vertex reduces to the old heavy comparison; any merger is already an old
merger at the appropriate boundary. This proves both directions. QED.

The compression is therefore an exact recognition device, not a claim
that the alternating insertion string is a transitive chain or a module.
Its internal matching can be complicated and is counted by (5).

## 4. The owner predicate is independent of the bulk matching

Every added endpoint lies alone between two consecutive heavy vertices.
It is therefore fixed in every good order. All old heavy vertices also
keep their order. Candidate good orders of a completed record correspond
bijectively to those of its primitive kernel: they permute only the same
old light gaps. Each completion adds exactly s to every such order's
backward cost, since its s matching edges have unit weights and every
other new comparison is forward. The residual kernel budget is j=k-s.

It remains to check the **lexicographic** owner convention; a bijection
of candidate orders alone would not suffice. We retain exactly the
inherited tuple order `(q,delta,u,ew,gaps,edges)`.

* The first three coordinates are unchanged by insertion: r gains s,
  t gains 2s, and u is unchanged.
* If q=0, every endpoint weight is one, so endpoint-weight words agree.
  If q=1, `delta+u=1`. When u=0 all light gaps are singletons; there is
  no local permutation. When u=1 there is at most one nonsingleton light
  gap. For two compared records with the same q and delta, the endpoint
  count in that one gap agrees. Inserting bulk unit endpoints outside it
  therefore preserves comparison of the endpoint-weight words, including
  the position of the possible single weight-two endpoint.
* Each kernel gap word has the same ordered heavy symbols. Bulk insertion
  replaces its i-th heavy symbol by a fixed positive number of copies.
  The heavy symbol is larger than every exceptional weight. This common
  expansion preserves lexicographic comparison: at the first difference,
  either two smaller symbols differ, or a heavy symbol is compared with
  a smaller symbol. Previously matched heavy runs expand identically.
* If endpoint-weight and gap words both agree, they specify the same
  endpoint positions in the full order. The old defect-edge labels are
  carried by one common increasing embedding. The added matching edges
  are identical in both candidate records. Removing these common edges
  from the sorted comparison leaves precisely the old comparison, and
  the increasing embedding preserves it.

Thus the relative order of **every pair** of good candidate records is
unchanged. A record is its canonical owner iff its primitive kernel is
its canonical owner at residual budget j. This proves the needed bulk
invariance, including cross-family changes in endpoint roles.

Distinct completed owner records cannot represent the same weighted core,
by the inherited ownership theorem. Their primitive extraction is unique.
Therefore (5) counts actual different classes for each accepted kernel,
without an isomorphism divisor or a global deduplication table.

## 5. Four gap symbols suffice

Use gap symbols `0,1,2,L`, where L denotes a light weight in `3,...,k`;
H denotes a heavy weight. Endpoints have weights one or two. A good D3
record has q<=1 and every defect endpoint weight is at most q+1<=2.
Consequently an L vertex can never be a defect endpoint in any good order.

Every comparison incident to such a vertex is therefore forward in every
good order. Its relative position to every other vertex is fixed. Replacing
its value by any other value at least three changes neither the set of
good orders nor the pair-free predicate. It also preserves the owner
comparison: endpoint words do not contain these vertices, and gap words
contain the same ordered list of L vertices. When two L entries align
before the first differing position, they are the same vertex. At the
first difference an L can only be compared with 0,1,2 or H, whose order
relative to it is fixed. Two different L values cannot decide that first
difference.

Each accepted symbolic kernel with ell occurrences of L therefore has
exactly `(k-2)^ell` independent light-weight realizations. Different
realizations are different owner records. Define Cj(y) to count primitive,
pair-free, self-owned symbolic kernels of budget j, with y per L gap.
This makes Cj a finite polynomial with nonnegative integer coefficients
and degree at most u<=3. To compute a kernel, use weight three for L and
a distinguished larger prohibitive value for H; H remains a separate
symbol even when j is small. In the implementation H=100 and j<=6.

Equations (4) and (5), together with this weight multiplicity, now prove
(1) for every k>=2. At k=2 the specialization y=0 excludes L entirely.

## 6. Degree, leading coefficient, and the recurrence

Direct factorial cancellation gives

    binom(2k-3,2j-3)(2k-2j-1)!!/(2k-3)!!
        = [2^(j-2)/(2j-3)!] (k-2)_(j-2).

This proves (2) and degree at most seven. To determine the highest term,
only the y^3 coefficient of C6 can contribute. Here u=3, q=delta=0:
the defect graph is a matching with six edges on twelve active endpoints.
All six edges must meet the endpoints beside the three exceptional gaps.
Those gaps must be internal and have disjoint endpoint pairs, and every
one of the six marked endpoints must be matched outside that marked set.
There are

    binom(9,3) * 6! = 60480                                (8)

such kernels: choose three disjoint adjacent pairs in the twelve-position
path, then match those six marked positions bijectively to the other six.

Every such kernel is pair free. Its three exceptional light clusters are
transitive, and all six backward edges go between distinct intrinsic
heavy gaps. The fixed crossing cost already equals the budget six;
any nontrivial internal reordering adds positive cost. Thus its good order
is unique, so every kernel in (8) is accepted. Conversely activation of
all twelve endpoints forces exactly this configuration. Therefore

    [y^3]C6=60480,
    [x^7]P = (2^4/9!)60480 = 8/3.

The eighth forward difference of the degree-seven polynomial in (2)
vanishes. Multiplying that identity by `(2k+13)!!` gives (3). This is
the source of the recurrence; no extrapolation from N3(1),N3(2),N3(3)
enters its proof.

## 7. Exact small kernel coefficients and controls

The modest finite computation enumerates only primitive symbolic kernels
at j=2 and j=3, rather than entire strata at large k. It finds

    C2(y)=756+792y+270y^2+30y^3,
    C3(y)=13380+12974y+4140y^2+435y^3.                       (9)

The explicit universes contain 2244 and 34698 primitive records before
the pair-free and owner filters. Every test uses the inherited exact
homogeneous-pair criterion and all good orders in the intrinsic light
gaps. The filtered counts in (9) imply

    N3(2)=C2(0)=756,
    N3(3)=3 C2(1)+C3(1)=36473,

agreeing with the independent earlier full capped-class census.

Additional controls choose accepted and rejected primitive records from
each available `(q,delta,u,L-count)` category. They enumerate every
completion by one or two additional matching edges, preserving actual
endpoint locations and partners: all 7480 tests retain the owner decision.
They also replace L entries independently by 3,4,7; all 469 substitution
tests retain it. No removable Python assertions or external dependencies
are used.
The leading-degree mechanism is checked on 168 explicit twelve-endpoint
port kernels, and 100 exact basis differences check the recurrence algebra.
Normal and `-O` runs pass with byte-identical 501-byte output in
[`owners_check.out`](owners_check.out); their combined replay took about
ten seconds on this workspace.

Reproduce from the repository root:

```text
python 05-knowledge/results/rotation_owners_spectrum_20261020e/owners_check.py
python -O 05-knowledge/results/rotation_owners_spectrum_20261020e/owners_check.py
```

The decisive surviving limitation is numerical, not a missing unbounded
interface: C4,C5,C6 still need evaluation. The external-port kernels remain
necessary; replacing them by local gap isomorphism types would identify
the inherited distinct-matching hostile and invalidate the sum.
