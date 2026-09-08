# Equal three-cycles and the complete H4 type (3)(3) obstruction

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
September 8, 2026. The two-cycle labelling shortcut is REFUTED. The
repaired conclusion is uniform in ambient degree.

## 1. Theorem and inheritance

Let a,b,c,d be permutations of a finite set. Each has exactly two
disjoint three-cycles, with all other letters fixed. Suppose

    aba=bab,  bcb=cbc,  cdcdc=dcdcd,
    [a,c]=[a,d]=[b,d]=1.                               (1)

Then

    a=b=c=d.                                           (2)

In particular their generated action is never transitive: the common
permutation has two separate nontrivial orbits of length three, besides
any fixed letters. There is no assumption of full retention, an Euler
identity, a faithful action, an involution quotient, or a degree cutoff.

The [actual mixed-cusp braid supplier](planar_jc48_sep08_mixed_cusp_braid.md)
pays a positive-meridian generating tuple with the H4 chain 3,3,5.
Its labels (z,y,w,x) become (a,b,c,d) in (1). Since the actual connected
cover has transitive monodromy, its hypothetical whole-support Keller
map cannot have positive meridians of type (3)(3), in any mapping degree.
The same marked consequence transports through the separately proved
connected intrinsic good family of
[the (5,3,3) cusp boundary](planar_jc48_sep08_three_cusp_boundary.md).
This uses the paid generating tuple, not merely a matching list of
singularities. It does not assert that the actual complement equals
the Artin group, nor that an abstract action realizes a Keller map.

The closest proved mechanism is the local odd-braid run argument in
[H4 single-cycle support](planar_jc48_sep08_h4_single_cycle.md). The
parallel [mixed (3)(2) analysis](planar_jc48_sep08_h4_mixed32.md) explains
why local factor information must be retained, but its theorem is not
a dependency here. Its cycle factors have different lengths and are
canonically extracted by powers. That operation cannot distinguish the
two equal three-cycles in the present problem.

The hostile below already occurs on six labels. The corrected near miss
is an assumed globally consistent ordering of the two equal-length
cycles. The least-used sidecar is the action of the fourth generator on
the **joint orbits of an ordinary adjacent pair**. The live concepts are
unordered cycle intersections, exact cyclic order, joint-orbit size,
centralizer semiregularity, and actual positive-meridian generation.
No tournament or arbitrary ambient-degree census is introduced.

Root assigned the type (3)(3) frontier after the (3)(2) result. The
bounded initial probe classified all pairs, found the factor-labelling
hostile, and led to the centralizer repair below. There is no external
group-classification or priority claim.

## 2. Minimal hostiles to a two-triple block interpretation

On six labels put

    sigma=(1 2 3)(4 5 6),
    tau3=(1 2 4)(3 5 6),
    tau5=(1 2 4)(3 6 5).                               (3)

The pair (sigma,tau3) satisfies braid3 and not braid5. Its generated
group has order twelve and is transitive on the six labels. The pair
(sigma,tau5) satisfies braid5 and not braid3; its generated group has
order sixty and is also transitive. These order statements are literal
finite group closures, not imported identifications with named groups.

Both pairs have the same canonical unordered intersection matrix

    [1 2]
    [2 1].

Neither pair preserves a partition into two triples. There are precisely
ten unordered two-triple partitions on six labels; the source tests all
ten for each pair and finds none invariant under both generators.
Consequently there is no way to rescue a common two-triple block system
merely by swapping the equal cycles' names. Six is the smallest possible
ambient cardinality for type (3)(3), so these block hostiles are minimal.

The first failed implication is local: a cycle of one generator need
not be an invariant factor for an adjacent generator. The relative
orientation of the second three-cycle changes braid3 to braid5 while
leaving the support intersection matrix unchanged. That matrix alone
therefore cannot decide an odd braid relation.

Neither pair is asserted to extend to an H4 tuple. The full graph
relations, rather than arbitrary cycle labelling, repair the obstruction.

## 3. Complete two-generator classification

Fix sigma=(1 2 3)(4 5 6) on twelve labels. For a partner tau, form the
two-by-two matrix of intersections of their three-cycle supports. Since
the cycles are unordered, its canonical representative is the
lexicographically least four-entry row tuple after independently swapping
the two rows and the two columns. The full permutations, including
cycle orientation, are retained when testing relations.

There are exactly

    (1/2)*binom(12,3)*binom(9,3)*2^2
       = binom(12,6)*40
       = 36,960                                         (4)

partners of type (3)(3). The source chooses two disjoint three-element
supports as an unordered pair and both cyclic orientations on each.
There is no preliminary support, braid, transitivity, parity, Euler or
retention filter.

Here is the complete matrix table. A displayed tuple lists the two
matrix rows in order.

| Relation | Canonical cells | Partners | Nontrivial joint orbit sizes |
|---|---|---:|---|
| braid3 | (0,2,2,0) | 270 | 4,4 |
| braid3 | (0,2,3,0) | 36 | 3,4 |
| braid3 | (0,3,3,0) | 1 | 3,3 |
| braid3 | (1,1,1,1) | 270 | 8 |
| braid3 | (1,2,2,1) | 9 | 6 |
| braid5 | (0,1,1,0) | 3,240 | 5,5 |
| braid5 | (0,1,3,0) | 180 | 3,5 |
| braid5 | (0,3,3,0) | 1 | 3,3 |
| braid5 | (1,2,2,1) | 18 | 6 |
| commute | (0,0,0,0) | 40 | 3,3,3,3 |
| commute | (0,0,0,3) | 160 | 3,3,3 |
| commute | (0,3,3,0) | 4 | 3,3 |

Thus braid3 has 586 partners, braid5 has 3,439 and commutation has 204.
The only ordinary partner with joint orbit sizes 3,3 is tau=sigma.

“Joint orbits” means orbits of the group generated by the two actual
permutations; jointly fixed singleton labels are omitted from the size
list. The source computes them by closure under both actions. Equivalently
their nontrivial components are read from the connected cycle-support
incidence graph with the private labels on each cycle retained. No
ordering of the two cycles is assumed.

### Why this finite universe proves an all-ambient pair theorem

Any two permutations of this type move at most twelve labels in their
union. Restrict to that union, remove jointly fixed padding, and add
jointly fixed labels if necessary to reach twelve. A simultaneous
relabeling puts the first generator in the fixed form. Its partner is
one of all 36,960 permutations in (4). Exact words, nontrivial joint
orbits and unordered intersection matrices are unchanged by this process.

Accordingly every **non-diagonal ordinary pair** in arbitrary ambient
degree has one of these four joint orbit types:

    (3,4), (4,4), (6), (8).                              (5)

This is the only finite classification needed by the global proof.
A general four-generator tuple is not assumed to fit into twelve
labels. Each ordinary pair separately has the stated intrinsic orbit
type in the original ambient set, and the proof retains the actual
action of the other generators on those orbits.

### Independent centralizer-orbit check of the full universe

The centralizer of sigma is

    (C3 wreath S2) x S6.

It permits independent rotations of the two cycles, their interchange,
and arbitrary permutations of the six fixed labels. These operations
generate the whole centralizer: a commuting permutation preserves
cycle length, maps each three-cycle to one of the two three-cycles with
a cyclic phase, and permutes the fixed points freely.

The source closes each orbit under these explicit generators, beginning
with the lexicographically least unused raw partner, and verifies every
image remains in the complete raw universe. There are **61** complete
centralizer orbits. Exactly **5, 4, 6** of them satisfy braid3, braid5
and commutation, respectively. Predicate, canonical matrix and joint
orbit-size invariance are checked on every orbit element. Weighting by
orbit sizes recovers every raw matrix count in the table.

In the commuting table, identical unordered support matrices can have
different relative cycle orientations: the shared-one-block row splits
into two centralizer orbits, and the shared-two-block row into three.
Those distinctions are kept in the full orbit bank; they are not removed
by an unjustified ordered-cycle gauge. Independent direct letter action
also checks every odd word equality against the permutation-composition
engine.

## 4. Three elementary structural lemmas

### 4.1 Odd runs on an individual three-cycle

For an odd braid of length 2r+1, with rightmost-first composition,
every run x,tau x,...,tau^r x in a moved tau-cycle meets supp(sigma).
Otherwise sigma fixes all its entries and the relation

    (sigma tau)^r sigma = tau(sigma tau)^r

evaluated at x gives tau^r x=tau^(r+1)x, contradicting that x is moved.
Repeated positions cause no exception.

Thus an ordinary braid pair of type (3)(3) has **at least two**
intersection labels in each individual three-cycle of either generator.
In a three-cycle, two outside labels would be adjacent cyclically.
Its total overlap is at least four and its support enlargement is at
most two. A fifth braid pair requires at least one intersection label
in each individual three-cycle. The latter is not the ordinary
two-label bound.

### 4.2 Commuting equal cycles are genuine whole blocks

Suppose a and c commute and both have type (3)(3). Conjugation by c
permutes the two nontrivial a-orbits. Since c has order three, that
permutation of two objects is trivial. Hence c preserves both a-cycle
supports individually. Its restriction to either three-element support
commutes with a three-cycle, and is either identity or a power of that
cycle.

Consequently every c three-cycle is either exactly one of a's whole
three-cycle supports or is disjoint from supp(a). This fact concerns a
**commuting** pair only. The six-label ordinary and fifth hostiles show
why it cannot be extended to an adjacent odd braid pair.

### 4.3 An order-three centralizer on small joint orbits

Let Gamma be generated by elements of order three, let O be a transitive
Gamma-orbit, and let d have order three and commute with Gamma. Assume
d preserves O setwise.

The action of d on O is semiregular unless it is trivial. Indeed, one
fixed point would imply every point is fixed, by transitivity and
commutation. A nontrivial d would therefore partition O into cycles
of length three. If |O| is four or eight this is impossible, so d fixes
O pointwise.

If |O|=6, a nontrivial d would have exactly two three-cycles there.
Every order-three generator of Gamma commutes with d and permutes these
two cycles; its induced permutation must be trivial. Both d-cycle
supports would be Gamma-invariant, contradicting transitivity on six
labels. Thus d also fixes a six-element transitive orbit pointwise.

The setwise preservation hypothesis will be paid below. This lemma
does not assume d fixes the jointly fixed labels outside O, nor claim
that every centralizer of a transitive group is trivial.

## 5. The uniform H4 proof

Put A=supp(a), B=supp(b), C=supp(c), D=supp(d), keeping D here as a
support set rather than an ambient degree.

First suppose a=b. The relation [a,c]=1 and the ordinary relation
bcb=cbc imply c=a. Then [a,d]=1 and the fifth relation cdcdc=dcdcd
imply d=a. In any group, commuting elements satisfying an odd braid
must be equal, since the commuting word powers cancel. This proves
the desired conclusion in the diagonal case.

Suppose for contradiction that a!=b. By the complete ordinary pair
classification, the nontrivial joint orbits of Gamma=<a,b> have sizes
in (5). There are at most two such orbits. Since d commutes with a and
b, it maps Gamma-orbits equivariantly to Gamma-orbits of the same size.
Its order is three, so it cannot interchange a pair of nontrivial
orbits. In (3,4), the two sizes differ anyway. A nontrivial orbit also
cannot map to one of the jointly fixed singleton orbits. Thus d
preserves every nontrivial orbit individually, paying Lemma 4.3's
setwise hypothesis.

If the orbit type is (4,4), (6) or (8), Lemma 4.3 makes d fix all of A
pointwise. On the other hand, each c-cycle meets B in at least two
labels by the ordinary b,c relation. If both c-cycles were disjoint
from A, their total intersection with B would be at most |B minus A|,
which is at most two by the ordinary a,b relation, against the required
total of at least four. Thus one c-cycle meets A. By [a,c]=1 and
Lemma 4.2, that cycle is a whole a-cycle and is fixed pointwise by d.
It is therefore disjoint from supp(d), contradicting the fifth c,d
relation on that individual c-cycle.

It remains to consider orbit type (3,4). Here the supports of a,b have
union size seven and each has size six, so |B minus A|=1. Every c-cycle
must meet B in at least two labels; hence no c-cycle can be disjoint
from A. Lemma 4.2 makes both c-cycles exactly the two a-cycle supports.
One lies in the four-element joint orbit. Lemma 4.3 makes d fix that
orbit pointwise, so this c-cycle again misses supp(d), contradicting
the fifth relation. The other common three-element orbit can carry a
nontrivial d action; the proof has not silently excluded it.

All non-diagonal ordinary pair types are impossible. Thus a=b and the
first paragraph yields a=b=c=d. This proves (2) uniformly.

There is no separate cycle-labelling compatibility assumption. The
common blocks used in the proof arise only from a paid commuting
relation; joint-orbit centralization then supplies the incompatible
fixed cycle. No full-tuple census, actual Euler ledger, finite reflection
group or power-map claim is used.

## 6. Consequence and exact stopping boundary

Taking all four generators equal to any (123)(456), with arbitrary fixed
padding, satisfies every H4 relation. This is the sharp abstract equality
control. Its two three-element orbits remain separate, so the action
is never transitive, even with no additional fixed labels.

For the actual marked mixed-cusp supplier, this excludes meridian cycle
type (3)(3) in every degree. The conclusion only concerns this type.
It does not settle (4)(2), more than two three-cycles, general mixed
orders, or JC(2). In particular, the centralizer proof would change with
three or more equal cycles: an order-three element could then permute
three joint orbits or three cycle blocks. That is an explicit failure
boundary of the present argument.

The source-to-target map first retains two actual adjacent meridians,
then their complete unordered cycle cells and joint orbits. It preserves
the exact odd relation and original support sets. The cell matrix alone
loses relative cyclic order, as (3) demonstrates. The needed sidecar is
the full joint action, and the map from the surrounding generator to its
centralizer on that action. The final proof uses actual orbit sizes
together with order three; it does not try to label factors globally.

The next bounded test at another cycle type should therefore retain
joint-orbit centralizers before attempting a canonical factorization.
The mixed (3)(2) route uses unequal lengths; the present route uses
equal lengths but a nontrivial small-orbit centralizer obstruction.
Neither mechanism justifies an arbitrary higher-degree census.

## 7. Reproduction and audit status

Source: [h4_mixed33.py](../../04-computation/planar_jc48_sep08_h4_mixed33.py).
Frozen output: [h4_mixed33.out](planar_jc48_sep08_h4_mixed33.out).

Run from the worktree root:

    python3 04-computation/planar_jc48_sep08_h4_mixed33.py
    python3 -O 04-computation/planar_jc48_sep08_h4_mixed33.py

The source is standalone standard-library Python. All gates use explicit
exceptions and remain active under optimization. The universe is exactly
all 36,960 partners in (4); the complete raw relation table is independently
checked by direct rightmost-first letter actions and by the 61-orbit
centralizer decomposition. The named six-letter hostiles retain all ten
unordered two-triple partitions and their complete generated groups.
Equal-generator controls are retained on 6,7,12 and25 labels.

Normal and optimized replays pass **448,965 always-active gates** and
are byte-identical, 3,111 bytes. The large gate count includes every
centralizer edge and every raw orbit-invariance check; it is not a
search over a range of tuple degrees.

Frozen pins:

| Artifact | Bytes | SHA256 |
|---|---:|---|
| Source | 10,487 | cc7da93f7326a3b3488dd48c63f52a22cb9c8b9a3d14810f9042beb04988f25c |
| Output | 3,111 | fad18d6a25dd81fed7286a8beb49881624aa4687f1da2d068e91f915cb446dde |

The semantic raw-pair hash is
db623f00fb368d3b9aea7fc8201bf3f62e4ee5db10f6eacdabc50dfffb0cc632.
The complete centralizer-orbit hash is
46983c0f3a5b015b5be6edf22893b76c68cdc0f919f0c0138d0a56950c8b85b4.

Three_ray_geometry independently challenged and accepted the analytic
centralizer/global proof, conditional only on the now-frozen complete
ordinary-pair inventory. The [independent complete audit](planar_jc48_sep08_h4_mixed33_audit.md)
accepts every support/centralizer step, the ambient-unbounded conclusion
and both frozen replays. A separate union-support-first/direct-letter
census reproduces all36960 partners, all61 centralizer orbits and the
complete admitted tables. All earlier artifacts are unchanged.

