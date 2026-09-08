# Independent audit: the uniform H4 type (3)(3) equality theorem

**Status: PROVED / INDEPENDENT ANALYTIC, SOURCE AND FULL REPLAY AUDIT PASS.**
The theorem is accepted for every finite ambient set: if all four H4 generators
have exactly two disjoint three-cycles and satisfy the stated chain 3,3,5 and
complementary commutations, they are all equal. Consequently their generated
action is never transitive. This is an abstract permutation theorem with a
separately paid actual-positive-meridian consumer. No full-tuple degree cutoff,
full-retention assumption, finite reflection quotient or Euler argument is used.

Auditor: `three_ray_geometry`, independently of `certificate_audit`, the producer.
I read the complete [h4_mixed33 proof](planar_jc48_sep08_h4_mixed33.md) and its
standalone source. Both independent normal and optimized replays reproduce the
frozen **448,965-gate, 3,111-byte** output exactly. A separate implementation
reconstructs every raw partner, table and centralizer orbit using a different
universe construction, direct letter predicates and union-find. No mathematical,
source or prose correction was required. This closes the earlier conditional
analytic acceptance: the finite pair inventory is now independently paid.

## 1. Exact theorem and the finite-to-unbounded bridge

Write the H4 relations as

    aba=bab, bcb=cbc, cdcdc=dcdcd,
    [a,c]=[a,d]=[b,d]=1.

Each generator has order three and support size six, but its two three-cycle
factors have no intrinsic order. The proof correctly avoids choosing a global
labelling of equal factors. Arbitrary additional fixed labels are permitted,
and the tuple need not be transitive in the abstract theorem.

Any pair of such permutations moves at most twelve labels in its support
union. Restricting to this union and padding with jointly fixed labels preserves
all exact words, the nontrivial joint orbits and the unordered cycle-support
intersection matrix. Simultaneous conjugation can put the first permutation in
the fixed form `(123)(456)` on twelve labels. Thus the second belongs to the
complete conjugacy class of size

    C(12,6) * [C(6,3)/2] * 2^2 = 924*40 = 36,960.

This proves the ambient-unbounded pair reduction. It does not put the full
four-generator tuple on twelve labels. In the later argument the other
permutations continue to act on the original ambient set, including the
jointly fixed labels outside the pair's support.

I checked that the producer includes every partner before imposing any relation.
There is no unrecorded support-overlap, transitivity, parity, retention or Euler
filter. Its unordered block choice has both cyclic orientations on each block.
The support matrix is only an observable; every relation is tested on the full
oriented permutations. Joint-orbit sizes refer to the actual generated action,
with jointly fixed singleton orbits omitted.

## 2. Independent complete pair and centralizer reconstruction

My independent construction first chooses a six-element union support, then
chooses the first three-element block to contain that support's least label.
The other block is its complement. Choosing both cyclic orientations on each
block yields exactly 36,960 distinct permutations. This differs from the
producer's successive choices of unordered disjoint triples from twelve labels.

For each partner p and fixed permutation s I evaluated the following direct
letter expressions at all twelve labels, without importing or calling the
producer's word/composition routines:

    braid3:  s[p[s[i]]] = p[s[p[i]]]
    braid5:  s[p[s[p[s[i]]]]] = p[s[p[s[p[i]]]]]
    commute: s[p[i]] = p[s[i]].

A union-find on the undirected edges of both cycle actions independently
recovers the nontrivial joint components. I computed unordered matrices by
minimizing over the two possible row orders and column orders, preserving all
orientations in the word test. The complete result is:

| Predicate | Canonical cells | Count | Nontrivial joint orbit sizes |
|---|---|---:|---|
| braid3 | (0,2,2,0) | 270 | 4,4 |
| braid3 | (0,2,3,0) | 36 | 3,4 |
| braid3 | (0,3,3,0) | 1 | 3,3 |
| braid3 | (1,1,1,1) | 270 | 8 |
| braid3 | (1,2,2,1) | 9 | 6 |
| braid5 | (0,1,1,0) | 3240 | 5,5 |
| braid5 | (0,1,3,0) | 180 | 3,5 |
| braid5 | (0,3,3,0) | 1 | 3,3 |
| braid5 | (1,2,2,1) | 18 | 6 |
| commute | (0,0,0,0) | 40 | 3,3,3,3 |
| commute | (0,0,0,3) | 160 | 3,3,3 |
| commute | (0,3,3,0) | 4 | 3,3 |

The predicate totals are respectively 586, 3,439 and 204. The sole ordinary
partner with joint orbit sizes (3,3) is exactly s. Therefore a non-diagonal
ordinary pair has one of precisely four joint-orbit types:

    (3,4), (4,4), (6), (8).

I also reconstructed the full centralizer decomposition independently. A
permutation commuting with s may rotate its two nontrivial cycles, exchange
them, and permute the six fixed labels arbitrarily, giving
`(C3 wreath S2) x S6`. Instead of the producer's eight generators, my traversal
uses four: rotation of the first triple, exchange of the triples, a six-cycle
on the fixed labels, and one adjacent transposition of fixed labels. These
generate the same full centralizer: conjugation by the block exchange supplies
the second triple rotation, and the last two generate S6.

For each raw p and generator g I formed its conjugate directly by the letter
rule `q[g[i]]=g[p[i]]` and merged the two indices with union-find. Every result
lies in the raw class. There are exactly 61 resulting orbits. All three
predicates, the unordered matrix and joint sizes are constant on each orbit.
Exactly 5,4,6 orbits satisfy braid3, braid5 and commutation, respectively;
weighting by their sizes reproduces every entry of the raw table above.
This independent reconstruction pays both the complete universe and the
centralizer gauge, including the relative orientations that split the
commuting cells into several orbits.

## 3. The three structural lemmas

For an odd braid of length 2r+1, put the rightmost-first action in the form

    (sigma*tau)^r*sigma = tau*(sigma*tau)^r.

If sigma fixes every label of a run `x,tau*x,...,tau^r*x` in a moved tau-cycle,
the two sides at x are `tau^r*x` and `tau^(r+1)*x`, a contradiction. This
argument allows repeated entries, so a whole small cycle outside the other
support is not an exception. In a three-cycle the ordinary r=1 relation forces
at least two labels to meet the other support, whereas the fifth r=2 relation
requires at least one. The proof never imports the ordinary two-label estimate
at the fifth edge. For an ordinary pair, summing over its two disjoint cycles
also gives support overlap at least four and enlargement at most two.

For commuting a,c, conjugation by c permutes the two nontrivial a-orbits. Since
c has order three, this permutation of two objects is trivial. Its restriction
to either three-element block is a power of that block's three-cycle. Hence
every c three-cycle is either an entire a-cycle support or disjoint from the
whole support of a. The complement is invariant as well. This is asserted only
for a commuting pair; the named ordinary/fifth hostiles refute its extension
to adjacent braid pairs.

For the centralizer lemma, suppose Gamma is generated by order-three elements
and is transitive on a set O preserved by an order-three permutation d
centralizing Gamma. If d fixes one point, commutation and transitivity force
it to fix all points. Otherwise its action is semiregular of order three.
This cannot happen on four or eight points. On six points a nontrivial d would
have two three-cycle orbits. Each order-three generator of Gamma preserves
these two orbits individually, since its induced action on two objects is
trivial. They would be proper Gamma-invariant subsets, contradicting
transitivity. Thus d also fixes a six-element joint orbit pointwise.

This lemma neither assumes d belongs to Gamma nor says that every centralizer
of a transitive group is trivial. Its setwise preservation hypothesis is paid
explicitly in the next step.

## 4. Uniform global implication

If a=b, the commuting relation for a,c and the ordinary relation for b,c force
c=a. The commuting relation for a,d and the fifth relation for c,d then force
d=a. Commuting elements satisfying an odd braid are equal by cancellation;
this does not depend on their cycle type.

Assume instead that a and b differ. Put A=supp(a), B=supp(b), and let
Gamma=<a,b>. The complete pair theorem gives at most two nontrivial joint
orbits, of types (3,4), (4,4), (6), or (8). Since d commutes with both a and b,
it sends a Gamma-orbit equivariantly to one of the same cardinality. It cannot
send a nontrivial orbit to a singleton. Its order is three, so it cannot
exchange the two nontrivial orbits; in type (3,4) their different sizes already
prevent exchange. Thus d preserves each nontrivial orbit, justifying the
centralizer lemma. Jointly fixed labels may still be moved among themselves
by d, and the proof correctly leaves them available.

In types (4,4), (6) and (8), d fixes all of A pointwise. Each individual c-cycle
must meet B in at least two labels, by the ordinary b,c edge. If both c-cycles
were disjoint from A, they would require at least four distinct labels in
`B minus A`, although the ordinary a,b edge gives at most two. Therefore one
c-cycle meets A. The commuting a,c lemma makes it an entire a-cycle. It is
fixed pointwise by d, contradicting the requirement that each c-cycle meet
supp(d) at the fifth edge.

In type (3,4), the union A union B has size seven, so `|B minus A|=1`.
Neither c-cycle can be disjoint from A, because each needs two labels of B.
Both are consequently the two complete a-cycle supports. One lies in the
four-element Gamma-orbit, which d fixes pointwise; this again contradicts the
fifth run condition. The common three-element orbit is allowed to carry a
nontrivial d-action and is not erroneously discarded.

All non-diagonal cases are impossible. Hence a=b and then a=b=c=d. The argument
is uniform in ambient cardinality and uses no full-tuple census. A general
tuple could initially have support union larger than twelve; none of the
steps rules this out by assumption. Equality is concluded from the intrinsic
pair orbit theorem and the other generators' original actions.

## 5. Minimal hostiles, equality and actual scope

I independently recomputed both named six-label controls:

    sigma=(123)(456),
    tau3=(124)(356),
    tau5=(124)(365).

The first pair satisfies braid3 and fails braid5; the second satisfies braid5
and fails braid3. Their generated transformation groups have respectively
orders 12 and 60 and both are transitive on six labels. My group closures
included inverse generators explicitly. All ten unordered partitions of six
labels into two triples were tested for each pair; none is invariant under
both generators. Thus these are genuine failures of a two-block interpretation,
not merely a poor ordering of the displayed factors. Six is the smallest
possible support size for this type.

Both have the same canonical unordered intersection matrix (1,2,2,1). The
orientation change therefore proves that this matrix alone loses information
needed even to select which odd braid holds. Neither hostile is asserted to
extend to an H4 tuple. The exact all-equal tuples, including fixed padding,
provide the correct equality controls for the final theorem and retain two
separate three-element orbits.

For the actual application, the proved
[mixed_cusp_braid supplier](planar_jc48_sep08_mixed_cusp_braid.md) supplies four
positive meridians generating the actual affine complement, satisfying the
marked H4 relations. Its (z,y,w,x) labels map to (a,b,c,d) here. A connected
finite covering has transitive monodromy; if those generating meridians all
have type (3)(3), the theorem would make their generated image cyclic with two
nontrivial orbits, impossible. Additional relations in the actual complement
cannot restore transitivity. The independently proved
[three_cusp_boundary family](planar_jc48_sep08_three_cusp_boundary.md) transports
the relevant actual marked supplier on its connected intrinsic good locus.
This consequence does not follow from singularity types alone, and the primary
does not claim an isomorphism of the complement group with the Artin group.
No assumption about full retention, sheet intersections, Euler characteristic
or a faithful action is used in this type-specific consumer.

The stated failure boundary is material. Three or more equal three-cycles
could admit a nontrivial order-three permutation of three blocks or three joint
orbits. Mixed type (4)(2) is not handled. The two unequal-cycle factors in the
separate (3)(2) result are not a proved dependency of this argument.

## 6. Source and frozen acceptance

The standalone standard-library source was read in full. All gates call an
explicit check that raises under both normal and optimized execution. Its
composition convention, independent direct letter tests, cycle extraction,
finite forward orbit closure and generated-group closure are coherent:
forward closure suffices because permutations generate a finite group. Every
centralizer edge is retained and checked against the unfiltered raw class.
The raw ordering and lexicographic orbit representative ordering make the
semantic hashes deterministic. No control substitutes a full-tuple finite
search for the unbounded structural proof.

Independent commands:

    python3 04-computation/planar_jc48_sep08_h4_mixed33.py
    python3 -O 04-computation/planar_jc48_sep08_h4_mixed33.py

Both pass 448,965 always-active gates and reproduce the same frozen 3,111-byte
output. In addition to those replays, the separate direct-letter/union-find
construction described in Section2 recovers all table and orbit data; it imports
no producer functions. Source, output and primary were not edited by this
auditor.

| Artifact | Bytes | SHA256 |
|---|---:|---|
| `04-computation/planar_jc48_sep08_h4_mixed33.py` | 10487 | `cc7da93f7326a3b3488dd48c63f52a22cb9c8b9a3d14810f9042beb04988f25c` |
| `05-knowledge/results/planar_jc48_sep08_h4_mixed33.out` | 3111 | `fad18d6a25dd81fed7286a8beb49881624aa4687f1da2d068e91f915cb446dde` |
| Primary proof at acceptance, before status promotion | 17315 | `006d1bbf916b00c13ab4f0c6bf7b876d098de5c7e1f3ae9f37c003049398ecf2` |

The frozen semantic raw-pair hash is
`db623f00fb368d3b9aea7fc8201bf3f62e4ee5db10f6eacdabc50dfffb0cc632`.
The complete centralizer-orbit hash is
`46983c0f3a5b015b5be6edf22893b76c68cdc0f919f0c0138d0a56950c8b85b4`.
Both are reproduced by both independent source replays.

This completes the final analytic, exact-inventory and replay audit. No
correction remains. Root owns subsequent primary promotion and integration;
those status changes do not alter the accepted source/output pins.
