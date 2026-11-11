# Fixed branching defect has finitely many transitive-inflation families

**PROVED IN NOTE; independently reviewed.** This strengthens the regular
language theorem in [PROOF](PROOF.md), whose positive defect identity and
finite residual-kernel bound are inherited here. No source-return closure
or uniform bound as the defect grows is claimed.

## 1. Exceptional singleton vertices are bounded

Let (T,r) be strong, I=N^-(r), d>=2, R=T-({r} union I), and
b=F_r(2)-H(T). Put delta=sum_{a<c in I}c_ac(T-r), so b>=delta.
The residual strong SCC list is bounded in terms of b. Its singleton
SCCs have intrinsic linear order. Their contact letters are subsets of I,
recording which inports they beat. Call the full subset the ordinary
letter; all other singleton contacts are exceptional.

The positive tagged-capacity theorem
[THM-4167](../../../../01-canon/theorems/THM-4167-tournament-exposure-capacity-deletion-support-moment-and-parity-holonomy.md)
assigns weight2 to each tagged odd cycle by itself. The following are
distinct atoms in the sum delta; none uses the marked root.

If d>=3 and a singleton x has mixed contacts, choose three inports, two
on one side of x and one on the other. As in
[the minimum-defect proof](../metagraph/MINIMAL_DEFECT.md), they give a
four-vertex directed path with inport endpoints and x internal. Close it
by the tagged two-edge ear to obtain a five-cycle. Distinct x give
different cycles because each such cycle has exactly one residual vertex.

If a singleton x loses to all I, strongness supplies a path from x to I.
Before this path first meets I, it must reach a vertex y of a LATER
residual SCC with y->c for some c in I: x cannot return to r, and every
residual SCC before x points toward x. Choose any a in I other than c.
The path a->x->y->c closes by its tagged ear to a five-cycle. Select one
y,c deterministically for each x. Different zero-contact x give distinct
cycles: the other residual vertex y beats an inport and therefore is not
a zero-contact vertex. These atoms also differ from the one-residual-vertex
mixed-contact atoms. Consequently, for d>=3,

    number of exceptional singleton vertices <= floor(b/2).       (1)

For d2 there are two mixed letters,10 and01. Two singleton vertices x<y
with the SAME mixed letter give a->x->y->c for the corresponding ordered
pair of distinct inports. Its tagged five-cycle has exactly that pair of
mixed residual vertices. Different pairs and the two different contact
types give distinct atoms. They are also distinct from the zero-contact
cycles, each of which has a zero-contact vertex. If m0,m10,m01 count these
three types, then

    m0 + binom(m10,2) + binom(m01,2) <= floor(b/2).                (2)

In particular all three numbers are bounded in terms of b. A convenient
nonsharp total bound is

    E_b=floor(b/2)+2 floor((1+sqrt(1+4b))/2).                     (3)

The argument counts atoms, not disjoint cycles in one packing. The same
auxiliary ear can occur in different summands of the capacity; demanding
the cycles be disjoint would incorrectly destroy this bound.

## 2. Only finitely many skeletons remain

Retain r, I, every nontrivial residual SCC and every exceptional singleton.
Between these residual objects retain a slot for a run of ordinary
all-one singletons. The residual SCC bound and(1)/(2) leave finitely many
possible skeletons. Every source has a unique such encoding before the
finite simultaneous isomorphism quotient on its inports and cores.

Split each slot into the cases empty or nonempty. Within a nonempty slot,
all vertices form a transitive module: they have identical contacts to
every vertex outside that slot. Changing its positive length is exactly
transitive inflation of one unmarked residual vertex. Strongness is
unchanged under such nonempty inflation, since the quotient has at least
the root and two inports and is strong exactly when the inflation is.

Every positive summand of the defect identity is nondecreasing under
this inflation. To inject a path or ordered two-path cover into the
one-vertex-longer block, insert a new smallest clone immediately before
the old smallest clone in its path. Any predecessor is external to the
block, because no later block vertex points to the smallest one. It has
the same incidence to the new clone. Specified starts/ends are in I, so
none is replaced by a clone. Deleting the new clone recovers the object.
Thus b is coordinatewise nondecreasing in all positive slot lengths.
This assertion need not hold when an EMPTY slot is first populated;
the separate zero/nonzero case split is essential.

## 3. Finite boxes, rational series, and eventual polynomials

The positive capped automaton of PROOF section3 recognizes exact defect b
and strongness for each skeleton. Let N exceed its number of states.
Fix the other slot lengths, and vary one positive length. Its DFA state
is eventually periodic by length N. But the actual integer defect is
nondecreasing, so the set of lengths giving exact b is an interval.
An eventually periodic interval has constant membership from N onward.
The same uniform N works for every prefix state and every suffix.

Therefore replacing any slot length >=N by N leaves acceptance unchanged,
independently of the other slots. Enumerate the finite grid of lengths
1,...,N, where N represents every length >=N. Accepted grid cells give
disjoint Cartesian boxes: each coordinate is either a specified integer
less than N or the infinite interval[N,infinity). This is an effective
complete classification, even though its finite automaton can be huge.

The finite isomorphism group acts on the skeleton and contact data and
preserves every ordered slot. It cannot exchange repeated strong SCCs
at distinct residual positions. Choosing one representative of each
finite skeleton orbit therefore leaves the same boxes; no length-dependent
permutation quotient is hidden. A box with s free lengths contributes
x^a/(1-x)^s to the ordinary ROOTED class series.

**Conclusion:** for every fixed b, all strong
rooted tournaments of defect b with at least two inports are a finite
disjoint union of explicitly computable transitive-inflation boxes.
Their order-generating series has denominator a power of1-x, and their
number at order n is eventually a polynomial in n. This improves mere
rationality: exponential growth and nontrivial eventual periods are absent.
The number of slots is at most one plus the bounded exceptional-singleton
count plus the bounded number of nontrivial residual SCCs.

The rooted automorphism group is constant throughout each positive-length
skeleton. Every automorphism fixes the intrinsic residual SCC order, its
singleton positions, and each rigid transitive slot. Its remaining finite
action on I and the nontrivial cores is independent of slot lengths.
Consequently the same eventual-polynomial conclusion holds after selecting
any specified rooted automorphism-group type. Labelled counts follow by
multiplying each ordinary rooted class by n!/|Aut(T,r)|; those counts are
not themselves asserted to be polynomial. If the root label is fixed in
advance, the factor is(n-1)!/|Aut(T,r)| instead.

The boundary d1 remains outside the theorem. Its inherited endpoint-budget
funnel grammar has genuinely different growth. No finite bank uniform in b,
nor a bounded-defect representative for every Hamilton count, follows.

## 4. Concrete controls and connection to the reset

At b2 and b4 the old pentagon and sixteen neutral tails have one free
ordinary prefix, giving eventual constants. At b5 the branching stratum
is finite. At b6 the residual triangle opens two independent slots and
the complete count becomes5n+2 for n>=8. These supply the first growth
transition predicted by this theorem, not its proof.

The uniform reset's target has defect24 and therefore belongs to one of
these finitely many boxes. Its particular slot is a transitive prefix,
and its residual core is S5a. Iterating its gadget increases defect to268,
3120,...; those steps leave each fixed-defect classification. This marks
precisely why a complete family dictionary at every fixed budget does
not yet supply an induction at unbounded Hamilton count.

The independent proof audit is recorded in
[the metagraph referee note](../metagraph/audit_templates.md).
