# Canonical augmentation between anchor languages

**PROVED, using the inherited proved fixed-anchor owner
construction:** a finite-state canonical augmentation language for every
fixed feedback rank; length-preserving weighted transfers between successive
anchor ranks; a parent forest whose components have fixed maximum transitive
size and finite height. **PROVED elementary consequence:** unrestricted
A000568 is not P-recursive. **FINITE-EXACT:** all classes through order seven,
including the complete maximum-transitive-size-three column. An efficient
uniform unrestricted counting recurrence remains **OPEN**.

## Inheritance and the operation being added

The closest mechanism is the fixed-a canonical regular language in
[zero_heavy.md, section5](../zero_heavy_ports_20261020g/zero_heavy.md).
Its canonical hostile is an ordered join of directed triangles, with
unbounded feedback vertex number. The corrected near miss is counting
accepting witnesses rather than distinct canonical words. A second near
miss, exposed below by the regular pentagon, is assuming that an ordinary
canonical child record projects to an ordinary canonical parent record.

The least-used sidecar is the exact identity

    tau(T) = |T| - alpha(T),

where tau is minimum feedback **vertex** number and alpha is maximum
transitive subtournament size. These are not weighted feedback-arc budget
or heavy-coordinate deficiency. The current
[degree-local enumeration](../degree_local_tournament_enumeration_20261021/session22-enumeration/README.md)
already accelerates unrestricted Burnside arithmetic without changing the
future-response quotient. This note adds a different structural transition;
no speed comparison with that engine is claimed.

The concept board is canonical words, rank-changing augmentation,
transitive remainders, exact clone disagreements, finite columns, and
growth barriers. The anchor is canonical augmentation; the niche is the
finite alpha=3 column; the wildcard is the contrast between rational
fixed-rank rows and a non-P-recursive unrestricted diagonal.

## 1. Definitions and the cross-rank regularity test

An a-anchor record is (F,w), with ordered anchor tournament F and a word
over {0,1}^a. The word lists the transitive remainder in its unique order;
bit i is one exactly when anchor i points to that vertex. Let L_a be the
inherited regular language of lexicographically least records, one per
class admitting a feedback set of size a. A fixed encoding orders F before w.

We also need the regular predicate tau<=b on a-records, for fixed a,b.
The inherited re-anchoring construction extends without an equal-size
assumption. Retain r old anchors, select b-r old base vertices, and insert
the a-r displaced anchors into the new transitive base. Guess the finite
anchor labels, selected letters and insertion order. Check transitivity
against each surviving input letter. Selected input positions are exempt
from base tests. Comparisons with selected future anchors use their guessed
letters and whether their positions have been passed. Every possible
b-anchor set is represented in this way.

At any time the difference between consumed and emitted letters is bounded
by max(a,b); its final value is b-a. Thus the relation is finite-state,
with a bounded queue when output constraints or lexical comparisons are
required. Existence of such a record is regular by projection. In
particular tau=a is a regular predicate on a-records: remove the tau<=a-1
language. For a=0 every represented tournament is transitive.

This paragraph is an existence construction of a finite automaton. A
determinized general automaton has not been materialized in the experiment.

## 2. The canonical augmentation language

Fix a>=0. An augmentation record is an (a+1)-anchor record with the final
anchor x distinguished as the added vertex. Define pi by deleting x from
the header and deleting its bit from every word letter. The underlying
tournament pi(F,w) is literally T-x in those coordinates. Word length is
unchanged.

Let R_(a+1) contain exactly the augmentation records satisfying:

1. tau(T)=a+1;
2. the projected parent record pi(F,w) belongs to L_a.

Both are regular predicates, so R_(a+1) is regular. The parent automatically
has exact rank a: its record gives tau(T-x)<=a, while
tau(T)<=tau(T-x)+1 gives the reverse inequality.

Every class of rank a+1 has at least one R_(a+1) record. Take a minimum
feedback set A, choose x in A, and canonically re-anchor T-x at a anchors.
Appending x to those parent anchors gives the required child record.

Now retain only the lexicographically least R_(a+1) record of each child
isomorphism class. Call this language G_(a+1). To recognize its complement
inside R_(a+1), use the inherited bounded-delay re-anchoring transducer,
requiring its output to be in R_(a+1) and smaller than its input. The
existence-of-smaller-competitor language is regular. Determinization and
relative complementation therefore show that G_(a+1) is regular.

**Conclusion.** G_(a+1) has exactly one record per class of exact feedback
vertex number a+1, and its parent projection is a canonical exact-rank-a
record. Multiple deletions, anchor permutations and automorphisms that
produce the same record are not counted repeatedly. This is why the
construction counts deterministic accepted words, not witness paths.

The restricted augmentation owner can differ from the ordinary owner L_(a+1).
That distinction is essential, not a different name for the same choice.

## 3. A weighted finite-state transition, not just an existence relation

Choose a deterministic automaton for G_(a+1), with finite header branches.
Fix a parent header F. Each extension header chooses its a incident
arc bits, so there are 2^a possibilities. Incorporate their allowed initial
states in an integer row vector i_F. For a parent word letter b in {0,1}^a,
its two possible child letters are b and b+2^a. Let B_b be the sum of their
two deterministic transition matrices, and let f be the accepting vector.
Then

    h_a(F,b_1...b_l) = i_F B_(b_1) ... B_(b_l) f             (1)

counts exactly those child classes of rank a+1 whose canonical augmentation
parent is that record. Header choices remain distinct until their accepted
child words are counted. Since G_(a+1) has one word per child, the sum
introduces no isomorphism multiplicity.

**Recanonicalization between steps is required.** A child word in G_(a+1)
is generally not its ordinary L_(a+1) owner. Before applying h_(a+1),
transport it to that owner using the inherited fixed-rank re-anchoring
relation. Restricted to G_(a+1) on the input and L_(a+1) on the output,
this is a bijection at the canonical-record level. The bounded-delay
relation can be synchronized by a finite queue and recognized as a regular
language of paired words. Determinize and count distinct paired words,
not nondeterministic re-anchoring witnesses. This normalization is part of
the transition machinery; omitting it is exactly the pentagon failure in
section5. The physical label permutation comes from the selected-anchor
and displaced-anchor metadata when literal ports must also be transported.

Writing D_a(n) for the number of order-n classes of exact rank a, we obtain

    D_(a+1)(n+1) = sum_(T: |T|=n, tau(T)=a) h_a(owner_a(T)). (2)

For fixed a the right side can be evaluated by a finite product automaton
and matrix powers; no list of all parent cores at each growing order is
needed in that construction. The matrices themselves can be extremely
large, and their construction is not shown to be efficient in a. Formula
(2) is a structural weighted-language recurrence, not a scalar closed
formula for the unrestricted sequence.

At a=0 there is one parent TT_n, and the inherited one-anchor formula gives

    h_0(TT_n) = 2^n - n - binom(n-1,2) - 1,       n>=1.   (3)

The first values are 0,1,3,8,20,47. At rank1 and parent order4, the three
parents have respectively zero, one and two canonically assigned rank2
children. Thus order and rank alone do not give a common per-parent factor.

## 4. The unrestricted parent forest and its finite columns

For every nontransitive isomorphism class T, define p(T) to be the class
represented by pi of its canonical augmentation record. Each edge T->p(T)
decreases order and tau by one. Consequently alpha=|T|-tau is preserved.
After exactly tau(T) steps the chain reaches TT_(alpha(T)). This makes a
genuine parent forest on **all** tournament isomorphism classes, including
the isolated empty and one-vertex roots.

Every tournament on 2^l vertices has a transitive (l+1)-subtournament:
one vertex has an in- or outneighbourhood of size at least 2^(l-1), and
induction extends a transitive l-subtournament by that vertex. Therefore

    alpha(T)=l  implies  |T|<=2^l-1,

and the tree rooted at TT_l has height at most 2^l-1-l. Each component is
finite. Every class occurs in exactly one component, so, with b_(l,d)
counting its depth-d vertices,

    A000568(n) = sum_(l=ceil(log2(n+1)))^n b_(l,n-l),       n>=1. (4)

Equation (1) supplies the exact rank-raising transfer used along each
column. This is an unrestricted decomposition with finite columns and
rational fixed-rank rows. It does not bound the uniform number of states,
make enumeration polynomial, or yield all-order black-core connectivity.

The parent map need not preserve strongness or any marked Hamiltonian path.
Its alphabet retains the complete tournament, but a particular path must
be carried as an additional labelled permutation before discussing its
black or blue transport. A class-count forest is not a metagraph detour.

## 5. Exact clone changes and the first ownership hostile

For surviving vertices u,v set

    D_T(u,v) = {z outside {u,v}: [u->z] != [v->z]}.

After deleting x, the pair {u,v} is homogeneous exactly when
D_T(u,v) is contained in {x}. It is newly homogeneous exactly when
D_T(u,v)={x}. Existing homogeneous pairs persist. Recompute the connected
components of the resulting homogeneous-pair graph to obtain the maximal
transitive clone chains. This retains every merger; no chain is contracted
before its size and incidences are recorded. The smallest merger is C3
with one vertex deleted, leaving the two-vertex transitive clone chain.

The fixed-anchor three-case clone predicate (base/base, anchor/anchor,
anchor/base) can be imposed on G_(a+1) or its parent projection separately.
Pair-freeness of the child does not imply pair-freeness of its parent.

The regular pentagon, with bit rows

    (6,12,24,17,3),

has rank2. In the inherited integer header/letter convention its ordinary
two-anchor owner is

    (0; 3,1,0).

Deleting the last anchor projects to the one-anchor word110, which is not
canonical: the same parent has canonical word100. The canonical
augmentation instead uses

    (1; 3,2,0) -> (0; 1,0,0).

This is the smallest failure in the complete census through order7. The
first failed implication is 'canonical child implies canonical projected
parent'. The repaired construction is the restricted owner language G.
The missing information is the interaction between the parent projection
and the child anchor gauge, not a missing scalar path count.

## 6. A complete nine-node microcosm: alpha=3

The bound n<=7 makes the following entire column finite. The exact census
recovers nine classes, with order counts1,3,3,1,1. Labels below are local
receipt IDs; the receipt retains every adjacency row.

|Node|Parent|Hamilton paths|Transitive clone-chain sizes|
|---|---|---:|---|
|3.1|root|1|3|
|4.1|3.1|3|1,1,1,1|
|4.2|3.1|3|1,1,1,1|
|4.3|3.1|5|1,1,2|
|5.1|4.1|13|1,1,1,1,1|
|5.2|4.3|15|1,1,1,1,1|
|5.3|4.3|15|1,1,1,1,1|
|6.1|5.1|45|1,1,1,1,1,1|
|7.1|6.1|189|1,1,1,1,1,1,1|

The endpoint at order7 is the Paley heptagon. Its TT4-free uniqueness is
already present in [shaved tournaments](../shaved_tournaments_unavoidable_cores_20261001.md),
and is independently recovered here, not claimed as a new classification.
The new coordinate is the specified canonical parent forest. Its H=13
vertex is a concrete link to the inherited one-anchor H-spectrum boundary:
rank2 is needed there. The H-values label the nodes; no monotonicity theorem
or graph-minor map is inferred from this particular chain.

## 7. A rigorous restriction on the requested scalar recurrence

Let t_n=A000568(n). There are 2^(n(n-1)/2) labelled tournaments, with at most
n! labellings per isomorphism class. Hence

    t_n >= 2^(n(n-1)/2)/n! = exp(Omega(n^2)).              (5)

This elementary lower bound is already used in the incoming
[fixed-degree rationality review](../fixed_degree_metagraph_rationality_20261021/session22-fixed-degree-review/REVIEW.md)
and [Cartesian fibres](../session56_rmq_lca_connections_20261025/cartesian-fibres/PROOF.md).

Suppose t_n satisfied a nonzero fixed-order linear recurrence with polynomial
coefficients. Solve for its highest-shift term once that leading polynomial
has no further integer zeros. The absolute value of each coefficient ratio
is at most C n^K for some constants. The maximum absolute value among a
sliding window of terms then grows by at most C' n^K per step. Iterating
gives an upper bound exp(O(n log n)), contradicting (5).

Thus t_n is **not P-recursive**; in particular its ordinary generating
series is neither rational nor D-finite (extracting coefficients from a
polynomial-coefficient differential equation gives a P-recurrence).
This is an elementary consequence
of the standard counting bound, with no literature-priority claim. It rules
out a fixed-order polynomial-coefficient linear scalar recurrence. It does
not rule out nonlinear recurrences, recurrences with faster-growing
coefficients, growing-order recurrences, or useful growing-state engines.

This specifies what the fixed-a rational machines cannot become when a is
allowed to grow. A bounded-state, bounded-alphabet one-letter-per-vertex
owner language would already fail the weaker exponential growth bound.

## 8. Exact audit scope and remaining target

Run `python 05-knowledge/results/ports_owners_f1_20261025/anchor_check.py`.
The imported class generator is identified by normalized SHA256 in the
receipt. Everything specific to the new owner projection, feedback rank,
augmentation, clone predicate and forest is local code. The class universe
has sizes1,1,1,2,4,12,56,456 at orders0..7. Feedback ranks are checked both by
maximum transitive subsets and by hitting all directed triangles. All525
positive-rank classes receive exact parents, and every root chain has the
claimed length and constant alpha. The alpha=3 HP counts are checked by
literal permutations and an independent subset DP.

Normal and optimized Python use the same explicit-exception guards and
must produce the same receipt. These finite controls do not materialize or
validate a general determinized automaton; the finite-state result rests on
the proof in sections1–3 and its inherited owner construction.

The next useful object is a compact realization of G_2 with its weighted
parent response, rather than fitting a recurrence to the first rank2
counts. A further quotient must preserve every future attachment response,
the parent owner test and the clone disagreement sets. The unrestricted
problem permits a growing state space; the goal is to understand and
compress that growth, within the scalar-recursion boundary just proved.
