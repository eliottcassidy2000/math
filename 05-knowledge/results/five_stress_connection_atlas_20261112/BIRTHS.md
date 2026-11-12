# A finite certificate for ninety-eight infinite seed families

**PROVED for the explicitly listed finite template bank, with FINITE-EXACT
certificate inputs.** Independently inflate the original stable unselected
ports of every stable-clone birth through order8. These98 infinite families
are pairwise disjoint as marked tournament isomorphism types. Their generating
function is

    G(z)=(z³+4z⁵+8z⁶+20z⁷+46z⁸)/(1−z)
         +(5z⁷+14z⁸)/(1−z)².                         (1)

Consequently this union has exactly **19n−49 members at every n>=8**.
Its first nonzero counts, starting at n3, are

    1,1,5,13,38,103,122,141,160,179,198,... .

This is an all-order enumeration of a specified collection, not a formula
for every stable seed. Larger births and growth at newly emergent ports are
outside its definition. Chan–Pak and universal odd-coefficient coverage
remain OPEN.

## 1. Canonical compression retains the missing coordinate

Let (T,O) be a marked tournament. Join two unselected vertices when their
orientations to every other vertex agree. Each nontrivial connected component
of this relation is a transitive interval, and its adjacent vertices are
exactly its twin pairs.

Here is a direct structural proof. A vertex can have at most one twin which
it beats: if it beats both a and b, the twin relation with a forces a→b,
whereas the relation with b forces b→a. Similarly it has at most one twin
which beats it. Directed twin edges cannot form a cycle: the relations
u→v→w imply u→w. Every component is therefore a directed path; propagating
the twin relation makes its induced tournament transitive and makes every
outside vertex see the entire component uniformly. Thus the components
are intrinsic maximal unselected transitive intervals. Selected vertices
are retained as singleton components.

Contract these components, retaining a positive integer length at each
quotient vertex. The quotient has no unselected twins. Indeed two such
quotient twins would form one transitive interval, with a twin edge between
the boundary vertices of their original intervals. This contradicts maximality.
Conversely, expansion of this weighted quotient reconstructs the marked
tournament. An isomorphism of expanded tournaments must preserve these
intrinsic components, and hence is precisely a marked quotient isomorphism
preserving the lengths.

This is ordinary twin/module compression with an essential retained weight;
no priority claim is made for that general method. Forgetting the lengths
is specifically unsafe here:28 of the98 births already contain unselected
twins at ports which are not stable.

## 2. Why a finite check proves disjointness at arbitrary lengths

At a stable-clone birth, each stable port lies in a singleton component.
Otherwise an adjacent twin supplies an eligible contraction, contrary to
the definition. In the weighted quotient, an original stable port therefore
has variable length m>=1; every other component has its fixed birth length.
No new twin components can appear under expansion, by the preceding lemma.

For two templates, enumerate the finitely many marked isomorphisms between
their unweighted quotients. For each such alignment, their cones intersect
if and only if every coordinate fixed on both sides has the same length.
This criterion is sufficient as well as necessary: choose the fixed length
where one exists and choose1 otherwise. There is no length cutoff in this
test. For one template, quotient automorphisms likewise detect possible
duplicate length vectors.

The exact certificate [birth_cones.json](birth_cones.json) contains:

-98 templates, all with full joint-defect tensor identically1;
-88 distinct marked unweighted quotients, each with trivial automorphism group;
-10 remaining pairs sharing a quotient;
-an explicit incompatible fixed-length coordinate for every one of those pairs.

Therefore every length vector represents exactly one member of its template
family, and no two template families overlap. In particular, this proof
does not assume the general stable-port hypothesis or all-order confluence.

## 3. Counting and the origin of nineteen

The birth orders and original numbers of independent ports are

| Birth order | One port | Two ports |
|---|---:|---:|
|3|1|0|
|5|4|0|
|6|8|0|
|7|20|5|
|8|46|14|

A birth of order b with k ports contributes z^b/(1−z)^k, because each
positive length contributes its excess over1 to the total vertex count.
Disjointness permits summing these terms, giving(1). For n>=8, the79 one-port
families each contribute1. The five two-port births at7 contribute5(n−6),
and the fourteen at8 contribute14(n−7). Their sum is19n−49.

Thus19 is the number of two-port birth families in this finite bank.
The intercept49 follows from their birth sizes. This does not identify
either number with the nineteen-state king-strip automaton or the older49
stable-pair orbits; those objects have different definitions and maps.

The independently constructed literal graphs agree with(1) through order10.
Their complete marked-isomorphism key sets agree with the entire stable
bank through8 and with the separate descendant enumeration at9. These
finite controls supplement, rather than replace, the no-cutoff proof.

## 4. The complete nine-vertex descendant layer

First, the entire infinite operation graph of the specified bank is already
classified. Allow only inflation/contraction of an original stable port.
It is the disjoint union of **79 directed rays and19 directed integer
quadrants**: the coordinates are the excess clone lengths. Birth orders
grade these components. Distinct coordinate operations commute, and the
injectivity/disjointness proof prevents hidden identifications between them.

At each n>=9 this operation graph therefore has

    vertices                         19n−49,
    contraction edges to order n−1   38n−215,
    states with two child types      19n−166,
    states with one child type       117.

Indeed a two-port birth at order b has n−b+1 states and2(n−b) outgoing
contraction edges. Its n−b−1 interior states have two children when n>=b+2;
the two boundary states each have one. Add the79 one-port rays and the
five/fourteen quadrant origins at7/8. The exact generating identities and
sequences through64 are in the certificate. All finite-bank contraction
counts through8 and the following order9 layer agree with these formulas.
These are ordinary edges between marked types, not a count of every labelled
twin pair. Newly emergent ports could add operations at larger orders and
are not included in this classified graph.

An exact small-number intersection is worth typing correctly: at n33 the
edge formula gives1039. It counts contraction edges in this specified
family. It is neither a Hamiltonian-path count nor a seed coefficient, and
does not construct a seed with c=1039.

The [order9 receipt](descendants_nine.json.gz) enumerates every strong marked
unit seed of order9 admitting an eligible stable-twin contraction. Completeness
follows from the [contraction lemma](SEEDS.md): its contraction is a strong
order8 unit seed with a stable port, hence belongs to the complete order8 bank.
Contracting a transitive interval preserves strong connectivity here because
every path projects to a directed walk and the quotient has at least two
vertices. Reversing the move generates the original graph.

The190 order8 port presentations yield122 distinct marked descendants.
Their numbers of stable ports are

| Stable ports |2|3|4|5|7|
|---|---:|---:|---:|---:|---:|
| Descendants |46|48|23|4|1|

All complete stable sets are acyclic and jointly stable by full tensors.
Every contraction retains the other stable ports in this finite layer.
There are127 edges to order8 marked types:117 descendants have one child
type and five have two. All terminal birth types are unique. Their birth
orders have multiplicities3:1,5:4,6:8,7:35,8:74.

This does not enumerate the order9 births themselves. It isolates the exact
part of an order9 search that can be generated without a fresh unrestricted
tournament census. The five two-child cases are finite joining branches;
general confluence and the absence of newly stable ports at all larger
lengths remain OPEN.

## 5. What remains for coverage

The structural coordinates are now explicit: birth quotient, fixed component
lengths, original stable-port variables, marking, and the full endpoint
tensor. The size-counting series is rational, but the Hamilton slope U is
an exponential polynomial in the lengths. A rational count of graph classes
does not make the odd parts of those U values surjective.

The useful next changes are operations which create a new birth, change a
marking, or create and certify a new stable port. Each needs a complete mixed
defect test and a retained endpoint interface. These are precise remaining
obligations, not evidence that every odd coefficient has already been reached.

## 6. General finite-bank rationality

The compression argument extends beyond this bank. Take any finite collection
of marked tournaments, with a specified set of unselected vertices in each
allowed independent transitive inflation. The ordinary generating function
counting the distinct marked isomorphism types in their union is rational,
and its coefficients are eventually quasipolynomial. Stability is unnecessary
for this counting statement; it is an additional requirement if the counted
objects are to be unit seeds.

To prove it, group templates by their intrinsic twin quotient Q. In each
component, no inflatable vertex means a fixed positive length. One or more
inflatable vertices mean every integer length at least its original total
length is possible. Thus a template image is a Cartesian set of fixed
coordinates and lower-bounded coordinates. Include its finitely many images
under Aut(Q,O), obtaining a finite invariant union of such sets.

Burnside's lemma counts the orbits of this union. For any automorphism g,
a fixed weight vector is constant on each cycle of g. Inclusion–exclusion
over the finite union reduces the fixed-point count to intersections. On
each cycle these impose either an inconsistent condition, a fixed length f,
or a lower bound L. The respective generating factors are0, z^(d f), and
z^(d L)/(1−z^d), where d is the cycle length. Multiply over cycles, sum the
finite intersections, and average over the group. All poles are roots of
unity, proving the eventual quasipolynomial statement.

This is a direct application of classical orbit counting and rational
generating functions, not a claim of external novelty. It gives a reusable
algorithm for finite clone banks even when their cones overlap or quotient
symmetries return. In the98-template bank every quotient automorphism group
is trivial and every overlap is excluded, so it reduces to(1). It also states
exactly which coordinate was needed to make infinitely many graph sizes
manageable: the positive length of every intrinsic unselected twin block.

## 7. A sharp scope obstruction for the global programme

No finite bank of fixed tournaments, using only independent transitive clone
inflations, can cover all tournament isomorphism types at every size.
If the bank has B templates and at most K vertices per template, there are
at most B times the sum of the positive composition counts with at most K
parts, hence O(n^(K−1)) possible expanded types at order n. Identifications
only decrease that number. On the other hand,

    A000568(n) >= 2^(n(n−1)/2)/n!,

because each unlabelled type has at most n! labelled representatives.
The latter grows faster than every polynomial. This elementary comparison
proves the obstruction without any asymptotic estimate for automorphism groups.

Thus a classification of all tournaments needs an unbounded supply of new
templates or operations beyond fixed-bank transitive cloning. It says nothing
against a finite-bank classification of some much smaller seed subclass, and
it does not show that the Hamiltonian counts of a finite bank miss a specified
odd integer. Those are different universes and different observers.

## 8. A weaker condition would imply global confluence

**Subsequent resolution (12 November): the universal hypothesis below is
REFUTED.** An order-eleven strong unit seed has two eligible blocks and
two distinct terminal order-ten births. The conditional implication and
the finite older receipts remain valid. [Explicit defect and contraction
diagram](../joint_stability_forbidden_20261112/README.md#3-j2-contractible-blocks-can-fail-to-commute).

**CONDITIONAL theorem.** Suppose that, in every strong marked unit seed,
any two individually stable ports lying in different eligible twin blocks
are jointly stable. Then every sequence of eligible stable-twin contractions
ends at the same marked birth type. This is weaker than requiring the full
stable-port set to be jointly stable and transitive.

Proof. The system terminates because each contraction removes one vertex.
Consider two available contractions. If their pairs lie in the same maximal
unselected twin interval, either contraction shortens the same transitive
interval by one and yields the same marked type. This includes overlapping
pairs. If the pairs lie in different intervals, choose a stable member of
each. Their assumed joint stability and the joint contraction lemma imply
that after either contraction the other pair remains eligible. Performing
both gives the same marked graph, independently of their order. Thus every
pair of first steps joins.

Induct on the number of vertices. Every smaller descendant has a unique
terminal type by the induction hypothesis; the common descendant just
constructed identifies the terminal types of any two first steps. Hence
the original graph also has a unique terminal type. This is the elementary
termination-plus-local-confluence argument, with its actual graph guards
spelled out. ∎

This does not assume that all pairwise stable ports are simultaneously
stable. A genuinely three-variable defect could obstruct full independent
growth while leaving the pairwise condition, and therefore unique reduction
roots, intact. Conversely, a failed joint pair among two contractible blocks
is the precise missing condition for this proof; it is not by itself a
counterexample to confluence. The full stable sets through8 and the complete
order9 descendant layer satisfy the hypothesis, but its universal statement
is now REFUTED by the linked continuation. This separates the birth-classification obligation from the
stronger unrestricted composition obligation.
