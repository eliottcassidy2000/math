# Central switches need a declared consumer

Status: PROVED elementary contextual criterion, with literal finite inputs;
not an all-target construction. Generic cover substitution and central fusion
transport are inherited mechanisms. The concrete contrast below is the result.

## 1. Literal sources and three different information levels

Let C_(2s+1) have vertices modulo2s+1 and arcs i→i+j for1≤j≤s.
Add root r which loses exactly to A and beats every other old vertex.
For C5 the two incoming sets A={0,1} and A={0,2} give X,Y with
H(X)=41,H(Y)=43. Both have the same rooted fusion interface
(h,Q,U,V,d,e)=(15,120,78,57,6,9).
Thus their fusion difference is2E14, not a change in an endpoint channel.
The digits lane records the source construction and inheritance of the
central-switch mechanism from Session229 and the five-module bridge.

Write p_j(T) for the number of unordered covers of all vertices by j
nonempty directed paths. Direct exact subset recurrence gives

    p(X)=(0,41,167,186,81,15,1),
    p(Y)=(0,43,169,186,81,15,1).

In particular their full cover-polynomial difference is2z+2z².
Rooted endpoint equality does not mean equality of this full cover profile.

## 2. Exact contextual amplification theorem

For any nonempty tournament B with marked vertex v, let B[v←G] denote
ordinary substitution: replace v by G with every old incoming/outgoing arc
copied uniformly to every vertex of G. Then

    H(B[v←Y])−H(B[v←X])=2H(B)+4a₂(B,v),                 (1)

where a₂(B,v) is the following explicit nonnegative integer. Delete v,
arrange all remaining labelled vertices in a linear order, and place two
identical nonempty G-run placeholders, with at least one old vertex between
them. Count those arrangements for which all old adjacencies are arcs of B
and every old/placeholder adjacency obeys the original v incidence.

Consequently the difference is at least2H(B). It is exactly2 IF AND ONLY IF
B is transitive. In particular a strong outer context with at least3
vertices amplifies this switch by at least6. This contrasts with fusion
contexts of root-deletion count1, which can be strong and preserve the+2.

Proof. A Hamilton path through B[v←G] breaks G into j maximal runs.
The runs are an unordered j-path cover of G, with j! possible orders.
The outside arrangement defines a_j(B,v) as above with j placeholders.
This is a bijection, so

    H(B[v←G])=Σ_(j≥1) a_j(B,v) j! p_j(G).

For j=1 contraction is a Hamilton path of B, so a₁=H(B).
All a_j≥0. Applying the displayed exact cover difference proves(1).
If the difference is2, H(B)=1, hence B is transitive by the inherited
unique-Hamilton-path characterization. Conversely in a transitive B no
path can leave the substituted block and later return: all vertices before
v dominate G, and G dominates all vertices after v. Thus a₂=0 and H(B)=1.
This also covers the singleton outer tournament. QED.

The general run-cover expansion is standard/inherited in the repository's
cover-profile contextual-equivalence work. No novelty is claimed for that
expansion or for the H=1 characterization. The sharp distinction between the
consumers of these concrete central switches is what is established here.

## 3. Root-central order need not survive another consumer

For C7, let X7 use incoming{0,2}, and Y7 use incoming{0,1}.
Their counts are559 and561 and their rooted endpoint interfaces agree.
Their full cover differences, however, are

    p(Y7)−p(X7)=2z−10z²−12z³−2z⁴.

Take outer B=C3 and replace one vertex by G. Every path has either one G run
or two. In the latter case the other two vertices form the unique directed
bridge out of and back into G. Therefore

    H(C3[G,1,1])=3p₁(G)+2p₂(G).

The Y7−X7 difference is6−20=−14, reversing the scalar ordering.
This is a hostile example against using an exact fusion-central switch as an
operation-independent count derivative. It does not invalidate fusion
centrality, which belongs to its specified consumer algebra.

## 4. An infinite consumer family, not just one hostile triangle

Let B_m=C3[{v},TT_m,{w}], m≥1, with v→TT_m→w→v,
and substitute G at v. A path uses at most two G-runs, since w is the
only return port from the outside. The one-run coefficient is
H(B_m)=2^m+1. For two runs the path must be

    G-run1, TT_S, w, G-run2, TT_(complement S),

where S is any nonempty subset of the m ordered transitive vertices.
Each subset is internally ordered uniquely. Hence a₂=2^m−1 and

    H(C3[G,TT_m,1])=(2^m+1)p₁(G)+2(2^m−1)p₂(G).        (2)

The empty S would make the first G-run directly precede w, which is
forbidden; this is why the coefficient is2^m−1, not2^m. At m=0 the outer graph is the ordinal join w→G and the count is H(G);
formula(2) does NOT apply, because its displayed one-run coefficient would
incorrectly be2. This zero-length boundary is separate.

Thus the C5 switch has exact difference6·2^m−2, whereas the C7 switch
has exact difference22−18·2^m. The latter reverses sign for EVERY m≥1,
while both have the original+2 difference at m=0. This supplies an all-order
classification of the consumer's response, not an extrapolation from the
triangle experiment. Its two response modes are the constant and2^m; the
root-central step can become an exponentially large reversed step as the
consumer grows. No inversion of that consumer is physically asserted.

## 5. Consequences for constructing a corridor

A proposed return map must carry its consumer type as well as its numeric
step. Ordinary substitution cannot be the silent replacement for a unit-
deletion fusion context: the C5 step is taxed by H(B). For C7 the sign can
even reverse. Conversely C5's nonnegative whole-cover difference gives a
context-uniform ordering theorem and a lower bound, stronger than its raw
count comparison but weaker than a gapless corridor.

A productive next question is to design multiple marks so that one consumer
acts centrally while another resets its count baseline without amplifying
the step. That reset must be a literal graph operation with a proved cover
response; no quotient by H(B) is allowed. This note supplies a cheap decisive
consumer test for such proposed resets, not their existence.

## Inspiration and scope

The seeded older Smith-selection-rule reflection distinguished which
transformations can legally act on the full carrier from transformations of a
shadow. Here no topological theorem is imported: the exact source/target map
is the run decomposition, the lost datum is p₂ and higher covers, and the
small triangle is the decisive probe. Old score/spectral regression work
reinforces the separate point that excellent scalar prediction is not a
contextual equality. These analogies motivated a concrete test, not a claim
of equivalence between the subjects.

Run `python3 -B 05-knowledge/results/session276_corridors_20261112/selector/check.py`.
The all-context statement rests on the bijection, not on the small replay.
