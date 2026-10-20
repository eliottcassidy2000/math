# Uniform obstruction to equivariant type bijections

Status: PROVED by the argument below, independently read and accepted by the session parent; additional audit remains a publication gate. This resolves both assertions of HYP-9172(a), for every n >= 5. It makes no assertion about an arbitrary bijection of unlabelled types or about HYP-9172(b)'s requested combinatorial pairing.

## Main result

For every n >= 5, neither of the following exists:

1. An S_n-equivariant map from labelled tournaments to labelled Royle-sign-even graphs that induces a bijection on isomorphism types.
2. An S_n-equivariant map from labelled tournament switching classes to labelled graphs that are both Royle-sign-even and Eulerian, inducing a bijection on isomorphism types.

Here a graph is sign-even if every automorphism reverses an even number of its edges relative to a fixed reference orientation. Eulerian separately means every degree is even. A type-bijective equivariant map is allowed to collapse labelled objects; the result is therefore stronger than nonexistence of an equivariant labelled bijection. The bound is sharp: existence for n <= 4 is inherited from THM-4531, with independent finite controls supplied here through n = 6.

## Inheritance and the missing connection

The closest mechanism is THM-4531, `01-canon/theorems/THM-4531-twisted-mallows-sloane-principle-and-no-natural-bijection-for-switching-classes.md`, Theorem G and Lemma N1. Its old block certificates covered bounded orders and depended for further progress on prime-additive representations. The hostile at n = 5 was competition between order-3 and order-5 symmetries for the empty graph.

The previously underused sidecar is the product of stabilizer orders in the Walsh orthogonality identity (THM-4721, `01-canon/theorems/THM-4721-twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md`, Theorem W). Its quotient version is reproduced below. The recent `05-knowledge/results/prime_local_walsh_bijection_20261020/PRIME_VALUATION_THEOREM.md` and `QUOTIENT_COROLLARY.md` prove stronger one-prime statements, but the present proof needs only the elementary product identity, not their matching constructions.

The key connection is:

- Source: a hypothetical equivariant map inducing a type bijection.
- Target: paired stabilizer orders.
- Preserved predicate: divisibility, forced by full stabilizer containment.
- Sidecar: equality of total odd-prime valuations across the two type sets.
- Consequence: every paired odd-prime valuation must be equal, simultaneously.
- Hostile: the empty graph would require a source stabilizer carrying full ambient Sylow-3 and Sylow-5 subgroups. A local 3-cycle obstruction makes this impossible.

## Lemma 1: monotonicity plus an odd-product identity

Let A and B be finite G-sets with equally many orbits. Suppose, for every odd prime p,

    sum_[a] v_p(|G_a|) = sum_[b] v_p(|G_b|).

If f:A -> B is G-equivariant and induces a bijection of orbits, then

    v_p(|G_a|) = v_p(|G_f(a)|)

for every a and every odd p.

Proof. Equivariance gives G_a <= G_f(a), so every difference on the right minus left is nonnegative. Summing over the orbit bijection gives zero. Each difference is therefore zero. This conclusion applies to all odd primes for the same hypothetical map; it does not assert that independently chosen one-prime bijections can be combined. QED.

## Lemma 2: the required product identity for every invariant reversal quotient

Let V=F_2^{E(K_n)}, let U<=V be S_n-invariant, and let A be the set of tournament classes modulo reversal along U. Put c=|A|=2^(binomial(n,2)-dim U). Let B be the sign-even graphs X in U-perp. Then A and B have equally many S_n-orbits and

    product_[a] |Stab(a)| / product_[X] |Aut(X)|

is an integral power of 2, possibly with a negative exponent. In particular, the hypothesis of Lemma 1 holds.

Proof. Work over R=Z[1/2]. Choose a reference tournament. For X in U-perp, the Walsh character phi_X on tournament reversal differences descends to A. These characters form an R-basis because the inverse Fourier matrix is c^-1 times its transpose. The S_n-action permutes them with the graph edge-reversal sign cocycle.

For an odd graph, invariance forces its coefficient to equal its negative, hence to vanish in R. For a sign-even graph, transport of one coefficient around its orbit is consistent, yielding one signed orbit sum Phi_X. Thus these Phi_X, one for each sign-even graph type in U-perp, form an R-basis of the invariant functions on A. Orbit indicators of A form another R-basis. Their integral evaluation matrix M, with M[a,X]=Phi_X(a), is consequently square and invertible over R.

Write D_A=diag(|Stab(a)|) and D_B=diag(|Aut(X)|), one entry per type. Character orthogonality grouped by S_n-orbits gives

    M^t D_A^-1 M = c D_B^-1,

and hence

    M D_B M^t = c D_A.

Taking determinants yields

    det(M)^2 product_[X] |Aut(X)| = c^N product_[a] |Stab(a)|.

Both c and det(M) are units of R, so their odd-prime valuations vanish. This proves the asserted equality of all odd-prime valuation sums. QED.

For U=0, B is all sign-even graphs. For U the cut space, U-perp is the Eulerian graph space, so B consists precisely of sign-even Eulerian graphs.

## Lemma 3: local oddness on the source

A tournament has no nonidentity involutory automorphism: if an involution interchanges u and v, it reverses the unique arc between them. Thus its automorphism group has odd order.

A tournament switching class has no nonidentity involutory automorphism fixing a vertex v. Indeed, every switching class contains a unique representative in which v is a source. To construct it, switch precisely those other vertices whose arcs point into v. Uniqueness follows because the cut between two such representatives must have no edge incident with v, hence is empty. Relabelling by an automorphism of the class fixing v preserves this unique representative. Such an involution would therefore be a tournament automorphism, impossible by the preceding paragraph.

In particular, for n >= 5, neither source stabilizer can contain a double transposition: a double transposition fixes at least one vertex.

## Lemma 4: 3-cycle supports and the 5-adic gap

Let n >= 5 and H<=S_n contain no double transposition. If

    v_3(|H|)=v_3(n!),

then, with k=floor(n/3),

    v_5(|H|) <= v_5(k!) < v_5(n!).

Proof. First, supports of distinct single 3-cycles in H are equal or disjoint. Suppose two supports overlap in two points. After replacing one cycle by its inverse and relabelling their four points, the two cycles can be written (1 2 3) and (1 2 4); their product is (1 3)(2 4), a forbidden double transposition. If two supports overlap in one point, conjugate the first cycle by the second: its support and the original support overlap in two points, reducing to the preceding case. This also proves the assertion at n=5, without any appeal to group classification.

The 3-adic equality says a Sylow-3 subgroup of H is a Sylow-3 subgroup of S_n. Every Sylow-3 subgroup of S_n contains k disjoint single 3-cycles: the subgroup generated by any k such cycles extends to a Sylow-3 subgroup, and all ambient Sylow subgroups are conjugate. Their k disjoint supports leave r=n-3k in {0,1,2} vertices. Since all single-3-cycle supports in H are pairwise disjoint, these k supports are the entire support collection.

Conjugation by H permutes that collection. Therefore H preserves the partition into k triples and the residual r-set, giving

    H <= (S_3 wr S_k) x S_r.

Its containing group's order is 6^k k! r!, whose 5-adic valuation is v_5(k!). Finally v_5(k!)<v_5(n!) for n>=5: for n=5,6,7 the factor 5 lies between k+1 and n; for n>=8 the interval from k+1 through n has at least five consecutive integers and therefore contains a multiple of 5. QED.

Only the primes 3 and 5 are used. There is no prime-additive conjecture, Feit-Thompson theorem, Jordan theorem, or Hall-subgroup classification in this proof.

## Proof of the main result

Apply Lemma 2 with U=0 or with U the cut space. Suppose a map as in the theorem exists. The empty graph belongs to the target and has stabilizer S_n. Its unique preimage type has a representative a with H=Stab(a). By Lemma 1,

    v_3(|H|)=v_3(n!),    v_5(|H|)=v_5(n!).

Lemma 3 excludes double transpositions from H, and Lemma 4 makes these two equalities incompatible. This contradiction proves both assertions for every n>=5. QED.

## Corollaries: a joint-prime obstruction and an explicit Hall family

For either pair of type sets in the main theorem, and every n>=5, no arbitrary type bijection can preserve BOTH of these predicates:

- its stabilizer contains a Sylow-3 subgroup of ambient S_n;
- its stabilizer contains a Sylow-5 subgroup of ambient S_n.

The source has no type satisfying both, by Lemmas 3 and 4, while the target has at least the empty graph satisfying both. This statement concerns **full ambient Sylow subgroups**. It does not claim that a large-order tournament cannot have some C_3 and some C_5 symmetries simultaneously.

There is also a uniform, explicitly described Hall-deficient family. Let A_p and B_p be the source and target type subsets carrying a full ambient Sylow-p subgroup. Then |A_p|=|B_p| for each odd p. This is the inherited one-prime refinement, but follows directly from Lemma 2's matrix as follows. Its inverse satisfies

    (M^-1)[X,a] = |Aut(X)| M[a,X] / (c |Stab(a)|).

Since all inverse entries lie in Z[1/2], a p-unit entry M[a,X] implies v_p(|Stab(a)|)<=v_p(|Aut(X)|). The determinant is a p-unit, so its expansion contains a permutation all of whose entries are p-units. Along that permutation, Lemma 2's equality of total valuations forces equality at every pair. In particular the numbers at the ambient maximal valuation v_p(n!) agree.

Now set S=A_3 union A_5. Its two parts are disjoint by Lemma 4. Any equivariant image of S lies inside B_3 union B_5, whereas

    |S| - |B_3 union B_5| = |B_3 intersection B_5| >= 1.

Thus S is a Hall violator for THM-4531's stabilizer-containment graph. In fact the matching deficit is at least the number of target types simultaneously carrying full ambient Sylow-3 and Sylow-5 subgroups. The finite controls check this bound at orders five and six, in both settings.

## Exact scope and failure boundary

- PROVED: HYP-9172(a), both directions from the stated source to its twisted Mallows-Sloane partner, at every n>=5.
- The word “directions” here refers to the tournament and switching-class source cases. Reverse maps from sign-even graphs to tournaments were already impossible by the fixed empty graph; reverse switching maps were already excluded in THM-4531.
- A map merely between unlabelled type sets has no stabilizer-containment requirement and is not obstructed; such bijections exist by the count identity.
- The two separate one-prime bijections remain valid; the impossible step is simultaneous monotone preservation forced by a labelled equivariant type bijection.
- At n<=4, the factor 5 is absent and the contradiction disappears. At n=4 a double transposition has no fixed point and may stabilize a switching class, so the source lemma's boundary is also explicit.
- The proof does not resolve the combinatorial-pairing demand in HYP-9172(b).

## Reproduction

`python3 -B exact_controls.py --max-n 6` independently enumerates all labelled tournaments, switching classes, and graphs at every order through six, computes sign-evenness from actual stabilizers, compares total odd-prime valuations, checks the stabilizer-containment matching, and probes every overlapping pair of single 3-cycles through seven vertices. It also checks the 5-adic gap through 10,000 as a finite control of the all-order elementary inequality.


The complete run reports tournament type counts 1,1,1,2,4,12,56 and switching-class type counts 1,1,1,1,2,2,6 at orders 0 through 6. Matching deficits are zero through order four, then 1,2 for tournaments and 1,1 for switching classes. Every stated finite check passes. These computations are controls, not a replacement for the all-order proof.
