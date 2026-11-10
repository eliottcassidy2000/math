# Exact cyclic orbit costs and the sharp sevenfold floor

Status: PROVED from inherited OCF, with finite input independently reproduced.
No external novelty claim. Chan–Pak completeness remains OPEN.

## Inheritance and scope

The mechanisms are THM-002 (`01-canon/theorems/THM-002-ocf.md`), the
orbit-product inequality in THM-4812
(`01-canon/theorems/THM-4812-presentation-fibres-7-and-21-are-absent-at-every-order-orbit-product-lemma-vertex-transitive-lemma-closed-window.md`),
and THM-4987 (`01-canon/theorems/THM-4987-subgroup-path-orbit-gaps-and-small-symmetry-cost.md`).
The eight seven-circulant census is in Session190,
`05-knowledge/results/session190_constructive_links_20261104/groups/session190-groups/PROOF.md` section3,
with audit `groups/session190-groups-audit/AUDIT.md` in that session packet.
The census and the hostile Q7/P7 pair are inherited,
not discoveries of this session. THM-4987 already records
(H,Aut order)=(175,7) and (189,21). The all-order statement below assembles
these inputs into the exact cost for a specified cyclic action, including
orbit multiplicity and the equality boundary. The optional low-q corollary
uses the inherited, computer-assisted THM-4812(c); the main theorem does not.

Concept board: OCF positive weights; cyclic orbit sizes; Hamilton-path orbit
counts; equality as absence of mixed cycles; full versus specified symmetry.
The hostile is P7: dividing by its full group21 gives9, but dividing by its
cyclic subgroup7 gives27. Forgetting the denominator destroys the conclusion.
The least-used sidecar is the whole orbit-size multiset, not just group order.

## Definitions

A tournament has one directed edge on each unordered pair of distinct vertices.
H(T) counts directed Hamilton paths as ordered vertex sequences. Let
A=C_m act faithfully by automorphisms. All tournament automorphism orders are
odd: an involution would swap two vertices and reverse their edge, and Cauchy's
theorem excludes even group order. Each orbit has odd size d dividing m, and
its induced tournament is a Cayley tournament on Z/d: identify an orbit with
the quotient of A by a point stabilizer. A skew connection set S picks exactly
one element from each pair {j,-j}; x points to y iff y-x belongs to S.

Set mu_1=1, and for odd d>1 define mu_d as the minimum H among the
2^((d-1)/2) skew connection sets on Z/d. This definition is exact, not an
asymptotic estimate. A nontrivial vertex-transitive tournament is strong:
its strongly connected components have a canonical linear order, so a
transitive automorphism action cannot have more than one such component.

## Theorem 1: arbitrary cyclic actions, sharp for each orbit profile

Let the orbit sizes of A=C_m on T be d_1,...,d_t (including fixed vertices).
Then

    H(T) >= product_i H(T[O_i]) >= product_i mu_(d_i),
    q_A(T)=H(T)/m >= (product_i mu_(d_i))/m.

Moreover, equality in the final bound holds if and only if the orbit blocks
are an ordinal sum and every block attains its cyclic minimum mu_(d_i).
An ordinal sum means all edges between earlier and later blocks point forward.

For every multiset of odd sizes with lcm(d_1,...,d_t)=m, the bound is attained:
choose minimizing Cayley tournaments and ordinally concatenate them. The
simultaneous translation generator has order m and the specified orbits.
Thus this is an exact best bound for the prescribed profile, not merely a
universal lower bound whose equality may be impossible.

### Proof

OCF expresses H(T) as the sum of 2^k over collections of k vertex-disjoint
directed odd cycles. Collections wholly internal to the blocks have total
weight product_i H(T[O_i]). All other terms are nonnegative, proving the
first bound and showing equality iff no directed odd cycle meets two blocks.
The second bound follows from the definition of mu_d.

For the equality characterization, no mixed odd cycle implies no mixed
triangle. If a vertex v outside a strong block B has both an out-neighbor
and an in-neighbor in B, take a directed path in B from an out-neighbor to an
in-neighbor. At the first change of relation to v, its consecutive vertices
with v form a mixed directed triangle. Hence v sees B uniformly. Applying
this to both strong blocks shows every pair of blocks has uniform direction;
singletons cause no difficulty. The quotient has no directed triangle, so it
is transitive. Conversely an ordinal sum has no cycle crossing blocks, and
its H is the product of the block counts. Finally, the order of a permutation
is the lcm of its orbit lengths, and simultaneous translations supply the
claimed faithful cyclic action. QED.

## Theorem 2: prime actions and the sevenfold price

If C_p acts nontrivially with s nontrivial orbits, all these orbits have size p,
so

    H(T)/p >= mu_p^s/p.

Equality is exactly an ordinal sum of s minimizing p-circulants and any number
of singleton fixed points, with the prescribed action transitive on each
p-block. Exact finite input gives mu_3=3, mu_5=15, and mu_7=175.
For p=7, all minimizers are isomorphic to Q7=Cay(Z/7,{1,2,3}). Therefore

    H(T)/7 >= 175^s/7 >=25,

and H(T)/7=25 iff T is an ordinal sum of one Q7 and fixed singletons.
For fixed s, equality in the first inequality is precisely an ordinal sum of
s copies of Q7 and fixed singletons. The diagonal C7 action on these blocks
is intended: the full automorphism group may be much larger when s>1.

The finite universe has eight skew connection sets: six affine images of
{1,2,3} have H=175; {1,2,4} and its negative have H=189. The attached checker
uses subset DP and independently enumerates all7! ordered sequences for each
of these eight tournaments. Analogous independent checks cover p=3,5.
An optional reproduction for all32 eleven-circulants gives mu_11=92411,
already known in the repository (see THM-132-z11-ocf-alpha1-not-monotone.md).
No claim about general prime minimizers is made.

## Corollary: small path-orbit count for an arbitrary cyclic subgroup

Assume the inherited THM-4812(c) vertex-transitive classification. If A is
cyclic and q_A<23, then |A| divides15. Thus only 1,C3,C5,C15 occur as possible
cyclic subgroup types in this range; each occurs, using ordinal combinations
of K1,C3,R5. This is a classification of possible groups, not of tournaments
or attainable q values.

Proof. Let O_i be the full-G orbits, where G=Aut(T). Put G_i=Aut(T[O_i]).
The faithful restriction embeds A into product_i A_i, where A_i is its cyclic
projection into G_i. OCF therefore gives

    q_A >= product_i H(T[O_i])/|A_i| >= product_i H(T[O_i])/|G_i|.

Every individual factor on the far right is an integer at least1, by freeness
on ordered Hamilton paths. THM-4812(c) forces every full-G orbit to induce
K1,C3,R5, or P7. Their groups have orders1,3,5,21 respectively. Aut(P7) is
C7 semidirect C3 with nontrivial action, so its cyclic subgroups have order
at most7. A P7 block would contribute at least189/7=27 to the middle product,
while every other factor is at least1, contradicting q_A<23. Thus only
K1,C3,R5 remain, and A embeds in a product of C3 and C5 groups, of exponent15.
A cyclic subgroup consequently has order dividing15. QED.

The group structure of P7 can also be read from its explicit affine
transformations x -> a*x+b, a in {1,2,4}: this gives21 automorphisms, the
inherited full order; all elements outside translations have order3.

## Verification and boundary

Run from repository root:

    python3 05-knowledge/results/session238_planar_chanpak_20261110/cyclic_cost/check.py

`results.json` records every skew set for p=3,5,7,11 and its exact count.
Positive controls include singleton padding and two ordinal Q7 blocks;
the hostile cross-block triangle strictly raises H over the internal product.
The p=11 computation uses only DP; no second independent11! enumeration is
claimed. The theorem is all-order because OCF and the orbit action are
all-order; finite data are used only to evaluate specified constants.

Map: a specified cyclic action -> induced orbit Cayley tournaments -> positive
OCF subfamilies. Preserved target: lower bound and exact equality for H/|A|.
Lost data: mixed cycle contributions, recovered by the equality sidecar.
The C7 bound alone does NOT say q_A>=25 whenever7 divides a larger cyclic
|A|; that would change denominators. The low-q corollary uses a separate proof.
Neither result proves Chan–Pak completeness or a universal constructive bank.

Verification gates use explicit exceptions and survive Python optimization.
Passing `--inject-failure` deliberately corrupts the expected sevenfold
minimum; both ordinary Python and `python3 -O` must exit nonzero.
