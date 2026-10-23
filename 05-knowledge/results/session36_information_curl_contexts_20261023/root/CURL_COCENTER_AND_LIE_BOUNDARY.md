# Curl is a trace quotient, not a compositional algebra quotient

Status: PROVED elementary algebraic interpretation and explicit obstruction. Standard names HH_0, free Lie algebra and PBW are identified as such, not new theories. The new bridge is the actual tournament-curl classification.

## 1. The ordered-join algebra is free

Let V be the rational vector space with basis the isomorphism types of nonempty strong tournaments, grading each letter by its vertex count. The rational algebra with tournament types as basis and ordered join as multiplication is exactly the tensor/free associative algebra T(V): unique SCC decomposition is its word representation. The empty tournament is its unit. This is an algebra of formal linear combinations, not a new multiplication on individual vertices.

The older `05-knowledge/results/tournament_even_bijection_20261011c.md` already used the weighted SCC alphabet and its Lyndon factorization to build an order-dependent tournament/even-graph type codec. We recover that mechanism, rather than presenting the alphabet or Lyndon count as new.

## 2. The curl quotient is its cocenter

For any associative algebra A, let [A,A]_lin be the linear span of commutators uv−vu, not the two-sided ideal they generate. Its quotient is the cocenter, also called zeroth Hochschild homology HH_0(A). Standard definition: MIT18.706, Lecture13, https://ocw.mit.edu/courses/18-706-noncommutative-algebra-spring-2023/mit18_706_s23_lec13.pdf .

For a word algebra, uv−vu identifies a word with a cyclic rotation. Conversely every cyclic rotation has that form, and bilinearity expresses any commutator as a sum of such differences. Hence A/[A,A]_lin has a basis of cyclic words. Together with the Session36 lattice theorem, its positive-degree basis is EXACTLY the isometry classes of tournament augmentation lattices equipped with their metric and alternating form.

This is more than equality of dimensions: the basis map sends a tournament's SCC word to its actual metric-curl class, and its kernel is precisely the commutator subspace. The word's initial gap is the data discarded by this map.

## 3. A sharp compositional failure at eight vertices

Let a be the singleton strong type and b the directed triangle. The tournaments ab and ba have the same metric-curl class: their strong-component words are rotations. Append the same tournament ab to both. The results are

abab,    baab.

These are not cyclic rotations: the former alternates its two letter types, while the latter has adjacent equal letters in its cyclic word. Their metric-curl classes differ, although both have H=9 and the same complete path-cover polynomial.

Thus the curl quotient is NOT stable under ordered join. One cannot multiply its elements by concatenating arbitrary representatives. This is precisely why [A,A]_lin is not an ideal. The tensor algebra's commutative algebra quotient would identify ALL permutations of the letters, a strictly coarser object than cyclic rotation.

The example has eight vertices and is minimal for this join-congruence failure. Below eight vertices every possible SCC-type multiset has only one cyclic ordering: if there is at most one nontrivial strong component this is immediate; two nontrivial components must both have size3, possibly accompanied by one singleton, and the size3 strong type is unique; any third nontrivial component would require at least9 vertices. The only other two-component possibility of sizes3 and4 has just two letters, whose two orders are rotations. Therefore no equal-curl representatives can yield different cyclic order after a common join at smaller total order.

This supplies a concrete restoration sidecar: an oriented SCC necklace plus a chosen initial gap recovers its linear SCC word. A marked vertex alone need not select this gap, since a gap can occur on either side of its component; the companion lattice note gives a marked vertex WITH its score as a sufficient numerical selector. Repeated isomorphic blocks can identify equivalent gap choices at type level.

## 4. The Lie-algebra connection, and what it does NOT preserve

The free Lie algebra on V embeds in T(V) with commutator bracket. Its universal enveloping algebra is T(V). PBW yields the Hilbert identity

1/(1−S(z)) = ∏_(n≥1)(1−z^n)^(−ell_n),

where ell_n is its degree-n dimension. Standard primary references: Reutenauer's PBW chapter, https://reutenauer.math.uqam.ca/wp-content/uploads/2024/04/Reutenauer1986_Chapter_TheoremOfPoincare-Birkhoff-Wit.pdf , and Melançon–Reutenauer, 'Lyndon Words, Free Algebras and Shuffles', Canadian Journal of Mathematics41(1989),577–591, https://doi.org/10.4153/CJM-1989-025-2 . Equivalently a Lyndon basis has one vector for each primitive cyclic word after choosing an alphabet order.

The inherited even-graph Euler identity therefore identifies ell_n with the number of connected automorphism-even graph types. The new curl interpretation identifies the corresponding primitive cyclic-word COUNT with primitive curl types. This is a valid representation-theoretic interpretation of these integers, not a resemblance between dimensions.

But the obvious map from the free Lie algebra to the cocenter kills EVERY Lie bracket. In letter length at least2 every Lie monomial is a bracket, so its image is zero. The surviving one-letter generators are independent in the cocenter. Thus the image of FreeLie(V) in HH_0(T(V)) is exactly V, not the full primitive-necklace span.

Consequently matching the Lyndon-basis index set to primitive curl classes is a vector-space basis correspondence, NOT the natural projection, not a Lie embedding into curl space, and not a recovered tournament Lie bracket. Dimension matching alone misses precisely this distinction.

## 5. Hamiltonian observables are characters, but not faithful traces

H extends linearly to an algebra character T(V)→Q because it multiplies under ordered join. The full path-colour polynomial similarly gives a character into a commutative polynomial ring (the inherited THM-3166 normalization). Both kill commutators and hence factor linearly through the cocenter. They actually forget more than the cocenter: all reorderings of the same component multiset have equal character values, and strong tournaments can also have equal path-cover profiles.

The hierarchy is therefore operational, not a ladder of increasingly grand names:

- the SCC word supports ordered composition;
- the cyclic word supports trace-like curl observations but needs a starting gap for composition;
- a commutative path-cover character supports the particular Hamiltonian module-response task but loses further structure;
- the free Lie primitives describe formal noncommutative interaction, yet those brackets vanish under every such scalar character.

This is an explicit instance of the user's information-flow viewpoint. It proves which 'order' is lost and which operation notices the loss. It makes no claim about Yang–Mills, vertex-algebra realization, or a natural graph bijection.
