# Equal path-color data with prescribed recursive symmetry

Use Session 310's directed tournaments A (seven vertices) and B (eight vertices), with equal path-color polynomials Q_B(t)=t Q_A(t). Write D=A→TT1 and use a star for taking the converse. Reversing every path in every slot shows that taking the converse preserves Q. Thus B,B*,D,D* have identical Q, and all four are pairwise nonisomorphic. Indeed B is strong and not self-converse (verified from the given outgoing rows), while D has its two strong components in sizes 7,1, and D* has them in sizes 1,7. The two strong seeds A and B are rigid; this is also checked from their rows. Consequently the four eight-vertex leaves are rigid.

## 1. Grow a supply of indistinguishable rigid colors

For a cyclic quotient C₃, the three substituted fibres are its maximal proper modules, regardless of whether the inner tournaments are strong. To see this, a module meeting two fibres must contain all vertices of the third: any outside vertex there has opposite relations to those two fibres. It then meets all three, and any outside vertex in a partly included fibre again has opposite relations to the other two. Thus any module meeting two fibres is the whole graph.

Choose any ordered triple of the four leaf types which is not constant. Inflate a directed triangle by that triple. The result is strong. Any automorphism induces a cyclic permutation of the three fibres and must preserve their isomorphism types, so the nonconstant triple forces this action to be trivial. The leaves are rigid, so the result is rigid. Two such triples give the same class exactly when they differ by cyclic rotation. All these tournaments have the same full Q because uniform substitution with equal slot profiles preserves the polynomial (see Session 310): for each global slot retain the itinerary of runs through the replaced module; across all slots its fillings are precisely an ordered path cover of the module by the required total number of runs.

If q_h is the number of available rigid colors at depth h, nonconstant cyclic triples give exactly

    q_(h+1)=(q_h³−q_h)/3,   q_0=4.

This yields 4,20,2660,… colors at respective orders 8,24,72,…; the positive-depth colors are all strong. Arbitrarily many strong rigid colors therefore share the same full Q at the same order. The recurrence counts a chosen recursively built family, not all members of its Q-fibre.

## 2. Write the symmetry of an arbitrary quotient into the construction

Let R be any nonempty tournament, and choose a depth h≥1 with at least |R| available colors. Replace each vertex of R by one of the equal-order strong rigid colors, with inter-fibre arcs following R. All choices of colors have the same full path-color polynomial by uniform substitution.

Every automorphism permutes the declared fibres. For a proof, let C be any module of the size of one fibre. If C contains a vertex outside a fibre U and meets U partially, choose an outside-of-C vertex x in U. Uniform incidence to the vertex of C outside U forces x to see every vertex of C∩U in the same direction, and the direction is the same for all x in U\C. This gives a uniform directed cut inside U, impossible by strongness. Hence a module of size |U| is a whole fibre: it either lies inside one, where size decides it, or is a union of whole fibres. Apply this to the image of a declared fibre under an automorphism.

The induced group is exactly the stabilizer of the coloring in Aut(R). A color-preserving quotient automorphism lifts by the unique isomorphism between rigid copies of each color, and an automorphism in the kernel fixes each rigid fibre pointwise. Thus

    Aut(inflated R) ≅ {g∈Aut(R): color(gv)=color(v) for every v}.

A monochromatic coloring has the full group of R; an injective coloring is rigid. Intermediate colorings realize every stabilizer of a partition into individually named color classes. This does not claim every abstract subgroup of Aut(R) is a color stabilizer. When R is strong the inflated tournaments are strong, so the full path-color polynomial fails to distinguish trivial from nontrivial global symmetry even among strong tournaments. The common Hamiltonian count already satisfies any divisibility forced by the group of the monochromatic realization; its rigid partner need not display that group.

For an explicit unbounded separation, take R to be the height-d iterated cyclic ternary tournament. Its automorphism order is 3^((3^d−1)/2): each node of the ternary tree can independently rotate its children, and the maximal proper modules above show these are all the automorphisms. Choosing enough colors and applying the construction yields same-order, equal-Q strong pairs with one automorphism group of that order and the other trivial. This is group erasure at unchanged *entire polynomial*, not just a comparison of Hamiltonian totals.

## Interpretation

The recursion uses two genuinely different observables: Q certifies that every uniform exterior count is preserved; the isomorphism type of a leaf controls the allowed permutation of the surrounding tree. This is the reverse of using equal scalar counts to infer local symmetry. It does not produce new prescribed counts, or a guarantee that selective endpoint ports survive a color replacement.

### Every finite odd-order group occurs alongside a rigid strong twin

By [Moon's realization theorem](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/CF0B648313DD8BE07D726932193D9044/S0008414X00035598a.pdf/tournaments-with-a-given-automorphism-group.pdf), every finite group G of odd order is Aut(S) for a tournament S. If S is not strong, embed it as the size-|S| fibre of C₃[S,TT1,TT2]. This is strong, its three fibres are canonical by the argument above, and the two rigid transitive fibres of distinct sizes prevent any cyclic permutation of the fibres. Thus the embedded tournament has group G as well. Apply the monochromatic and injective constructions to that strong quotient. It follows that **for every finite odd-order group G there exist two same-order strong tournaments with the same full Q, one with automorphism group G and one rigid**. This combines the classical representation theorem with the explicit equal-Q color construction; it is not a new proof of Moon's theorem.

## 3. Let reversal symmetry grow without ordinary automorphisms

Use a different equal-Q leaf alphabet of size nine per leaf:

    X=TT1→A→TT1,   Y=B→TT1,   Z=TT1→B*.

All three have Q=t Q_B, and are rigid and pairwise nonisomorphic. The converse fixes X and exchanges Y,Z. Use nonconstant cyclic triples as before. At every depth the resulting family is closed under converse: on triples it reverses their cyclic order and takes the converse of each child. Let q be the number of current colors and f the number of self-converse ones. Then the numbers at the next depth are

    q'=(q³−q)/3,        f'=f(q−1),        (q,f)=(3,1) initially.

For the second formula a nonconstant triple has a free orbit under the three rotations. Its orbit is fixed by converse exactly when any one of its representatives is fixed by one of the three twisted reflections (converse followed by a rotation); the choice of reflection is unique per representative. Each such reflection has one fixed position and exchanges the other two, so it fixes f q triples, of which f are constant. Dividing the total over the three reflections by three gives f(q−1). Thus q begins 3,8,168, while f begins 1,2,14,2338. At every positive depth these are strong rigid tournaments; the count f measures order reversal, not a nontrivial automorphism of the directed tournament. The two symmetries can therefore be varied separately at exactly the same full Q.

### Any finite list of odd groups can share one full polynomial

Let G₁,…,Gₖ be finite odd-order groups and choose tournaments S_i realizing them by Moon's theorem. If k≥2, first form the ordinal join of S_1,…,S_k and retain their separate vertex intervals; every chosen color assignment may use disjoint palettes on the different intervals. Automorphisms respecting the colors therefore preserve these intervals, even if S_i are not individually canonical. Set the chosen interval monochromatic and all other intervals injectively colored, with all colors in use distinct across intervals. The color stabilizer is exactly G_i, by independent automorphisms within each ordinal interval. For a rigid member use all distinct colors. A common strong quotient can be obtained by surrounding the entire ordinal join with two rigid transitive fibres in the cyclic construction above. Give those two added fibres distinct vertex colors to force them fixed pointwise. The equal-Q inflation therefore yields a single polynomial fibre at one order, represented by strong tournaments with each of the prescribed groups and by a rigid member. The case k=1 is already covered above. This is a finite-list theorem, not a fixed finite-order family with infinitely many different groups.
