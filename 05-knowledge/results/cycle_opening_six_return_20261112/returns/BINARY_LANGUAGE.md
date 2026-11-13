# A binary indegree-one sublanguage, with its Fibonacci restriction

Status: **PROVED** restricted-language corollary of the inherited adaptive-ancestry theorem, with exact three-coordinate response laws. **VERIFIED** all binary words through ten births. This is not a new general unique-ancestry theorem, a new ordinary tournament census, or a Chan-Pak proof.

## Recovery before derivation

The earlier [ANCESTRY, A1-A3 and A5](../../adaptive_metagraph_reformulation_20261103/ANCESTRY.md) proves a stronger general mechanism: an adaptive birth from a strong parent has a unique legal inverse; ordinary histories are root-word orbits; and automorphisms fix all chronological newborns and the used original root vertices. It already gives (N-1)!/6 ordinary descendants from a triangle at order N>=4. The independent [ancestry audit](../../adaptive_metagraph_reformulation_20261103/ancestry_audit.md) verifies the scope and orbit formula.

The local two-indegree-one peeling argument in [FIBRE_GRAPH](FIBRE_GRAPH.md) is a particularly cheap recognition rule inside this inherited ancestry framework. Its contribution is the sharply specified neutral macro fibers and their closed response coordinates, not a stronger general descendant-count lower bound.

The earlier [classical unique-cycle grammar](../../operation_closed_synthesis_20261029/THEORY.md#ii-a-classical-grammar-with-a-new-finite-observer-experiment), attributed there to Shelor/Garey, has ordinary count F_(2N-6) at order N>=4; [TWO_ENTRY section2](../../two_entry_funnels_20261103/TWO_ENTRY.md) explicitly preserves that attribution. The Fibonacci restriction below has a different variable and a smaller language. The two counts must not be conflated.

## Two physical markings of the same birth

Keep the state v=(x,b,d), with H=x+d, b=F-H=d-s, and marked root indegree one. Operation L appends a vertex z beaten only by the old marked root r and marks z. Define K to perform the same graph birth but keep the old r marked. The exact maps are

    L(x,b,d)=(x+2b,d-b,x+d),
    K(x,b,d)=(2x+2b,0,d).

The L law is already proved in NEUTRAL_LANGUAGE. For K, every Hamilton path starting at old r must use z next, so there are d such paths. Deleting r leaves z as a global source over the old deletion, also with d paths. Hence the marked defect is zero; the total count is the same birth count 2x+2b+d, proving K's formula. Both marked roots have indegree one, and strongness is preserved.

Chronologically K then L is exactly the former two-vertex operation T: both newborns are beaten only by the same old r, the second beats the first, and the second is marked. Matrixwise T=L*K. This is a physical identity, not only equality of response coordinates.

## Free rooted words and exact ordinary classes

Fix a strong rooted source (T,r) of order n. After any birth there are exactly two indegree-one vertices: its old port p and its newborn z. L marks z; K marks p. Every later birth chooses its port from this marked pair.

Canonical peeling recovers every chronological newborn v_i and the port p_i used at its birth. The choice of marking immediately after birth i is recovered from the port p_(i+1) of the next birth: it was L if p_(i+1)=v_i, and K if p_(i+1)=p_i. At the final birth the distinguished output root resolves the last letter in the same way. Thus rooted L/K words are recovered uniquely from the rooted output. Rooted isomorphism classes at exactly c births are therefore exactly 2^c, not merely at least that many.

Forgetting the final root forgets exactly the final letter: L and K have the same physical last birth. All preceding letters remain recoverable from chronological ports. Thus at every positive c the exact number of ordinary output classes is

    B_c=2^(c-1).

There are no other identifications: the same intrinsic chronology is respected by every unrooted isomorphism. The two final root choices are distinct root orbits because every output automorphism fixes the chronological newborn and the preceding port individually. For c=0 there is the single original ordinary source. The selected ordinary generating function is

    x^n * (1-x)/(1-2x),

and the rooted-word generating function is x^n/(1-2x). For c>=1, the automorphism group of each ordinary output is Aut(T,r), by the inherited ancestry theorem or the direct extension argument in FIBRE_GRAPH. Consequently each has (n+c)!/|Aut(T,r)| labelings.

These are exact subfamily enumerations within already classified adaptive ancestry. They do not improve the earlier unrestricted triangle-descendant lower bound (N-1)!/6, nor supply a new A000568 term. Their value is that this binary sublanguage has a closed three-coordinate Hamiltonian response compiler.

## Complete count distribution from the response compiler

Let V_c be the integer multiset of three-coordinate states after all binary words of length c. Then

    V_0={v0},
    V_(c+1)=L_*V_c + K_*V_c,

where each pushforward retains multiplicity. For c>=1, the number of ordinary classes with full count h is exactly one half the total multiplicity of states (x,b,d) satisfying x+d=h. The factor two is the proven two final markings of every graph; it is not an assumption about distinct numerical states. The receipt computes this full distribution through c=10.

This representation is three-dimensional, not a claim of only three possible numerical states. The empty, L, and LL count observations recover all three coordinates as proved in NEUTRAL_LANGUAGE, so scalar H alone cannot replace this continuation interface.

## The L/T Fibonacci sublanguage

Under T=KL, words in L,T expand to binary words with no consecutive K and with final letter L (unless the word is empty). Conversely every such binary word has a unique parsing into tiles L and KL. Thus the L/T language is free on the fixed rooted strong source and has exactly

    C_c=F_(c+1)

ordinary outputs with c added vertices, with C_0=C_1=1 and C_c=C_(c-1)+C_(c-2). Its ordinary generating function is x^n/(1-x-x^2). The analogous state multiset recurrence is W_c=L_*W_(c-1)+T_*W_(c-2), with the evident empty/negative-cost conventions, and no factor two is needed because every expanded representative ends in L.

Give the initial marked root birth-clock time zero; this clock is not its literal vertex label in a general seed. Equivalently, the chronological port gaps are tiles (1) for L and (1,2) for T. Gap one chooses the preceding marked/new vertex; gap two returns to the port immediately before the first birth of T. The graph reconstruction through ports, rather than gaps alone stripped of their initial root, proves uniqueness.

The full binary language is strictly larger beginning at c=3: its ordinary counts are1,2,4,..., whereas the L/T restriction has1,2,3,... at c=1,2,3. Neither is the complete classical unique-Hamilton-cycle grammar. For the triangle source, the classical order-five count is3, whereas these restricted languages have2.

## The literal five-core boundary

Take the directed triangle rooted at vertex0. After one birth the underlying graph is the strong four-vertex core, with H=5. At two births, ordinary representatives LL and KL=T have H=9 and are isomorphic respectively to the two distinct five-vertex H9 cores recorded in the earlier five-core packet. The receipt gives the exact row tuples, marked states, and independently checked isomorphisms to S5a and S5b. Thus the binary marking choice has a direct graph-level realization in the existing five-core construction. The singleton fifth core belongs to the larger inherited core picture; it is not a strong starting source for the present theorem.

## Recovered book branch and Fermat candidates

On the rooted triangle, v0=(2,0,1), so direct iteration gives

    K^c(v0)=(2^(c+1),0,1),
    H=2^(c+1)+1,  order=c+3.

This is exactly the earlier transitive book family C3[1,1,TT_(c+1)], not a new prime family. An explicit isomorphism sends the two singleton blocks to old vertices2 and0, and orders the transitive block as newest newborn, ..., first newborn, old vertex1. The marked old vertex0 is the singleton receiving the arc from old vertex2. All arcs and the marked indegree-one root agree with the cyclic wrapper.

If c is odd, H=1+(2^((c+1)/2))^2, an exact near-square candidate. If 2^(c+1)+1 is prime, c+1 must be a power of two: an odd divisor d>1 of c+1 would factor it as U^d+1. Thus this factorization leaves only the Fermat-number indices c=2^j-1, at tournament order2^j+2. This is only a necessary primality restriction; no new Fermat prime or infinitude assertion follows. For j>=1 these are n^2+1 candidates; the j=0 value3 is outside that integer-square-plus-one form.

The transitive unit branch of the [strong six-step reservoir construction](../dynamics/PROOF.md) has module order 2m(|A|+|B|), hence is the slice c+1=2m(|A|+|B|) of this same rooted book family. This connection is through an explicit graph isomorphism. The exact period-six behavior still lives on the retained resource-and-cut carrier; the unmarked book graph alone does not acquire a six-cycle from this arithmetic identification.

## Verification and remaining target obligation

    python -B 05-knowledge/results/cycle_opening_six_return_20261112/returns/binary_language.py

The certificate checks all2046 nonempty binary rooted words through ten births, intrinsic decoding before and after reversed relabeling, the exact two-to-one final-mark quotient, every L/T word in the same cost range, Fibonacci coefficients, complete H-distributions, both H9-core identifications, and the all-K closed form plus explicit book isomorphism through c=16. Full subset-path counting independently checks all62 rooted binary outputs through five births. Later controls use the already proved matrices and exact adjacency/ancestry reconstruction, not exponential Hamiltonian DP.

Unique decoding gives an inverse from a generated graph to its word. It does not give an inverse from an arbitrary desired integer to a generated graph. The affine target-return problem therefore remains open despite the exact language classification.
