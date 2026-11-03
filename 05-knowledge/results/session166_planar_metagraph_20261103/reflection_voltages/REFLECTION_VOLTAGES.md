# Reflection voltages classify the converse-sensitive metagraph patterns

**Status: PROVED, all n ≥ 3; FINITE-EXACT controls through n=7.**

There is an exact structural classification of the information lost when the original path-line metagraph is folded by tournament converse. Its black part is a two-sheet graph over the folded nodes, with the sheets identified at self-converse nodes. Once the reflection pairing of individual black lines is retained, the remaining data are binary cycle signs on the subgraph of non-self-converse nodes. Its blue part lies entirely in the fixed locus. This is a classification of equivariant graph structure, not a classification of tournament classes, a formula for all adjacency entries, or a proof that arbitrary abstract signed graphs are tournament-realizable.

The general language of signed graphs, switching, and double covers is standard. What is specialized here is the actual S151 marked-path line graph, including its loops, individual line identities, fixed vertices, and blue lines. No literature-priority claim is made.

## 1. Inherited structure and the exact object retained

Use S151's convention. A tiling t on the nonconsecutive pairs of the fixed path 0→1→⋯→n−1 specifies a tournament T(t). Let

- F complement every chord bit;
- S reflect chord positions (i,j)↦(n−1−j,n−1−i);
- c(t) be the tournament isomorphism class;
- J(C)=C^op be converse on classes.

The tilings are the **incidence flags** of the original multigraph. The fixed-point-free involution F joins the two flags of one individual line, and c attaches a flag to a vertex. A loop still has two distinct flags. Reflection satisfies

    FS=SF,             cS=Jc.

A line is blue exactly when S fixes its flags. Otherwise it is black. The chord (0,n−1) is fixed by reflection for every n≥3. Hence S(t)=F(t) is impossible: its bit at that chord would have to equal its complement. Consequently:

1. every S-fixed line is blue, and S fixes both of its flags;
2. S acts freely on black lines, including black loops;
3. no line has its two flags exchanged by S.

These facts, the black-pair theorem, and the fiber parities are inherited from S151 and its cited predecessors. They are not claimed as new results here.

Write V₀ for self-converse (SC) nodes of the folded graph and V₁ for its non-self-converse (NS) nodes. Above u∈V₀ there is one original class; above u∈V₁ there are two, denoted u⁺ and u⁻ after choosing an arbitrary sheet label.

Let B retain every individual blue line. Define Q to have the folded vertices and **one edge per S-orbit of black lines**. Its incidence flags are the S-orbits of black flags. F induces a free flag involution on this quotient, by the preceding observation. Thus Q is an ordinary multigraph with loops, not a graph with half-edges. The ordinary folded black metagraph has twice the number of lines at every Q edge; Q also remembers which two lines form a reflection pair. A multiplicity matrix alone does not retain these identities.

## 2. Complete reconstruction at the equivariant-graph level

Let Q₁ be the subgraph of Q induced by V₁, retaining every parallel edge and loop. For an edge e=uv in Q₁, choose either one of its two original black lines and put

    σ(e)=+1 if that line joins equal sheet labels,
    σ(e)=−1 if it joins opposite sheet labels.

Choosing the reflected line changes both sheet labels, so σ(e) is well-defined. There is no sign attached to an edge meeting V₀: its fixed endpoint carries no sheet choice.

**Reconstruction theorem.** The original colored path-line multigraph, together with its converse involution, is recovered up to equivariant incidence isomorphism from B, Q, the SC/NS vertex types, and σ. More generally, every finite abstract datum of this form produces exactly such an equivariant multigraph with no edge inversion. It need not produce a tournament metagraph.

**Construction.** Keep one vertex at u∈V₀ and two vertices u⁺,u⁻ at u∈V₁. Keep B unchanged. Replace each edge of Q as follows:

| Quotient edge | Two original black lines |
|---|---|
| Distinct SC vertices u,v | Two parallel u−v lines |
| Loop at SC vertex u | Two loops at u |
| SC vertex u and NS vertex v | One u−v⁺ and one u−v⁻ line |
| Distinct NS vertices u,v, σ=+1 | u⁺−v⁺ and u⁻−v⁻ |
| Distinct NS vertices u,v, σ=−1 | u⁺−v⁻ and u⁻−v⁺ |
| Loop at NS vertex u, σ=+1 | One loop at u⁺ and one at u⁻ |
| Loop at NS vertex u, σ=−1 | Two parallel u⁺−u⁻ lines |

Reflection exchanges the two displayed black lines, exchanges ± at NS vertices, and fixes the entire blue graph. Each line retains its two flags, including in the two loop rows.

These are seven sheet presentations, but only six local switching types: the two presentations between distinct NS nodes become isomorphic after changing the sheet label at one endpoint. The two NS-loop types remain different.

**Proof.** Choose one original flag over each chosen quotient flag. Its reflected flag determines the other lift. At an NS endpoint, the class of that flag fixes its sheet. At an SC endpoint, both lifted flags attach to the same class. F pairs the chosen lift with either the same or opposite sheet at another NS endpoint; these are exactly the two values of σ. All listed choices exhaust the possibilities because black edges have trivial S-stabilizer and there are no S-inverted edges. The construction reverses the quotient operation on vertices, flags and edges. ∎

This recovers the incidence graph, not the original bit-vector names of tilings or a tournament representative at each vertex. Those require the corresponding original data.

## 3. Global classification by cycle holonomy

Changing the sheet label at each NS vertex by g(u)∈{±1} replaces the signature by

    σ'(uv)=g(u)σ(uv)g(v).                                  (1)

This is switching, or a binary gauge transformation. For a loop the two factors g(u) cancel, so its sign is unchanged. Two lifts over a **fixed labelled quotient**, preserving each quotient flag and edge identity and the blue graph, are equivariantly isomorphic exactly when their signatures are related by (1).

Let Q₁ have v vertices, e edges and c connected components, including isolated NS vertices. Put β=e−v+c, with β=0 for the empty graph. Then the number of isomorphism classes of such abstract lifts is exactly

    2^β.                                                   (2)

To prove this, choose a spanning forest. Starting at one root in each component, switch recursively to make every forest edge positive. The remaining β edges have arbitrary signs. Their signs are the products of the original signs around the corresponding fundamental cycles, including one-edge loops and two-edge parallel cycles. Every further switching that preserves positive forest edges is constant on each component and changes no edge sign. Thus these β signs are complete and independent invariants.

Equivalently, the signatures modulo switching are H¹(Q₁; Z/2). This is the elementary graph cohomology identity C¹/δC⁰, proved by the forest argument above. If the quotient itself is allowed to be relabelled, one must further take the orbit under its type-preserving, incidence-preserving automorphisms. Formula (2) is deliberately for a fixed labelled quotient.

There is also a precise categorical count. If a is the number of black quotient edges with both endpoints in V₀, loops included, then the equivariant automorphism group over that fixed quotient is

    (Z/2)^(c+a).                                           (3)

The c factors exchange sheets throughout individual components of Q₁, extending uniquely across edges with an NS endpoint. Each SC–SC black pair may independently be exchanged, giving a further factor. Blue flags and all quotient flags are fixed under this definition, so there are no extra loop-flag flips or permutations of distinct parallel quotient edges. This specifies the isomorphism classes and their isotropy groups, rather than discarding the symmetry data of the lift.

Geometrically, the black graph is an ordinary two-sheet cover away from the SC vertices, with the two fibers identified at each SC vertex. Blue edges lie in the fixed locus. Thus binary monodromy lives on the NS part; it does not continue as a local system through a vertex where the sheets have been identified.

## 4. Recurring cycles and the exact role of a blue step

An NS-only black walk has a unique lift from either chosen starting sheet. Its final sheet is the starting sheet multiplied by the product of its edge signs. Consequently a simple cycle of Q₁ lifts to:

- two cycles of the same length if its sign product is +1;
- one cycle of twice the length if its sign product is −1.

This includes loops and parallel-edge cycles under the multigraph convention. For a connected component of Q containing an SC vertex, its entire black lift is connected: paths to that fixed vertex connect the two sheets. A component with no SC vertex has two lifted components if every cycle is positive, and one otherwise. These are classifications of the lifted black graph, not of a dynamical system on tournaments.

**Blue-word theorem.** Fix any finite word w in the colors blue and black that contains at least one blue step. For any original classes C,D, the number of individual-line walks with color word w from C to D equals the number from C to D^op. It also equals the count from C^op to D.

**Proof.** Choose the first blue step. Its ending vertex is fixed by converse. Apply S to the suffix after that step, preserving the prefix. This is an involutive bijection of walks, preserves the color word, and changes the terminal class from D to D^op. If the suffix is empty, D is already SC and the identity works. Reflecting the whole walk gives the equality at the initial class. Loops and parallel lines are retained, so this is an incidence-level bijection, not merely a statement about unweighted adjacency. ∎

A blue step therefore erases the distinction between the two converse sheets for subsequent walk counts. This statement concerns counted graph walks. It does not say that the unmarked graph supplies a deterministic controller, or that a chosen arithmetic program forgets its state.

## 5. The exact representation-theoretic split

Let A_B,A_K be the original symmetric incidence matrices, with a loop contributing 2 on the diagonal. Let M_B,M_K be their folded incidence matrices, summing all original entries over each ordered pair of converse orbits. Let D be the diagonal matrix with D_uu=1 on V₀ and D_uu=2 on V₁. Finally let A_Q be the ordinary incidence adjacency of Q and A_σ the signed adjacency of Q₁, where a signed loop contributes 2σ(e) on the diagonal.

Then

    M_K=2A_Q.                                              (4)

The converse permutation J commutes with both original matrices. The invariant subspace has orthonormal basis consisting of each SC vertex and (u⁺+u⁻)/√2 for each NS node. The anti-invariant subspace has basis (u⁺−u⁻)/√2. Directly checking the seven rows of the reconstruction table gives, for all scalars x,y,

    xA_B+yA_K  ≃
      D^(−1/2)(xM_B+2yA_Q)D^(−1/2)  ⊕  yA_σ.             (5)

Here ≃ is orthogonal similarity. In particular the blue matrix is zero on the anti-invariant subspace, and black converse-sensitive modes are exactly the signed adjacency modes of Q₁. The ordinary merged matrix by itself is not the invariant restriction: the orbit-size normalization D is necessary.

An equivalent rational basis uses the orbit indicators and sheet differences. In that basis the invariant block is D^−1(xM_B+2yA_Q), which is an integer matrix at x=y=1. Thus the characteristic polynomial factors exactly into invariant and signed parts. This is standard C₂ representation theory applied to this particular incidence object.

For ℓ≥0, putting A=A_B+A_K gives the character identity

    tr(A_σ^ℓ)
      = 1/2 [tr(A^ℓ) − tr(JA^ℓ)].                          (6)

The right side compares ordinary closed walks to walks ending at the converse of their start. Loop traversals use the incidence convention above. More generally any colored-word matrix containing A_B annihilates the anti-invariant part; entrywise it has equal columns at D and D^op and equal rows at C and C^op. This recovers the blue-word theorem algebraically. The only converse-sensitive words are all-black words, whose anti-invariant action is A_σ^ℓ.

## 6. Necessary physical constraints and the remaining gap

S151 proves that a class fiber has odd size H(C)/|Aut(C)| and that an SC fiber has an odd number of S-fixed flags. Translating those inherited facts into the lift model gives

    deg_Q(u) is odd for every NS node u,
    deg_B(u) is odd for every SC node u,
    original black degree at an SC node = 2 deg_Q(u).        (7)

Loops count twice. With m=(n−1)(n−2)/2 and d=⌊(n−1)²/4⌋, the inherited global counts become

    |E(B)| = 2^(d−1),
    |E(Q)| = (2^m−2^d)/4.                                 (8)

These are necessary constraints on a physical tournament quotient. They do not constrain the signature beyond the available quotient. Changing a cycle sign preserves all vertex degrees, the folded colored graph, the global counts, and all these parity tests. This exhibits exactly why those tests cannot recover the original equivariant graph.

An especially visible lost bit is the sign of an NS loop. Its positive lift has two original loops; its negative lift has two original lines between converse classes. Both fold to the same two black loops at the merged node. Both loop signs occur in the actual n=6 graph. A proposed all-order reconstruction that retains only merged multiplicities necessarily loses this information.

What remains unresolved here is the physical realization problem: determining the base Q, the blue graph B, their particular signature, and tournament representatives without an exact tournament/path enumeration or another structural theorem. Nothing here asserts that the other 2^β−1 abstract signatures over an actual quotient occur as metagraphs of another tournament order. The object at each order n is fixed; the abstract family describes information the quotient forgets.

## 7. Exact finite controls

Run `python3 check.py` and `python3 -O check.py`. Their standard output agrees byte-for-byte with `results.json`. The checker uses no third-party packages and no Python `assert` statements. It enumerates every fixed-path tiling and exact tournament isomorphism orbit through n=7, constructs reflection orbits of individual black lines, reconstructs every original colored incidence entry, and verifies both matrix intertwining identities using integers and rational numbers rather than floating point. It checks every colored word through length four for n≤6, including all endpoint equalities for words containing blue and the all-black character trace identity. It checks all possible gauges for n≤5 and independently enumerates all signatures and gauges for small abstract multigraphs with loops, parallel edges, and disconnected components.

| n | Original classes | SC nodes | NS nodes | Blue lines | Black pair-orbits | NS edges | NS components | NS cycle rank |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 2 | 2 | 0 | 1 | 0 | 0 | 0 | 0 |
| 4 | 4 | 2 | 1 | 2 | 1 | 0 | 1 | 0 |
| 5 | 12 | 8 | 2 | 8 | 12 | 0 | 2 | 0 |
| 6 | 56 | 12 | 22 | 32 | 240 | 157 | 1 | 136 |
| 7 | 456 | 88 | 184 | 256 | 8064 | 5104 | 1 | 4921 |

Thus the first nonzero cycle rank in this finite range is 136 at n=6. This is a finite onset, not a theorem that ranks are monotone or that the NS graph remains connected at all larger orders. The fixed labelled n=6 quotient admits 2^136 abstract switching classes; n=7 admits 2^4921. No physical realizability is inferred from these counts.

At n=6, the black quotient has 3 positive NS loops and 9 negative NS loops. At n=7 these counts are 13 and 35. For explicit n=6 witnesses, tiling 11 and its complement 1012 both have class 22, while their reflections 776 and 247 both have its converse class 18: this is an untwisted loop-orbit. Tiling 13 and its complement 1010 have classes 26 and 20, while reflections 584 and 439 have classes 20 and 26: this is a twisted loop-orbit. Exact witnesses at n=7 are also recorded in `results.json`. Counts of signs on nonloop edges depend on the sheet convention; the cycle rank and loop signs do not. The number of negative normalized fundamental cycles also depends on the chosen spanning forest and is not claimed as a canonical invariant.

## Sources and provenance

- S151, [original/merged structure](https://github.com/eliottcassidy2000/math/blob/6b24421039fd2ffa820e9b56d107e458d75498ed/05-knowledge/results/session151_spectrum_metagraph_20261103/metagraph/STRUCTURE.md): precise F,S line convention, fixed-coordinate obstruction, black pairing, fiber parity, all-order line counts, and warnings about original versus merged loops. The implementation follows its exact canonical-code convention.
- S151, [colored incidence](https://github.com/eliottcassidy2000/math/blob/6b24421039fd2ffa820e9b56d107e458d75498ed/05-knowledge/results/session151_spectrum_metagraph_20261103/incidence/INCIDENCE.md): path-orbit and filtered double-coset multiplicities, diagonal factor, and folded incidence matrices.
- S160, [actual metagraph hub](https://github.com/eliottcassidy2000/math/blob/6b24421039fd2ffa820e9b56d107e458d75498ed/05-knowledge/results/session160_ordered_programs_20261103/metagraph/NOTE.md): actual line walks, the need to retain marked paths and selected line identities, and the distinction between graph structure and externally supplied programs.
- S165, [FI interfaces](https://github.com/eliottcassidy2000/math/blob/6b24421039fd2ffa820e9b56d107e458d75498ed/05-knowledge/results/session165_geometry_representation_20261103/fi_interfaces/FI_INTERFACES.md) and [guard geometry](https://github.com/eliottcassidy2000/math/blob/6b24421039fd2ffa820e9b56d107e458d75498ed/05-knowledge/results/session165_geometry_representation_20261103/guard_adjunctions/NOTE.md): methodological inspiration for retaining the correct category and the exact information lost by a quotient. The present graph lift is a different object from those FI interfaces and chamber banks.
- Thomas Zaslavsky, [Signed Graphs and Geometry](https://people.math.binghamton.edu/zaslav/Tpapers/sggm.pdf), and [Petersen Signed Graphs, §15](https://people.math.binghamton.edu/zaslav/Tpapers/petersen.pdf): standard signed-graph switching and double-cover background. The proofs above are self-contained, and no priority for the general covering or switching machinery is claimed.
