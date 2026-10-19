# A three-free-color release dichotomy and tight-kernel certificates

19 October 2026. Status: proved structural results, plus two finite-exact arithmetic audits. These statements concern reconfiguration of already graceful trees on the fixed labels 0,...,m. They do not preserve abstract tree shape or prove the Graceful Tree Conjecture.

## Result in brief

Suppose a nonempty graceful completion fiber has exactly three free difference colors. Release a fixed edge e whose color becomes connected after the remaining fixed forest is contracted. Then:

- If e is not a bridge of the graph formed by the fixed edges and all edges of the three old free colors, the enlarged fiber is two-color connected.
- If e is a bridge of that graph, releasing it preserves the original exchange components exactly, with additional realizations of the released color inside each component.

Color 1 always satisfies the connected-color hypothesis. Thus this is an exact criterion for whether releasing a fixed color-1 edge repairs a disconnected three-free-color fiber. A new 15-label example shows that even a nonforced color-1 release can fail. A general tight-kernel argument supplies the proof and extends the sufficient criterion to more free colors.

## 1. Conventions and the classical input

Let E_d consist of every edge {i,i+d} on 0,...,m. A fixed forest F contains one specified edge of each fixed color. A completion selects one edge of each remaining color so that all selected and fixed edges form a spanning tree. A two-color move changes at most two selected edges, retaining their colors. All actual edges remain distinct after contraction, including parallel edges.

For a forest G write M_G for the graphic matroid after contracting G. For a set I of colors, put

    q_G(I) = rank(G union (union over d in I of E_d)) - rank(G).

For one color write q_G(d). Connectivity of a color means that its contact graph on the components of G is connected, equivalently its rank is the full contracted rank.

We use the following classical theorem, with attribution rather than a novelty claim. D. Fon-Der-Flaass, “Local equivalence of transversals in matroids,” Electronic Journal of Combinatorics 3(1), R24 (1996), Theorem 1, proves that independent transversals of color classes A_1,...,A_s in an arbitrary matroid are connected by one-coordinate replacements whenever

    rank(union over i in I of A_i) > |I|

for every nonempty I. The paper also gives diameter at most s^2-s+1. Its statement allows arbitrary matroids and subsets, hence applies to finite graphic multigraphs with distinguishable parallel realizations. Primary source: https://www.combinatorics.org/ojs/index.php/eljc/article/download/v3i1r24/pdf/ , pages 2–3; DOI https://doi.org/10.37236/1248 .

Only this strict-rank connectivity theorem is imported. The contraction reduction and arithmetic release consequences below are proved explicitly.

## 2. The maximal tight-color subset

Let A_1,...,A_r be classes in a finite matroid M with at least one independent transversal. For I contained in {1,...,r}, write A_I for their union. Call I tight when

    rank(A_I)=|I|.

The empty set is tight. Existence of a transversal implies rank(A_I)>=|I| for every I.

### Lemma 1: tight-kernel factorization

There is a unique largest tight subset K, the union of all tight subsets. Put k=|K| and H=cl_M(A_K). Then:

1. Every independent transversal restricted to K is a basis of H.
2. The remaining coordinates form an independent transversal of the classes A_i in M/H.
3. These two choices are independent: any such basis transversal on K combines with any such quotient transversal outside K.
4. For every nonempty J outside K,

       rank_(M/H)(A_J)>|J|.

Thus the outside-K coordinate graph is connected by singles, by Fon-Der-Flaass.

Proof. If I and J are tight, rank submodularity and the transversal inequalities give

    |I|+|J| >= rank(A_(I union J))+rank(A_(I intersection J))
              >= |I union J|+|I intersection J| = |I|+|J|.

(The first inequality follows even for overlapping color classes, since A_(I intersection J) is contained in A_I intersection A_J.) Equality forces I union J to be tight. Finite union gives K. Since an independent transversal uses k independent elements of A_K, those elements span H. The standard contraction rank identity now proves items 2 and 3. For item 4,

    rank_(M/H)(A_J)=rank(A_(K union J))-k>=|J|.

Equality would make K union J a larger tight set, contrary to maximality. QED.

This is an exact reduction of local obstruction, not merely a sufficient rank test. For any move width w>=1, the components of the independent-transversal graph correspond exactly to those of the basis-transversal graph on K with width w. Projection of any full move is a legal kernel move. Conversely, a kernel path lifts with all outside coordinates fixed, and outside coordinates can then be adjusted by singles. A kernel path preserves the flat H at every step.

### Corollary 1.1: small tight kernels

If |K|<=2, any two independent transversals are joined using singles and closure-preserving changes of at most two coordinates. First change the K coordinates directly to their target values; this changes at most two coordinates and preserves H, hence preserves the closure of the entire chosen set. Then use singles in M/H with the K coordinates fixed.

A larger kernel also works if it has a chain

    empty=K_0 subset K_1 subset ... subset K_t=K

of tight subsets with |K_j minus K_(j-1)|<=2. Replace successive blocks by their target values in increasing j. At step j the prefix remains a basis of the same flat cl(A_(K_j)), so the later coordinates remain independent and the closure of the full selection is unchanged. Afterwards adjust the quotient coordinates by singles. This is a sufficient tight-chain certificate; no necessity is asserted.

## 3. Adding a connected buffer

Suppose M has rank r+1, the old classes A_1,...,A_r have an independent transversal, and a further class B has rank r+1. Every old transversal P extends to a basis P union {b} with b in B, because B is not contained in cl(P).

### Lemma 2: lifting partial moves

Suppose old transversals are connected using singles and closure-preserving changes of at most two old coordinates. Then basis transversals with the additional B coordinate are two-color connected.

Proof. For an old single change P to P', choose any b' in B outside cl(P'). Change the old coordinate and, if needed, the B coordinate simultaneously. This changes at most two colors and has a valid endpoint. There is no requirement that the same buffer edge be valid for both P and P'. For a closure-preserving two-coordinate change retain the current buffer, which remains outside the same flat. At the end change the buffer alone to its target if necessary. QED.

### Theorem 2: a general one-edge-release certificate

Start with any nonempty arithmetic completion fiber with r old free colors C. Let e in F have color b, put G=F minus {e}, and suppose q_G(b)=r+1. Form the maximal tight subset K of the old colors using q_G(I)=|I|.

If |K|<=2, or more generally K has a tight chain whose increments have size at most two, then releasing only e makes the whole enlarged fiber two-color connected.

Proof. In M_G the old colors have an independent transversal supplied by any original completion, and the released color is a full-rank buffer. Apply Lemma 1, its corollary, and Lemma 2. Every other originally fixed edge remains unchanged, and paths whose endpoints contain e restore it automatically. QED.

A convenient sufficient condition for |K|<=2 is

    q_G(I)>=|I|+1 for every old-color subset I with |I|>=3.

No rank-staircase ordering is required. Testing all these subsets can be exponential in r; the statement is a structural certificate, not a claimed efficient universal release algorithm.

## 4. The exact three-color dichotomy

### Theorem 3

Let C have exactly three colors, suppose the original completion fiber is nonempty, and let e in F have color b. Define

    G=F minus {e},
    Q=F union (union over d in C of E_d).

Assume q_G(b)=4. Then:

(a) If e is not a bridge of Q, the enlarged fiber with only G fixed is two-color connected, with diameter at most 8.

(b) If e is a bridge of Q, the components of the enlarged two-color exchange graph correspond bijectively to those of the original fiber. In particular a disconnected original fiber remains disconnected. This half does not require q_G(b)=4.

Proof of (a). The original nonempty fiber implies that Q is connected. Since e is not a bridge, Q minus {e}=G union E_C is connected. The forest G has five components, so q_G(C)=4. Consequently the full set of three old colors is not tight; the maximal tight subset has size at most two. Theorem 2 applies.

For the bound let k=|K| and s=3-k. Change K in at most one closure-preserving move when k>0. The outside coordinates require at most s^2-s+1 singles by the cited theorem; lift these through B and finally correct B in at most one step. For k=0,1,2 these bounds are respectively 8,5,3. QED.

Proof of (b). The graph Q minus {e} has exactly two connected components U and V. All G edges and all old-color edges lie inside these two sides. Since G has five components, the old-color union has contracted rank three and spans the fixed flat H consisting of the contacts internal to U and V.

Every valid old-color selection has three independent edges and is a basis of H. Such a selection is valid in an original completion exactly when it is a basis of H, because e crosses U,V. In the enlarged fiber the B coordinate must likewise be an edge crossing U,V; every such actual b-edge works with every old basis selection. Thus the enlarged state set is the direct product of the original old-color state set and the nonempty set of crossing b-edges.

A full move projects to a legal old move changing at most two old colors. If it changes the buffer and one old coordinate, its projection changes at most one old coordinate; if it changes two old coordinates, the buffer stays fixed. Conversely, every old two-color path lifts while retaining any chosen crossing buffer edge, and all buffer realizations are connected by singles. The component correspondence follows. QED.

### Corollary 3.1: exact color-1 test

If 1 is fixed, its full edge set is the path 0--1--...--m. Its contact graph is connected after contracting any forest, so q_G(1)=4 automatically. For a disconnected three-free-color fiber, releasing its color-1 edge repairs the fiber if and only if that edge is not a bridge of Q.

In particular, if the graph E_C itself is connected, releasing color 1 always repairs such a fiber. This conclusion applies uniformly to every nonempty choice of fixed edges with those three free colors.

## 5. Two arithmetic audits

### 5.1 The prior 13-label obstruction

The fixed edges and free colors are recorded in [the first checkpoint](../dynamics_labeling_entry_20261019/GRACEFUL_RELEASE_AND_MATCHING_SUPPORT.md). There C={6,9,11} on 0,...,12. The graph E_C itself is connected. Release 3 has free-color ranks (4,2,2,2), so no rank staircase exists, but Theorem 3(a) now proves its connectivity without enumerating completions. The verifier independently confirms its 24 states form one component. Release 10 has no full-rank buffer and remains outside this theorem, although the earlier exact census proves it connected.

### 5.2 A nonforced color-1 release that fails

Use labels 0,...,14, free colors C={8,11,13}, and fixed edges

    1: 6--7      2: 1--3       3: 10--13     4: 8--12
    5: 0--5      6: 3--9       7: 6--13      9: 4--13
    10: 1--11    12: 0--12     14: 0--14.

The four fixed components are

    A={0,5,8,12,14}, B={1,3,9,11}, C0={2}, D={4,6,7,10,13}.

The nonloop free contact supports, including multiplicities, are

    color 8:  A--D (3 realizations), C0--D (1)
    color 11: A--B (3 realizations), C0--D (1)
    color 13: A--B (1 realization), A--D (1).

Their contact union is the tree B--A--D--C0. Its color-to-contact incidence graph is a chordless six-cycle. Hence the matching reduction from the first checkpoint gives exactly two exchange components:

    color types (8,11,13)=(C0D,AB,AD): 3 states;
    color types (8,11,13)=(AD,C0D,AB): 3 states.

Vertex 7 is incident to no old-free-color edge and has only the fixed edge 6--7 in Q. That edge is therefore a bridge. After releasing it, the two crossing color-1 choices are 6--7 and 7--8. Thus the enlarged fiber has 12 states in two components of size 6, as Theorem 3(b) predicts. The released color has full rank four and is not globally forced; it has fourteen possible edges before feasibility is imposed.

This disproves the claim that an arbitrary connected-buffer release, or specifically release of color 1, always repairs a three-color obstruction. It does not disprove the existence of some successful one-edge release: for this example, release 3 has 26 states in one component and is certified by Theorem 3(a), again with ranks (4,2,2,2). The verifier records the complete release census.

## 6. Scope of what remains open

The general existence question is still open here: must every disconnected nonempty integer-difference completion fiber have some single nonforced fixed edge whose release makes the entire enlarged fiber connected?

The new result resolves the outcome of every proposed connected-buffer release when there are exactly three old free colors. For larger numbers, a proper tight kernel can still contain a genuine three-color or larger obstruction, so full rank of the total old-color union alone is not enough. The tight-kernel or tight-chain certificate identifies the additional condition that the proof needs.

All proof claims above are independent of the exploratory 22,200-case support census. The published verifier checks exact rank calculations, tight-subset data, component counts, and the negative witness using every permitted integer-difference edge. No absence of a sampled counterexample is used as proof.
