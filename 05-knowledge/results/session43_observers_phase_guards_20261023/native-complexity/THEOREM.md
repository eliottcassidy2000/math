# NP-complete native template realization, including a pinned regular-host variant

Status: proved reduction for the precisely defined decision problems below. The construction was proposed in the Session43 root synthesis and is verified and completed here. No literature-priority claim is made. This result concerns restricted observer families, NOT existence of an ordinary Hamilton path in a tournament.

## 1. The problem and the apparent paradox

Input: an actual tournament T, a supplied partition into directed C3 cells with cyclic labels0,1,2, and a spanning vertex template W containing each vertex exactly once. The template need not already be a directed path. Choose one phase a_i in F3 per cell, and replace every occurrence of local label u in cell i by u+a_i. Ask whether the resulting literal spanning order is a directed Hamilton path.

Call this **C3-template realization**. The independent cell rotations are internal automorphisms, not presumed automorphisms of the complete host. The choice is shared across ALL occurrences of the cell. An indexed-owner pin may restrict a cell phase to one value.

**Theorem A.** C3-template realization is NP-complete, even when the supplied C3 partition is equitable and every cross block is circulant or uniform. It remains NP-complete with one prescribed original-cell owner pinned.

**Theorem B.** A pinned version is NP-complete even on regular tournaments with a supplied equitable partition into directed C3 cells and one fixed singleton. The pins fix one actual owner in every cell of a supplied genuine suffix; only the other cells have free phases.

Neither statement contradicts the existence of a Hamilton path in every tournament. The constraint is membership in a specified shared-chart family. Our negative finite examples explicitly have ordinary Hamilton paths outside that family.

Membership in NP is immediate: the phase tuple is a polynomial-size certificate, and the resulting consecutive arcs and pins are checked in polynomial time.

## 2. Starting problem and connected inputs

We use the classical NP-completeness of graph3-colourability. A primary reference is Garey, Johnson and Stockmeyer, *Some simplified NP-complete graph problems*, Theoretical Computer Science1 (1976),237–267, DOI https://doi.org/10.1016/0304-3975(76)90059-1 . An accessible scan is https://www.acs.uwinnipeg.ca/ychen2/journalpapers/Garey-Johnson-2maxsat.pdf . This is the only externally inherited complexity premise.

We may assume the input simple graph G is connected and has at least one edge. To reduce general inputs to this case, connect one representative from each connected component in a chain. If all original components are3-colourable, permuting the colours independently in each successive component makes the new bridge valid; conversely restricting any colouring recovers all original colourings. Trivial small inputs may be handled separately. Thus this connected restriction remains NP-hard without requiring a separate cited theorem.

## 3. Bounded-degree equality copies

For every incidence (v,e) of G create a separate variable. For each original vertex v, join its incidence copies in a chain, each edge imposing equality of the two F3 phases. For each original edge e=uv, join its u and v copies with a phase-inequality constraint.

The resulting constraint graph L is simple and connected. Every variable has one original-edge incidence and at most two chain incidences, hence degree at most3. Its solutions correspond bijectively to proper3-colourings of G: equality chains identify all copies, and inequality enforces the original edges.

Write constraints as phase differences belonging to E:

    equality: E={0}; inequality: E={1,2}.

Both sets are invariant under negation, so traversal orientation does not change these constraints.

## 4. Eulerization without changing feasibility

Pair the odd-degree vertices of L arbitrarily. For each pair u,v, add a fresh vertex z and the two-edge path u-z-v. Both new constraints are unrestricted, E=F3. No parallel edge is introduced. The former odd vertices now have even degree at most4; each fresh vertex has degree2. The graph remains connected and every degree is even.

Choose one original copy root. Add two fresh leaves s,t joined to root by unrestricted edges. Now s and t are the only odd-degree vertices. All degrees are at most6, and an Euler trail from s to t exists.

A non-endpoint vertex of degree d occurs d/2 times in the Euler trail's vertex sequence; s,t each occur once. Consequently every vertex occurs r<=3 times. No constraint has changed its meaning; the fresh vertices carry arbitrary phases.

## 5. Compile the Euler trail into an actual tournament template

Replace each constraint variable by a directed C3. If a variable appears r times along the Euler trail, replace its first occurrence by a consecutive run of length4-r and its other occurrences by one vertex each. Assign local labels0,1,2 chronologically. Each cell supplies exactly three template vertices, with no repetitions. Every within-run step is local u->u+1, so it remains a valid C3 arc after every phase shift.

Every cross-cell boundary corresponds to one Euler-traversed constraint edge. The constraint graph is simple and every edge is traversed once; hence no cell pair receives incompatible block definitions from different boundaries.

Suppose that boundary goes from local label u in cell i to local label v in cell j, and its phase-difference palette is E. Define the whole actual cross block by

    (i,x) -> (j,y) iff y-x belongs to D := E+(v-u).

Then its actual shifted boundary arc is valid exactly when

    (v+a_j)-(u+a_i) in E+(v-u),

which is equivalent to a_j-a_i in E. Orient the remaining vertex pairs in the opposite direction. On all cell pairs not used by the constraint graph, choose an arbitrary uniform direction.

This produces an actual tournament, not a partial orientation. Its C3 partition is equitable because every defined cross block is circulant or uniform. The template is a literal spanning order. It becomes a Hamilton path precisely when every original equality/inequality constraint is satisfied.

There are O(|E(G)|) cells. Eulerization, padding, and the explicitly written tournament adjacency matrix take polynomial time. This proves Theorem A's NP-hardness.

The construction also preserves global phase translation: adding the same value to every phase leaves all differences unchanged. Pinning one original cell's actual indexed owner to phase0 therefore does not change feasibility, proving the pinned assertion.

### Exact lift multiplicity

If k fresh parity-repair vertices were added, the k repair phases and the two leaf phases are unconstrained. Thus the number of accepted template paths is exactly

    3^(k+2) * number of proper3-colourings of G.

Every phase tuple gives a distinct literal spanning order. Pinning one original copy to0 divides this count by3. This is a precise fibre statement, not an additional complexity classification.

## 6. Regularization with a fixed genuine suffix

Let T be the N-vertex tournament just constructed, and let d_i be its vertex scores. Form D(T) on two copies A,B of T and a fresh center c. Both A and B have internal tournament T. Put all A->c and all c->B. Across the two copies set

    A_i -> B_j iff j->i in T, for i!=j;
    B_i -> A_i on the diagonal.

Then

    score(A_i)=d_i+1+(N-1-d_i)=N,
    score(B_i)=d_i+(N-d_i)=N,
    score(c)=N.

So D(T) is regular of order2N+1.

Find one genuine Hamilton path Q in B by the standard insertion construction: when inserting a new vertex, put it immediately before the first existing vertex it beats, or at the end if it beats none. The preceding vertex, if present, then beats the new one, so directedness is preserved.

Use the full template

    [the original A-template], c, [Q in B].

Pin one actual indexed owner in every B cell, fixing each of its rotations to0, and fix the singleton c. All A-to-c and c-to-B boundary arcs are automatic. The suffix Q is already genuine. Therefore the full template is realizable if and only if the original A-template is realizable, with exactly the same number of free-phase solutions.

The partition consists of the inherited C3 cells in A and B plus {c}. It is equitable: original inter-cell blocks are circulant/uniform; cross-copy blocks between different cells are their reversed-transpose forms; the matched-cell cross-copy block is the reverse C3 adjacency with a zero diagonal, also row- and column-regular; and center attachments are uniform. This proves Theorem B with all types and pins retained.

No claim is made here for a regular-host variant without these pins or without the singleton center. Nor are the allowed phase tuples claimed to be complete-host automorphisms.

## 7. Positive structural regime: bounded interaction treewidth

The hardness identifies a limit, not a reason to abandon native compilation. For any supplied C3 template, each cross-cell consecutive pair produces a binary relation on its two shared phases. Multiple boundaries between the same cells intersect their relations; within-cell steps and unary pins are checked exactly. This is a finite-domain CSP on the cell-interaction graph.

Given a tree decomposition of width w, ordinary bag-state dynamic programming with3^(w+1) assignments per bag solves feasibility and exact counting in polynomial time for fixed w. Separators retain the SAME actual phase values on both sides. Factors can be assigned to a bag containing their endpoints, and forget steps sum over the forgotten phase. Additional finite guard states may be included at their explicit cost, but no arbitrary electrical guard is assumed additive or bounded-memory.

Thus small interaction treewidth gives a genuine positive native compiler. A fixed spanning template has at most three occurrences per cell, but this only bounds the interaction degree by6, not its treewidth. The reduction exploits precisely that distinction.

## Checks

The checker covers all connected labelled graphs on2,3,4 vertices, twelve seeded connected5-vertex graphs, and K5:56 source cases total. It independently rebuilds each native binary relation from the actual tournament and template, solves those relations using generic domain propagation/backtracking, and compares with direct source colouring counts and the exact lift factor. Small instances additionally enumerate all phase tuples literally.

It also constructs an ordinary Hamilton path in every host, including the empty-family K4/K5 cases, and verifies every regularized vertex degree, the equitable partition, actual suffix pins, and source/regular-template equivalence on deterministic phase controls. The proofs supply the all-size conclusions; these checks guard the implementation and normalization.
