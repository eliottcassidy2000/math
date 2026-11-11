# Focused review: bipartite domination below225

**Verdict: PASS for the stated new all-order bounded-range theorem and boundary-state formulas.** No material mathematical correction required. This review inspects `constructors/BIPARTITE_COVERAGE.md`, `bipartite_boundary.py`, and `bipartite_spectrum.cpp`. It does not rerun the large adjacency-matrix census, reopen inherited even/odd-graph theorems, or certify an unbounded four-exception conjecture.

## 1. The arbitrary-order cutoff is valid

The rooted-tree states(A,B,C) record exactly the information needed at an attachment edge. The formulas A(a+b+c), B(a+b)+Ca, Cb distinguish selection, external domination, and the unresolved root without dropping a case. Cutting an edge incident to the root produces two smaller rooted trees. Conversely attachment of any two such trees gives a tree. All size splits therefore generate every rooted tree, and merging only identical triples is lawful for all future attachments.

The computed order10 minimum225 is a finite exact input. It extends to every connected graph of order10 because adding edges only enlarges the collection of dominating sets. Every larger connected graph admits a connected induced ten-vertex prefix, obtained for example by an initial segment of a connected vertex exploration. Adding each subsequent vertex preserves an injection of old dominating sets by selecting the new vertex. Thus the lower bound extends to all larger orders; no unproved monotonicity of tree minima is needed.

The same argument extends the computed order7 minimum45, proving the smaller cutoff used for43.

## 2. The enumeration and component accounting match the theorem

Every nontrivial connected bipartite component below the225 count threshold has at most9 vertices and two nonempty shores. Swapping shore names permits the smaller size p≤floor(n/2). Enumerating all0/1 p-by-(n−p) matrices therefore covers it, despite retaining many duplicate colourings and isomorphisms. The C++ code counts dominating sets by literal closed-neighbourhood unions; the code range and integer types suffice for its stated n≤9 universe.

Disjoint component counts multiply. Isolates contribute1 and cannot evade the bound or alter a target. The closure code processes positive odd targets in increasing order; every factor in a nontrivial factorization is smaller than its product, so the recursive multiplicative test is complete. Allowing already disconnected census entries only adds valid generators. The claimed all-order spectrum below225 is therefore a legitimate finite-base theorem, conditional only on the explicitly provided exact census receipts, not an inference from missing entries alone.

In particular, primeness of43 leaves only one nontrivial component. Its order is at most6 by the45 cutoff, and the finite absence there proves arbitrary-order absence. The21 product case is handled separately; an excluded factor7 cannot rescue it.

Minor wording suggestion only: the sentence “every positive value in that range is already present” is clearest as “every attained positive odd value in that range already has a witness on at most9 vertices.” The surrounding text already supplies the intended restriction.

## 3. Boundary arithmetic is correctly guarded

- **Leaf:** selecting the old root permits either leaf state; not selecting it forces the leaf. This gives(2A,B+C,0). For a nonisolated root, adding a fixed neighbour injects C-sets into B-sets because a C-set contains no neighbour. Oddness then sharpens d_new≤2d to d_new≤2d−1. The isolated-root case correctly gives3d instead.
- **False twins:** the three classes in equation(4) distinguish a nonempty clone selection with an externally hit clone neighbourhood, forced selection of all clones when it is not hit, and no clone selected. They are exhaustive and yield the stated possibly negative affine offset. The complete-bipartite +1 ladder is lawful but fixes its coefficient; it is not a free digit system.
- **Shore hub:** the selecting-hub contribution is F_R. With no hub selected, the forbidden old dominating set with empty R is uniquely all of L, and it exists exactly when epsilon_R=1. The F_L,F_R updates follow from the same empty/nonempty distinction. Nonempty shores are necessary and explicitly assumed.
- **Two-shore batches:** the forced-all selection when the opposite shore has no selected vertex explains the−1 corrections. Equation(8) agrees with iterating the one-hub recurrence. It cannot be obtained by multiplying a one-hub increment blindly. The single-batch formula likewise retains its own boundary constant.

The final separation between bipartite support coverage and a zero-defect tournament adapter is explicit. Nothing here claims that every bipartite witness transports its count unchanged to a tournament.

## 4. Tiny independent controls

An independent Prüfer-sequence enumeration of every labelled tree through order5, with every possible root and direct subset domination counting, gives the following distinct-state counts and minima:

    n2:1 state, minimum3;
    n3:2 states, minimum5;
    n4:4 states, minimum9;
    n5:8 states, minimum15.

These agree with the beginning of the compositional table without reusing its recurrence. As a forced-all boundary hostile, two isolated old vertices on opposite shores plus one new hub on each shore form P4 and have9 dominating sets by direct subset enumeration; equation(8) gives4+1+1+3=9.

No large matrix enumeration was repeated. The verdict covers the mathematical cutoff, complete-state argument, product accounting, operation formulas, and declared computational universe. It does not establish a global spectrum beyond225 or a general tournament adapter.
