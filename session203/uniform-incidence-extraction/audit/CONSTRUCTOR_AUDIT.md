# Independent audit of the uniform-incidence constructor obstructions

Status: **PASS in the precisely stated constructor classes.** No eventual Hamiltonian-path spectrum threshold is proved.

## Frozen mathematical inputs

1. `session203-complement-gadget/PROOF.md`, SHA256 `f2d1a2b8c9463ac37df5c3242aca79edc5e0670a849d66577580a424724bd57c`.
2. `session203-structural-induction/PROOF.md`, SHA256 `06a16ce5044b629bb162b99ea5d46c70a7b6ffe73c8a37d2afc1f9fc1df25fd7`.
3. `session203-structural-induction/FINITE_MENU.md`, SHA256 `7d56566c7649fb047abda3de716b8f581904ffbe5c796d6815e44d250738e7e1`.

The first file's earlier hash was replaced only to correct its inherited citation from `geometry/PROOF.md` to the actual published Session192 `weighted_cuts/PROOF.md`. The auditor read that published capacity proof through origin/main and verified the correction. Its mathematical text and verifier were unchanged. The finite-menu file later gained reproduction commands and explicit non-application to the separate twelve-/forty-form occupancy conjectures; those additions were reviewed and do not change its theorem or constants.

## 1. Single uniform-incidence constructor, arbitrary core size

The proof's dichotomy is exhaustive.

- With zero sensitive core vertices, the graph and count are independent of the input bits, so cannot equal a nonconstant affine function of the index.
- With exactly one sensitive vertex, removing it leaves a fixed tournament. The physical weighted-cut formula gives every digit-pair coefficient exactly the negative exposed-path capacity, for either common sensitive polarity. The inherited lower bound c_ij>=2H(U-{i,j})>=2 is valid, including the empty remainder. Therefore a nonzero mixed difference defeats any affine function on the full Boolean cube.
- With at least two sensitive vertices, the words L,a,M,c,R have four distinct boundary digits when L,R are nonempty and |M|>=2. Each boundary condition has probability1/2 independently. Their count is exactly `3^m-(3+m/2)2^m+2m+3`. The induced-path insertion injection retains that lower bound in an arbitrarily large core.

At m=13 the expected count floor is 94,783, strictly above 73,732, the average of 6u+1 over 8193<=u<=16384. The minus target has still smaller average 73,730. This proves the claimed block obstruction; it does not say any individual count in that block is absent.

**Independent controls.** `independent_uniform_gadget.py` imports no author code. Literal permutations verify the cut identity on all 75 labeled parents through order 4 and all 1,098 cuts, and verify all 410 positive capacities. It also checks all 480 two-apex graphs with 2..5 digits under both sensitive polarities and both apex arcs, and independently counts the three-bin assignments. Normal and `-O` outputs agree. The all-order theorem rests on the displayed proof, not these controls.

## 2. Fixed graft retaining E(u): pumping and exact residual coefficient

The module-pumping argument is valid: two vertex-disjoint outside A->B arcs give a three-module-visit partial path, deterministic insertion retains those separations, and replacing its three module vertices by nonempty increasing subsequences gives `3^t-3*2^t+3` distinct paths. This contradicts affine growth along the specified binary contexts.

Both nonuniform incidence types have the stated disjoint cross-edge witnesses in the words10(1^t) and(1^t)0; their indices are respectively 3*2^t and2^(t+1)-1. The second pump10^t forces B->A and B->r. Hence deletion of c is exactly B->P->S with S=A union{r}.

The inherited ordered-root-amalgamation law applies twice, including B empty. Its non-singleton endpoint requirement is satisfied by the private-chain factor and by the right factor containing r. Expansion gives the displayed full formula, with positive g coefficient H(B)H(S). All remaining fixed-length terms are bitwise affine. The adjacent-bit mixed difference of g is−2, so the quadrilateral is exactly−2H(B)H(S), never zero.

**Independent controls.** `independent_fixed_graft.py` directly counts 180 literal grafts from 45 independently generated contexts, with B orders 0..3 and S orders 1..4, through total order 12. It checks every unit-source coordinate, the full formula including its constant term, and the exact quadrilateral. Normal and `-O` outputs agree. No author verifier is imported.

## 3. A finite target-dependent menu, including length-dependent cores

The shared pinned subcube is essential and correctly constructed. Prefix 0011 supplies distinct first(M),first(R) of the required polarities; suffix 0011 supplies distinct last(L),last(M). Unused prefix vertices go in L, unused suffix vertices in R. All t middle vertices may be assigned independently to three increasing bins. Thus every middle bitword has at least 3^t paths for every constructor with at least two sensitive vertices, regardless of its size or other arcs. The pinned bits themselves do not depend on the chosen constructor.

When this exceeds the affine target, only zero- and one-sensitive constructors can succeed on that subcube. The2q-by-[q(2q−1)+1] grid uses `(q+1)(2q−1)` free coordinates. Each column has at least q equal-color row pairs; there are exactly q^2(2q−1) possible(row pair,color) labels. Pigeonhole gives a monochromatic rectangle in two distinct columns. For a one-sensitive constructor the rectangle difference is the sum of strictly negative pair coefficients across two nonempty disjoint increments. A zero-sensitive constructor cannot hit two distinct affine targets. Both cases contradict success at all four corners.

The quantitative bound is correct: for6u±1, growth holds for m>=28 since3^20=3,486,784,401>3,221,225,473=6*2^29+1. The free-coordinate requirement is m>=8+(q+1)(2q−1). Thus any menu of at most 3 constructors fails at m=28, and a covering menu for sufficiently large m must satisfy the strict lower bound `q>(sqrt(8m−55)-1)/4`.

**Independent controls.** `independent_menu_anchors.py` separately constructs and validates 12,440 injected paths across 248 literal two-apex graphs, including every free middle word through length 4 and every polarity/apex-arc choice. It exhausts 291,305 color-class histograms for q<=8, confirming the sharp per-column pair floor, and checks 1,000 seeded colored grids for q<=10. It checks the exact m=28 constants. Normal and `-O` receipts agree. The graph injection and pigeonhole proof, not the finite grid tests, supply the all-order quantifiers.

## Scope safeguards

These are obstructions to one fixed incidence table or a bounded menu on a specified digit encoding. They permit positional incidences, other encodings, unbounded menus, and other graph reconstructions. They do not contradict or disprove the separate twelve-template central-arithmetic occupancy conjecture: that family has two freely selected correlated source indices and is not the same direct-affine digit constructor. They imply neither a missing unrestricted odd count nor an eventual threshold. The physical capacities are inherited, and no literature-priority claim is made.

No original proof or shared dirty-checkout file was edited by the auditor. Independent scripts, outputs and frozen input hashes are retained in this directory.
