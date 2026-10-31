# Independent audit: rooted forests, CG words, and the G32 origin

**PASS.** The reviewed explicit bijections, weighted refinement, boundary
conventions, and Schur/Green responses are correct. The origin-model review
is in the separate companion origin/numeration packet. All-order proofs are distinguished
from finite controls below. The work credits classical results and does
not claim historical novelty or a new Hamiltonicity theorem.

## 1. Actual objects, not merely matching counts

For a rooted path component of k vertices with root in position r, the code
`D^(r-1) M D^(k-r)` has k-1 dominoes and one monomer. A separator monomer
between components gives an odd tiling of 2m-1 cells. Conversely, an odd
tiling has an odd number of monomers. Alternating root-marker and separator
roles uniquely recovers the ordered blocks and root positions. This proves
both directions for every m>=1, including the singleton boundary.

The word transducer reads LSD first while processing pairs of cells from
right to left. Its Free and Locked cases exhaust the possibilities for
those cells in any tiling: a domino inside the pair; two monomers; a monomer
and left-crossing domino; or, if the right cell is already occupied, either
a monomer or a left-crossing domino at the remaining cell. These are exactly
the five table rows. A left-crossing domino accounts for the next Locked
state. The final one-cell boundary is consequently forced, with no gaps,
overlap or protrusion. Reversing this local decision at every cell pair
provides a genuine inverse, not an argument from equal cardinalities.

The rejected transition is precisely a second 2 without an intervening 0;
1's preserve the lock. Therefore the language is exactly the words avoiding
`2 1* 2`. This forbidden-factor language is reversal invariant, justifying
LSD processing of displayed MSB-first words.

**Direction is nevertheless a required state coordinate.** Reversal preserves
acceptance, not the final automaton state. Word `201` ends Free when scanned
MSB first, but Locked in the forest's LSD-first scan. Its tiling is `DMMMD`;
reversed word `102` has tiling `MDMDM` and the opposite final state. A forest
pending-domino state must not be identified with an MSB append-digit or
Beatty eligibility state without specifying the reading direction.

The fixed-length CG integer interval follows independently from the free/
locked suffix intervals: length k has free range `[0,A_k-1]` and locked
range `[0,A_k-A_(k-1)-1]`, with `A_-1=0`. At the next digit the legal blocks
are disjoint adjacent intervals because `A_k=3A_(k-1)-A_(k-2)`. This proves
rank and uniqueness, rather than inferring arithmetic order from forest
counts. Thus length L corresponds to `P_(L+1)`, not `P_L`.

## 2. Preserved statistic, matrix, coefficients and boundaries

Each nonzero digit creates exactly one domino; each domino is one forest
edge. Digit 2 creates one edge, so **digit sum is not edge count**. A tiling
of 2L+1 cells with k dominoes has 2L+1-k tiles total; selecting its k domino
positions gives exactly `binomial(2L+1-k,k)` codes. This proves the coefficient
formula for every L and k, independently of the finite histograms.

The Free/Locked transition matrix is `[[1+t,t],[1,t]]`: the two Free-to-Free
choices have weights 1 and t, and every crossing/within-pair domino has
weight t. Starting row `(1,0)` and final column `(1,1)^T` give the polynomial.
Cayley-Hamilton gives recurrence coefficient 1+2t and determinant t². At
t=1 the matrix is the square of the Fibonacci transfer matrix, proving
`f_m(1)=F_(2m)` with the correct initial values.

For completeness the general weighted matrix-forest identity can be checked
without extrapolating any finite census. Let B be an oriented incidence
matrix. Cauchy-Binet expands `det(I+t BB^T)` as a sum over column subsets F
of `t^|F| det(B_F^T B_F)`. Cyclic subsets have dependent columns and contribute
zero. For a tree on k vertices, every row-deleted square incidence minor is
±1 (successively remove a leaf), so another Cauchy-Binet application gives
its Gram determinant k. Forest components separate into blocks, giving
`prod(component sizes)`, exactly the number of root choices. This proves
`det(I+t L)=sum_rooted_forests t^(number of edges)` for all finite simple
graphs. This classical result is credited to the matrix-forest literature.

The formal recurrence continuation has f0(t)=0; the actual empty graph has
one rooted forest and determinant 1. The empty word maps to `P1`, not the
empty graph. Prepending a leading zero prefixes two monomers, hence one
rooted singleton component and shifts the old path labels by one. It
preserves integer value while changing forest order. Fixed-length numeric
successor only reorders forests of the same path; it is not vertex arrival.

The final curated operation corollary is also valid. Right-appending 0 is
the initial Free-to-Free LSD step and creates exactly two rightmost monomers,
so it appends a rooted singleton. Its value is the shifted weight sum, equal
to `ceil(beta*n)` by the companion Beatty theorem (and tested here with exact
integer square roots). In contrast, leading-zero padding prepends a singleton
and preserves value. Canonical positive words have no first singleton;
with exactly one edge their first component must be `(2,1)` or `(2,2)` and
all later components are singletons, precisely the words `10^k` and `20^k`.
The recurrence preserves the quadratic expression Q, giving values 1 and 4
on these two subclasses. No converse classification of all norm solutions
is being asserted. Individual forest rank and total path-forest count remain
separate quantities; the illustrative ranks 1,3,8 occur on path orders 2,3,4.

## 3. Schur and named-port Green identities

For `A=I+L(G)`, positivity follows from
`x^T A x = sum x_v² + sum_(uv)(x_u-x_v)²`, so A and the enlarged diagonal
matrix M are invertible. Appending a vertex with distinct old neighbor set S
has exact block matrix `[[M,-b],[-b^T,1+s]]`, where `M=A+diag(b)` and b is
the indicator of S. Block elimination gives
`det(M)(1+s-b^T M^-1 b)`, equivalently the division-free adjugate formula.
This covers the empty neighbor set, isolated old vertices and disconnected
old graphs; connectedness is not a hidden hypothesis.

For U selecting the named ports and `C=U^T A^-1 U`, determinant lemma and
Woodbury give `det(M)=det(A)det(I+C)` and
`U^T M^-1 U=C(I+C)^-1`. Since `s=1^T1`, the remaining factor becomes
`1+1^T(I+C)^-1 1`, proving the stated Green identity. The zero-port case
uses the usual determinant of the 0-by-0 matrix, namely 1. For one port,
the identity becomes `2det(A)+adj(A)_(vv)`.

For P3 the displayed inverse is correct: endpoint diagonal 5/8 and middle
diagonal 4/8, with determinant 8. Thus a leaf at an endpoint yields 21 and
a leaf at the middle yields 20. The marked port, not just arrival degree
or old count, is necessary. The named-port inverse submatrix is sufficient
for the stated one-vertex forest-counting query, not for future queries at
forgotten ports and certainly not for Hamiltonian continuation.

## 4. Replays and independently constructed controls

The producer's normal/-O replays are byte-identical, with **945,738 checks
per run**, covering 265,720 ternary words of lengths 0..11, 75,024 odd
tilings and 75,024 rooted path forests, all 1,100 labelled simple graphs
of orders 0..5 and all 1,099 one-vertex arrivals from old orders 0..4.
The curated producer test includes three additional 201/102 direction checks,
plus the right-singleton/Beatty operation and canonical one-edge controls.
Current exact check totals are bound in the publication manifest and receipts.

The separate `independent_audit.py` imports no producer implementation.
It reconstructs the code from a completed cell covering, enumerates path
forests from compositions and root products, and enumerates general rooted
forests as parent functions with only self-rooted directed cycles allowed.
That method is different from the producer's edge-mask/union-find count.
Determinants use rational Gaussian elimination and Green entries use signed
cofactors, not the producer's Bareiss determinant and linear solver.
Its independent graph and arrival universes are likewise unfiltered.

The original combined audit also freshly enumerates every degree-two edge
set for G31, G32 and G33. Completeness follows because degree propagation
only forces or excludes an edge when the residual degree equation forces
it, and each branch enumerates every subset that completes one vertex to
degree two. These branches are disjoint and exhaust all completions. Global
connectedness is examined only after the degree equations are satisfied.
This is finite-exact evidence, not an all-order Hamiltonicity argument.

Independent rational determinants confirm the five actual square-sum graph
counts in the producer table, including `R(G12)=9240=8*21*55`. Negative
boundaries include digit-sum confusion, the true-empty/formal-zero mismatch,
MSB/LSD terminal-state confusion, and the false scalar-only growth rule.
All checks use explicit exceptions and remain active under Python -O.

## 5. Credit, verification limits and access

Primary sources independently inspected on 2026-10-31:

- [Burns, Chung-Graham and Zeckendorf representations, Theorem 2](https://arxiv.org/html/2502.18870v1): classical Chung-Graham grammar and numeration.
- [Chebotarev–Shamis, Matrix-Forest Theorems](https://arxiv.org/html/math/0602575v2): weighted forest determinant and inverse interpretation.
- [Chebotarev, Spanning Forests and the Golden Ratio, Theorem 1](https://arxiv.org/html/math/0612792v4): path counts F_(2m), with earlier fan-graph history noted there.

These sources support the classical attributions; the local proofs above
verify the specified maps and boundaries. No Hamiltonicity result follows
from a determinant or CG digit normalization. The G32 threshold requires
the separate finite exclusions and CITED Gerbicz all-order theorem.
There was no tool access denial. No source-worktree or remote edits were
made by this lane. Publication is owned by the parent session.
