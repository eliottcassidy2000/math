# Session85: exact missing-bit budgets at LCA and separator interfaces

**Status:** PROVED elementary linear-algebra statements; FINITE-EXACT supplementary controls. No literature-priority, nonlinear-guard, arbitrary-tournament entropy, or global arithmetic-dynamics claim. All logarithms and entropies below are base two. A “bit” is an independent binary dimension under the explicitly uniform affine law, not merely a named Boolean test.

## 1. What is actually inherited

Sources pinned at repository commit `83b601158076669edc8729ce461f150da1c3877e`:

- `05-knowledge/results/session56_rmq_lca_connections_20261025/tournament-hierarchy/REPORT.md`: the canonical strong-module tree of a tournament, with original leaf labels, named child ports, and the quotient tournament at each internal node. For leaves u,v, their arc is read from the quotient at their LCA. Tree shape alone loses the quotient orientation. The packet has prime order-five star-tree hostiles with H=9,13,15.
- `05-knowledge/results/session56_rmq_lca_connections_20261025/correlated-gluing/THEOREM.md`, with its `05-knowledge/results/session56_rmq_lca_connections_20261025/audit-gluing-debt/REPORT.md`: actual shared choice identities are glued once; correct compatible union is not independent resampling. Even-parity on three bits has full ab and bc projections, so their perfectly separator-consistent join invents one extra bit. Running intersection proves exact *conjunction of retained local constraints*, not recovery of a constraint discarded by marginalization.
- `05-knowledge/results/session56_rmq_lca_connections_20261025/debt-range/RESULT.md`: scalar guarded summaries do not permit shared owners to be spent twice. Its CRT/mod-30 and shared-port hostiles remain valid; affine rank below is not a replacement for nonlinear capacity constraints.
- `05-knowledge/results/session56_rmq_lca_connections_20261025/hierarchy-covers/REPORT.md`: same-label cluster families need exactly the chromatic number of their crossing graph many laminar hierarchies.
- `05-knowledge/results/fixed_twists_conductors_reverse_potential_20261025/lca/LCA_AND_PETERSEN_MERGES.md`: ten two-subsets of five labels have crossing graph L(K5), hence need five hierarchies; Petersen is their compatibility graph. Five here is a hierarchy-cover count, not automatically a bit count.

A bounded exact-term check of `01-canon/MISTAKES.md` found no Session56/LCA correction to these sources; the independent gluing audit accepts their stated scopes. They are inherited mechanisms, not discoveries of this session.

Closest mechanism: exact correlated join. Canonical hostile: even-parity projection. Corrected near miss: counting named guards as independent information. Least-used useful sidecar: the *space of constraints supported on a separator*, rather than only a separator's variable count. Live board: LCA ownership / affine row spaces / correlated guards / named ports / Fano duality / laminar obstruction.

## 2. Scope-loss theorem: exactly how many global bits were discarded

Let C=x0+K be a nonempty affine subset of F2^E, with W=K-perp its space of homogeneous parity checks. Thus C consists of x satisfying w·x=w·x0 for every w in W. Let B1,...,Bt be arbitrary coordinate scopes; they need not be laminar or cover E. Write V_B for vectors supported in B and define

    W_B = W intersect V_B,
    W_local = sum_i W_Bi,
    J = {x in F2^E : x|Bi belongs to projection_Bi(C) for all i}.

**Theorem (exact scope information loss).**

    J = x0 + W_local-perp,
    d = dim W - dim W_local >= 0,
    |J| / |C| = 2^d.

Consequently the projections reconstruct C exactly iff d=0. Exactly one missing binary linear constraint means d=1. A basis of W/W_local gives d parity guards that repair the join, and fewer than d parity equations cannot cut J down to C. Recovering the entire lost syndrome also requires at least d binary-valued observations because it distinguishes 2^d affine cosets. This is NOT a lower bound for arbitrary nonlinear membership tests: the single Boolean test “x belongs to C” can recognize C but is not a parity measurement or a full syndrome decoder.

**Proof.** A parity on B is constant on projection_B(C) exactly when its zero extension annihilates K, so the annihilator of that projection is W intersect V_B. A finite-dimensional affine subspace is determined by all its constant parity checks. Joining the projections therefore imposes exactly the sum W_local, with the consistent right-hand sides inherited from x0. Rank-nullity gives the ratio. Each independent added quotient check halves the current affine set. Within J, the W/W_local syndromes take all 2^d values uniformly, proving the repair lower bound. This is a representation-relative recovery statement, not a claim that arbitrary external nonlinear knowledge is absent.

Here d counts constraints lost by the projection representation, not the entropy of H(X) under the true prior uniform on C: that true syndrome is already fixed. The relaxed join J invents d syndrome degrees of freedom. Its entropy is dim C+d, whereas the original code has only dim C free coordinates.

This computes *which* constraint is missing, not only a general warning about observers. Gaussian elimination first finds each supported intersection W_B, then a basis for their sum, then explicit representatives of the quotient. If a supplied guard has right-hand side incompatible with C or the other guards, the feasible set is empty; augmented rank must be tested before quoting a free-bit count.

**Conditional targets.** For X uniform on C and any observed linear functionals spanning O, a target q·X is determined by those observations iff q is in W+O. Otherwise it is exactly balanced on every nonempty observation fibre. For several target rows Q, their remaining joint entropy is

    dim((W+O+span Q)/(W+O)).

Proof: on an observation fibre the translation space is K intersect O-perp, whose annihilator is W+O. A nonzero linear functional on that fibre has equal zero and one fibres. This is not an independence assertion for arbitrary nonlinear guards.

## 3. A junction tree turns the loss into an exact separator invoice

Suppose the scopes B_v form a tree with **running intersection**: bags containing each coordinate form a connected subtree. Let S_e be the intersection of the two endpoint bags, and assume the bags cover E. Then

    dim W_local = sum_v dim W_Bv - sum_e dim W_Se,                 (1)

and therefore

    d = sum_v dim projection_Bv(C)
        - sum_e dim projection_Se(C) - dim C.                    (2)

**Proof of (1).** Delete a leaf bag B with separator S. Its coordinates shared with all remaining bags are exactly S. Hence its local row space intersects the sum of the remaining local row spaces exactly in W_S: containment follows from support, and the reverse containment holds because its neighboring bag contains S. Apply dim(U+V)=dim U+dim V-dim(U intersect V) and induct. For (2), use dim projection_B(C)=|B|-dim W_B and the running-intersection identity sum|B|-sum|S|=|E|.

For two bags A,B, formula (2) is exactly

    d = I(X_(A\B); X_(B\A) | X_(A intersect B)).

The entropies equal affine dimensions because X is uniform. Thus a missing global parity is a *conditional dependence across the separator*, even when all projected separator assignments agree. This identifies the old Session56 parity hostile quantitatively: for C={a+b+c=0}, bags ab and bc have dimensions 2,2; separator b has dimension 1; C has dimension 2; d=2+2-1-2=1.

A precise algebraic chain complex is available without a loose topological analogy:

    0 -> direct_sum_edges W_Se -> direct_sum_bags W_Bv
      -> W_local -> 0.

The first map inserts a separator check into both incident bags (signs coincide over F2), and the second sums all checks in ambient F2^E. Exactness follows by the same leaf removal: the separator map is injective and every cancelling bag tuple is generated by these edge duplications. Hence duplicate shared-port constraints are subtracted once, and no first homology remains for this *specified tree complex*. The genuinely lost information is the separate quotient W/W_local. No claim identifies it with relative homology of an unspecified space.

## 4. Affine separator messages are small and algorithmic

The preceding theorem starts with a global C. There is also an exact construction from local affine equations, with no pre-known C. Assume the bags of the running-intersection tree cover every variable whose assignments are being counted. Put every equation into a containing bag, once, and root the tree. Coverage is essential: an unconstrained bit omitted from every bag would otherwise contribute an uncounted factor of two. A child-to-parent message is either empty or

    (an affine support A_v on separator S_v, an integer k_v),

meaning every supported separator assignment has exactly 2^k_v extensions into that subtree outside the separator.

At v, conjoin its owned equations and the child message supports on the coordinates of B_v; call the resulting affine set G_v. If empty, propagate empty. Otherwise project it onto S_v and set

    A_v = projection_Sv(G_v),
    k_v = sum_children k_child + dim G_v - dim A_v.                (3)

At the root S is empty, so the final count is 2^k_root (or zero if empty).

**Proof.** Every fibre of a linear projection of a nonempty affine space has size 2^(dim G_v-dim A_v). Running intersection makes different descendant branches share no unexposed variables; hence, conditional on the bag assignment, their extensions combine bijectively. Multiply extension counts and induct. Constant fibre size is essential; it fails for general nonlinear relations.

Gaussian elimination computes supports and fibre exponents. A separator of size s can store its affine support with at most s independent equations, rather than explicitly listing up to 2^s assignments; storage is polynomial in s. Computation is polynomial matrix work in bag width and equation-input size, not a universal claim of constant-time LCA queries or a benefit over global elimination for every instance. Ordinary junction-tree tables for arbitrary constraints can still be exponential in separator size. Treewidth chooses useful bags; it does not certify losslessness of projections of an unknown global relation.

## 5. LCA rank ownership: a guard is charged only for new rank

Fix a rooted leaf hierarchy on E. Assign each nonzero parity row to the LCA of its support; singleton rows belong to their leaf. Keep each row's actual coefficients and right-hand side. At a node v, let R_v be the span of assigned rows and U_v the sum of the already accumulated child row spaces. Child supports are disjoint. Define

    delta_v = dim(U_v+R_v)-dim U_v
            = dim R_v-dim(R_v intersect U_v).

Then

    rank(all rows) = sum_v delta_v,
    free bits = |E|-sum_v delta_v                               (4)

for a consistent affine system. Proof: dim U_v is the sum of child dimensions; telescope the recurrence from leaves to root. Zero rows require the ordinary consistency test 0=b, rather than an LCA.

Thus the LCA is where a cross-child guard is *tested for new information*, not a guarantee that it contributes one bit. A single parent parity either contributes zero rank (already forced, or contradicts the old right-hand side) or one rank (halves the feasible space). Repeating a shared guard twice contributes one bit, not two. Relations among several guards can cause additional savings even when no pair is a duplicate.

Tournament transfer is conditional, not universal: the inherited strong-module LCA tells where arc data live, but arbitrary prime quotient choices need not form an affine family or be uniformly distributed. Formula (4) applies when a specified tournament family/guard representation really is affine in identified arc coordinates. It cannot infer Hamilton counts from the hierarchy or turn an arbitrary blackness predicate into a parity equation.

## 6. Fano/Hamming is the sharp bridge, with the dual kept straight

Let H have the seven nonzero vectors of F2^3 as columns (integer labels 1,...,7), and let C=ker H. Then dim C=4, while W=row(H) has dimension 3. Every nonzero word in W has weight four, and its support is the complement of one Fano line. Therefore:

1. Every projection of C to at most three coordinates is full. In particular, retaining **all seven Fano-line projections** retains no parity check at all: W_local=0, J=F2^7, and d=3. Even seven perfectly consistent local views leave three lost syndrome bits, not one.
2. Retaining the seven complementary-four projections recovers W and hence C. Each four-scope projection has one parity check. But the seven named checks jointly have **rank three**, so all seven zero guards pass on a uniform unconstrained seven-bit vector with probability 1/8, not 1/128.
3. Among triples of distinct checks, exactly seven triples have rank two and 28 have rank three. A geometrical Fano line of check labels is dependent; any other triple is a basis. Three guards need not mean three bits.
4. C has 15 nonzero vectors: seven line words of weight three, seven complementary-four words of weight four, and the all-one word of weight seven. Its 14 support-minimal nonzero vectors are the matroid circuits; the all-one word is not a circuit.

The column-matroid circuits are ambiguity directions in ker H. Recovery checks belong to the rowspace/dual cocircuits. For this particular binary self-orthogonal presentation the seven weight-four supports occur on both sides, but that coincidence must not erase the primal/dual distinction.

**Why ordinary graph cycles are insufficient.** The seven columns represent a simple rank-three matroid. Any simple graphic matroid of rank three has at most six elements: for graph components with vertex counts v_i, rank=sum(v_i-1)=3, while edges<=sum binom(v_i,2)<=binom(4,2)=6. Hence this seven-element matroid is not graphic. Its nullity four is still a legitimate linear information count, but not the cycle rank of a graph carrying these seven elements and their exact circuit dependencies. Matroid circuits, not merely graph cycles, are the correct carrier.

**Why one laminar tree cannot simply hold every Fano test.** Distinct Fano lines cross; so do their four-element complements. Each seven-scope family has crossing graph K7 and requires seven same-label hierarchies under Session56's node-cluster model. This does not forbid a tree decomposition with enlarged bags, nor does it imply seven information bits. The counts answer different questions: hierarchy slots = 7, independent parity checks = 3, free Hamming coordinates = 4, nonzero codewords = 15. The inherited K5/Petersen five-hierarchy result is another geometry-of-scopes count, not an unexplained fifth free bit. For Fano, three basis parity measurements reconstruct all seven named parity values algebraically, although literal representation of all seven four-scopes as same-label tree nodes needs seven trees. This changes the allowed interface from literal scope nodes to syndrome-valued messages; it is not a five-to-three compression theorem for the distinct K5 pair-relation problem.

## 7. Hostile controls and reproducibility

- **Hidden parity:** even parity on abc, bags ab/bc, one lost bit despite a valid tree.
- **Correlated guards:** seven Fano complementary-four checks impose only three independent conditions.
- **Higher dependency with zero pairwise overlap:** one-dimensional row spaces generated by a+b, b+c, a+c have pairwise intersection zero but total rank two. Pairwise duplicate accounting alone misses a three-way relation.
- **Inconsistent duplicate:** imposing a+b=0 and a+b=1 gives no solutions; “one free bit” would be wrong if augmented rank were omitted.
- **Nonlinear exclusion:** an inequality, capacity allocation, Hamilton-path condition, or arbitrary guard need not have power-of-two fibres. Do not apply (1)--(4) without an affine representation and a specified probability law.

Run from the flat published packet directory (the development source `lca/check.py` is published as `check_lca_gluing.py`):

    python -B check_lca_gluing.py > lca-normal.json
    python -B -O check_lca_gluing.py > lca-optimized.json
    cmp lca-normal.json lca-optimized.json

The portable verifier uses explicit exceptions and only the Python standard library. It compares scope annihilators against actual projection/join enumeration in 11,820 exhaustive two-bag affine cases (all subspaces through dimension four, all covering bag pairs, two offsets); 1,200 deterministic arbitrary-scope cases; 1,200 LCA row-ownership cases; and 1,500 affine join-tree systems, with a genuine elimination-message implementation compared against enumerative separator messages, direct global assignments, and augmented global rank. The exact Fano controls and hostiles are separate. Normal and optimized receipts agree byte for byte: **127,511 explicit checks**, PASS. Enumeration is a supplementary implementation audit; the displayed all-size proofs do not rely on these finite bounds.

**Stopping boundary:** the linear interface and its lost-bit budget are closed here. Applying it to a named arithmetic or tournament guard family requires its actual affine matrix or a proved affine parameterization; neither the tree skeleton nor a numerical coincidence supplies that premise.
