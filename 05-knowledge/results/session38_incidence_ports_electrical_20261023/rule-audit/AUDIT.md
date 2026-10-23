# Independent rule audit: transport, incidence lifts, dyadic constraints

23 October 2026. Status: **PASS**. No Git state was changed by this audit.

This report distinguishes algebraic proof verification from independent finite controls. It does not independently authenticate the large order-eight ledger or recompute the order-fifteen path counts; those have separate auditors.

## 1. Projected transport uniqueness

The three equivalent conditions in `session38-root/PROJECTED_TRANSPORT_UNIQUENESS.md` are correct for the stated transportation model:

1. Every feasible integer table has the same prescribed linear observation.
2. Every positive-capacity residual directed cycle has zero signed observation cost.
3. Every residual strongly connected component admits compatible vector potentials.

The integral hypotheses are sufficient: every present finite residual capacity is a positive integer, so a simple residual cycle admits a unit augmentation. Differences of feasible tables give residual integral circulations and decompose into directed cycles. A directed root-to-vertex path in an SCC gives the potential, with well-definedness proved by appending the same return path to two competing paths. Inter-SCC arcs need no potential constraint.

Parallel channels must remain distinct. A forward and reverse traversal of the same cell cancels both table change and cost, whereas a residual two-cycle using different channels can matter. A merged loop is not a same-side residual self-loop: the two transport sides remain disjoint, and its observation coordinate may represent a loop. A black-orbit observation may legitimately carry coefficient two, conditional on the separately stated reflection-compatible channel model.

The six-cycle example has exactly the claimed unique table. All six permitted vertices have degree two and all margins are positive, so degree-one peeling stalls. This does not show that every capacity or tight-cut peeling implementation stalls; the root note correctly limits the claim.

Independent exhaustive controls use K_(2,2) with a second labelled channel at one endpoint pair. Every one of the 243 capacity vectors in {0,1,2}^5 is tested, including every feasible base table and every margin fiber. Six observations include the identity, merged parallel channels, total mass, the zero map, a signed scalar, and a mixed literal-black/literal-blue projection. The verifier compares brute-force fiber projections, explicit simple-cycle enumeration, and separately computed SCC potentials. Every nonzero-cycle witness is actually augmented and its bounds, margins, and changed projection are checked.

Results: 7,776 base tables, 46,656 table/projection comparisons, 11,298 verified changed-projection cycle witnesses, and 13,962 instances with more than one table but only one projection. A two-channel black-loop example retains two tables while projecting to two literal black lines in both cases. The unbounded six-cycle is checked separately.

These controls test the abstract model, not whether a proposed black-orbit channel system is genuinely reflection-compatible or whether every relaxed table is tournament-realizable.

## 2. Incidence collision propagation

The formulas in `session38-root/INCIDENCE_COLLISION_PROPAGATION.md` are correct. The key recoverable datum is c(u), the sum of second coordinates in the outgoing E row of u. Thus the global multiset of c values is retained even without explicit root tags.

Sources create one new row consisting of (c(v),0); sinks append (0,0) to every old row and create an empty row. Both preserve Hamiltonian path counts because a new source must be first and a new sink must be last.

For T -> b -> P_r -> T, independently counting triangles gives:

- c'(u)=c(u)+r for u in T; c'(b)=nr; c'(p)=n for p in P_r.
- Internal T edge incidences are unchanged.
- Edge incidences on T->b, b->P_r, and P_r->T are respectively r, n, and 1.
- Internal P_r edge incidences are zero.

These determine the new E solely from the old E and r. The composition is strong for every nonempty T and r>=1, regardless of whether T is strong.

A direct module-word count, independent of the inherited substitution kernel, verifies the Hamiltonian formula. With one T run, the possible words are T b P, b P T, P T b, and P T b P. They contribute [3+(2^r-2)]H. With two T runs, the words are T b P T and P T b P T. They contribute [1+(2^r-2)]F2, where F2 counts ordered spanning covers by two nonempty directed paths. Hence

H(C3[T,{b},P_r])=(2^r+1)H(T)+(2^r-1)F2(T).

This also covers singleton T and r=1, with impossible two-run cases contributing zero. Every labelled tournament on one through four vertices was tested for r=1,2,3: all 225 strong-lift signature and Hamiltonian formulas pass, as do all 150 source/sink extensions. Counts are computed by subset endpoint dynamic programming in this independent checker.

The separate primary circulant audit supplies differences Delta H=12,330 and Delta F2=-100,740. Substitution gives 113,070-88,410*2^r, which is strictly negative for every r>=1. There is no exceptional positive r for these seeds. This audit checks that arithmetic and its use, but does not claim an additional independent order-fifteen path-count calculation.

## 3. Fixed-alphabet counting boundary

`session38-root/FIXED_ALPHABET_BOUNDARY.md` is correct as written. There are at most (n+1)^L row types and at most (n+1)^(Ln) ordered lists of n row types, an upper bound also valid for outer multisets. Fixing an oriented spanning Hamiltonian cycle yields 2^(n(n-1)/2-n) labelled strong tournaments. Division by n! gives a valid lower bound on strong isomorphism classes.

The displayed inequality (n-3)/2 > (L+1)log_2(n+1) is sufficient. For L=12, its dyadic proof works starting at n=256: at k=8 the lower bound on the left is 126.5 while the upper bound on the right is 117, and the gap increases. This is only an explicit sufficient threshold; it does not assert minimality, different Hamiltonian counts, or failure of authenticated inventory reconstruction.

## 4. Dyadic tower and switch detection

The all-size arguments in `session38-electrical/PROOF.md` check algebraically.

For symmetric nonsingular integral A, use qualifying integer lifts from M_k={z:Az divisible by 2^k}. Replacing a lift by another with the same binary reduction changes the quotient pairing by an even integer, because symmetry transfers A to the other qualifying lift. This proves that q_k on V_k=M_k mod 2 is well-defined.

For y in M_k, write b=Ay/2^k. Being in the radical means b dot z is even for every z in M_k. Perfect duality on (Z/2^k)^r makes this equivalent to solvability of A t=-2^(k-1)b mod 2^k. The corrected lift y+2t then belongs to M_(k+1). Thus rad(q_k)=V_(k+1). The Smith-factor dimension formula follows via a unimodular change of domain. No symmetric diagonalization over the dyadics is needed.

The warning against simply adding 2^k times a binary correction to an arbitrary prior lift is essential and correct: lower lift bits may need alteration.

For a disjoint switch ab,cd -> ac,bd, the stated integer Laplacian update has the correct sign. Writing x=e_a-e_d and y=e_b-e_c gives Delta A=2(xy^T+yx^T). On the common blue kernel, the half-form change is alpha tensor beta plus beta tensor alpha, which is zero exactly when the two binary functionals are dependent and otherwise has rank two. The common unreduced blue kernel is always valid; identifying full-component quotients requires the component-partition caveat already present in the note.

Independent exhaustive matrix controls check all 108 nonsingular symmetric 2x2 matrices with entries in [-2,2], at k=1,2,3. All qualifying lifts modulo 2^(k+1) are enumerated. The 324 stages verify lift independence, radical equality with the next layer, and Smith dimensions from the 2x2 determinant/gcd formula. All 1,024 binary symmetric 4x4 matrices and all 24 labelled distinct-endpoint choices are also tested: 24,576 switch-form comparisons, of which 1,848 have rank two. The latter tests the functional statement on arbitrary symmetric kernels, hence includes the graph-kernel case.

The current identity and nonnegative mismatch-energy equivalence are algebraically correct conditional on the separately authenticated target data. In particular, computing the target current from the candidate itself would erase the test, and zero short-signature energy does not imply realizability.

## Reproduction

Run:

    python session38-rule-audit/check_rules.py
    python -O session38-rule-audit/check_rules.py

The frozen `receipt.json` is compared by default. `--write` intentionally regenerates it. Checks use explicit exceptions, not removable Python assert statements. There are no producer-checker imports. Ordinary and optimized replays both pass.
