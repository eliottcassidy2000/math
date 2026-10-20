# Odd contraction, twin-class cycle indices, and unbounded defect weights

20 October 2026. **PROVED: elementary local formulas and counterexamples, freshly audited with an independent literal-packing implementation.** This note takes the marked-triple carrier in the incoming THM-4723 seriously. It derives a compressed formula and tests the proposed split/join cancellation beyond the initial census. The least-order boundary and additional failures of the proposed sign partition are proved in `SHARP_DEFECT_BOUNDARY.md`. No unrestricted odd-group equality is asserted, and no new canon identifier is assigned by this staged report.

## 1. Inheritance and the exact carrier

The incoming `01-canon/theorems/THM-4723-sylow-profile-form-of-hyp-9340-and-the-defect-reduction-for-the-three-cycle-type.md` defines, for a graph X,

    b(X) = sum_(j>=1) (-1)^(j+1) times the number of Aut(X)-orbits
           of j disjoint independent twin triples whose contracted
           edge-sign is trivial on their setwise stabilizer.

Write e(X)=1 if X is sign-even and has an independent twin triple, and 0 otherwise. The incoming reduction identifies sum_X b(X) with the number of tournament types admitting a single 3-cycle symmetry. The inherited prime-local Walsh theorem already proves sum_X(b(X)-e(X))=0. The issue here is a structural cancellation, not whether that global equality remains open.

The proposed rule splits a large independent twin class into two triples and a remainder, joins the two triples, and tries to pair negative and positive defects. Initial examples through eight vertices had coefficients only -1 and +1. That unit-weight premise fails.

## 2. Odd independent-module contraction preserves the sign character

**Lemma.** Let a graph X be partitioned into independent modules B_i of odd cardinality, with every pair of blocks either completely joined or completely nonadjacent. Let Q be its quotient graph. If g is an automorphism preserving this block partition and inducing sigma on Q, then

    sign_X(g) = sign_Q(sigma).

**Proof.** Orient each interblock edge according to a chosen total ordering of the blocks. There are no intrablock edges. An edge between quotient vertices i,j reverses direction under sigma exactly when all |B_i||B_j| edges above it reverse direction under g. That multiplicity is odd. Summing reversal counts modulo two proves the identity. Internal permutations within a block contribute nothing. Since the edge-reversal character on automorphisms is independent of the reference orientation, this computes the intrinsic sign. ∎

For disjoint marked independent twin triples, take the marked triples and all unmarked singleton vertices as the partition. Every block has size 3 or 1, so the lemma applies. Consequently the contracted sign in THM-4723 is exactly the original graph sign restricted to the stabilizer of the marked packing.

This does not validate the rejected fixed-point-count correction formula in the incoming script's old comment. It supplies a different direct proof. Contracted sign cannot be replaced by an unrelated permutation sign, but in this odd-block setting it equals the original edge-reversal character for the stated reason.

## 3. Compress all triple packings to a weighted false-twin quotient

Let C_1,...,C_s be the maximal classes of vertices with identical open neighborhoods. They are independent false-twin classes, of sizes c_i. The quotient Q has one vertex per class and an edge for each complete join. Let

    L = Aut(Q,c),          k_i = floor(c_i/3).

The full graph automorphism group is the semidirect product of the within-class symmetric groups with L. Indeed every automorphism preserves equality of open neighborhoods; conversely, every weight-preserving quotient automorphism lifts.

A packing of disjoint twin triples is, up to the within-class symmetric groups, completely described by

    a=(a_1,...,a_s),       0<=a_i<=k_i,

where a_i is the number of marked triples in class i. Two packings are Aut(X)-equivalent exactly when these vectors are L-equivalent. The stabilizer projects onto L_a. The internal kernel has trivial graph sign. After contraction, class i has c_i-2a_i vertices, of the same parity as c_i.

The graph sign of a quotient automorphism sigma is

    chi(sigma) = (-1)^(sum_(i<j, ij in E(Q), sigma(i)>sigma(j)) c_i c_j).

It depends only on the parity-weighted quotient and is a character of L. By the preceding lemma, the packing contributes to b precisely when chi is trivial on L_a.

Let E_0(X)=1 if X is sign-even, and 0 otherwise. For a cycle C of sigma on the quotient vertices, its class size and hence k_C are constant around C. Define

    F(k,l) = sum_(a=0)^k (-1)^(l a)
           = k+1                     if l is even,
             1                       if l is odd and k is even,
             0                       if l is odd and k is odd.

**Theorem.** For every finite graph X,

    b(X) = E_0(X) - D(X),
    D(X) = (1/|L|) sum_(sigma in L) chi(sigma)
                         product_(cycles C of sigma) F(k_C,|C|).

If X has a twin triple, then b(X)-e(X)=-D(X). If it has none, b(X)=e(X)=0.

**Proof.** Include also the empty packing a=0. The alternating, sign-trivial orbit sum over all a is D. To verify this, signed Burnside gives for each orbit the average of chi over its stabilizer, which is 1 when that restriction is trivial and 0 otherwise. The weight (-1)^(sum a_i) is constant on each L-orbit. Interchange the finite sums. A vector fixed by sigma is constant on each quotient cycle C; its contribution factor is F(k_C,|C|), proving the displayed formula.

The empty packing contributes E_0. Removing it and reversing the alternating sign gives b=E_0-D. If some k_i>0, e=E_0. If all k_i=0, every product F is 1 and D=E_0, proving the final boundary. ∎

The rational-looking average D is an integer because it is the signed orbit sum just derived. The formula replaces enumeration of every triple packing and the full factorial automorphism kernel by a weighted quotient calculation. It is a cycle-index evaluation with one parity-dependent factor per quotient cycle.

## 4. The first explicit failure of unit-weight split/join cancellation

Take X=K_(6,6). It is sign-even: within-part automorphisms preserve the orientation across the bipartition and a part swap reverses 36 edges. Its quotient is a single edge with two weights 6, so L=C2 and both k_i=2. The identity contributes 1 and the part swap contributes 3. Thus

    D(X)=(1+3)/2=2,
    b(X)=-1,       e(X)=1,       b(X)-e(X)=-2.

This can also be checked without Burnside. The nonempty packing-vector orbits are

    (0,1), (0,2), (1,1), (1,2), (2,2),

with alternating contributions +1,-1,-1,+1,-1, totaling -1. Every contracted stabilizer sign is trivial because both contracted part sizes remain even.

Splitting either six-vertex part into two triples and joining them gives Y=K_(3,3,6). It is sign-odd: swapping its two three-vertex parts reverses nine edges. In its weighted quotient, L=C2 exchanges the two weight-three vertices, chi takes values +1,-1, and k=(1,1,2). The identity contribution vanishes, while the swap contributes -2. Therefore

    D(Y)=-1,
    b(Y)=1,        e(Y)=0,        b(Y)-e(Y)=1.

Hence this proposed split/join pair has weights -2 and +1, not opposite weights. Even if the operation is made into an involution on these two underlying graph types, it does not cancel the integer defect sum. The initial unit-weight census cannot justify an unweighted cancellation argument at all orders.

## 5. An unbounded family of multiplicity obstructions

For every integer r>=1 let X_r=K_(6r,6r). Then k_1=k_2=2r, chi is trivial, and

    D(X_r)=(1+(2r+1))/2=r+1,
    b(X_r)=-r,      e(X_r)=1,      b(X_r)-e(X_r)=-(r+1).

Now perform the proposed split/join on one part. For r=1 the result is K_(3,3,6). For r>1 its weighted quotient has vertices of sizes 3,3,6r-6,6r, with edges between the first two, from each of the first two to the last, and from the third to the last. Its only nontrivial weighted quotient automorphism swaps the first two vertices, of odd graph sign. The identity factor vanishes, and the swap factor is -2: the other two k-values are even. Thus every resulting Y_r has

    D(Y_r)=-1,      b(Y_r)-e(Y_r)=1.

The mismatch grows without bound. A structural cancellation must carry multiplicity information or use a different grouping. This is analogous to the separately proved branch-integrality boundary, where the 27 unadorned even-Aut types at n=8 cannot be paired despite an integral weighted total. The carriers and weights are different; neither result proves that every possible decorated pairing is impossible.

## 6. A closed formula for balanced complete bipartite graphs

For X=K_(c,c), c>=1, let k=floor(c/3) and s_k=1 if k is even and 0 otherwise. Then

    D(X) = (s_k + (-1)^c (k+1))/2,
    b(X) = [c even] - D(X).

For c<3 this correctly gives b=0. For c>=3 the defect is -D. In particular K_(9,9) is sign-odd with b=2, showing that positive defect weights need not be one either.

This formula is a useful exact control of any proposed decorated split/join construction: it must account for the full multiplicities, not only the sign of each graph type.

## 7. Scope and next target

The new local results are the odd-module sign identity, the weighted quotient formula, and the explicit unbounded failure of the naive unit-weight cancellation. The global equality sum(b-e)=0 remains guaranteed by the inherited p-group theorem. A successful combinatorial proof would have to pair marked orbit data, account for stabilizer weights, or identify a more general cancellation complex. Neither the branch-integrality proof nor the mere existence of an unlabelled matching supplies that construction.

The incoming THM-4723's reduction is credited as the source of b. Its finite data through n=8 are compatible with these formulas. The stronger unweighted all-order pairing suggestion requires the refinement demonstrated above. No unrelated incoming result or historical certificate is retracted here.

## 8. Reproduction and fresh independent audit

All programs use Python's standard library and explicit exception checks that
remain active under `python3 -O`.

    python3 verify_twin_quotient.py > verify.out
    cp verification.json verification.normal.json
    python3 -O verify_twin_quotient.py > verify.optimized.out
    cmp verification.json verification.normal.json
    cmp verify.out verify.optimized.out
    python3 independent_packing_audit.py > independent_packing_audit.json
    python3 -O independent_packing_audit.py > independent_packing_audit.optimized.json
    cmp independent_packing_audit.json independent_packing_audit.optimized.json
    python3 verify_sharp_defects.py > sharp_defects.json
    python3 -O verify_sharp_defects.py > sharp_defects.optimized.json
    cmp sharp_defects.json sharp_defects.optimized.json

The original verifier checks 8,468 weighted false-twin-free quotients, 142
explicit odd-module contractions, balanced complete bipartite sizes 1 through
60, and split/join parameters 1 through 20, with 25,986 runtime checks.
It compares the cycle-index average with direct packing-vector orbit sums.

The fresh independent verifier imports no original code. It exhausts every
graph type on one through six vertices using the full symmetric groups,
enumerates actual labelled disjoint twin-triple packings, tests their actual
automorphism stabilizers, forms their contracted graphs, and compares the
contracted signs and alternating orbit sums with a separately implemented
quotient formula. It also checks all 32 false-twin-free labelled four-vertex
quotients with weights (3,3,3,3), the only nontrivial finite subcase in the sharp
order-twelve classification. Its 994 runtime checks pass normally and with
optimization. This is a distinct algorithmic check of the packing carrier,
not merely another evaluation of the same quotient formula.

The sharp-boundary verifier exhausts 16,798 positive-weighted, false-twin-free
labelled quotients of orders one through four and total weight at most twelve,
then checks both unbounded odd-sign families for r=1 through 40. It performs
77,005 runtime checks. Its scope is explicitly smaller than all graphs through
twelve: the universal boundary is proved in the companion argument, which
reduces its only finite classification to the five four-vertex quotient types.

`SHARP_DEFECT_BOUNDARY.md` proves the all-graph bound through eleven vertices,
classifies its first failures at twelve, and specifies the exact remaining
decorated-cancellation target. That proof is not inferred from a graph census
through six. `SHA256SUMS.txt` pins the publication files.
