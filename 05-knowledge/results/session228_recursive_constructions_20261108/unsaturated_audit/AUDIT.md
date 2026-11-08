# Independent audit: unsaturated six-vertex module fibers

Status: PASS, including all-order conclusions below. Inputs are frozen in inputs/frozen.json. The verifier imports no producer code and uses explicit exceptions rather than assertions. Its independent local transfers are computed from subset endpoint dynamic programming, rather than the producer's permutation/cut enumeration. Normal and optimized receipts must agree.

## Definitions and literal evidence

Adjacency row bits denote outgoing arcs. Root out-neighbors in the module are masks 17 and 40. Singleton bit 1 means root -> singleton. Earlier blocks dominate later blocks. The construction is singleton1, B_e0, singleton1, ..., B_e(t-1), singleton0, together with the root. It has 7t+2 vertices and requires t>=1.

The frozen block has H=19, ordered two-path-cover count Q=256, alpha=(21,20,29,35,73,59). Both masks have alpha sum94. The block is strong. Independent literal graph subset DP verifies all six graphs for t=1,2, comparing Hamiltonian count, root-ending count E, rooted-start cover U, deletion count d, and strong connectivity. Binary transfers verify every word through t=10.

The four-vertex negative control has H=5,Q=32,alpha=(5,6,9,7); all16 subset sums are distinct. This is a negative control for that representative only, not a classification assertion about arbitrary larger blocks.

## Exact algebra and all-order formulas

State order is (pre-root last category, post-root begun), with last category empty=0,in=1,out=2; the second coordinate runs 0,1. Column vectors are propagated chronologically from the right. Let S_b be singleton transfer, M_b module transfer, F=M_1-M_0, A=S_1 M_0, v=e_0, seed=S_1 v, and ell=(1,1,1,1,0,0). Thus D0(t)=ell S_0 M_0 A^(t-1) seed.

The verifier checks exactly:

* (A-512I)(A-19I)^2 seed=0.
* F(A-512I)(A-19I)=0 as an entire matrix.
* K=span(e_2-e_4,e_3-e_5) is invariant under both M_b and S_b.
* F kills K and its image is contained in K (the latter is visible directly in the frozen matrix receipt).

The first identity and the first three baseline values prove, for all t>=1,

D0(t)=[19950*512^t+(13862t-2289)*19^t]/15979.

The second identity and first two defect values prove

F A^i seed=(0,0,a_i,b_i,-a_i,-b_i),
a_i=-2*19^i,
b_i=-(8400*512^i+474*19^i)/493.

The exact transfers act on K coordinates (a,b) as follows: S0 maps to (0,b), S1 to (0,a+b), and M0 maps to (0,8a+19b). M1 has the same restriction because F kills K. Therefore a last-position defect contributes b_i; a nonlast defect contributes 19^(t-1-i)(a_i+b_i). Every term containing two defects vanishes. Consequently

D(e)=D0(t)+sum_i e_i Delta_i,

Delta_i=-(8400*512^i*19^(t-1-i)+1460*19^(t-1))/493 for i<t-1,
Delta_(t-1)=-(8400*512^(t-1)+474*19^(t-1))/493.

These are all-order identities supported by annihilator certificates and initial data, not an extrapolation of finite tests.

## Sidecar and strongness

E=d=19^t. For E, every block must precede root and the final singleton points into root. For deletion, Hamiltonian paths through an ordered sum traverse every block in order, giving the product of H.

For U, in any root-started path plus complementary path, let j be the first nonempty block on the root-started path. This gives U=d+sum_j alpha_j product_(k<j) H_k product_(k>j) Q_k. Here alpha is1 on singleton1,0 on singleton0,94 on both module choices; Q is2 on a singleton and256 on a module. Summing the geometric series yields

U(t)=(700*512^t-207*19^t)/493.

This also proves the sidecar is independent of every bit. Strongness holds for all t because root reaches the first singleton, that singleton reaches all subsequent vertices, every vertex reaches the final singleton, and the final singleton reaches root. Strongness of individual modules is independently checked but is not needed for this global argument.

## Distinctness and gcd

Write w_i=-Delta_i>0. For consecutive nonterminal weights, w_(i+1)>20 w_i because after dividing the numerator difference by19^(t-1) one obtains 8400(512/19-20)(512/19)^i-19*1460>0, whose minimum i=0 is positive. For the final ratio, the same expression replaces -27740 by 474-20*1460=-28726, still positive at i=0. Hence every weight exceeds the sum of all preceding weights. All2^t subset sums, and therefore scalars, are distinct.

For t>=2, w0=20*19^(t-1). The last weight is2 modulo4 and nonzero modulo5 and19, since493 is invertible modulo all three moduli,8400 is nonzero modulo19, and474 is nonzero modulo5. Thus gcd(w0,wlast)=2. D0 is odd from its formula with odd denominator15979. The gcd of all scalar values is therefore1. For t=1 the values653 and635 have gcd1.

## General separated-defect boundary

For arbitrary intervening block transfers, the only needed universal property is preservation of K. On K every nonempty block transfer maps (a,b) to (0,c*a+H*b), where c is the number of Hamiltonian paths in that block whose first vertex is a root out-neighbor; this follows by cancellation of the two previous endpoint categories whenever the new pre-root path is nonempty. Thus F R F=0 for any product R of such intervening block transfers. This proves additive separation of multiple defect insertions in an otherwise fixed chronological product. It does not assert arbitrary blocks preserve the same sidecar or strongness, nor that all defect positions have nonzero weights.

No prime coverage, interval coverage, universal representation theorem, or asymptotically optimal tournament order is claimed. The proven family provides2^t distinct strong-tournament scalars in one physical sidecar fiber, with fiber gcd1.

## Audit of the stronger universal separator lemma

The producer's frozen PROOF.md (SHA256 4d36a631d5f40e1bc14a04eb9b136d36242583f8d436633c65b0481f9aa7c5dd) proves a broader statement than the specific K argument above. Our independent verifier now checks all eight elementary endpoint-category maps on the basis of V={z0=0,z2+z4=0,z1+z3+z5=0}. Every map preserves V; both singleton maps send V into W={support on1,3,5; sum zero}. Each of four basis defects satisfying equal H,Q,alpha maps all columns into V and kills W. Every physical block preserves W (indeed acts by H), by cancellation of nonempty pre-root paths and the H empty-pre-root-path contributions. Thus any product of intervening physical blocks containing a singleton at any position sends a right defect's image to W, where the left defect kills it. The verifier also confirms that empty separation can fail for the universal defect basis. This validates the stated arbitrary-intervening-block generalization, without claiming fixed sidecars for arbitrary inserted blocks.

The producer's canonical verify.py is independent of its exploratory scripts, uses no SymPy, and has SHA256 a466966300be92a32b347d5fe1dd095b2f42c1bbad385e04341a6496a772640c. Those exploratory scripts are not part of this audit or its dependencies. The independent audit additionally confirms internal B strongness, t=1 baseline653 and alternate635, and the worst t=2 final-weight comparison. Exact divisibility is checked in finite regression, and for all orders follows from equality with the integral transfer expressions proved by recurrence, rather than an assumption about the displayed rational formulas.

Final negative-control extension: the independent verifier also exhausts all64 labeled four-vertex tournaments, finds24 with H=5,Q=32, and confirms injective root-mask alpha subset sums for each. This is an exhaustive order-four check only.
