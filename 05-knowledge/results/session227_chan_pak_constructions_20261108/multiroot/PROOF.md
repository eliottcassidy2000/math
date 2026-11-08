# Set-root grafts: positive cover expansion, strong realizations, and a finite-menu boundary

**Status: PROVED, independently audited PASS.** Exact finite controls are in `verify.py` and `CHECK.json`. Chan–Pak completeness remains **OPEN**. This packet claims no literature priority and no improved unrestricted scalar coverage interval.

## Inheritance and the question

The inherited single-root operation is Session225's graft: retain nonempty tournaments A and B, put B -> a and A-a -> B. Session226 gives its strongness iff and four-coordinate context equivalence. The broader path-cover substitution language is inherited from **THM-3121-path-cover-walk-content-substitution-kernel.md**; multiplicativity of the path-colour polynomial under ordered sums is **THM-3166-falling-factorial-order-join-path-colour-transform.md**. Neither general transform is claimed as new here.

The question is whether varying and reordering blocks inside an attached gadget escapes the arithmetic thinness of fixed words. The new scoped conclusions are the sharp positive set-root expansion, its exact connectivity guard, the resulting finite-menu boundary (including arbitrary outer ordered sums), and an explicit order-faithful strong collision family. The general expansion is an elementary specialization of the inherited cover calculus, not an independent replacement for it.

The concept board is: (1) actual graph construction versus scalar recurrences; (2) interface size versus cover degree; (3) ordered component words versus their commutative colour statistics; (4) strong-component recoverability; (5) additive logarithmic cost after outer ordered sums. The hostile cases are empty/all root sets, singleton sources, and permutations of nonisomorphic strong gadget blocks.

## 1. The operation and exact positive coefficients

Let A and B be nonempty finite tournaments and R be any subset of V(A), including the empty and full subsets. Put C=V(A) minus R. Define G(A,R;B) by retaining both induced tournaments and orienting every cross-edge as

    B -> R,       C -> B.

For a nonempty tournament B let F_B(l) be its number of ordered l-tuples of **nonempty** vertex-disjoint directed paths covering B. A path is its directed vertex sequence. Thus F_B(1)=H(B), and F_B(l)=0 for l>|B|.

For each permutation pi=(a_1,...,a_n) of V(A), regard the n+1 slots before, between, and after its vertices as potential nonempty B-runs. A slot is eligible exactly as follows:

- the initial slot is eligible if a_1 belongs to R;
- the final slot is eligible if a_n belongs to C;
- an internal slot (a_i,a_(i+1)) is eligible if a_i belongs to C and a_(i+1) belongs to R.

Every backward A-adjacency, meaning a_(i+1)->a_i, must be separated by a B-run. If any such slot is ineligible, discard pi. Otherwise let b(pi) be its number of backward adjacencies and o(pi) the number of eligible optional slots, including eligible boundary slots and eligible forward internal slots. Define, for l>=1,

    c_l(A,R) = sum over compatible pi of binom(o(pi), l-b(pi)),

where binomial coefficients outside 0<=l-b<=o are zero.

**Theorem 1 (positive expansion and sharp degree).**

    H(G(A,R;B)) = sum_(l=1)^L c_l(A,R) F_B(l),
    L = min(|R|, |C|)+1.

The coefficients are nonnegative integers independent of B. If R is empty or all of A, the formula reduces to H(A)H(B).

**Proof.** Remove B from a Hamiltonian path of G. The resulting order of the A-vertices is a permutation pi. Its nonempty B-runs occupy distinct slots. Each occupied internal slot requires C->B->R; the two boundary conditions are exactly those above. An unoccupied internal slot must already be an arc of A. Thus all backward slots must be occupied, and any remaining eligible slots may be chosen independently. If exactly l are occupied, assigning an ordered spanning l-path cover of B reconstructs the Hamiltonian path uniquely. Conversely every counted choice gives a path. This proves the formula with all l.

For the degree bound, there are l runs of B and at least l-1 transitions B->A and at least l-1 transitions A->B. The first kind enters distinct vertices of R, and the second leaves distinct vertices of C. Hence l-1<=min(|R|,|C|). This proves the stated truncation. If R is empty or full, only one boundary run is possible and the graph is an ordered sum. ∎

**Sharpness at every interface size.** Let A be the transitive tournament whose ordered vertex classes are R followed by C, with r=|R|>=1 and c=|C|>=1. Let m=min(r,c), M=max(r,c). Then

    c_(m+1)(A,R) >= (m!)^2 S(M,m) > 0,

where S(M,m) is the number of partitions of an M-element set into m nonempty unlabelled blocks. Count first those paths with m+1 B-runs that begin and end in B; the intervening A-runs all start in R and end in C. If r<=c, the r roots occur once each as the ordered starts of r runs, in r! orders; the c C-vertices are distributed surjectively over the runs, in r!S(c,r) ways, and occur increasingly inside each run. The case c<=r is dual. When r=c=m, maximal run count forces both boundary B-runs, so equality holds and the coefficient is (m!)^2. Unequal shores can have additional maximal-run paths lacking one boundary B-run; only the lower bound is asserted there. Thus the degree cannot be reduced uniformly.

This is a positive combinatorial expansion; signed coefficients introduced below are only a change of basis.

## 2. Exact strongness and the exponential order bound

**Theorem 2 (strongness iff).** G(A,R;B) is strongly connected if and only if the first strongly connected component of A meets R and the last strongly connected component of A meets C. This is independent of the internal arcs of the nonempty gadget B.

Equivalently:

1. every vertex of A can be reached within A from some vertex of R;
2. every vertex of A can reach some vertex of C within A.

**Proof.** If the graft is strong, a path from B into A first enters R; after every further visit to B it re-enters R. Thus every A-vertex must already be reachable in A from R. Similarly, a path from A to B first leaves A through C, so every A-vertex must reach C within A. Since tournament SCCs are totally ordered, these two conditions are exactly the first/last-SCC conditions.

Conversely, each A-vertex reaches C and then every vertex of B. Each B-vertex reaches R and then every vertex of A. Any B-vertex therefore reaches any other via A. This proves all ordered-pair reachability, without a connectivity assumption on B. ∎

For R={a}, the second condition is automatic once |A|>=2 and a reaches all of A. This recovers exactly Session226's criterion. Empty/full R and a singleton A never satisfy the condition.

**Lemma (a positive two-run coefficient).** Under the strongness condition, c_2(A,R)>=1.

**Proof.** A has a Hamiltonian path starting in R and ending in C. If A has one SCC, it is a nontrivial strong tournament, and a Hamiltonian cycle (Camion's inherited theorem) has an arc from C to R; break that arc. If A has at least two SCCs, choose a Hamiltonian path in the first SCC starting at a vertex of R and a Hamiltonian path in the last SCC ending at a vertex of C, and concatenate them with Hamiltonian paths in all middle SCCs. These endpoint paths exist by breaking Hamiltonian cycles in nontrivial strong components; singleton components are immediate. The path obtained supplies one term of c_2 by occupying the initial and final B-slots. ∎

For any nonempty B of order N,

    F_B(2) = sum_(empty != S != V(B)) H(B[S]) H(B-S) >= 2^N-2,

using existence of a Hamiltonian path in every tournament. Moreover c_1=H(G(A,R;K1))>=3, because the singleton gadget gives a strong tournament. Consequently every strong-compatible source satisfies the sharp all-order bound

    H(G(A,R;B)) >= 2^|B|+1.                         (1)

It holds even when B is not strong. Equality holds for A a directed edge, R its source, and B transitive, at every gadget order. The source guard is essential: for R empty and transitive A,B, H=1 at arbitrarily large gadget order.

## 3. Scalar diagonalization and explicit strong families

Define Q_B(j) to count ordered j-bin spanning path covers of B, with bins allowed to be empty. There is one empty path per empty bin. Selecting the nonempty bins gives

    Q_B(j) = sum_l binom(j,l) F_B(l),
    F_B(l) = sum_(j=0)^l (-1)^(l-j) binom(l,j) Q_B(j).

Since B is nonempty, Q_B(0)=0. Theorem 1 becomes

    H(G(A,R;B)) = sum_(j=1)^L lambda_j(A,R) Q_B(j),
    lambda_j = sum_(l=j)^L (-1)^(l-j) binom(l,j)c_l.     (2)

The all-order ordered-sum identity

    Q_(B1=>...=>Bt)(j) = product_i Q_Bi(j)               (3)

is inherited from THM-3166: intersect each colour's path with each component and concatenate in component order. It is restated solely to apply it to (2).

For B=TT_N, Q_B(j)=j^N. Thus every strong-compatible (A,R) gives an explicitly constructed infinite strong family with

    H(G(A,R;TT_N)) = sum_(j=1)^L lambda_j j^N, N>=1.    (4)

The recurrence with characteristic product_(j=1)^L (x-j) follows from this **physical graph identity**, not from assuming a scalar recurrence is realizable. For A=TT_(2m), R its first m vertices, the coefficient of (m+1)^N is exactly (m!)^2, so every integer exponential base >=2 occurs as the dominant base in this construction.

**An explicit three-mode strong family.** Take A=TT_4 on labels 0,1,2,3 and R={0,2}. Slot enumeration gives (c_1,c_2,c_3)=(9,10,2), hence (lambda_1,lambda_2,lambda_3)=(-5,4,2). For B an ordered sum of m>=1 directed triangles, the inherited polynomial Q_C3(t)=t^3+2t gives

    H(G(A,{0,2};C3=>...=>C3)) = -5*3^m+4*12^m+2*33^m.

These are strong tournaments of order 3m+4, starting with counts 99,2709,78651. The graph construction and Theorem 1, not the signs or an extrapolated recurrence, establish realizability for every m. This example was supplied independently by the audit and is verified literally through m=3.

More generally, any word of nonempty gadget blocks B_i produces the physical count (2)-(3). This count is invariant under every permutation of the block word. The internal rooted structure of the resulting graph need not be invariant.

### 3.1 Exact minimal rational linear gadget summary

Fix r>=1 and allow all source pairs with min(|R|,|A-R|)<=r. Equations (2)-(3) show that the r+1 quantities

    (Q_B(1),...,Q_B(r+1))

are sufficient to compute every such graft's Hamilton count by a rational linear functional. They are also necessary in dimension among summaries whose every context response is a rational linear functional. This is about the **unrooted gadget side** of a one-layer scalar observation, not Session226's four-state rooted-source observer.

For a physical rank witness, first take the source edge with its source in R. Its row in (Q1,Q2) is (1,1). Next take TT_3 with R its source. Its row is (-1,3); these two rows have determinant 4. For m=2,...,r add the source TT_(2m) with R its first m vertices. Its coefficients vanish above Q_(m+1), and the coefficient of Q_(m+1) is exactly (m!)^2 by equal-shore sharpness. Thus these r+1 context rows are linearly independent. Every one is strong-compatible.

Evaluate these contexts on the physical gadgets B=TT_n, n=1,...,r+1. The gadget-coordinate matrix is (j^n), with 1<=j,n<=r+1, an invertible Vandermonde matrix after factoring one j from each row. The context-response matrix therefore has rank r+1. Any model assigning each B a vector in Q^d with every allowed context response linear factors this matrix through Q^d, so d>=r+1. This proves exact minimality. The argument does not exclude arbitrary nonlinear encodings or give a physical inverse for an arbitrary moment vector.

## 4. A finite-menu boundary, even after arbitrary outer ordered sums

Fix a finite bank of pairs (A,R), all satisfying Theorem 2. Fix a menu of k>=1 nonempty tournaments B_1,...,B_k. A first-stage construction chooses any nonempty finite word in this menu as its ordered-sum gadget and grafts it onto one source from the bank.

**Theorem 3.** The number of distinct first-stage Hamilton counts at most X is O((log(X+2))^k). Even if one takes arbitrary finite ordered sums of first-stage outputs, the number of distinct resulting counts at most X is

    exp(O((log X)^(k/(k+1)))) = X^(o(1)).               (5)

For the first-stage bound take X>=3; below 3 there are no first-stage values. The constants depend on the fixed source bank and gadget menu. In particular, neither this language nor its arbitrary outer ordered-sum closure realizes all sufficiently large odd integers.

**First-stage proof.** If a gadget word has t letters, |B|>=t. By (1), an output at most X has t<=log_2(X-1). By (2)-(3), the count depends on the word only through its k multiplicities. There are at most binom(floor(log_2(X-1))+k,k) such nonnegative multiplicity vectors, including the empty vector as an upper bound. Multiply by the fixed number of sources. This proves the first assertion.

**Outer-sum proof.** Let S be the set of distinct first-stage counts. Every first-stage graph is strong and has at least three vertices, so its count is at least 3. Let N(t)=#{g in S: log g<=t}. The preceding result gives N(t)<=C t^k for all t>=0 after adjusting C (N is zero below log3).

Every outer ordered sum has count a product of elements of S. Different orders give the same product, so the number of distinct products is at most the number of multisets of elements of S with the corresponding logarithmic cost. For s>0 the multiset generating product has

    Z(s) = product_(g in S) (1-g^(-s))^(-1),
    log Z(s) = s integral_0^infinity N(t)/(exp(st)-1) dt
             <= C s^(-k) integral_0^infinity u^k/(exp(u)-1) du
             = C k! (sum_(r>=1) r^(-k-1)) s^(-k).

All exchanges involve nonnegative terms. The integral is finite for k>=1, establishing convergence as well as the bound. The number of multisets of cost at most log X is at most X^s Z(s), hence at most exp(s log X+C' s^(-k)). Taking s=(log X)^(-1/(k+1)) proves (5). The empty ordered sum may add the value 1 and does not change the bound. ∎

**Scope.** This rules out one finite-menu construction architecture, with arbitrary word lengths and arbitrary outer SCC assembly. It does not bound recursive grafts whose newly produced graphs become new sources or gadget letters, an infinite gadget alphabet, freely changing root sets/sources, or general tournaments. In particular it does not strengthen the Session226 finite-line no-go into an unrestricted obstruction. It identifies the precise information loss here: the full block order disappears from the scalar because the native operation is the commutative product (3).

## 5. An explicit family of rooted-graft-indistinguishable strong outputs

There is also a constructive consequence of the same loss. Let

    C3 = [2,4,1],      F4 = [8,1,3,6],

where rows are outgoing-neighbour masks. Both are strong. Their (H,Q(2)) pairs are (3,12) and (5,32). These follow either by direct small enumeration or from Q_C3(t)=t^3+2t and Q_F4(t)=t^4+4t^2.

Fix k,l>=0. Let W range over all words with k letters C3 and l letters F4. Put B_W equal to the ordered sum of W followed by one fixed final C3. Graft B_W onto the source-rooted edge a->c, and mark vertex 0 of the final C3 as the new root. There are binom(k+l,k) such words.

Write h=3^k5^l and q=12^k32^l. Every resulting tournament has order 3k+4l+5, is strong, and has

    H = 3h+12q.                                        (6)

More strongly, in Session225's four-coordinate convention x=(D,E,U,d), every such rooted output has the same state

    x = (3h+12q, 3q, 6h+14q-1_(k+l=0), h+4q).        (7)

Here D=H, E counts Hamiltonian paths ending at the root, U counts ordered two-bin covers with the second nonempty and starting at the root, and d=H after root deletion.

**Proof of the state formula.** Let P be the ordered-sum prefix W. Its H and Q(2) are h and q, even for the empty prefix, where both are 1. For B=P=>C3 with its root in the final triangle, the gadget statistics of Session225 are

    H_B=3h, Q_B=12q, s_B=1_(P empty), e_B=h,
    alpha_B=3h, delta_B=3q, T_B=7q,
    H(B-r)=h, Q_(B-r)(2)=4q.

For example, in a rooted three-bin cover starting at r, the prefix cannot use that bin and has q ways to occupy the other two bins, while C3 contributes 7. The corresponding rooted two-bin formulas force the entire prefix into the ordinary bin (alpha), or allow both bins in the prefix (delta). These are direct componentwise path-bin bijections. Substitution in Session225's proved graft matrix, acting on the physical edge state (1,0,2,1), gives (7). The first coordinate is also the independent source-edge identity H(G)=H(B)+Q_B(2), from c_1=3,c_2=1 in Theorem 1.

**The graphs are pairwise nonisomorphic.** In every output, a is the unique vertex of outdegree 1 and c is the unique vertex of indegree 1. Indeed a only beats c, while c only loses to a. Every vertex of B belongs to a strong component of B of order at least three, hence has at least one internal in-neighbour and out-neighbour; after grafting it has both indegree and outdegree at least 2. Any isomorphism must therefore fix a and c and restrict to an isomorphism of the induced B. Tournament SCCs have a unique order, so B determines its word of C3/F4 components. Distinct W give distinct B and hence nonisomorphic outputs. This proof does not assume the mark is preserved; the outputs are distinct even as unrooted tournaments.

Thus (7) supplies binom(k+l,k) nonisomorphic strong physical realizers that **no common finite future rooted-graft context, using the fixed marking convention, can distinguish by its Hamilton count**, by the inherited Session226 context equivalence. This is a constructive collision theorem, not scalar coverage. Arbitrary rerooting, graph-dependent choice of different contexts, and other graph tests are outside the indistinguishability statement. At k=l=m, it gives binom(2m,m) nonisomorphic strong realizers of the same count 3*15^m+12*384^m at the same order 7m+5. The state quotient loses exponentially many graph types, even inside strong tournaments. This gives a concrete hostile to any proposed physical inverse rule that silently reconstructs a unique tournament from the four-state vector.

## 6. Reproduction and precise finite universe

Run in this packet's directory:

    python3 -I -B verify.py > replay.json
    cmp replay.json CHECK.json
    python3 -I -B -O verify.py > replay.optimized.json
    cmp replay.optimized.json CHECK.json

Only the Python standard library is used, and all checks are explicit exceptions rather than assertions. The exact universe contains all 75 labelled tournaments of orders 1 through 4; every subset R is considered on every source. All source/gadget pairs are checked except those with both orders 4, with no connectivity filter. Direct Hamilton subset DP checks the positive cover law, the spectral transform, and the strongness iff on the literal graft. Permutation-cut enumeration checks gadget covers; an independent subset convolution checks the empty-bin transform. All degree-sharpness cases 1<=|R|,|C|<=4 are checked, and the three-mode family is checked on literal graphs for m=1,2,3. Physical response matrices check the exact minimal-summary ranks through r=4. Separate literal checks verify the strong collision family for (k,l)=(0,0),(1,0),(0,1),(1,1),(2,1), including all corresponding words and their complete four-state vectors. The receipt includes the literal witnesses.

The all-order quantifiers come from the proofs above. The finite checks do not supply a Chan–Pak cutoff, an arbitrary-state realization result, or an unrestricted nonexistence theorem.
