# Positive path-cut lift and sparse lower ear bands

Status: PROVED elementary consequences of the inherited all-order OCF;
FINITE-EXACT audits described below. No global spectrum completeness claim.

## Inheritance and scope

The closest mechanisms are THM-4097 (`01-canon/theorems/THM-4097-order-nine-strong-ear-spectrum-solid-interval-and-lane-extension.md`), its exact cut-plus-field carrier, and THM-4114 (`01-canon/theorems/THM-4114-ocf-mobius-positivity-tropical-defect-layer-and-opposite-ear-cut-curvature.md`), its OCF positivity and ear submodularity. THM-4118 (`01-canon/theorems/THM-4118-ear-response-lattice-and-stateful-unit-component-intervals.md`) proves that arithmetic saturation does not certify an interval. The canonical hostile is strong order-five code 8, with saturated gcd but missing ear values. THM-1900 (`01-canon/theorems/THM-1900-the-insertion-response-calculus.md`) contains the corrected all-order neutrality theorem; we recover that theorem by a different proof, not as a new result.

The new material is a canonical **positive integer network lift**, its invisible circulation coordinate, exact mincut optimization, a chain/compatible-split theorem for small responses, and a quantitative obstruction to an individual parent's solid image. These are not another mean/max recurrence or a restatement of the scalar response gcd.

Concept board: (1) full strong spectrum / cross-parent overlap; (2) ear quadratic / positive network lift; (3) low response band / compatible splits; (4) scalar quotient / invisible circulation; (5) singleton cuts / minimal hostile. No claimed LRC connection.

## 1. Canonical positive lift

Let T be any tournament on V, n>=2. Put H(empty)=1 and H=H(T). Adjoin x with x->i iff i in S and write F(S)=H(T+x_S). For distinct i,j let

    R_ij = sum_P H(T[V\V(P)]),

where P ranges over simple directed paths in T starting at i and ending at j with an **even number of vertices**, at least two. Set R_ii=0. Then R is a nonnegative integer matrix and

    F(S) = H + 2 c_R(S),
    c_R(S) = sum_(i in S,j notin S) R_ij.                 (1)

Proof: split the OCF collections into those avoiding x and those using x. At most one cycle of a vertex-disjoint collection can use x. Deleting x from that odd directed cycle gives precisely one even-cardinality directed path P, beginning in S and ending outside S. Conversely each such P closes through x to one odd directed cycle. The remaining collection is any odd-cycle collection on T-V(P), whose total weight is H(T-V(P)); the cycle through x contributes the factor 2. The cycle has a distinguished vertex x, so there is no rotation multiplicity. This proves (1).

Compare coefficients with THM-4097's w and h:

    w_ij = R_ij+R_ji,
    h_i = sum_j R_ij - sum_j R_ji,
    A_i = 2 sum_j R_ij.                                  (2)

There is one directed one-edge path on each unordered pair, giving

    w_ij >= H(T-{i,j}) >= 1.                              (3)

The omitted terms are exactly paths with at least four vertices. Thus equality in the first inequality holds iff there is no directed even-cardinality simple path of at least four vertices between i and j in either direction. In particular, pair positivity is strict, regardless of strongness.

Immediate consequences, with their scope separated:

* F(S)>=H (inherited monotonicity, now termwise).
* F(S)-H >= 2 sum_(i in S,j outside S,i->j) H(T-{i,j}).
* F(S)=H iff no arc goes from S to its complement: the forward implication follows from the one-edge contributions; for the converse no directed path can leave S. These are exactly the SCC suffixes, recovering corrected THM-1900.
* For strong T every nonconstant S has c_R(S)>=1, and its child is strong.

## 2. The circulation quotient

Two real matrices R,R' with zero diagonal have the same directed-cut function iff D=R'-R satisfies

    D_ij=-D_ji and sum_j D_ij=0 for every i.               (4)

Proof: singleton cuts recover row sums; singleton/pair cuts recover R_ij+R_ji. These are necessary. Conversely a skew-symmetric zero-divergence matrix sums to zero across every cut, since its total row sum over S is zero and internal skew terms cancel. Thus the additional coordinate in the positive lift is a circulation, invisible to every one-vertex ear response. The admissible nonnegative lifts form a bounded circulation polytope with fixed unordered capacities w and prescribed divergence h, containing the canonical integer lift; no claim is made that every admissible lift is realized by a tournament path kernel.

Minimal example: the two labelled cyclic orientations of C3 have R equal to their adjacency matrices. They differ by a nonzero directed triangle circulation but have the same full ear response: H=3, every nonconstant ear equals 5. A cut carrier cannot recover this orientation, even when it retains all labelled cuts.

## 3. Exact mincut and terminal chain geometry

Let D consist of all nonempty proper subsets and let

    lambda = min_(S in D) c_R(S),  m = H+2 lambda.         (5)

Then m is exactly the smallest one-vertex ear response; for strong T it is the smallest strong-ear response from this parent. Given the kernel R, lambda can be found with 2(n-1) ordinary directed minimum-cut computations: fix r and take the minimum of all (r,t) and (t,r) terminal mincuts. Every nonconstant cut separates r from at least one t, proving equality. This is polynomial only **after R is supplied**; computing R is not asserted polynomial-time in the tournament order.

For A,B subsets,

    F(A)+F(B)-F(A intersect B)-F(A union B)
      = 2 sum_(i in A\B,j in B\A) w_ij.                  (6)

Thus the ear response is strictly submodular whenever A and B are incomparable, with a pair-deletion lower bound on its uncrossing cost.

Fix distinct terminals s,t. Among all cuts with s in S and t outside S, the minimum-response cuts form a chain under inclusion, so there are at most n-1 such labelled cut states. Indeed intersection and union remain feasible. Two incomparable minimizing cuts would, by (6), give a sum strictly greater than twice the minimum, a contradiction.

This does not mean the unrestricted minimum is a singleton cut. The minimal strong hostile is order six, LSB-first lexicographic tournament code 83, adjacency masks [38,8,2,5,15,30]. It has H=17 and unique minimizing cut S=14={1,2,3}, with F(S)=27. All singleton and complement-singleton cuts give at least 29. The independent auditor has checked these values by literal permutation enumeration.

## 4. Compatible low cuts and a forced-hole theorem

Let

    rho = min_(i<j) H(T-{i,j}),                            (7)

which is a positive odd integer. Consider the low band of cut states

    L = {S in D : F(S)<=m+rho-1}.                         (8)

**Theorem.** The unordered splits {S,V\S}, for S in L, are pairwise compatible: no two have all four intersection corners nonempty. Consequently

    |L| <= 4n-6.                                         (9)

Proof: if all four corners of A,B are nonempty, their intersection and union belong to D. By (6) and (3),

    (F(A)-m)+(F(B)-m) >= 2rho,

but membership in L bounds the left side by 2rho-2. Thus the splits are compatible. Fix a root r and represent each split by its side not containing r. Distinct representatives form a laminar family of nonempty subsets of an (n-1)-element ground set, so there are at most 2(n-1)-1=2n-3 of them. Each split has at most two orientations in L, giving (9). Complementary cuts need NOT have equal responses; the count only uses an upper bound of two orientations.

For n>=4 the entire numeric band [m,m+rho-1] lies inside the response hull. Choose distinct a,b,c and leave a fourth vertex outside {a,b,c}. For A={a,b},B={b,c}, all four sets in (6) are nonconstant. If M=max_D F, then

    2(M-m) >= 2w_ac >= 2rho.

Since M-m is even and rho is odd, M-m>=rho+1. Therefore the bottom band contains at least

    max(0,(rho+1)/2-(4n-6))                               (10)

**missing odd integers strictly inside the ear-response hull**. In particular rho>=8n-11 certifies that this single parent has a nonsolid ear image. This theorem does not constrain the union across parents: other parents may fill every one of these holes, as the existing strong spectrum intervals require.

The same argument for fixed s,t gives a chain for cuts within rho-1 of the terminal minimum and hence at most n-1 states in that terminal low band. The global hole statement (10) does not require a terminal range assertion.

## 5. Symmetry compresses the bottom band to a divisibility chain

Suppose T is vertex-transitive, n>=3, and let G be a transitive subgroup of its automorphism group. Its order n is odd, since all outdegrees are equal to (n-1)/2. The canonical R, w and h are G-invariant. Thus h is constant and its zero sum forces h=0. In particular

    F(S)=F(V\S).                                         (11)

Represent each low-band split by its smaller side A; n is odd, so this side is unique. The resulting family is G-invariant and laminar. Indeed two overlapping incomparable small sides have proper union, as both sizes are at most (n-1)/2, and would contradict split compatibility.

Every such A is a block of imprimitivity for G. Its equal-sized translates must be disjoint or identical: otherwise their intersection, both differences, and the complement of their union would all be nonempty, contradicting compatible splits. The translates cover V by transitivity, hence form a partition and |A| divides n.

For a fixed size d there is just one orbit of low small sides: the orbit of any one partitions V, and any second small side of size d meets one member of that partition; laminarity and equal size force equality. Thus all small sides of the same size have the same response. Equation (11) gives that response to both orientations.

The partitions belonging to different occurring sizes refine one another. If d<e, a d-block and an e-block that meet are nested, so each e-block is a disjoint union of d-blocks. Consequently the occurring sizes satisfy

    d_1 | d_2 | ... | d_r | n,
    d_(i+1)/d_i >= 3,
    d_r <= n/3,                                         (12)

because all sizes and n are odd. Hence 3^r<=n.

**Vertex-transitive low-spectrum theorem.** The low band contains at most

    floor(log_3 n)                                       (13)

distinct response values. For n>=5 it therefore contains at least

    max(0, (rho+1)/2 - floor(log_3 n))                    (14)

missing odd integers inside the ear-response hull. The logarithm bound is to be computed by integer powers of three, not floating rounding. The theorem bounds values more sharply than states; each block partition can have many cut states but just one value.

If G acts primitively, every small side is a singleton. Thus all 2n singleton/complement-singleton cuts are precisely the low-band states, and the only low-band response is m. For n>=5 there are exactly

    (rho-1)/2                                           (15)

missing odd integers in the band. This applies in particular to every vertex-transitive tournament of prime order, since a block size divides the degree. No converse or claim that a composite-degree action is primitive is made.

For cyclic T9, (14) forces at least 51 of its 53 bottom-band odd levels to be missing; the exact audit finds 52. For cyclic T11, primitive symmetry gives exactly 1062 missing levels. These are one-parent holes, not holes in the global tournament spectrum.

### Divisor-cut compiler for circulant parents

Suppose T is circulant on the additive group Z/nZ. A block for the regular cyclic action is a coset of a subgroup. For completeness, translate a block A to contain 0. For each a in A, the translates A and A+a meet at a, so they coincide. Thus A is closed under translations by its own elements and, being finite, is a subgroup. Every minimum-response small side is such a block. Translation and complement equality show that

    min_(empty != S != V) F(S)
      = min_(d|n, d<n) F({0,n/d,2n/d,...,(d-1)n/d}).      (16)

Hence the exact minimum over 2^n-2 cuts requires only tau(n)-1 representative subgroup cuts after the circulant parent is supplied. This reduces the number of queried child counts; it does not make Hamiltonian-path counting polynomial-time. Every minimum cut is a subgroup coset or the complement of one, but not every subgroup cut need minimize.

For a primitive transitive parent, all singleton responses equal their cardinality-one average from THM-4123, and the primitive theorem makes this the global minimum:

    m = H + 2W/n = ((2n-1)H+F_1(T))/n,
    W = ((n-1)H+F_1(T))/2.                               (17)

This is valid for any primitive transitive action, not only a cyclic action. The right side is automatically an odd integer under these hypotheses. It does not hold for arbitrary strong parents, as code83 demonstrates.

Here F_1(T) counts permutations of the old vertices with exactly one backward adjacency. The normalization 2W=(n-1)H+F_1 is exactly THM-4111, equation (6), and THM-4123, equations (4)--(5), in `01-canon/theorems/THM-4111-uniform-ear-average-and-recursive-selected-bank-growth.md` and `01-canon/theorems/THM-4123-balanced-cardinality-ear-average-growth-and-johnson-layer-lattice.md`.

## 6. Connection contract and limits

Source: tournament plus OCF packing weights. Target: a positive directed capacity network. Map: sum even-cardinality paths weighted by complement H. Preserved: every oriented ear response, exact min response, all response cuts, and neutrality. Destroyed on passing from R to (w,h): circulation. Extra sidecar: a particular positive lift (canonical R retains path combinatorics). Cheapest tests: C3 circulation; code83 balanced-minimum hostile; literal child enumeration through order five.

Positive capacity does not imply interval coverage, and the new theorem explains a concrete opposite pressure near the minimum. The full H-spectrum conjecture (all positive odd values except 7 and 21), all-order cross-parent interval overlap, and strong witnesses beyond the inherited finite prefix remain OPEN. No new all-order completeness, global maximum, or general tournament realizability theorem is asserted.

## 7. Audits

Primary scripts and outputs are in this directory. Exact scan universe: all labelled tournaments through order five for response identities and direct child DP; all strong labelled tournaments through order six for lower-band geometry and singleton-minimum hostiles. Larger cyclic controls use the orientation i->i+d modulo n, 1<=d<=(n-1)/2. Independent auditor additionally uses literal permutation counting on hostile and cyclic small controls.

The complete strong singleton-minimum census is: zero failures among 2,24,544 strong labelled parents at n=3,4,5; 720 failures among 22,320 at n=6. This is a finite-exact boundary, not a classification beyond order six.

The primary exact referee reports PASS on 1,098 labelled parents through order five, 33,864 literal child-DP comparisons, all 22,320 strong labelled order-six parents, and cyclic controls. Normal and optimized streams byte-match. Its semantic SHA256 is

    bff3e227aef79c6562ca6d5487286048716de1f82dbaf62af99e25b9baa8fdec

The larger cyclic controls are:

| order | H(T) | rho | minimum ear | maximum ear | low states | low values | actual missing odd levels |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 7 | 175 | 9 | 409 | 659 | 14 | 1 | 4 |
| 9 | 3267 | 105 | 8053 | 15459 | 18 | 1 | 52 |
| 11 | 93027 | 2125 | 236205 | 528411 | 22 | 1 | 1062 |

All cyclic ears through order nine are directly checked by a separate child DP; at order eleven both extremal children are directly checked. The kernel itself is compiled by the even-path/OCF formula. An independent agent has accepted the proof, literally enumerated the code83 child, and independently checked the cyclic9/11 H and rho values; its own artifacts are outside this subdirectory.

Reproduction from this directory:

    python3 -B probe.py > probe.out
    python3 -B audit.py > audit.out
    python3 -B -O audit.py > audit_optimized.out
    cmp audit.out audit_optimized.out

The file `symmetry_audit.py` separately verifies the block/divisibility and primitive controls without importing the primary referee.
