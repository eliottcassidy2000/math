# Tournament attachment/deletion: exact Weyl algebra and refresh spectrum

## Status and scope

**Recovered internal inheritance (2026-10-29):** the repository already proves the source/deletion Weyl pair, integral source-free induced-count basis, complete deletion cokernels, and all finite-field Hasse ranks in `PADDING_FROBENIUS_PROOF.md`, §5 (2026-10-19; content blob `9f035d2324589187ee3fb4bc3ff1b3131eb4a917`, retrieved and read directly). Repository link: [current repository path](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/tournament_hasse_frobenius_20261019/PADDING_FROBENIUS_PROOF.md). These portions below are independently rechecked inherited results, not new discoveries. That prior note does not state the normalized uniform-refresh spectrum or its explicit cyclic-triangle eigenfunction.

The identities below are proved for every degree over Q. An independent exhaustive checker verifies them through degree 5 using its own permutation-orbit canonicalizer and exact SymPy rational matrices, with identical successful results under ordinary Python and `python -O`. This is an attachment/deletion operator on actual unlabeled tournaments, **not** the Hamilton-path-word metagraph and not a gauge change identifying different tournaments. No source vertex or Hamiltonian marking is imposed. Multiplicities are essential.

Write V_n for the rational vector space with basis the isomorphism classes of n-vertex tournaments, including V_0=Q. Let a_n=dim V_n. U_n attaches a new named vertex with each of its 2^n orientation masks and then forgets all labels; D_n sums the n vertex-deleted cards, including repetitions. Define E_n=2^{-n}U_n and D_0=0.

## Commutator and harmonic decomposition

Deleting the newly attached vertex contributes 2^n copies of the input. Deleting an old vertex instead commutes with attachment except for the edge to that deleted vertex, whose two possible orientations have disappeared. Hence

    D_(n+1) U_n = 2 U_(n-1) D_n + 2^n I,
    D_(n+1) E_n - E_(n-1) D_n = I.

The second term is zero at n=0. These maps are independent of the chosen labeled representative: changing representatives merely permutes the masks and deletion choices.

Use E^k,D^k for the correctly graded compositions. From [D,E]=I,

    D E^r = E^r D + r E^(r-1).

E is injective: if Ev=0 and m is maximal with D^m v nonzero, then D^(m+1)Ev=(m+1)D^m v, a contradiction. D is surjective, for the finite degree-dependent right inverse

    R(v) = sum_(k>=0) (-1)^k E^(k+1)D^k(v)/(k+1)!

satisfies DR(v)=v by telescoping. Therefore H_n=ker D_n has dimension h_n=a_n-a_(n-1), setting a_-1=0. The finite operator

    P_n = sum_(k=0)^n (-1)^k E^k D^k/k!

has DP_n=0 and restricts to the identity on H_n. It also satisfies P_n E_(n-1)=0. Thus V_n=H_n direct-sum E_(n-1)V_(n-1); iteration gives the canonical decomposition

    V_n = direct-sum_(j=0)^n E^(n-j) H_j.

“Canonical” describes the summands and projections, not a preferred basis inside H_j. Harmonics are signed class combinations, **not prime tournaments**. Indeed for n>0 the sum of the coordinates of D_n v is n times the sum of the coordinates of v, so every nonzero vector in H_n has both positive and negative coefficients. H_3 is spanned by [directed triangle]-[transitive triple].

## Uniform vertex-star refresh on isomorphism classes

For n>=1 define the column-stochastic transition matrix

    K_n = E_(n-1) D_n / n.

It first deletes a uniformly chosen vertex, then attaches a vertex with independent fair orientations. Equivalently, on any labeled representative choose a vertex uniformly and resample every incident edge independently; then forget labels. The independent checker compares these two constructions entry by entry.

For h in H_j, E D(E^(n-j)h)=(n-j)E^(n-j)h. Consequently K_n is diagonalizable and its full spectrum is

    eigenvalue 1-j/n with multiplicity a_j-a_(j-1), 0<=j<=n.

Zero multiplicities are omitted. Since h_0=1,h_1=h_2=0,h_3=1, the spectral gap is exactly 3/n for every n>=3. (At n=1,2 the state space is a singleton.) All eigenvalues are nonnegative; this also gives the absolute spectral gap. This gap is for the unlabeled quotient, not the labeled state space, whose additional nonsymmetric modes can be slower.

A concrete backward eigenfunction attaining the gap is

    F(T) = number of cyclic triples in T - binomial(n,3)/4.

A triple not meeting the refreshed vertex is unchanged. A triple meeting it becomes cyclic with probability 1/4, regardless of its surviving edge. Thus E[F(T_next)|T]=(1-3/n)F(T). This explicitly connects the first harmonic degree to the first tournament intransitivity statistic.

## Stationarity, automorphisms, and coherence

For classes T of size n and S of size n+1, write U(T,S) and D(S,T) for the integer multiplicities. Counting pointed labeled realizations gives

    U(T,S) |Aut(S)| = D(S,T) |Aut(T)|.

There is **no extra factor n+1** in this identity. A way to count is to choose a vertex of S whose deletion gives T, then an isomorphism of that deletion onto the fixed representative of T, and quotient by Aut(S). The action is free.

Uniform measure on labeled tournaments pushes forward to

    pi_n(T) = n! / (2^binomial(n,2) |Aut(T)|).

The pointed identity implies the adjacent-level coherence relation

    pi_n(T) E_n(S,T) = pi_(n+1)(S) D_(n+1)(T,S)/(n+1).

Hence pi_n is stationary for K_n and K_n is reversible with respect to pi_n. It is generally not uniform on isomorphism classes. Irreducibility also follows directly by successively refreshing all labeled vertices with target edge choices, an event of positive probability; the self-loop probability is positive.

## What this does and does not solve about A000568

This proves monotonicity a_n>=a_(n-1) and interprets every difference as a harmonic dimension. It does not give a scalar recurrence computing a_n without independent information about h_n. The spectrum packages the entire preceding class-count sequence; it does not remove the enumeration problem. No efficient class generator follows merely from this decomposition.

## A second Weyl pair with integer coefficients and a positive counting interpretation

Define delta_n[T] to delete the universal source if T has one, and to give zero otherwise. A universal source means a vertex beating every other vertex, so it is unique when present (including the single vertex in degree 1). Set delta_0=0. With the original unnormalized attachment U,

    delta_(n+1) U_n - U_(n-1) delta_n = I.

Indeed a newly attached universal source occurs for exactly one mask and its removal returns T. If an old vertex is the universal source after attachment, it was already the universal source of T and must beat the new vertex; all remaining incident orientations are unrestricted. Their contributions are exactly U_(n-1)delta_n[T]. These are disjoint exhaustive possibilities. This is the inherited source/deletion Weyl relation in an automorphism-weighted basis, not a new tournament Weyl algebra. If A_n is diagonal with entries |Aut(T)|, then A_(n+1)U_n A_n^(-1)=D_(n+1)^T and A_(n-1)delta_n A_n^(-1)=delta_n: adding a universal source leaves the automorphism group unchanged. Thus the present pair is rationally conjugate to the old function pair (h,d), where h=delta and d=D^T. The conjugacy need not be unimodular. This is a Weyl relation with integer matrices, independently checked through degree 4 using banks through degree 5.

Unlike ker D, ker delta is literally spanned by source-free tournament classes. Adding a universal source gives a bijection from classes in degree n-1 to classes with a universal source in degree n. Consequently the source-free class count is exactly a_n-a_(n-1), including the convention that the empty tournament is source-free. The checker independently obtains [1,0,0,1,2,8] through degree 5. This is a positive combinatorial interpretation of the harmonic dimensions; it does not identify the actual vectors in ker delta with those in ker D.

Iteratively stripping universal sources leaves a unique source-free core, possibly empty, preceded by a uniquely determined initial transitive source chain. This does not assert that the core is strongly connected or substitution-prime: a source-free tournament can have a nontrivial initial strongly connected component.

Let S denote deterministic addition of a universal source. Source stripping gives an integral basis decomposition by the S-lifts of source-free classes. It must not be confused with the rational Weyl decomposition using U-lifts. Even for this integer Weyl pair, U-lifts do not automatically form an integral lattice basis: U^2[empty]=2[transitive pair]. Factorial-index issues remain, and no integral Fock-space isomorphism is claimed.

## Integral and finite-field boundary: extra torsion already at degree 5

The preceding splitting is over Q, not an isomorphism of integral lattices. If epsilon denotes coordinate sum, then

    epsilon_(n-1) D_n = n epsilon_n.

Thus augmentation modulo n induces a surjection coker_Z(D_n) onto Z/n. Rational full rank makes this cokernel finite, but D_n is not integrally surjective for n>1. For example D_2=[2] and D_3=[3,3], so they vanish modulo 2 and 3 respectively. The factorials in P and the powers of 2 in E cannot simply be ignored in positive characteristic.

Exact Smith diagonals from the independently constructed integer deletion matrices are:

| n | Nonzero Smith diagonal | Integral cokernel |
|---|---|---|
| 1 | (1) | 0 |
| 2 | (2) | Z/2 |
| 3 | (3) | Z/3 |
| 4 | (1,4) | Z/4 |
| 5 | (1,1,1,10) | Z/10 |

Degree 5 therefore has extra 2-torsion beyond the augmentation quotient Z/5. There is a transparent parity certificate: the sum of cyclic-triangle counts over all four-vertex deletion cards of a five-vertex tournament is twice the original cyclic-triangle count. Modulo 2 this defines a nonzero left-kernel vector for D_5. In the checker's degree-four canonical-code order [0,2,4,5], that vector is [0,1,0,1]. It detects the two classes consisting of a directed triangle plus a universal source or sink; the other classes have zero or two cyclic triples. The code independently checks this vector by exhaustive binary left-kernel enumeration. These finite computations independently confirm the **already proved all-degree classification** in the inherited §5, rather than motivate a conjecture. Transposing its unimodular function-basis diagonalization gives

    coker_Z D_n = direct-sum_(r=0)^(n-1) (Z/(n-r)Z)^(h_r),
    rank_Fp D_n = sum_(r=0)^(n-1) h_r * indicator[p does not divide n-r].

The displayed cyclic factors need not already be Smith invariant factors: prime valuations must be sorted separately and recombined. For n=5 the inherited diagonal is (5,2,1,1), giving Smith (1,1,1,10). The checker compares each independently computed Smith form to this all-degree formula through n=5. More generally, unordered deletion of m vertices from n has diagonal binomial(n-r,m), repeated h_r for 0<=r<=n-m; this is the transpose of the inherited Hasse theorem, with no factorial division after reduction.

The local parity certificate is an instance of the general induced-pattern double count. For a fixed j-vertex tournament H, let c_H(T) count induced j-subsets of type H. Every occurrence survives exactly n-j vertex deletions, so

    sum_(v in V(T)) c_H(T-v) = (n-j)c_H(T).

If p divides n-j, the card-feature row c_H modulo p annihilates D_n. Only nonzero, linearly independent such rows give rank-loss or torsion constraints; independence of arbitrary pattern features is not asserted. Degree 5 and the cyclic triangle give the explicit extra factor 2, consistent with the complete inherited classification.

## Exact trace identities, not scalar enumeration closure

The spectrum telescopes to

    n tr(K_n) = sum_(j=0)^(n-1) a_j,
    a_n = (n+1) tr(K_(n+1)) - n tr(K_n).

The first formula is independently checked through degree 5. The second follows by subtraction for n>=1 (or n=0 with the conventional zero product 0 tr(K_0)=0). These are exact identities, but computing the new transition trace still requires structural information; this is not a closed scalar recurrence for A000568.

## Classical attribution

The algebraic mechanism is standard Weyl/Heisenberg commutation, with the differential-poset precedent [Richard P. Stanley, *Differential posets*, JAMS 1 (1988), 919–961](https://math.mit.edu/~rstan/pubs/pubfiles/77.pdf). Its extension to paired graded-graph operators appears in [Sergey Fomin, *Duality of Graded Graphs*, J. Algebraic Combinatorics 3 (1994), 357–404](https://www.kurims.kyoto-u.ac.jp/EMIS/journals/JACO/Volume3_4/bx4c3m5xlcdmtgu3.fulltext.pdf). The normalized tournament operators have rational weighted edges and are not asserted to be an ordinary unweighted differential poset. The tournament construction and consequences are verified here; no literature-priority claim is made.

## Reproduction

From this packet directory:

    python3 updown/check_updown.py
    python3 -O updown/check_updown.py

Requires SymPy (tested with 1.14.0). All checks use explicit exceptions rather than Python assert, so optimization cannot erase them. Enumerated counts: 1,1,1,2,4,12; harmonic dimensions: 1,0,0,1,2,8. Checks include the marked automorphism identity, full rank, normalized commutators through n=4 (requiring level 5), exact projectors, direct sums, all transition entries against independent labeled refresh, stationary measures, detailed balance, cyclic-triangle eigenfunction, exact spectra through n=5, Smith forms through n=5, the degree-five parity obstruction, the trace identity, integer source-deletion commutators, and source-free class counts.

## Compressed Smith atlas without enumerating classes

`compressed_smith.py` realizes the inherited all-degree diagonal classification using only the A000568 count prefix from the vendored Session87 `degree_bucket_dp` implementation. It does **not** generate tournaments and does not expand a diagonal of length a_(n-1).

For each prime p<=n, accumulate the histogram of v_p(n-r), weighted by h_r, for 0<=r<n. Sorting these valuation buckets gives a run-length encoding of the p-adic valuations of the Smith invariant factors. Merge the cumulative-position breakpoints of all prime runs. On each resulting interval, all prime valuations are constant; multiplying p raised to these valuations gives its invariant factor, and the interval length is its multiplicity. Adjacent identical factors are merged.

Correctness follows from the inherited diagonal-to-Smith rule: sort the valuations independently for each prime and recombine at matching positions. The compressed representation computes exactly that operation. A prime p contributes at most floor(log_p n) internal breakpoints, since its valuation lies between 0 and floor(log_p n). Consequently the output has at most

    1 + sum_(p<=n, p prime) floor(log_p n)

runs, even when its rank has thousands of digits. All ranks, multiplicities, breakpoint positions, and factors are exact integers. This is a bound on the number of runs, not a claim that arbitrary-precision arithmetic or the count-prefix computation has constant cost. The algorithm is an implementation/compression of the recovered theorem, not a new torsion classification.

The selected atlas is:

| n | Decimal digits of rank a_(n-1) | Actual factor runs | Proven run bound |
|---|---:|---:|---:|
| 5 | 1 | 2 | 5 |
| 10 | 6 | 6 | 8 |
| 25 | 60 | 14 | 15 |
| 50 | 292 | 23 | 24 |
| 100 | 1305 | 36 | 36 |
| 200 | 5558 | 60 | 61 |

`smith_atlas_200.json` contains every exact invariant factor and its multiplicity, the full ranks, and the exact prime factorization of each finite cokernel's order. It also pins the count-prefix contents by SHA-256. Prime exponents encode the group order compactly; constructing that astronomical integer is unnecessary.

Verification uses independent direct expansion and sorting of valuation arrays through n=10 (largest expanded rank 191536), the separately enumerated tournament deletion-matrix Smith forms through n=5, and rank conservation, invariant-factor divisibility, group-order valuations, and the proven interval bound for every selected large atlas entry. The ordinary and optimized runs produced byte-identical JSON:

    python3 updown/compressed_smith.py
    python3 -O updown/compressed_smith.py
