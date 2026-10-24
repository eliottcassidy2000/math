# Matrix-valued distance constraints for the merged tiling metagraph

Status: PROVED, all n >= 3, for the identities and necessary constraints below; FINITE-EXACT for the n=3..6 audits and the n=5 hostile. These are new derived constraints in this session, not a claim of novelty for the underlying Walsh or cube linear algebra. No sufficiency or unique-arrangement claim is made.

## Inheritance and scope

The closest proved mechanisms are THM-781, `01-canon/theorems/THM-781-hamiltonian-path-inverse-metagraph-fibre.md` (the exact HP/Aut fibres), THM-346, `01-canon/theorems/THM-346-tiling-quotient-bucket-balance.md` (all-layer balance and non-equitability warning), and THM-796, `01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md` (commuting complement/reflection, coloured incidence, and recursive loss). THM-830 and THM-851 already describe defect-labelled pair-groupoid composition, so that construction is not claimed again here. THM-848's Walsh needle stalk uses the **full arc cube** and one observable H; this note uses **all fibre indicators in the fixed-path cube**, a different carrier.

The hostile inherited from HYP-3809, `05-knowledge/hypotheses/HYP-3809-merged-metagraph-bucket-conjecture-atlas.md`, is a count/colour-preserving rewiring. The new constraints couple all distances at once and retain matrix positivity and rank, not merely degrees. The independent signed-count input is THM-1936, `01-canon/theorems/THM-1936-signed-redei-join-multiplicative.md`; its label-sign warning is respected below.

## 1. Definitions and loop conventions

Put m=binom(n-1,2), N=2^m and X={0,1}^m. Fix the directed path 0 -> 1 -> ... -> n-1. Coordinates are the nonconsecutive pairs (i,j), i+2 <= j, in lexicographic order; bit 1 means i -> j. This is the ascending relabelling of the standard descending-path convention.

Let pi:X->U send a tiling to its tournament isomorphism class modulo converse, and write F_u=pi^{-1}(u), M_u=|F_u|, D=diag(M_u). For 0 <= d <= m define

    W_d(u,v) = #{(x,y) in F_u x F_v : dist(x,y)=d}.

Set W_d=0 when d<0 or d>m, so the local-flip formula also covers m=1. These are **ordered half-edge counts**. Consequently W_0=D; for d>0 an internal unordered line contributes 2 to W_d(u,u); an unordered line between distinct nodes contributes 1 to W_d(u,v) and 1 to W_d(v,u). W_m is precisely the complement-line half-edge matrix of THM-796, with doubled loop diagonal. Every W_d is symmetric, nonnegative and integral, and

    W_d 1 = binom(m,d) M,               sum_d W_d = M M^T.             (1)

Reflection sigma permutes coordinates by (i,j)->(n-1-j,n-1-i), preserves every F_u and does not invert bits. Put f=floor((n-1)/2), q=(m-f)/2. Its coordinate permutation has f singletons and q two-cycles.

## 2. The matrix Walsh moment theorem

For a coordinate subset S define chi_S(x)=(-1)^{sum_{e in S}x_e}, and

    Z_j(u,S)=sum_{x in F_u}chi_S(x),        |S|=j.

Let

    K_j(d)=sum_a (-1)^a binom(d,a)binom(m-d,j-a),
    Q_j=sum_{d=0}^m K_j(d) W_d.                                  (2)

Then, for every j,

    Q_j=Z_j Z_j^T >= 0.                                          (3)

In particular the Q_j are **integer Gram matrices**, not just entrywise inequalities. Also

    Q_0=M M^T,
    Q_j 1=0 for j>0,
    sum_j Q_j=N D,
    W_d=(1/N)sum_j K_d(j) Q_j.                                  (4)

If

    F_j=[t^j](1+t)^f(1+t^2)^q,
    r_j=(binom(m,j)+F_j)/2,

then

    rank Q_j <= r_j,
    rank Q_j <= |U|-1 for j>0.                                  (5)

The reflection refinement is essential: rank <= binom(m,j) is weaker.

Proof. For fixed x,y, the sum of chi_S(x)chi_S(y) over |S|=j is K_j(dist(x,y)), by choosing how many coordinates of S meet the disagreement set. Summing over fibres proves (3). Orthogonality of the N Walsh characters proves (4). Since sigma preserves F_u, Z_j(u,S)=Z_j(u,sigma S); thus columns are equal on reflection orbits. There are r_j orbits, by Burnside: a fixed j-subset chooses arbitrary fixed coordinates and either both or neither member of each transposed pair, giving F_j. This proves (5). The row-null condition gives the second rank bound. All proofs are over R or Q; no equitability assumption enters.

Equivalently, E_j=(1/N)D^{-1/2}Q_jD^{-1/2} are positive semidefinite matrices summing to I. They are compressions of the cube's orthogonal Walsh-level projectors, but are generally **not projectors** and generally do not commute. Treating them as an association scheme on merged nodes would discard precisely the loss measured in Section 5.

## 3. The top moment and exact arithmetic

Call the following signed even-minus-odd count the **tiling parity charge**:

    z_u=sum_{x in F_u}(-1)^{|x|}.

Since there is one weight-m character,

    Q_m=z z^T,
    sum_d(-1)^d W_d(u,v)=z_u z_v.                               (6)

Let tau be the transitive node. Its fibre consists of the all-one tiling, so M_tau=1 and z_tau=(-1)^m. Therefore, if

    a_u=sum_d(-1)^d W_d(tau,u),

then a_u=z_tau z_u and the entire top matrix is forced by the transitive row:

    sum_d(-1)^d W_d(u,v)=a_u a_v.                               (7)

This is the **alternating-distance balance**, a pairwise arrangement identity determined by alternating distance counts from one node. Because m>=1, sum_u z_u=0. The total charge on positive nodes equals the absolute total on negative nodes.

THM-781 says an ordinary class fibre has odd size H(T)/|Aut(T)|. Thus a self-converse merged fibre is odd and z_u is odd. A non-self-converse fibre is a disjoint pair of odd ordinary fibres exchanged by sigma. Reflection preserves Hamming weight, so the two signed sums agree, each is odd, and z_u is twice an odd integer. Hence

    v_2(z_u)=0 for SC nodes, and 1 for NS nodes.                  (8)

In particular **no merged node has equally many even-weight and odd-weight tilings**, Q_m has rank exactly 1, and every entry of Q_m is nonzero. For every pair of merged nodes, the total number of even-distance ordered tiling pairs differs from the number of odd-distance pairs. If such node pairs are signed by which total is larger, positive pairs lie within the two parity-charge sign classes and negative pairs cross between them. Every triangle has positive sign product: this symmetric signed graph is balanced, not a tournament. Both sign classes are nonempty because all charges are nonzero and sum to zero. Quantitatively the even/odd totals are (M_u M_v +/- z_u z_v)/2. Its diagonal obeys

    Q_m(u,u)=1 mod 8 for SC,
    Q_m(u,u)=4 mod 32 for NS.                                   (9)

More generally,

    Q_j(u,u)=binom(m,j) mod 8 for SC,
    Q_j(u,u)=4 F_j mod 8 for NS.                                (10)

For SC, each Z_j(u,S) is odd, so its square is 1 modulo 8. For NS, every Z_j(u,S) is even. If sigma S=S, the two ordinary sheet sums are equal odd integers, so Z_j(u,S)=2 mod 4 and its square is 4 mod 8. If sigma S != S, equal columns occur in pairs and contribute 2 times an even square, divisible by 8. This proves (10); the unique top column gives the sharper modulus in (9).

## 4. Explicit signed-Rédei bridge

For a labelled tournament T put x_ij=+1 if i->j, and -1 otherwise, for i<j. Set

    g(T)=product_{i<j}x_ij,
    R(T)=sum_{p in HP(T)}sgn(p),
    J(T)=g(T)R(T)/|Aut(T)|.

The product g(T) and the signed Rédei count R(T) each change by the sign of a relabelling, so their product is relabelling-invariant. Converse multiplies both by (-1)^{binom(n,2)}, the second because reversing a Hamiltonian word has that sign. Thus J is a **signed, converse-invariant integer**, unlike the uncorrected pseudo-invariant R.

If u has representative T and epsilon_u=1 for SC or 2 for NS, then

    z_u=(-1)^m epsilon_u J(T),       a_u=epsilon_u J(T).           (11)

Proof. Normalizing a Hamiltonian word p to positions gives a tiling t with

    g(T_p)=sgn(p)g(T)=(-1)^{m-|t|}.

Thus (-1)^{|t|}=(-1)^m sgn(p)g(T). The automorphism group has odd order and hence consists of even permutations. Its free action on HP(T) makes sgn(p) constant on every orbit, so summing and dividing by |Aut(T)| proves the ordinary-class formula. Converse contributes the same tiling sign sum, giving epsilon_u. This also proves integrality of J directly.

This is a genuine preserved map: the signed HP-orbit enumerator equals a scalar read from all-distance line counts. It does **not** identify complement transport with a signed-Rédei eigenvector; that would require a separate invariant-subspace condition, generally absent.

## 5. The exact compression-loss matrix

Let A_d be distance-d adjacency on X, and let Y_{x,u}=1[x in F_u]. Then Y^T Y=D, W_d=Y^T A_d Y, and P=Y D^{-1}Y^T is orthogonal projection onto fibre-constant functions.

The cube intersection numbers are

    p_ab^c=binom(c,(a-b+c)/2)binom(m-c,(a+b-c)/2),

with the value zero unless both lower arguments are integral and in range. Define

    Gamma_ab=sum_c p_ab^c W_c-W_a D^{-1} W_b
            =Y^T A_a(I-P)A_bY.                                 (12)

The block matrix (Gamma_ab)_{a,b=0}^m is positive semidefinite. It is the Gram matrix of the blocks (I-P)A_bY. In particular

    mD+2W_2-W_1D^{-1}W_1 >= 0.                                 (13)

Equation (13) is the exact within-fibre covariance of local neighbour profiles. If h_x(v)=#{y in F_v:dist(x,y)=1}, then it equals

    sum_u sum_{x in F_u}(h_x-bar h_u)(h_x-bar h_u)^T,
    bar h_u(v)=W_1(u,v)/M_u.                                   (14)

Consequently equality holds if and only if every h_x is constant within its fibre, namely the merged partition is equitable for single-coordinate flips. This is not generally true. At n=3 equality holds. At n=4, in canonical-code order (0,2,4), masses are (1,2,5) and the matrix is

    [ 4/5, -4/5,  0 ]
    [-4/5, 14/5, -2 ]
    [  0,   -2,   2 ],

of rank 2. The equality boundary is thus sharp already at the first nontrivial merged fibre.

The intersection number follows by counting an intermediate y: among c differing coordinates of x,z it disagrees with x in (a-b+c)/2 places, and among the m-c agreeing coordinates in (a+b-c)/2 places. This gives A_a A_b=sum_c p_ab^c A_c and proves (12). The covariance expression is the entrywise action of I-P. It explains how commuting cube moves can acquire noncommuting compressed shadows, rather than incorrectly assuming lumpability.

## 6. A count-, parity-, and complement-preserving impossible arrangement

At n=5 use canonical merged codes

    (0,2,4,5,8,10,12,40,41,76)

with masses

    (1,2,10,1,9,9,11,5,13,3).

The code is the least full-arc bitmask over all vertex permutations and converse. Indices below refer to the displayed order, not explorer structural addresses.

Define a symmetric zero-diagonal matrix Delta by the four nonzero unordered pairs

    Delta(0,2)=2,  Delta(1,6)=2,
    Delta(0,1)=-2, Delta(2,6)=-2.

Replace W_1 by W_1+Delta and W_2 by W_2-Delta, leaving every other W_d unchanged. The affected upper-triangle counts are:

    pair   W1 old -> new     W2 old -> new
    0,1       2 -> 0            0 -> 2
    0,2       2 -> 4            6 -> 4
    1,6       0 -> 2            4 -> 2
    2,6       8 -> 6           26 -> 28.

All entries remain nonnegative. Every layer retains its complete row-sum vector, every diagonal, every entry parity, the total matrix sum M M^T, and consequently all scalar per-fibre within-distance distributions and their scalar Walsh diagonal tests. W_6, and hence the **entire original blue/black complement graph including line multiplicities and loops**, is unchanged. The layerwise cut-degree parity laws remain unchanged as well. This does not claim to preserve every finer tiling-colour or coordinate-labelled invariant.

Nevertheless the new Q_6=Q_6(old)-2Delta has on nodes 1 and 6 the principal submatrix

    [ 4, -10 ]
    [-10,  9 ],

whose determinant is -64. No partition of the six-dimensional cube, let alone the actual tournament partition, can realize it. In particular neither degrees plus parity, nor their supplementation by each fibre's scalar distance distribution, guarantees a coherent joint arrangement. This is a genuinely matrix-valued obstruction.

## 7. Exact audit and boundary ledger

The self-contained program `probe.py` independently canonicalizes the fixed-path tournaments, builds every W_d from ordered pairs, and checks each transformed Q_j against a direct Walsh Gram computation. It also checks reflection rank bounds, odd-fibre diagonal congruences, and the covariance identity by direct neighbour lists. No stored atlas is imported.

    n   tilings nodes   ranks Q_0,...,Q_m                 rank Gamma_11
    3       2      2    1,1                                      0
    4       8      3    1,2,2,1                                  2
    5      64     10    1,4,8,7,7,4,1                            9
    6    1024     34    1,6,16,24,27,26,26,17,14,5,1             33.

For n=6 the reflection caps are 1,6,25,64,110,132,110,64,25,6,1. These are upper bounds, not claimed attained. Exact rational ranks use SymPy DomainMatrix, not floating eigenvalues. `signed_redei_check.py` independently enumerates Hamiltonian words and automorphisms for the signed bridge. `hostile.py` produces and checks the displayed switch.

Reproduce from this result directory (Python 3 with NumPy and SymPy):

    python3 verify.py

This regenerates both size audits, the signed-Rédei checks, and the hostile witness. The standalone `probe.py`, `hostile.py`, and `signed_redei_check.py` provide narrower checks.

The n=6 output is `audit_n6.json` (the final verifier is consolidated in `verify.py`). All identities are proved for all n; the finite calculations are implementation controls and hostile witnesses, not induction evidence.

Loss ledger: pi forgets the literal tiling position and the HP/Aut orbit. W_d retains every two-fibre distance total but forgets which individual tilings meet. Q_j retains a Gram matrix but loses the labelled character columns and their compatible phases. Positivity, rank, integral-Gram, congruence, and covariance constraints are **necessary**, not proved sufficient for any cube partition or for tournament fibres. The precise next question is which additional compatible Fourier/Boolean-incidence constraints close that realizability gap. No map to LRC or any unrelated open problem is asserted.
