# Physical tournament profiles as positive lattice-path states

Session 197, 2026-11-05. **Exact all-order operator identities, positive integer states, physical continued fractions, and scoped obstructions.** The full Hamilton-path spectrum with exceptions `{7,21}` remains **OPEN**. No literature-priority claim is made.

## 1. Inheritance and the question answered

The [S195 colored path-cover algebra](https://github.com/eliottcassidy2000/math/blob/ee07a5aac9cdfd8d3834fa834df4fb357fcaf769/session195/pathcover-algebra/PROOF.md) gives an exact physical composition interface. The [S196 cloning proof](https://github.com/eliottcassidy2000/math/blob/cf3b44b9a442c3cff0f28cab924e73c4c00bb00c/session196/spectral-moves/PROOF.md) gives matching rank, binary affine clone rays, and limits of a finite cloning bank. Those mechanisms do not justify assuming that every positive polynomial change of coordinates preserves real-rootedness, total positivity, or physical realizability.

The exact inherited starting operator is Section 6 of [S193 profile substitution](https://github.com/eliottcassidy2000/math/blob/169792fb5eeb483e351112b5a09706e6d7bd89f0/05-knowledge/results/session193_triangulated_proof_20261105/profile_substitution/PROOF.md). It defines the labelled-bin polynomial and proves

    D f(x) = x^2 f(x) + x f(x+1) - x f(x-1).

The odd-cycle formula used below is inherited, not a new claim: see [S196 symmetry](https://github.com/eliottcassidy2000/math/blob/cf3b44b9a442c3cff0f28cab924e73c4c00bb00c/session196/symmetry/PROOF.md) and its [S194 weak-deletion source](https://github.com/eliottcassidy2000/math/blob/188bb29fa08b6dfd51ab94cc8e494d6d56f5e9ec/session194/weak-deletion/PROOF.md). We give a self-contained algebraic derivation to specify conventions. The parallel S197 wildcard work also uses its positive labelled-bin expansion.

The present bridge is an explicit change to two central-factorial bases in which **every physical tournament profile has nonnegative integer coordinates**, ordered joins have nonnegative integer linearization, and D becomes tridiagonal. A second actual wrapper has a positive four-band continuation and escapes the first arithmetic obstruction. These are unbounded state spaces, not a contradiction of the inherited infinite continuation-rank theorem. The coordinates do not retain a marked root or endpoint colors, so no claim is made that they close every S195 colored operation.

## 2. Bin profiles and an exact positive expansion

A path cover is an unordered set of nonempty vertex-disjoint directed paths covering all vertices. Let c_k(T) count covers with k paths. Put

    B_T(x) = sum_k c_k(T) (x)_k,
    (x)_k = x(x-1)...(x-k+1).

For nonnegative integer x, this counts covers placed in x labelled bins, allowing empty bins. Thus B_T(1)=H(T), B_T(0)=0 for nonempty T, B_empty=1, and B_Tr_n=x^n. Ordered joins obey B_(A->B)=B_A B_B, since each bin is its A-path followed by its B-path.

An odd directed cycle below has length at least three and is counted modulo rotation, not modulo reversal. A packing Gamma is an unordered set of vertex-disjoint such cycles. Then

    B_T(x) = sum_Gamma 2^|Gamma| x^[n - sum_(C in Gamma)(|C|-1)].       (1)

In particular B_T is monic, has nonnegative integer monomial coefficients, has parity n, and has no constant term for n>0. Also H(T) is odd: the empty packing contributes one and every other contribution is even.

Here is an algebraic proof, including the determinant signs. Work over the commutative squarefree ring with variables z_v satisfying z_v^2=0; let Z=diag(z_v), A be the 0/1 tournament adjacency matrix, and J the all-ones matrix. The polynomial

    P(z)=sum_(S subset V) H(T[S]) product_(v in S) z_v

is 1 plus the sum of all nonempty directed words. Repeated vertices vanish, so the matrix determinant lemma gives

    P = det(I-ZA+ZJ)/det(I-ZA)
      = product_v(1+z_v) det(I+ZA)/det(I-ZA).                         (2)

For the second equality use A+A^t=J-I, Z^2=0, and
(I+Z)(I+ZA^t)=I-ZA+ZJ; transposition and det(I+UV)=det(I+VU) identify the remaining determinants. All denominators have constant term one.

Taking formal logarithms in this finite nilpotent ring cancels even closed walks and gives

    P = exp(sum_v z_v + 2 sum_(odd directed cycles C) product_(v in C) z_v).

The coefficient of the full vertex monomial in P^d counts covers in d labelled bins. Expanding the last exponential to the d-th power proves (1) for every nonnegative integer d and hence as a polynomial identity. Distinct directed cycles on the same vertex set remain distinct objects. This derivation does not turn a signed determinant expansion into a physical subtraction operation.

## 3. A positive integer state for both parities

Write

    P_k(y)=product_(j=1)^k (y-j^2),       P_0=1,
    q_k(x)=x P_k(x^2),
    r_k(x)=x^2 P_k(x^2),                 k>=0.

Use q for odd-order nonempty tournaments and r for even-order nonempty tournaments.

**Theorem 1 (positive central-factorial state).** Every such tournament has a unique expansion

    B_T=sum_k u_k q_k       if |T| is odd,
    B_T=sum_k u_k r_k       if |T| is even,

where all u_k are nonnegative integers, u_0=H(T), and the highest coordinate is one. Ordered joins have nonnegative integer structure constants in these two bases.

**Proof.** The bases are monic and triangular in their respective positive-degree parity spaces. For either b=q or b=r,

    x^2 b_k = b_(k+1) + (k+1)^2 b_k.                                (3)

Starting with x=q_0 or x^2=r_0 shows that every positive-degree monomial of the appropriate parity expands with nonnegative integer coordinates. Apply (1). At x=1, b_0(1)=1 and b_k(1)=0 for k>0, giving u_0=H.

For products suppose i>=j. First factor out P_i(x^2), and multiply successively by x^2-s^2 for s=1,...,j. At every intermediate index k>=i>=s,

    (x^2-s^2)b_k=b_(k+1)+[(k+1)^2-s^2]b_k

has nonnegative integer coefficients. For q_i q_j the initial parity factor is r_i; for q_i r_j it is x^2 q_max(i,j); for r_i r_j it is x^2 r_max(i,j). Extra x^2 factors are handled by (3). This proves all three parity cases. QED.

Appending a universal source or sink is multiplication by x. Its positive transitions are

    X q_k=r_k,        X r_k=q_(k+1)+(k+1)^2 q_k.                     (4)

This theorem asserts that physical profiles lie inside a positive coordinate cone. It does **not** assert that an arbitrary nonnegative coordinate vector is a physical tournament profile.

## 4. The cyclic two-singleton wrapper is tridiagonal

Let W(T)=C3[T,K1,K1]. Concretely add a,b with a->T->b->a. The inherited full-profile identity is B_W=D B_T with D as in Section 1.

**Theorem 2 (two positive Jacobi operators).** With terms of negative index omitted,

    D q_k = q_(k+1) + (k^2+6k+3)q_k + 2k^2(2k+1)q_(k-1),          (5)
    D r_k = r_(k+1) + (k^2+6k+5)r_k + 2k(k+1)(2k+1)r_(k-1).      (6)

**Proof.** D preserves parity. For k>=1, Dq_k vanishes at 0 and at every integer +/-j with 1<=j<k, since the three input arguments j-1,j,j+1 are roots of q_k. Thus Dq_k is q_(k-1) times an even polynomial of degree four, leaving precisely the three displayed basis terms. The same reasoning applies to Dr_k, which also has a double zero at zero: r_k is even, and x[r_k(x+1)-r_k(x-1)] is divisible by x^2.

The leading coefficient is one. Comparing the next coefficient, using sum_(j<=k)j^2 for the second coefficient of each basis element, gives k^2+6k+3 or k^2+6k+5. Finally q_(k-1)(k)=(2k-1)! and q_k(k+1)=(2k+1)!, so evaluation at k gives the downward coefficient 2k^2(2k+1). Multiplication by the additional endpoint factors in r gives 2k(k+1)(2k+1). For k=0, direct evaluation gives D x=x^3+2x and D x^2=x^4+4x^2, agreeing with (5)-(6). QED.

Every continuation D^t B_T can therefore be evaluated as a nonnegative weighted lattice-path sum on heights 0,1,2,... . Upward steps have weight one, level steps at k have the middle coefficient a_k, and downward steps k->k-1 have the final coefficient b_k. The starting height k is weighted by u_k, and the terminal height is zero. This is an exact evaluation of a physically constructed graph, not a realization of arbitrary lattice paths as independent tournaments.

In particular, taking the zero-coordinate of (5)-(6) gives

    H(W(T)) = 3h+6u_1     for odd |T|,
    H(W(T)) = 5h+12u_1    for even |T|,                              (7)

where h=H(T). Missing u_1 means zero.

### Actual continued fractions and Hankel products

Let h_t^odd=H(W^t(K1)) and h_t^even=H(W^t(Tr_2)), including t=0. These are the vacuum states q_0 and r_0. Their formal ordinary generating functions are

    sum_(t>=0) h_t z^t
       = 1 / (1-a_0 z - b_1 z^2/(1-a_1 z - b_2 z^2/(1-a_2 z-...))),  (8)

with (a_k,b_k) from the corresponding line (5) or (6). This is the Jacobi, not the Stieltjes, continued-fraction normalization. All displayed weights are positive. First-return decomposition of a height-zero Motzkin excursion proves (8); each coefficient uses only finitely many heights.

For N>=1,

    det(h_(i+j))_(0<=i,j<N) = product_(k=1)^(N-1) b_k^(N-k) > 0.     (9)

For a short proof, symmetrize the tridiagonal matrix: its neighboring entries become sqrt(b_k). The moment matrix is the Gram matrix of e_0,J e_0,...,J^(N-1)e_0. Their coordinates are triangular, with diagonal products sqrt(b_1...b_k). Taking its squared determinant gives (9). This proof uses only finite matrix corners for each N, so no analytic moment convergence is needed.

The first sequences are

    odd:  1,3,15,123,1479,23955,496479,12725547,...
    even: 1,5,37,389,5557,103157,2399173,67980965,... .

For example the odd leading determinants start 1,6,1440,43545600. This is a physical positive continued fraction on a restricted family, not a universal moment theorem for arbitrary physical initial states.

## 5. A physical positive escape operator

The grammar generated from singleton tournaments by ordered joins and W misses every count coprime to three in residues 7 or 11 modulo 12.

Indeed the base count is one. An odd-order W output is divisible by three. An even-order W output coprime to three has residue 5h modulo 12 by (7); if h was 1 or 5 modulo 12 the result is again 1 or 5. A coprime-to-three product has only coprime-to-three factors, and {1,5} is closed under multiplication modulo 12. This proves the invariant for every expression tree, with both parities explicitly included.

Now replace one of the two exterior singleton blocks by a transitive edge:

    V(T)=C3[T,Tr_2,K1].

This is one transitive clone of an outer vertex of W(T), so it is an actual permitted graph change rather than a formal subtraction. In the S193 multi-port theorem take a transitive two-vertex core and both vertices as ports. Its bin polynomial is x^2; the summed first-bin start and end responses are 2x-1 and 2x+1, and the two-distinct-bin response is two. The exact operator is therefore

    E f=x[(x^2-2)f(x)-(2x-1)f(x-1)+(2x+1)f(x+1)].                  (10)

Equivalently, E=DX+XD-X^3-2X. This identity is algebraic, but the positive basis form is

    E q_k = r_(k+1)+(k^2+10k+5)r_k+2k(2k+1)^2 r_(k-1),             (11)

    E r_k = q_(k+2)+(2k^2+14k+13)q_(k+1)
            +(k+1)(k^3+19k^2+27k+11)q_k
            +2k^2(k+1)(2k+1)^2 q_(k-1).                           (12)

Substitution of (3)-(6) in E=DX+XD-X^3-2X proves these identities exactly. Every resulting coefficient is a nonnegative integer; at k=0 the last term is zero. In particular E changes order parity, as adding three vertices requires.

Taking the zero-coordinate yields

    H(V(T))=5h+18u_1     for odd |T|,
    H(V(T))=11h+36u_1    for even |T|.                              (13)

Consequently V(Tr_2) is a physical H11 escape from the previous modulo-12 invariant. In outgoing-bitmask notation its rows are

    [14,12,24,16,3].

Its parent W(Tr_2) has H5. No change to the internal input arcs is needed for this particular escape.

### The next obstruction is explicit, not hidden

The enlarged grammar generated by singletons, ordered joins, W, and V still misses H13 and H19. If a prime count is obtained by an ordered join, exactly one factor has that count and the others have count one; descend to that factor. Eventually a W or V operation must first produce the prime, since singleton leaves have count one. Its input h is a positive odd integer and u_1 is nonnegative.

For target 13 or 19 the odd W formula is excluded modulo three. In the even W and odd V formulas the bounds 5h<=target leave only h=1, or h=1,3 respectively; substitution into 5h+12u_1 and 5h+18u_1 yields neither target. In the even V formula 11h<=target forces h=1, and 11+36u_1 yields neither. This excludes both targets from the entire grammar, not merely small expression trees.

These are **grammar omissions**, not new tournament spectrum gaps. Both targets are known attainable outside this restricted language. The positive central-factorial cone and its first two physical operators do not resolve inverse targeting. Additional operations need fresh exact transfer and sign checks; positivity cannot manufacture them.

## 6. Three hostile checks against overextending the bridge

### 6.1 Arbitrary physical starting states need not give positive moment matrices

Start with the actual transitive tournament Tr_m. Directly from D,

    H(Tr_m)=1,
    H(W(Tr_m))=1+2^m,
    H(W^2(Tr_m))=2*3^m+5*2^m-1.

The first 2-by-2 Hankel determinant is

    2*3^m+3*2^m-2-4^m.

For m=4 this is -48, with physical counts 1,17,241. For m=5 it is -444, with counts 1,33,645. Thus positive coordinates and a positive transition matrix do not make every initial-state continuation a positive moment sequence. The vacuum restriction in (8)-(9) is essential. Both parity cases fail.

### 6.2 Even the simplest cloning ray eventually loses imaginary-rootedness

The S196 triangle cloning family is W(Tr_m), with H=2^m+1. Its exact bin profile is

    B_m(x)=x^(m+2)+x[(x+1)^m-(x-1)^m].                            (14)

All its zeros are purely imaginary (zero included) precisely for 1<=m<=5. For m=6, remove its x^2 factor and put y=x^2. The reduced polynomial is

    y^3+12y^2+40y+12,

whose discriminant is -8752. For m=7 the reduced polynomial is

    y^4+14y^3+70y^2+42y+2,

with discriminant -224218032. Hence both have nonreal y-roots and B_m has zeros off the imaginary axis.

For even m>=8, the reduced monic polynomial has degree d=m/2 and first coefficients a_1=2m, a_2=2 binom(m,3). If all its roots were real, positivity of its coefficients would force them to be negative, and Newton's inequality would require

    a_1^2 >= [2d/(d-1)]a_2,

which simplifies to m<=7. For odd m>=9, its degree is (m+1)/2 and the same inequality requires m^2-7m-2<=0, again impossible. For m=1,...,5 the reduced polynomials are y+2, y+4, y^2+6y+2, y^2+8y+8, and y^3+10y^2+20y+2. Their discriminants are positive when relevant (28,32,7092), and all real roots are negative. This proves the threshold, not merely a numerical root plot.

### 6.3 Odd-cycle packing polynomials are not universally matching polynomials

Let A_k be an ordered join of k directed triangles and let T_k=W(A_k). Define

    Psi_T(z)=sum_(odd-cycle packings Gamma) z^|Gamma|.

Every cycle not inside an individual triangle must use both added vertices a,b, so at most one such exterior cycle is selected. Within each triangle, an exterior path can select 0,1,2,3 vertices in 1,3,3,3 ways respectively. If zero vertices are selected, its own triangle may also be packed. The per-block polynomial is therefore

    1+z+3t+3t^2+3t^3.

An exterior cycle is odd exactly when the selected internal path has odd length. Extracting odd t-parity by evaluations at t=+1,-1 gives the exact all-order formula

    Psi_(T_k)(z)=(1+z)^k+(z/2)[(z+10)^k-(z-2)^k].                 (15)

This is a physical enumeration, despite the parity-extraction subtraction. At k=3,

    Psi=1+507z+147z^2+19z^3.

Its derivative has discriminant -29160 and positive leading coefficient, so the cubic is strictly increasing on the real line and has only one real zero. It cannot be the ordinary matching generating polynomial of a graph with nonnegative edge weights, whose zeros are real and negative. No universal matching-polynomial or real-rootedness transfer is available for these odd-cycle packings. Nevertheless Psi(2)=1755 exactly equals the physical Hamilton count.

## 7. Verification and scope

Run in this packet:

    python3 verify.py > normal.stdout
    python3 -O verify.py > optimized.stdout
    cmp normal.stdout optimized.stdout

The verifier uses the standard library and explicit exceptions, so its gates survive optimization. It independently computes physical Hamilton counts by vertex-subset endpoint DP and bin polynomials by an unordered subset-cover recurrence followed by the falling-factorial transform. Literal permutation-and-cut enumeration is a separate expected-cover control, not another use of that recurrence.

The finite universe is:

- all 1,099 labelled tournaments of orders 1 through 5: central-factorial coefficient positivity, bin parity, and Hamilton readout;
- all 75 labelled tournaments of orders 1 through 4: independent literal covers and both physical wrapper full-profile identities, 150 actual output graphs;
- all basis indices 0 through 20 in each displayed D/E identity, 84 exact polynomial controls;
- 196 basis products covering both input parities and indices 0 through 6;
- both vacuum families through five actual wrappers, and six leading Hankel determinants in each parity;
- actual negative-Hankel witnesses in both parities, actual clone graphs at m=6,7, the eleven-vertex odd-cycle-packing hostile, and the H11 escape graph.

Odd-cycle packing counts in the hostile are computed independently from rooted directed-cycle subset DPs and a disjoint-packing recurrence, not inferred from (15). The resulting JSON contains exact rows and coefficients. These checks support the proofs and fix normalization; they do not constitute a spectrum census. All files were created in separate Session 197 scratch, with no dirty-checkout edits or worktree operations.
