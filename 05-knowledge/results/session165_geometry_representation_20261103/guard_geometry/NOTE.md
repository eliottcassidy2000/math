# Guard-state geometry: a dense algebraic orbit with a thin physical closure

**Status.** Proved statements about the literal one-vertex C/D guard language of Session 164, exact finite controls, and an explicit strong source family outside that language. No tournament-spectrum conjecture is resolved. No claim of literature novelty is made.

## Source and conventions

The input is [S164's positive moving-guard calculus](https://github.com/eliottcassidy2000/math/blob/89416555c4d95c6fa20c9ebafca1a180ad3fcaf3/05-knowledge/results/session164_emergent_20261103/profile_algebras/GUARD_ALGEBRA.md), with its [canonical return classification](https://github.com/eliottcassidy2000/math/blob/89416555c4d95c6fa20c9ebafca1a180ad3fcaf3/05-knowledge/results/session164_emergent_20261103/profile_algebras/CANONICAL_RETURNS.md). This note studies the geometry and integral representation of those already proved physical operations. The response-polynomial discussion below is a finite identity-testing consequence, not a claim that an arbitrary polynomial is a realizable tournament response.

For a tournament X with guard q having unique in-neighbor r, use the column state

    v=(s,u,d),
    s = Hamilton paths of X starting at q,
    u = Hamilton paths of X not starting at q,
    d = Hamilton paths of X−q starting at r.

Thus H=s+u. Add a new vertex a with q→a and a→X−q. Move C marks a as the new guard; move D keeps q. The exact matrices are

          [0 1 0]           [1 0 1]
    C  =  [1 1 1],     D = [0 2 0].
          [1 0 0]           [0 0 0]

Each move costs one vertex and preserves strongness. A directed triangle has v0=(1,2,0). Words are written in **temporal order**: CD acts as the matrix DC. Statements involving matrix products, such as DCD below, use ordinary multiplication; the palindromic reset is the same in both conventions.

The distinction between three objects is essential:

1. the set of matrices made by positive words;
2. the integer algebra formed by sums and differences of those matrices;
3. the set of actual physical states reached from v0.

They are not interchangeable.

## 1. Three literal suffixes recover the state

For any guarded input, put

    h0=H(X), h1=H(CX), h2=H(C²X).

Then

    h0=s+u,
    h1=s+2u+d,
    h2=3s+3u+2d.

The observation matrix is

        [1 1 0]
    O = [1 2 1],       det(O)=2.
        [3 3 2]

Consequently

    s=(h0−2h1+h2)/2,
    u=(h0+2h1−h2)/2,
    d=(h2−3h0)/2.                                      (1)

This is an exact physical three-count test, using at most two extra vertices. In algebraic-geometric terms the observation map is a linear isomorphism after inverting 2; its constant Jacobian determinant is 2. There is no appeal to a general Jacobian conjecture or nonlinear invertibility principle.

For two C/D word matrices M,N, the following nine tests are equivalent to M=N:

    H(C^i ; M ; C^j applied to C3)
      = H(C^i ; N ; C^j applied to C3),    0≤i,j≤2.       (2)

Indeed these differences form O(M−N)V, up to transposing the array of scalar tests, where

        [1 2 3]
    V = [2 3 6],       det(V)=−2.
        [0 1 2]

Both matrices are invertible over Q. Every test in (2) is a literal strong tournament. Each tested word uses at most four additional C-vertices beyond the word itself. The same equivalence holds modulo every odd integer, because 2 is then a unit.

The scalar C-sequence already forces dimension three for any linear realization of this source language:

    det [[3,5,9],[5,9,17],[9,17,31]] = −4.                (3)

A realization of dimension at most two would factor this Hankel matrix through a space of dimension at most two. The displayed three-state realization is therefore minimal over Q and over every field of odd characteristic. Modulo 2 the scalar series is instead constantly 1, so that scalar language has a one-dimensional realization. These statements concern this restricted language, not arbitrary tournament operations or all path-cover profiles.

## 2. The full matrix algebra, with an exact characteristic-two defect

The physical reset is

    R=DCD=v0 (1,2,1).

Let l=(1,2,1). The row matrix with rows l,lC,lC² has determinant 2, while V above has determinant −2. Therefore the nine rank-one matrices

    C^i R C^j,       0≤i,j≤2,

are a basis of Mat_3(Q), and a basis of Mat_3(F) over any field F of odd characteristic. In particular the representation has no nonzero proper common invariant linear subspace over any such field. This does **not** say that every matrix is a physical positive word or that every finite-field state is reached.

There is also an exact integral version. Write a general matrix as

        [a b c]
    M = [e f g].
        [h i j]

Then Z[C,D] consists exactly of the integer matrices satisfying

    e=h,
    f=h+j,
    g=h+i,
    a=b+c+h+i+j                          modulo 2.     (4)

Thus its additive index in Mat_3(Z) is 16. Equivalently, modulo 2 it preserves the flag

    <(1,1,1)> ⊂ {(s,u,d):u=d} ⊂ F_2³,

and has equal induced scalars on the first and last one-dimensional quotients of that flag.

**Proof.** C and D have this flag property, and sums and products preserve it. The associated integral congruence lattice has index 2^4. The nine temporal words

    empty, C, D, CC, CD, DC, DD, CCD, CDD

have vectorized determinant 16. Their integer span is therefore already the entire congruence lattice, proving (4) and the index assertion. The determinant is an explicit nine-by-nine certificate reproduced by the checker. This proof also shows that all failure of the *algebra* to be the full integer matrix lattice is 2-primary. It makes no assertion that the positive semigroup lacks odd-prime reachability obstructions. ∎

## 3. Zariski density, with a bounded polynomial test bank

Let O be the set of actual C/D states obtained from v0. Then O is Zariski dense in affine three-space over Q: no nonzero polynomial in Q[s,u,d] vanishes on every state in O.

A stronger effective statement gives a practical certificate.

**Polynomial escape theorem.** If P(s,u,d) is a nonzero polynomial of total degree at most delta, one of the physical words

    D^m C D^n C^k,
    0≤m,k≤delta,  1≤n≤delta+1,                         (5)

has P(v)≠0. Every graph in this bank is strong, and has at most 3delta+5 vertices. Hence testing zero on the whole bank is equivalent to the polynomial identity P=0. The bank has at most (delta+1)^3 states, with possible duplications.

**Proof.** Put L={d=0}. Since det(C)=1, the planes L_k=C^k L are distinct. Otherwise C^t L=L for some t>0; but C^t v0 has positive d for every t>0, a contradiction. A polynomial of degree at most delta cannot vanish identically on delta+1 distinct planes unless it is zero: each plane's distinct defining linear form would divide it.

Choose k≤delta such that P restricted to L_k is nonzero. Direct calculation of the first three blocks of (5) gives

    D^n C D^m v0 = A_m (1,B_n,0),
    A_m=1+2^(m+1), B_n=2^n.

Therefore

    Q(A,B)=P(A C^k(1,B,0))

is nonzero and has degree at most delta in A and at most delta in B separately. The substitution is dominant onto L_k: over Q, it covers every point of that plane whose corresponding first coordinate is nonzero. The delta+1 values of A_m and the delta+1 values of B_n are distinct. Two applications of the univariate root bound show that Q cannot vanish on their whole grid. The word length is at most 3delta+2. ∎

For a **homogeneous** polynomial of degree at most delta, the smaller bank

    D^n C^k,       0≤n,k≤delta,                         (6)

suffices, using at most 2delta+3 vertices. On each plane only its projective slope must be varied, since scaling does not change whether a homogeneous polynomial vanishes. The same plane-factor and root-count proof applies.

This has a useful application to proposed response formulas. If two explicitly polynomial formulas of degree at most delta in this guard state agree on (5), they agree on every vector, and therefore on every physical source for which those formulas are valid. Conversely a spurious polynomial relation fitted to a short list of guards must fail in the stated finite bank. The assumption of polynomial dependence must be independently established; scalar counts alone do not establish it.

This is characteristic-zero algebraic density, not density among integers, complete modular reachability, or exact odd-count coverage. In particular it does not remove inequality constraints.

## 4. The projective closure is nevertheless thin

Every state in O has s,u>0, d≥0, and u≥s+d. In coordinates

    x=s/u, y=d/u,

the containing triangle is Delta={x≥0,y≥0,x+y≤1}, and the maps become

    C(x,y)=(1/(1+x+y), x/(1+x+y)),
    D(x,y)=((x+y)/2,0).                                 (7)

Let I={(x,0):0≤x≤1/2}. The seed lies in I, and every D sends Delta to I. It follows by taking the last D in a word that the entire projective orbit lies in

    union_(k≥0) C^k(I).                                 (8)

Each term is a projective line segment. Moreover these segments converge uniformly to

    p=(1/lambda,1/lambda²),
    lambda³−lambda²−lambda−1=0, lambda>1.                (9)

For completeness, C³ is strictly positive. The characteristic polynomial has one real root lambda>1 and a conjugate pair of modulus lambda^(−1/2). The positive left eigenvector (lambda(lambda−1),lambda,1) has a positive minimum pairing with the normalized compact segment I. Spectral decomposition therefore gives uniform convergence of the normalized images to (9).

Thus the projective closure of O is contained in the closed set

    K={p} union union_(k≥0) C^k(I).

This set has Hausdorff dimension at most one, hence empty interior in the projective plane. The Zariski closure is two-dimensional while the ordinary topological closure is thin. This is a concrete warning against reading algebraic density as geometric or arithmetic completeness.

## 5. A forbidden ratio interval and a strong family beyond the bank

For states with d=0, the restriction is sharper. Every projective state of O on that line satisfies

    s/u belongs to union_(m≥1) [3/2^(m+2), 1/2^m].       (10)

In particular no such state has

    1/4 < s/u < 3/8.                                    (11)

**Proof.** Begin on the reset segment with vector (a,b,0), where a,b>0 and b≥2a. After k≥1 consecutive C moves, put

    E_k=4(s_k+d_k)−3u_k.

Its first three values are a+b, 2(b−a), and 3a, all positive. By the characteristic polynomial of C, E_k satisfies the tribonacci recurrence; hence E_k>0 for all k≥1. The slack inequality also gives s_k+d_k≤u_k. A following D therefore has ratio in (3/8,1/2]. If that final D block has length m, each subsequent D halves the ratio, giving the interval in (10). Pure D words and the seed are right endpoints of such intervals. ∎

The obstruction is to the **marked state**, even up to multiplication by a positive scalar. It need not forbid a Hamilton count from another marked state.

Here is a literal all-order strong family inside the missing interval. For n≥6, start with the transitive order 0,1,…,n−1 and reverse exactly the two edges

    n−1 → 0,       n−1 → 2.

Call the result Y_n and mark q=0. The directed cycle 0→1→…→n−1→0 proves strongness, and q has unique in-neighbor n−1. Its exact state is

    (s,u,d)=(1+2^(n−4), 2^(n−2)−2, 0),
    H(Y_n)=5·2^(n−4)−1.                                 (12)

The counts begin 19,39,79,159,319,… . Writing t=2^(n−4)≥4 gives

    s/u=(t+1)/(4t−2),

which lies strictly between 1/4 and 3/8. Therefore none of these marked states is even projectively in O.

**Counting proof.** A Hamilton path can use at most one backward edge, since both have tail n−1. There is one completely increasing path. For a path using n−1→0, choose the vertices occurring before n−1 from {1,…,n−2}; each side is forced to increase. The only forbidden choices have maximum 2, namely {2} and {1,2}. This gives 2^(n−2)−2 paths. For a path using n−1→2, the initial increasing block must contain 0,1 and any subset of {3,…,n−2}, giving 2^(n−4) paths. Exactly these last paths and the completely increasing path start at 0. After deleting 0, a path starting at n−1 must next use 2 and then cannot reach 1, so d=0. These disjoint cases prove (12). ∎

There is a literal one-vertex extension Y_n→Y_(n+1): insert a new vertex immediately before the terminal vertex n−1 in the transitive presentation, retaining the same two backward targets. On this family it gives

    (s,u,0) → (2s−1, 2u+2, 0),
    H → 2H+1.                                          (13)

This operation is not asserted to satisfy (13) for arbitrary guarded tournaments. It is a source-family extension outside the C/D language. The known reset CD, valid on d=0, then turns any Y_n into the canonical state H(Y_n)v0 with two more vertices. This is a new source route to that canonical state, not a proof of a reusable multiplier H(Y_n) acting on every canonical input.

## Reproduction and boundaries

Run `python3 check.py` and `python3 -O check.py`. All checks use explicit exceptions rather than removable assertions. The self-contained checker verifies the transfer on all qualifying guarded labeled tournaments through order four; the observation and matrix-algebra certificates; the parity flag; finite polynomial-bank ranks; plane separation; all C/D words through a stated finite depth; and the literal family Y_n by independent subset/endpoint Hamilton-path DP. It records exact test counts and outputs as JSON.

The general proofs are the matrix factorizations, grid root bound, positivity/spectral argument, tribonacci inequality, and direct backward-edge path enumeration. Finite tests corroborate them rather than replacing them. The projective gap, the proper positive-word semigroup, characteristic-zero Zariski density, and the full odd-field *linear algebra* are compatible facts, not contradictory completeness claims.
