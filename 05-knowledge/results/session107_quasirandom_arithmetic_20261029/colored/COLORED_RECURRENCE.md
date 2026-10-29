# Named-color tournament recursion and the refresh spectrum

Session 107, 2026-10-29. **PROVED extension of the Session 105 attachment/deletion calculation, with exact finite checks.** This is an enumeration and Markov-operator bridge, not a new prime-distribution estimate, a closed scalar counting recurrence, or a literature-priority claim. All vector-space statements below are over Q. Colors have fixed names; isomorphisms preserve every color, and never permute the palette.

## 1. From the arithmetic character code to a colored A000568

The fixed-base character theorem in `../characters/CHARACTER_STRUCTURE.md` produces an affine feasible palette with q=2^d colors, where d is the squareclass matrix rank after fixing the sign of -1. The number of input ports need not equal d. Its fixed-pattern limiting law is uniform on q colors independently at each vertex and fair independent moving-vertex edges. This motivates the following finite combinatorial model for **every integer q>=1**.

Let V_n(q) have basis all named-q-colored tournament isomorphism classes of size n, let a_n(q)=dim V_n(q), and set a_0(q)=1 and a_-1(q)=0. Burnside gives

    a_n(q) = sum_(lambda partition of n, all parts odd)
                 2^e(lambda) q^ell(lambda) / z_lambda,
    e(lambda) = sum_i (lambda_i-1)/2 + sum_(i<j) gcd(lambda_i,lambda_j),
    z_lambda = product_(m>=1) m^(c_m) c_m!.

Here c_m counts parts equal to m. An even permutation cycle reverses an antipodal vertex pair under a suitable power and therefore fixes no tournament. For odd cycles all unordered-pair orbits permit two consistent orientations; the displayed e counts these orbits. A fixed coloring is constant on each vertex cycle, giving q^ell. Thus q=1 specializes exactly to the ordinary tournament sequence A000568. The independent incoming colored-Burnside analysis supplies the count-polynomial lane; the verifier here instead enumerates vertex-permutation orbits directly.

## 2. Attachment, deletion, and their exact commutator

Define U_n:V_n(q)->V_(n+1)(q) by attaching a new vertex in each of q colors and each of 2^n edge-orientation masks, then forgetting labels. Every choice contributes, including repeated isotypes. Define D_n by summing the n vertex-deleted cards, also with multiplicity, and D_0=0. The definitions do not depend on the representative, since a relabeling bijects the choices. Deleting the newly attached vertex contributes q 2^n copies of the input. For deletion of an old vertex, each remaining colored attachment has exactly two preimages, distinguished by its erased edge. Therefore

    D_(n+1) U_n = 2 U_(n-1) D_n + q 2^n I.

There is no extra q in the first term: the new color is retained on both sides. Set E_n=U_n/(q 2^n). Then

    D_(n+1) E_n - E_(n-1) D_n = I.                 (1)

Terms through degree -1 are zero. This is precisely the normalized Session 105 Weyl relation, now with the color factor in the normalization.

For completeness, (1) yields D E^r=E^r D+r E^(r-1). Local nilpotence of D implies injectivity of E: if Ev=0, take the largest m with D^m v nonzero and apply D^(m+1). Surjectivity of D follows from the finite right inverse

    R(v)=sum_(k>=0) (-1)^k E^(k+1)D^k(v)/(k+1)!.

The harmonic projector P_n=sum_(k=0)^n (-1)^k E^k D^k/k! has DP_n=0, is the identity on H_n=ker D_n, and kills E V_(n-1). Consequently

    h_n(q):=dim H_n=a_n(q)-a_(n-1)(q),
    V_n(q)=direct-sum_(j=0)^n E^(n-j) H_j.        (2)

No integral splitting or positive-characteristic extension is asserted here.

## 3. Color-and-edge vertex refresh

For n>=1 the column-stochastic matrix

    K_n=E_(n-1) D_n/n

chooses a vertex uniformly, deletes it, and inserts a uniformly colored vertex with fair independent incident edges. Equivalently it refreshes **both the color and all incident edges** of a uniformly selected labeled vertex and then forgets labels. Keeping its old color would be a different chain and would not have the following spectrum.

On E^(n-j)H_j, equation (1) gives E D=(n-j)I. Hence the full, diagonalizable spectrum is

    1-j/n, with multiplicity a_j(q)-a_(j-1)(q), 0<=j<=n.   (3)

For q>=2, h_1=q-1>0 and the spectral gap (also absolute gap) is exactly 1/n for every n>=1. Explicit backward eigenfunctions are the centered color counts

    f_c(T)=number of vertices of color c - n/q,
    K_n^T f_c=(1-1/n)f_c.

Indeed refreshing a vertex removes an expected N_c/n instances of c and adds one with probability 1/q. The q functions have the single relation sum_c f_c=0, and span a q-1 dimensional first-mode space. By contrast q=1 has h_1=h_2=0 and h_3=1, so the Session 105 uncolored gap is 3/n for n>=3; the smaller uncolored state spaces are singletons. The slower colored mode is a concrete effect of retaining observable named-color composition.

For classes T at level n and S at level n+1, pointed realization counting gives

    U_n(S,T) |Aut_color(S)| = D_(n+1)(T,S) |Aut_color(T)|.

There is no factor n+1 in this multiplicity identity. Uniform labeled colored tournaments push forward to

    pi_n(T)= n! / (q^n 2^binom(n,2) |Aut_color(T)|).       (4)

The pointed identity implies pi_n(T) E_n(S,T)=pi_(n+1)(S)D_(n+1)(T,S)/(n+1), and thus reversibility and stationarity of K_n. The chain is irreducible: successively refresh every labeled vertex with a target color and target incident edges, with positive probability; later refreshes can retain those same target edge choices. There are self-loops. Formula (4) is exactly the limiting probability of an unlabelled named-color isotype under the previously proved fixed-base prime-pattern law: there are n!/|Aut_color(T)| labeled realizations, each with limiting probability q^-n 2^-binom(n,2). This matching is a limiting-distribution bridge; it does not assert that prime-prefix samples themselves evolve by K_n.

## 4. The critical source-free distinction

A universal source, if present, is unique and may have any one of q colors. Deleting it and retaining its color is a bijection between source-bearing n-vertex classes and pairs (color, (n-1)-vertex class). Thus the number b_n(q) of source-free classes is

    b_0(q)=1,    b_n(q)=a_n(q)-q a_(n-1)(q), n>=1.        (5)

For q>1 this is **not** the harmonic dimension in (2). Already b_1=0 but h_1=q-1. The latter consists of signed color contrasts, not source-free objects. More generally h_n=b_n+(q-1)a_(n-1) for n>=1. This is a dimension identity, not a canonical identification of those bases.

There is still a positive counting interpretation of h_n after choosing one distinguished color c0: classes with a universal source of color c0 are in bijection with all (n-1)-vertex classes, so h_n counts the classes having **no universal source of color c0**. This choice breaks palette-permutation symmetry and does not identify a positive class basis with ker D. For n>0, summing coordinates of Dv gives n times the coordinate sum of v, so every nonzero rational harmonic has both positive and negative coefficients. At q=2,n=1 the one-dimensional harmonic is the difference of the two singleton colors; the distinguished-source-free class is instead the singleton of the other color.

Iteratively stripping the source gives a unique source-free core preceded by an ordered word of named colors on the stripped transitive chain. Consequently, as formal series,

    A_q(t)=B_q(t)/(1-q t).

This does not compute the unknown source-free core counts. Defining delta to delete the source when present and give zero otherwise yields the separate integer relation delta U-U delta=q I: new-source attachments contribute q, while an old surviving source contributes U delta. For q>1 its kernel is not merely the span of source-free classes: it also contains differences between source lifts of the same core in different colors, since delta forgets the removed color. Indeed delta is surjective and dim ker delta=a_n-a_(n-1)=h_n, whereas the source-free subspace has dimension b_n. At n=1 the difference of two singleton colors is the simplest witness. Neither this kernel nor its source-free subspace should be identified with ker D. If instead delta_c deletes only a universal source of one distinguished color c, then delta_c U-U delta_c=I, and ker delta_c really is the coordinate span of classes having no source of color c: its other columns biject onto the previous-level basis. This is consistent with the distinguished-color positive count above, not a color-symmetric source-free identification.

The trace identity from (3) also extends:

    n tr(K_n)=sum_(j=0)^(n-1) a_j(q).

Again the transition matrix already contains structural information: no fast enumeration algorithm or closed scalar recurrence has been extracted.

## 5. Exact verification and scope

From the Session 107 packet directory run:

    python3 colored/verify_colored.py > colored/verification.json
    python3 -O colored/verify_colored.py > colored/verification.optimized.json
    cmp colored/verification.json colored/verification.optimized.json

SymPy 1.14.0 was used. The standalone checker imports no earlier tournament banks, canonicalizer, or Burnside formula. It enumerates all q-colored labeled tournaments through n=4 for q=1 and q=2 and quotients only by vertex permutations. It checks both commutators on input degrees 0..3, deletion surjectivity and attachment injectivity through level 4, source stripping, the pointed automorphism relation, every entry of ED/n against independent direct star-and-color refresh, exact full spectra through n=4, color-count eigenfunctions, stationary probability, detailed balance, and traces. Checks use explicit exceptions and survive optimization. Both runs passed with byte-identical JSON.

| q | a_0,...,a_4 | h_0,...,h_4 | b_0,...,b_4 |
|---|---|---|---|
| 1 | 1,1,1,2,4 | 1,0,0,1,2 | 1,0,0,1,2 |
| 2 | 1,2,4,12,48 | 1,1,2,8,36 | 1,0,0,4,24 |

In particular K_4 for q=2 has eigenvalues 1,3/4,1/2,1/4,0 with multiplicities 1,1,2,8,36. Finite tests verify these instances, not all degrees; the all-degree claims rest on the proofs above and the inherited Session 105 rational Weyl argument. No publication action was taken by this lane.
