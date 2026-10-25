# Exact enumeration of valuation-balanced Gaussian gap triples

Status: PROVED, with independent FINITE-EXACT checks. These formulas count triples passing a necessary valuation filter, not Schur relations, magic hourglasses, or magic squares. No literature-priority claim is made. The prior prime-support packet is unchanged.

## 1. Object and arithmetic interpretation

For an exponent vector a=(a_1,...,a_k), with each a_i≥1, let B(a)=∏_{i=1}^k{−a_i,...,a_i}. Let P(a) be its nonzero vectors modulo simultaneous sign j~−j. In the Gaussian parametrization of square gaps at a center whose split-prime exponents are a_i, P(a) indexes the distinct positive admissible gaps. Inert primes and powers of 2 affect only a common square scale.

Call an unordered triple of distinct points in P(a) balanced if, at every coordinate i, the maximum of the three absolute coordinate values is attained at least twice. The predicate is independent of the sign chosen for each point. A signed three-gap relation necessarily satisfies this predicate by the unique-minimum valuation obstruction. Its converse is not asserted and is false in general.

When every a_i=1, P(a) is PG(k−1,3) and balance says that no coordinate is nonzero in exactly one of the three points. For k=0 we use the empty point set and count zero triples, without asserting a nonempty projective space.

## 2. General exponent-box theorem

Define

    L(a)=12a²+8a+1=(2a+1)(6a+1),
    U(a)=4a+1,
    V(a)=2a²+4a+1,
    W(a)=2a+1.

The number of balanced unordered triples of distinct points is exactly

    N(a_1,...,a_k)
      = [∏L(a_i) − 3∏U(a_i) − 6∏V(a_i)
         + 14∏W(a_i) − 6]/48.

Empty products are one, giving N(())=0. For one prime the numerator vanishes identically, agreeing with the one-prime valuation obstruction. Example values are N(1,2)=15, N(2,2)=54, N(3,3)=252, N(1,1,2)=492, N(1,2,3)=3315, and N(2,2,2)=5098.

### Proof: local row counts

Represent an ordered triple of signed vectors as a k-by-3 integer matrix. The rows are independent, and each row must have its largest absolute entry at least twice.

For a coordinate bounded by a:

- The all-zero row contributes one. At positive maximum m, exactly two entries have absolute value m in 3·2²·(2m−1) ways, and all three do in 2³ ways. Hence the row count is
  1+Σ_{m=1}^a[12(2m−1)+8]=L(a).
- If a specified column is zero, the other two entries have equal absolute value. There is one zero row and four choices at each positive magnitude, giving U(a).
- If the second column is a fixed sign times the first, write the first entry as x. The third entry z must satisfy |z|≤|x|. Thus there is one all-zero row and 2(2m+1) rows at |x|=m, giving V(a).
- If all three columns are specified signs times the first, the first entry is arbitrary, giving W(a). The same W(a) counts a fixed equal-up-to-sign pair with the third entry zero.

### Proof: remove zero and repeated sign-orbit columns

Initially there are ∏L matrices. A specified zero column leaves ∏U choices. Two zero columns force the third column zero as well, so inclusion-exclusion leaves

    A = ∏L − 3∏U + 2

matrices with all columns nonzero.

For a specified pair of columns, sign-orbit equality means equality up to one global sign, not coordinatewise independently chosen signs. For each of the two signs the local row count is V. Removing a zero third column subtracts ∏W; the first column being zero already forces every column zero and is included in that removal. The two signs are disjoint on nonzero columns. Thus each pair-equality event has size

    B = 2(∏V − ∏W).

The intersection of any two different pair-equality events forces all three columns equal up to global signs. Choose the two relative signs in four ways and a nonzero first column in ∏W−1 ways. All pairwise and triple intersections are the same locus, of size

    C = 4(∏W − 1).

Removing the three pair events by inclusion-exclusion gives

    A−3B+3C−C = ∏L−3∏U−6∏V+14∏W−6.

Finally every unordered triple of distinct nonzero sign-orbits has exactly 3!·2³=48 ordered signed presentations. The action is free: nonzero integer vectors are not fixed by sign reversal, and a permutation stabilizing the presentation up to signs would identify two sign-orbit columns. Division by 48 proves the theorem and its integrality.

The proof uses simultaneous sign-orbits in an integer box. For exponents larger than one, it does not reinterpret the box as a projective space over a larger field.

## 3. Squarefree specialization

At a=1 the four local counts are L=21, U=5, V=7, W=3. Therefore for every integer k≥0,

    N_k = [21^k − 6·7^k − 3·5^k + 14·3^k − 6]/48.

The first values, starting at k=0, are

    0, 0, 4, 150, 3736, 82860, 1771324, 37415490,
    787233136, 16542340440, 347461617844.

In the squarefree case the 21 local row words are precisely the vectors in F_3³ of Hamming weights 0, 2, or 3: 1+12+8=21. The local count 7 is a fixed pair-equality event with the remaining support contained in that pair's support. These are alphabet counts in this enumeration, not a Hamiltonian spectrum, tournament identification, or evidence of an unrelated 7/21 correspondence.

## 4. Rational generating function and minimal recurrence

For F(z)=Σ_{k≥0}N_k z^k,

    F(z) = (1/48)[1/(1−21z) − 6/(1−7z) − 3/(1−5z)
                          + 14/(1−3z) − 6/(1−z)]
         = (4z²+2z³−126z⁴)
           / [(1−z)(1−3z)(1−5z)(1−7z)(1−21z)].

The expanded denominator is

    1−37z+422z²−1982z³+3801z⁴−2205z⁵.

Thus, for k≥5,

    N_k = 37N_(k−1) − 422N_(k−2) + 1982N_(k−3)
          − 3801N_(k−4) + 2205N_(k−5),

with initial values N_0=0, N_1=0, N_2=4, N_3=150, N_4=3736.

The recurrence has minimal order five among homogeneous linear recurrences with constant rational (or real) coefficients, including eventual recurrences. Indeed the partial fractions have five distinct poles and all five residues are nonzero; no denominator factor cancels. Equivalently, a recurrence annihilating this sequence must vanish at the five distinct exponential bases 21,7,5,3,1, by Vandermonde independence of their exponential sequences.

## 5. Density among all triples

For k≥2 the number of all distinct unordered triples of points is

    T_k = binom((3^k−1)/2,3)
        = [(3^k−1)(3^k−3)(3^k−5)]/48.

Consequently

    N_k/T_k
      = [21^k−6·7^k−3·5^k+14·3^k−6]
        / [(3^k−1)(3^k−3)(3^k−5)]
      ~ (7/9)^k.

The balanced-candidate density tends to zero, while the absolute count grows like 21^k/48. No density claim about actual magic hourglasses follows, beyond their being contained in this necessary candidate set at a fixed center. No density is defined for k=0 or k=1, where there are no triples.

## 6. Exact checks and hostile controls

Run `python check.py` and `python -O check.py`. Both produce identical output; no checks depend on assertions.

Two independent finite paths are used:

1. Direct projective-point enumeration, normalized by first nonzero coordinate equal to 1, followed by every unordered distinct triple, for k=0,...,5. Counts are 0,0,4,150,3736,82860. The largest universe has 121 points and 287,980 triples.
2. Direct enumeration of matrices from the 21-word row alphabet for k=0,...,3, checking the nonzero-column count, fixed-pair equality count, all-pair intersection count, and final 48-to-one quotient separately.

General exponent boxes are independently enumerated for the empty vector, each singleton (a) with 1≤a≤7, and (1,2), (2,2), (3,3), (1,1,2), (1,2,3), (2,2,2). Local L,U,V,W counts are checked by exhaustive row enumeration for a=1,...,7. Exact polynomial arithmetic verifies the generating-function numerator and denominator, and the recurrence is checked through k=40.

A singleton-support row is explicitly rejected. The labels (1,0),(0,1),(1,1) pass the filter, but at center root e=65 their positive gaps are 4056,3000,3696, none the sum of the other two. All four triples at k=2 pass the filter, so the count N_2=4 is deliberately not presented as four hourglasses.

Artifacts: check.py, certificate.json, check.out, check-optimized.out. The finite computations validate the independently proved formula; they do not replace the all-exponent inclusion-exclusion proof.

## 7. Scope ledger

Source: the Gaussian index box together with the center's split-prime exponents. Target: triples satisfying coordinatewise repeated absolute maxima. Preserved: every actual signed three-gap relation passes the predicate; inert scaling and coordinate permutations do not change the count. Lost: the numerical Gaussian primes, coordinate phases, gap magnitudes, and the signed additive equation. Restoration sidecar: actual gap values, followed by an exact Schur test. Failure boundary: the e=65 three-gap hostile above.

The formulas are purely combinatorial counts of a necessary arithmetic filter. They do not assert a projective line gives an additive relation, do not enumerate all hourglasses, and do not settle any unrestricted magic-square question.

## 8. Primitive candidates and disjoint ancestry

For a triple at center exponent vector a, let b_i be the maximum of its three absolute i-th coordinates. The common divisor of the seven associated square roots is ∏p_i^(a_i−b_i), after inert/2 stripping. Dividing the roots by that common divisor changes the center to ∏p_i^b_i and leaves the Gaussian indices unchanged. Coordinates with b_i=0 are deleted.

This motivates a definition on candidate triples, whether or not their gaps actually satisfy a Schur relation: the triple is primitive if b_i=a_i in every coordinate. Together with balance, this says at least two absolute entries in coordinate i equal a_i. It means precisely that no common center factor can be divided from all associated roots, not that the numerical gaps are pairwise coprime.

Let M(a_1,...,a_k) count primitive balanced unordered distinct triples. Set M(())=0. For k≥1,

    M(a_1,...,a_k)
      = [∏(24a_i−4) − 3·4^k − 6∏(4a_i+2) + 14·2^k]/48.

Direct local proof: a row with repeated absolute maximum exactly a has 12(2a−1)+8=24a−4 choices. With one specified zero column, the two other entries must each have absolute value a, giving four choices; two zero columns are impossible. For a fixed global-sign equality between two columns, their common magnitude must be a and the third entry can be any integer from −a through a, giving 4a+2 choices. Removing a zero third column subtracts 2 choices per coordinate. Three columns equal up to two fixed signs give 2 choices per coordinate. Therefore the number of nonzero-column matrices is ∏(24a_i−4)−3·4^k; each pair event has size 2[∏(4a_i+2)−2^k]; and the common intersection has size 4·2^k. Inclusion-exclusion and division by 48 yield the formula. The restriction k≥1 is essential in these zero-column exclusions.

More fundamentally, every total candidate triple has exactly one maxima vector b with 0≤b_i≤a_i. Deleting its zero coordinates identifies it with a primitive candidate at b. Distinctness and nonzero columns survive because the deleted coordinates are zero in every column. Thus this is a disjoint, reversible combinatorial and arithmetic partition:

    N(a) = Σ_{0≤b_i≤a_i} M(b with all zero coordinates deleted).

The all-zero vector contributes M(())=0. This identity is not merely a formal equality of closed formulas: the maxima vector is the unique ancestry label, and the inverse map inserts the missing zero coordinates and rescales the roots. As an algebraic consistency check, taking a first difference in every exponent gives ΔL(a)=24a−4, ΔU(a)=4, ΔV(a)=4a+2, and ΔW(a)=2, recovering the primitive formula.

## 9. Squarefree primitive sequence, generating function, and densities

Write M_k=M(1,...,1), with M_0=0. For k≥1,

    M_k = [20^k − 6·6^k − 3·4^k + 14·2^k]/48.

The initial values are

    0, 0, 4, 138, 3160, 65640, 1327264, 26630688,
    533119360, 10665390720, 213325709824.

At k=0 the exponential expression would give 1/8, so the exceptional definition M_0=0 cannot be omitted. The exact hereditary partition becomes the binomial transform

    N_k = Σ_{j=0}^k binom(k,j) M_j.

For example, 150=138+3·4 and 3736=3160+4·138+6·4. At k=3 the 12 nonprimitive candidates are the four candidates on each of the three coordinate lines; no candidate triple lies in an intersection consisting of only one point, so these ancestry contributions do not overlap.

For G(z)=Σ_{k≥0}M_k z^k the corrected generating function is

    G(z) = (1/48)[1/(1−20z) − 6/(1−6z) − 3/(1−4z)
                              + 14/(1−2z) − 6]
         = (4z²+10z³−120z⁴)
           / [(1−2z)(1−4z)(1−6z)(1−20z)].

The denominator expands to 1−32z+284z²−928z³+960z⁴. Hence

    M_k = 32M_(k−1) − 284M_(k−2) + 928M_(k−3) − 960M_(k−4)

holds for k≥5. It is the minimal eventual homogeneous constant-coefficient recurrence, of order four: all four distinct exponential bases have nonzero coefficients. The word “eventual” and the start index are load-bearing. Applying the recurrence at k=4 with M_0=0 would predict 3280, whereas M_4=3160; the exact residual is −120, matching the numerator's z⁴ coefficient. An order-four recurrence starting at its order index is not being claimed.

The binomial transform also gives F(z)=G(z/(1−z))/(1−z), consistent with both rational functions.

Since M_k~20^k/48,

    M_k/T_k ~ (20/27)^k,
    M_k/N_k ~ (20/21)^k,

where T_k counts all projective triples and N_k counts all balanced candidates. The denominators in these density statements are different and are not interchangeable. In particular the second density says that an increasing fraction of balanced candidates descends from a proper subcenter. Neither density asserts that any candidate is an actual hourglass.

## 10. Additional primitive checks

The same independent direct enumerations in §6 also count primitive candidates by requiring the coordinate maxima to reach every exponent bound. They verify the primitive formulas on all stated exponent boxes, and the direct squarefree counts through k=5 are 0,0,4,138,3160,65640. The exact maxima partition is checked on every tested exponent vector; the binomial transform is checked through k=40. Exact polynomial arithmetic checks the primitive numerator and denominator, and the recurrence through k=40. The failed k=4 application is retained as an explicit hostile control with residual −120.

All formulas remain counts of necessary valuation candidates. The common-root ancestry is exact, but it does not supply the missing signed numerical gap relation.
