# Independent audit: tropical candidate-triple count

Status: PASS / PROVED (combinatorial formula, OGF, minimal recurrence, asymptotic). Independent finite-exact computations corroborate, but do not replace, the proof.

## Object and scope

For k>=1 let the point set be the nonzero ternary column vectors in {-1,0,1}^k modulo simultaneous sign, identified with PG(k-1,3). Count unordered triples of distinct points such that each coordinate has either zero, two, or three nonzero entries across the three points. This is exactly the squarefree maximum-twice allocation filter. Define the k=0 point set to be empty, without invoking a nonempty negative-dimensional projective space.

These are necessary candidates for arithmetic Schur triples of gaps, not actual Schur triples. No choice of split primes occurs in this count. Numerical gap values remain required. There is no Hamilton-spectrum consequence.

## Inclusion-exclusion audit

Temporarily order the three columns and choose their signs. A row is any element of {-1,0,1}^3 other than a weight-one vector. There are 1+3*4+8=21 rows.

First exclude zero columns. For a fixed zero column the other two entries must be both zero or both nonzero, giving 5 rows. For two specified zero columns, the third must also be zero. Inclusion-exclusion therefore gives

    A_k = 21^k - 3*5^k + 2.

This identity includes k=0: all named zero-column events coincide on the single empty matrix, leaving zero matrices.

Within these nonzero-column matrices, let E_ij mean that columns i and j represent the same projective point, that is, they are equal or negatives. For a fixed sign relation, if their common entry is zero then the remaining entry must be zero; otherwise it is arbitrary. There are 1+2*3=7 row choices. Requiring the third column nonzero removes 3^k choices; once the third is nonzero, the common column cannot be identically zero because the row filter would force the third to be zero too. Hence

    |E_ij| = 2*(7^k - 3^k).

The factor two is legitimate because the related columns are nonzero, so the two sign choices are disjoint.

An intersection of any two of E_12,E_13,E_23 is the same rank-one locus, in which all three columns represent the same nonzero projective point. Choose a nonzero signed first vector (3^k-1 choices), then either sign for each remaining column, giving 4*(3^k-1). The triple intersection is exactly that same locus. Thus adding the three pair intersections and subtracting the triple intersection restores 8*(3^k-1).

The resulting signed ordered count is

    21^k - 6*7^k - 3*5^k + 14*3^k - 6.

Each unordered triple of distinct projective points has exactly 3!*2^3=48 signed ordered lifts. The action is free: a fixed lift would force the permutation to fix each distinct projective column, then each sign to be positive because its column is nonzero. Consequently

    N_k = (21^k - 6*7^k - 3*5^k + 14*3^k - 6)/48.

There is no division-by-orbit-size error at k=0 or k=1 because the admissible set is empty. Directly, N_0=N_1=0. The prefix is

    0, 0, 4, 150, 3736, 82860, 1771324, 37415490, 787233136.

## OGF and minimal recurrence

Summing the five geometric sequences gives

    F(z) = sum_(k>=0) N_k z^k
         = [1/(1-21z) - 6/(1-7z) - 3/(1-5z)
            + 14/(1-3z) - 6/(1-z)]/48
         = (4z^2 + 2z^3 - 126z^4)
           / [(1-z)(1-3z)(1-5z)(1-7z)(1-21z)].

The denominator expands to

    1 - 37z + 422z^2 - 1982z^3 + 3801z^4 - 2205z^5.

Therefore, for k>=5,

    N_k = 37N_(k-1) - 422N_(k-2) + 1982N_(k-3)
          - 3801N_(k-4) + 2205N_(k-5),

with initial values N_0=0,N_1=0,N_2=4,N_3=150,N_4=3736.

All five partial-fraction coefficients are nonzero and the five poles are distinct. No denominator factor cancels. Thus order five is minimal for a homogeneous constant-coefficient linear recurrence over Q, including an eventual such recurrence. Equivalently, linear independence of the sequences 1^k,3^k,5^k,7^k,21^k forces any annihilating characteristic polynomial to vanish at all five bases. This minimality claim does not address polynomial-coefficient or inhomogeneous recurrences.

## Density

For k>=2, the total unordered distinct projective triples are

    T_k = binom((3^k-1)/2,3)
        = (3^k-1)(3^k-3)(3^k-5)/48.

The candidate density is exactly

    N_k/T_k = [21^k - 6*7^k - 3*5^k + 14*3^k - 6]
              /[(3^k-1)(3^k-3)(3^k-5)].

Since N_k~21^k/48 and T_k~27^k/48,

    N_k/T_k ~ (7/9)^k.

More precisely it is (7/9)^k*(1+O(3^(-k))). No density is assigned at k=0 or k=1, where there are no triples. This exponential sparsity is for combinatorial necessary candidates, not the probability that arithmetic gap addition holds.

## Independent exact checks

Reproduce:

    python check.py > normal.json
    python -O check.py > optimized.json

The independent implementation does not import the producer's code. It dynamically counts row matrices while tracking nine event bits: three zero-column events and six signed equality events. Each additional row intersects the persistent event mask with that row's event mask. Mask zero counts exactly the nonzero, projectively distinct columns. This checks the formula for every k=0..40, with exact arbitrary-precision integers.

A separate exhaustive enumeration of all 21^k signed row matrices for k=0..3 agrees. Polynomial convolution verifies the OGF numerator and recurrence against all 41 DP terms. Hostile controls verify that zero-column exclusion alone wrongly retains eight signed matrices at k=1, while the proper distinct-point count is zero. The explicit RuntimeError validity guard is tested and remains active under -O. Normal and optimized output files agree byte-for-byte.

No new literature or priority claim is made. All writes were confined to this audit directory.

## General exponent-box extension

Status: PASS / PROVED. Let a=(a_1,...,a_k) be nonnegative integer coordinate bounds, and let the index objects be nonzero vectors in product_i[-a_i,a_i] modulo **global sign only**. For actual retained center primes one normally uses a_i>=1; a zero bound is a redundant coordinate. Count unordered triples of distinct sign classes whose maximum absolute entry in each row is attained at least twice.

For one coordinate of bound a, the exact row counts are:

- L(a)=12a²+8a+1=(2a+1)(6a+1): all allowed rows. At positive maximum r, exactly two maximizers give 3*4*(2r-1) rows and three maximizers give 8; sum 24r-4 over r=1..a and include the zero row.
- U(a)=4a+1: one specified column zero. The other entries are both zero or have the same positive absolute value, with four sign choices per magnitude.
- V(a)=2a²+4a+1: one specified signed equality x=+y or x=-y. If x is zero then z must be zero. If |x|=r>0, there are two choices for x and 2r+1 choices for z, since |z|<=r. Sum 1+2*sum_(r=1)^a(2r+1).
- W(a)=2a+1: the same signed pair with third column zero; also the number of choices for a base coordinate when all three columns have specified signs relative to it.

Write products over coordinates as L,U,V,W. Excluding zero columns gives L-3U+2. Each signed-class equality event has size 2(V-W). Each intersection of two such pair events, and their triple intersection, has size 4(W-1). Thus the valid signed ordered count is

    L - 3U - 6V + 14W - 6,

and the same free group action gives

    N(a) = [product L(a_i) - 3 product U(a_i) - 6 product V(a_i)
            + 14 product W(a_i) - 6]/48.

The empty coordinate tuple gives zero, and a single coordinate of any bound gives zero: three distinct absolute index values have a unique maximum. Substitution in the polynomial formula verifies this identically. Setting every a_i=1 recovers the squarefree formula.

Important terminology guard: the pair events mean **equal up to sign**, not arbitrary proportionality. For exponent bounds greater than one, vectors (1,0) and (2,0) are distinct sign classes and can index distinct gaps, although they are proportional over Q. There is no projective-space claim for these general boxes.

Three independently specified routes, the general formula, row-event DP, and brute-force triples of canonical sign representatives, agree for boxes (), (1), (2), (3), (1,2), (2,2), (1,1,2), and (2,2,2). Their counts are respectively 0,0,0,0,15,54,492,5098. Both normal and -O runs agree. These remain necessary allocation-filter counts, not arithmetic Schur-relation counts.

## Primitive refinement

Status: PASS / PROVED with an essential k=0/recurrence boundary correction.

The exact common-root divisor audited in the companion prime-support result has p_i exponent a_i-max_r|j_i^(r)|. Therefore a triple at stripped center is primitive exactly when each coordinate maximum equals a_i. Combined with the tropical filter, at least two of the three entries in that row must have absolute value a_i. This is an arithmetic root-gcd statement about actual indexed gap triples, and a combinatorial condition for candidate index triples; it does not assert a Schur relation.

For positive coordinate bounds a_i and k>=1:

- Allowed primitive rows number 24a_i-4, the positive-max row count at r=a_i.
- With one column fixed zero, rows number 4.
- With a specified signed pair equal, rows number 4a_i+2.
- With that pair equal and the third zero, rows number 2.
- Two zero columns are impossible. The pair-equality intersections are the common sign-class locus, of size 4*2^k.

Consequently

    P(a) = [product_i(24a_i-4) - 3*4^k
            - 6 product_i(4a_i+2) + 14*2^k]/48,  k>=1.

Set P(empty)=0 separately. The displayed raw formula would give 1/8 at k=0, so blindly using it there is incorrect. For squarefree support, the valid all-k expression is

    P_k = [20^k - 6*6^k - 3*4^k + 14*2^k - 6 delta_(k,0)]/48.

Each tropical triple has a unique coordinate maximum vector b with 0<=b_i<=a_i. Delete its zero coordinates; it is then primitive at that reduced box. This gives the exact partition

    N(a) = sum_(0<=b_i<=a_i) P(nonzero coordinates of b).

No Gaussian prime changes or speculative closure is needed for this partition. In the squarefree case it is the binomial transform N_k=sum_(j=0)^k binom(k,j)P_j. Hereditary division of actual gaps by the square of the common root divisor realizes the same reduction and preserves any actual Schur relation, if one exists.

The squarefree primitive prefix is

    0, 0, 4, 138, 3160, 65640, 1327264, 26630688, ...

Its OGF is

    P(z) = (4z^2+10z^3-120z^4)
           /(1-32z+284z^2-928z^3+960z^4),

with denominator (1-20z)(1-6z)(1-4z)(1-2z). Thus

    P_k = 32P_(k-1)-284P_(k-2)+928P_(k-3)-960P_(k-4)

holds for k>=5. It does **not** hold at k=4: the residual there is -120. The numerator has the same degree as the denominator because of the k=0 correction. All four poles survive, so the minimal **eventual** homogeneous constant-coefficient recurrence order is four. If a convention requires an order-r recurrence to start at index r, this order-four recurrence does not qualify without the separately stated exceptional boundary; including the finitely supported zero-index correction gives an annihilator with an additional zero root.

Density denominators must be named:

- Among all tropical candidates, P_k/N_k~(20/21)^k.
- Among all unordered distinct squarefree-index triples, P_k/T_k~(20/27)^k.
- Among primitive unordered distinct index triples, the density is ~(10/13)^k. For clarity, that denominator is [26^k-9*8^k+23*2^k]/48 for k>=1: all rows are nonzero (26 choices), then remove zero columns and signed duplicate classes by the same inclusion-exclusion. Its leading base is 26, not 27. A density is defined only when this denominator is nonzero, hence k>=2.

The independent script compares primitive formula with direct canonical-vector triples for boxes (), (1), (2), (3), (1,2), (2,2), (1,1,1), (1,1,2), and (2,2,2). It checks the full maxima-vector partition for each, the squarefree binomial transform through k=40, and the primitive OGF convolution through k=40. It explicitly rejects starting the recurrence at k=4. Both normal and -O runs pass byte-identically.
