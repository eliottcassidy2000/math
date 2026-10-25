# Finite-field Parker classification via a six-factor slice

Status: PROVED using the cited Weil bound; algebra independently audited. Complete bounded census independently audited with a separate tuple-field implementation and Frobenius/gcd irreducibility checks.

## Scope and inheritance

A Parker field has no 3 by 3 magic square with nine distinct entries that are squares; zero is permitted. Work over the FIELD F_(p^r), not the ring Z/p^r Z. The standard three-parameter magic plane, normalization of nonzero square center to one, and center-zero constructions are inherited mechanisms (Cain, 2019, Sections 5–7). The hostile examples are characteristic two and F_27, which is Parker although Z/27 Z is not. Our niche is the collision arrangement, and our wildcard is the use of a conic slice to remove two character constraints before applying Weil.

Live concept board: magic plane / square-scaling quotient / conic slice / collision arrangement / squarefree character sums / extension-field versus residue-ring arithmetic. No implication about the integer problem is asserted.

## 1. Universal magic plane and collision directions

Every magic square over a commutative ring is uniquely of the form

    c+u      c-u-v     c+v
    c-u+v    c         c+u-v
    c-v      c+u+v     c-u.

Indeed, the sum of the middle row, middle column and two diagonals is 4T. The sum of all nine entries is 3T, and the former counts the center three extra times, so T=3c. Each opposite pair therefore sums to 2c. Set u=top-left minus c and v=top-right minus c; the remaining entries follow from the first row and opposite pairs. Direct addition proves the converse. No division was used.

The nine offsets are exactly su+tv with s,t in {-1,0,1}. In characteristic two opposite offsets coincide; thus every characteristic-two field is Parker, independently of whether all its elements are squares.

In odd characteristic, with v fixed and nonzero, distinctness is equivalent to

    u not in {0, v, -v, 2v, -2v, v/2, -v/2}.

To prove this, a collision means au+bv=0 with a,b in {-2,-1,0,1,2}, not both zero. If a=0 then v=0, excluded. Otherwise division gives exactly the listed possibilities. Conversely each listed ratio realizes a difference of two grid indices. Repeated listed values in small characteristic only reduce the forbidden set.

Projectively the forbidden directions are 0, infinity, ±1, ±2, ±1/2. There are four lines in characteristic three, six in characteristic five, and eight in characteristic at least seven. This explains some small-characteristic geometry, not the square predicate by itself.

## 2. Character-sum input

CITED standard Weil estimate: for odd q, quadratic character chi with chi(0)=0, and squarefree f of degree d>=1,

    |sum_x chi(f(x))| <= (d-1) sqrt(q).

A primary researcher-authored formulation is Theorem 4.3 in L. Babai's notes, https://people.cs.uchicago.edu/~laci/HANDOUTS/weil.pdf . A product of distinct affine linear factors is squarefree even in small characteristic; no restriction d < characteristic is needed here.

We also use the elementary exact identity, for a nonzero leading coefficient A and distinct roots r,s,

    sum_x chi(A(x-r)(x-s)) = -chi(A).

For proof, translate and scale to x(x-1). The equation y^2=x(x-1) is equivalent to (2x-1-2y)(2x-1+2y)=1, which has q-1 solutions. Subtracting q yields the identity.

## 3. Selecting a good conic slice

For q>25 odd, there exists v outside {0,±1,±2} such that 1+v and 1-v are nonzero squares.

Proof. Put

    W(v) = (1+chi(1+v))(1+chi(1-v))/4.

Its sum is (q-chi(-1))/4 by the quadratic identity and vanishing of linear character sums. Away from v=±1 it is exactly the indicator that both factors are nonzero squares. At the two roots W is at most 1/2 each. Hence the number of v making both nonzero squares is at least (q-5)/4. Removing the at most five values {0,±1,±2} leaves at least (q-25)/4>0. The estimate deliberately overcounts removals; in characteristic three the exceptional set is smaller.

## 4. Six-factor theorem

THEOREM (uniform effective eventual non-Parker bound). Every finite field of odd order q>=14400 is non-Parker. More generally, for each v selected as in Section 3, the number of u giving a distinct magic square whose entries are all nonzero squares is at least

    (q - 15 - 114 sqrt(q))/64 - 10.

Proof. The center is c=1, and 1±v are already nonzero squares. The remaining six entries, viewed as affine functions of u, are

    1+u, 1-u, 1+u+v, 1+u-v, 1-u+v, 1-u-v.

Their roots are ±1, ±(1+v), ±(1-v), which are pairwise distinct because v is outside {0,±1,±2} and the characteristic is odd. Define

    W(u) = 2^(-6) product_over_six_factors (1+chi(L(u))).

Expand and sum. The empty subset contributes q. All six singleton character sums vanish. Each of the 15 two-factor terms is at least -1 by the exact quadratic identity. For subset sizes 3,4,5,6, the total Weil error is at most

    [C(6,3)*2 + C(6,4)*3 + C(6,5)*4 + C(6,6)*5] sqrt(q)
    = (40+45+24+5) sqrt(q) = 114 sqrt(q).

Thus sum W >= (q-15-114 sqrt(q))/64. At each of the six distinct roots W is at most 1/2; removing roots costs at most 3. At all other points W is exactly the indicator of six nonzero square entries. Removing the collision set from Section 1 costs at most seven more. The displayed lower bound follows.

At q=14400, q-114 sqrt(q)-655 = 14400-13680-655 = 65>0. The expression increases for q>=14400, so at least one good u exists. All eight magic sums agree by Section 1. This proves the theorem.

The bound is conservative, not claimed optimal. It gives an O(q) field-operation constructive algorithm: enumerate v to find a conic slice, then enumerate u on that slice; for q>=14400 every admissible v succeeds. Square membership can be precomputed from x^2 in O(q) field operations. Below threshold failure on one slice is not an obstruction certificate.

## 5. Complete finite reduction and exact certificate protocol

For odd q<14400 the accompanying verifier enumerates every prime power, including extensions. It searches the center-zero normalized plane c=0,u=1, then the complete center-one plane c=1 with u,v constrained only by 1±u and 1±v being squares. These are necessary constraints, so they lose no solution.

Completeness of normalization: if c is a nonzero square, multiply all entries by c^(-1), itself a square, to get c=1. If c=0, distinctness forces the top-left entry u to be a nonzero square; multiply by u^(-1) to get u=1. Conversely square scaling preserves magic sums, square membership and distinctness. Hence absence in BOTH searches proves Parker.

Positive certificates include the field characteristic, extension degree, irreducible modulus in low-coefficient-first order, and the nine entries. An independent direct check verifies all nine entries are squares, pairwise distinct, and all eight sums agree. Negative cases are certified by complete enumeration, not a heuristic line search. Moduli are proved irreducible by trial division by every monic polynomial up to half their degree. Every reducible polynomial has such a factor.

Combining the above uniform theorem with a successful complete census supplies a full classification; the census's actual output and hash must be read before asserting that classification.

## 6. Literature status and boundaries

Cain's 2019 paper states conjectures, not a complete classification. Its conjecture on strict prime powers must be read with odd characteristic: characteristic two supplies infinitely many Parker fields. A newer author-maintained repository, https://github.com/jiji7879/magic_squares_fields , advertises a complete odd-characteristic classification and bounds 77 and 553736 in the two residue classes mod 4. The accompanying paper was located and inspected at commit 0677087714b7dd446095db3e7020c813e6eef7fb, dated 2026-09-29: https://github.com/jiji7879/magic_squares_fields/blob/0677087714b7dd446095db3e7020c813e6eef7fb/paper/main.tex . It explicitly gives the full sixteen-field odd classification and an eight-factor character bound. Its computational census has not been independently rerun here. Our proof and census stand independently and give a smaller sufficient bound than that paper; no priority claim about the classification is made.

Earlier primary sources include Giancarlo Labruna's 2018 thesis, https://digitalcommons.montclair.edu/etd/138/ , and Cain's https://arxiv.org/html/1908.03236v2 . Hengeveld, Labruna and Li, “Magic squares of squares over a finite field,” Contemporary Mathematics 773 (2021), 111–122, DOI https://doi.org/10.1090/conm/773/15536, has a primary institutional record at https://researchwith.montclair.edu/en/publications/magic-squares-of-squares-over-a-finite-field/ . Its abstract supplies infinite families, not the full classification; its full theorem content is not used here.

A field embedding F_q -> F_(q^m) transports every non-Parker witness. This is an implication along divisibility of extension degrees, not numerical monotonicity in q. Conversely, reduction of a ring witness modulo p can collapse distinct entries; no Parker obstruction propagates from F_p to Z/p^rZ without tracking collisions. The bounded interweaving of Parker and non-Parker orders is therefore not evidence against eventual odd-characteristic existence, nor evidence for an integer square.


## 7. Complete classification and center-type boundary

FINITE-EXACT: census.json enumerates all 1730 odd prime powers below 14400. The only negative fields have orders

    3, 5, 7, 9, 11, 13, 17, 19, 23, 25, 27, 31, 43, 47, 67, 243.

Together with Section 4 and the characteristic-two obstruction, this proves that the complete Parker fields are all characteristic-two fields and these sixteen odd fields. This reproduces the classification reported by Lai, with an independent proof reduction and code.

The smallest non-Parker field is F_29. The smallest field admitting a witness with NONZERO CENTER is F_59. These are different exact questions. Counts below are labeled normalized parameter counts, not dihedral equivalence classes and not counts of root choices:

| q | c=0,u=1 | c=1 |
|---|---:|---:|
|29|4|0|
|37|4|0|
|41|2|0|
|49|2|0|
|53|4|0|
|59|0|8|
|61|8|0|
|67|0|0|
|71|0|8|

For odd q, the total number of labeled distinct magic squares of square ENTRIES is ((q-1)/2) times the sum of these two counts. Each square-scaling orbit has exactly one representative in its respective normalization. The action is free because some entry is nonzero. A brute triple enumeration independently checks this identity for p=3,5,7,11,17,29.

The c=0 branch cannot be an integer nonnegative distinct-square witness: magic total zero forces every nonnegative entry zero. Furthermore, finite-field scaling does not preserve integer height, nonnegativity, carries or chosen roots. Neither the threshold 29 nor 59 predicts an integer solution height.

## 8. Two structural corollaries

PROVED (shorter zero-center threshold). If q is odd, q=1 mod 4 and q>=65, there is a zero-center distinct square witness. Set c=0,u=1 and seek t=v such that t-1,t,t+1 are nonzero squares. Since -1 is a square, all nine entries are then squares. The weighted product of the three quadratic-character indicators has total at least (q-3-2 sqrt(q))/8: the three linear terms vanish, each quadratic term is -1, and the cubic term has absolute value at most 2 sqrt(q). At the three roots the total spurious weight is at most 3/2. The remaining collision parameters are contained in {±2,±1/2}, costing at most four more. The number of valid t is at least (q-47-2 sqrt(q))/8, positive for q>=65. This is a sufficient bound, not the least non-Parker order.

PROVED (characteristic-three plane criterion). In characteristic three, a magic square has nine distinct entries exactly when u,v are F_3-linearly independent, and its set of entries is c+span_F3(u,v). Thus F_(3^r) is non-Parker iff its set of squares contains an affine F_3-plane of dimension two. Every such plane produces exactly 9*|GL_2(F_3)|=9*48=432 labeled magic squares: choose the center among its nine points, then an ordered basis for its translation plane. No two choices give the same labeled matrix, since center,u,v are recoverable.

PROVED (square-order elementary construction). Every odd square-order field F_(s²) with s>=9 is non-Parker. All elements of the subfield F_s are squares in F_(s²): for a nonzero subfield element, a^((s²-1)/2)=(a^(s-1))^((s+1)/2)=1. Select u=1 and v outside the at most seven forbidden collision values in F_s. As s>=9 such v exists, so the magic grid with c=0 has nine distinct subfield entries, all squares upstairs. This proof has no character-sum dependency.


## 9. Sharp characteristic-specific eventual behavior

The proved classification gives the following exact exponent sets (r>=1):

- p=2: all r are Parker.
- p=3: precisely r in {1,2,3,5} are Parker. In particular F81 is non-Parker but the larger F243 is Parker; every r>=6 is non-Parker.
- p=5: precisely r in {1,2} are Parker; every r>=3 is non-Parker.
- p in {7,11,13,17,19,23,31,43,47,67}: precisely r=1 is Parker.
- Every other odd prime p: no r is Parker.

Consequently all odd finite fields with q>243 are non-Parker, and the endpoint 243 cannot be lowered as a universal numerical cutoff. For prime fields only, every p>67 is non-Parker. These sharp classification consequences use the bounded census; the stand-alone analytic proof only supplies the sufficient threshold 14400.

## 10. Audit record

An independent implementation checked the full 1730-order universe, all 1714 positive witnesses, all 16 negative complete searches, all 71 small-count records, and every modulus using Frobenius/gcd irreducibility rather than this producer's trial-division algorithm. The external audit also checked collision and six-root equivalences over small prime and extension fields, including characteristics 3 and 5. The census hash is 37bfadf8dc7b223a82699e6c29e8207d58fab220acf8c8cce8b03e3130fa5150. The independently audited first nonzero-center field is 59.
