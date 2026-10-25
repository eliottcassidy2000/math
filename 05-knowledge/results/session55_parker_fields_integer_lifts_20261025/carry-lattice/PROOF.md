# Carry syndromes, integer height, and faithful modular lifting

Status: PROVED elementary reductions, with FINITE-EXACT controls. This packet does not produce a 3-by-3 magic square of nine distinct integer squares, prove that none exists, or establish an external novelty claim.

## 1. The integral magic lattice splits, in every characteristic

Index cells row-major from 0 to 8. Put

M(c,a,b) = (a, 3c-a-b, b, c-a+b, c, c+a-b, 2c-b, a+b-c, 2c-a).

For a vector v define E(v)=v-M(v4,v0,v2). Its coordinates 4,0,2 vanish. Write S(v) for the other six coordinates, in order 1,3,5,6,7,8.

**Theorem 1.** Over any commutative ring with identity, v has all eight row, column, and diagonal sums equal if and only if S(v)=0. The common sum is 3v4. Over the integers, S is a split surjection Z^9 to Z^6, and its kernel is the magic lattice M(Z^3). Thus the quotient is free: there is no hidden torsion or bad-prime exception in the linear magic equations themselves.

Proof. For a magic grid, add the three lines through the center consisting of its two diagonals and middle column. Their sum is the top row plus bottom row plus three copies of the center. Consequently 3T=2T+3c and T=3c without division. Opposite cells sum to 2c, and the top row and outside columns then give the displayed formula. Conversely all eight sums of M are 3c by direct addition. On the six non-pivot coordinate basis vectors, S is the identity. Finally every v is uniquely M(v4,v0,v2)+E(v). These are integer identities, so they remain valid modulo any modulus and in every field. In particular there are exactly m^3 magic arrays over Z/mZ, before square/distinctness restrictions.

The split linear lattice is not the nonlinear variety of square-root vectors. The quadratic constraints, collisions, and integer heights remain load-bearing.

## 2. The exact carry equation

Fix m≥2. Let r_i be canonical residues in [0,m), forming a modular magic square, and suppose a chosen integer root u_i satisfies u_i^2≡r_i (mod m). Define

h_i=(u_i^2-r_i)/m, and q=S(r)/m in Z^6.

**Theorem 2.** The chosen actual squares u_i^2 are an integer magic square if and only if

S(h)=-q.

Their distinctness and positivity must be checked separately. Proof: linearity gives S(u^2)=S(r)+mS(h)=m(q+S(h)). This tests precisely the lost coordinate: the root carries must cancel the residual-grid carry syndrome. It does not require that every carry vanish.

The six coordinates are an exact independent linear residual system. This is not a claim that six tests are optimal among arbitrary nonlinear tests.

A no-wrap magic residue grid has q=0, but the actual root squares still require S(h)=0. Replacing roots by least modular roots, or normalizing the center by a square unit, need not preserve that property or height.

## 3. Explicit hostile: no-wrap without integer squareness

For p=59, the canonical residue grid

28  3 17
 5 16 27
15 29  4

has nine distinct quadratic residues and all eight ordinary line sums equal to 48<59. In this grid exactly two entries, 4 and 16, are integer squares; the other seven are not. The certificate records explicit roots and their nonzero carry syndrome. Thus this satisfies the weak residue/no-wrap interpretation but not the integer-square target. No minimality claim is made in this packet; the separate bounded-search packet audits that question.

## 4. Height-faithful lifting and the minimum-modulus question

**Theorem 3.** Given nine distinct actual nonnegative integer squares, if their eight line sums are congruent modulo m and lie in one real interval of length strictly less than m, then all line sums are equal as integers. Conversely an integer magic square passes these tests for every m larger than its common line sum.

Proof: two unequal congruent integers differ by at least m. The endpoint version uses an interval with only one endpoint included, such as [0,m) or (0,m]. A closed interval [0,m] needs an additional exclusion of simultaneous 0 and m. For nine distinct nonnegative entries each line sum is positive, so the 'at most m' version is safe: line sums lie in (0,m].

**Corollary 3a.** Any nonnegative integer magic square of distinct squares has center e²>0, common sum 3e², and every entry strictly between 0 and 2e². Zero cannot occur: its opposite entry would be 2e², which cannot be an integer square for positive integer e. Therefore nonnegative and positive target conventions agree.

**Corollary 3b.** Conditional on existence, the least possible quotient-ring modulus for a no-wrap integer solution with sums at most the modulus is exactly 3e_min², where e_min is the least center root of any integer solution. Among prime moduli the least is the first prime strictly greater than 3e_min². Equality with a prime could only occur at e=1,p=3, impossible for nine distinct entries. Distinctness modulo these moduli is automatic because entries lie below 2e². Thus searching for the smallest faithful prime-field no-wrap solution is an equivalent height optimization of the original integer problem, not a weaker finite-field problem. Without existence there is no such least modulus.

For an integer solution, merely retaining distinctness under reduction requires modulus greater than its entry range, which can be much smaller than its common sum. Do not confuse the two thresholds.

## 5. Field order is not a modulus for integer height

Every unital ring map Z→F_(p^r) factors through F_p. If r>1, using p^r as a numerical modulus does not model this map. An extension-field entry outside F_p has no integer antecedent at all. A labeling of field elements by integers 0,...,p^r-1 is not an addition-preserving encoding.

Consequently increasing extension degree in characteristic two never makes reduction injective on nine distinct integers. Infinitely many Parker extension fields alone cannot exclude an integer square. In contrast, if unbounded integers m give Parker rings Z/mZ, an integer solution is impossible: all sufficiently large m preserve its nine distinct entries and equations.

## 6. Any fixed finite congruence bank misses integer carries

For any M≥2, let u_i=1+iM for i=1,...,9 in row-major order. The nine actual squares are positive and distinct and all equal 1 modulo M. They pass the magic equations modulo every divisor of M, but their integer rows have different sums. Thus any finite bank of congruence equalities alone admits false positives. This does not undermine a height-bounded sieve followed by exact verification; nor does it claim that this example passes local distinctness requirements.

For a valid modular square with distinct residues, adding independent multiples of M to chosen roots also preserves all modular equations and distinct residues, while the carry equation still needs checking. The missing datum is not another field name but the actual root-height/carry vector.

## Verification

Run `python3 -B verify.py` and `python3 -O -B verify.py` in this directory. Standard library only. The checker validates the integral splitting and equivalence on all 19,683 vectors in {-1,0,1}^9, verifies the primitive syndrome matrix, all parameterized arrays over moduli 2,3,4,5,7,9, the exact carry identity, the p59 hostile and five finite-sieve hostiles. Proofs establish all-size statements; finite tests are controls, not their replacement.

## Provenance and connections

User-requested inspiration: O. M. Cain, *Gaussian Integers, Rings, Finite Fields, and the Magic Square of Squares*, arXiv:1908.03236v2, https://arxiv.org/html/1908.03236v2 (2019). The displayed linear parametrization and modular/integer distinction are elementary independently derived here; no external priority claim. The packet supplies a loss ledger for passing from integer squares to modular square residues and links the repository's retained-coordinate principle to an explicit split exact sequence and six-coordinate carry equation.
