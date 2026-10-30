# The entire norm defect is one three-bit syndrome

**PROVED within the fixed seven-line XOR model.** This is a quantitative sequel to `NORM_SIGNS.md`; it does not test arbitrary bilinear algebras or claim the factoring engine tolerates incorrect identities.

Let z=(A,B,C,D,E,F,G) be any of the 128 Fano line orientations and let

    s(z)=(D+E+F+G+1, B+C+D+E+1, A+C+D+F+1) in F_2^3.

The norm composes iff s(z)=0. The 112 rectangles in the norm coefficient criterion are indexed by nonzero t, an unordered pair {a,a+t}, and an unordered pair {b,b+t}. Each rectangle defines a necessary affine parity equation on z.

## D1. An exact local defect count

Among these 112 equations:
- 70 are tautologies forced by normalized cyclic line multiplication;
- each of the seven nonzero linear functionals of the three syndrome bits occurs exactly six times.

This is a finite symbolic identity: use the seven lines listed in `NORM_SIGNS.md`, replace each product by its cyclic line bit, and collect identical equations. The checker independently enumerates all 112 equations and verifies precisely this multiplicity distribution, including right-hand sides.

Every nonzero vector in F_2^3 evaluates to1 under exactly four of its seven nonzero linear functionals. Consequently every invalid orientation violates **exactly 24** rectangles; all 16 valid orientations violate none. This is not merely an average across the 112 invalid tables.

A uniformly selected rectangle therefore detects every invalid promised-line table with probability `24/112=3/14`. Checking r independent uniform rectangles has one-sided failure probability `(11/14)^r`. On a violated rectangle, the two-term inputs in `NORM_SIGNS.md` give norm product4 but output norm0 or8. These are exact arithmetic checks, with no numerical tolerance. One may instead avoid the tautological rectangles: sampling uniformly from the42 nontrivial rectangles detects an invalid table with probability4/7.

The deterministic three syndrome equations already certify the whole table when its seven line bits are available. Random tests are useful only for a black-box multiplication interface; they are not claimed superior to reading all seven signs.

## D2. An exact mean-square energy gap

Put `Delta_z(x,y)=N(xy)-N(x)N(y)`. Let all coordinates of x and y be independent, centered, variance-one random variables, with x independent of y. For this expectation no Gaussian assumption and no fourth-moment assumption is needed: each mixed monomial uses distinct coordinates on either side.

Distinct mixed monomials `x_a x_c y_b y_e`, with a<c and b<e, are pairwise orthogonal in expectation. Each has mean square1. A violated rectangle contributes coefficient +4 or−4, and every other mixed coefficient is zero. Hence

    E[Delta_z(x,y)^2] = 16 · (# violated rectangles)
                     = 0 if s(z)=0, and384 otherwise.     (D1)

For independent uniform Rademacher vectors normalized by `1/sqrt(8)`, both input norms are exactly1 and

    E[(N(xy)-1)^2] = 384/8^4 = 3/32

for every invalid table. The normalization scales a biquadratic defect by1/64, hence its square by1/4096. A single dense test may hide cancellation; the nonzero mean-square gap is not a deterministic one-sample certificate.

## D3. Recovery and arithmetic meaning

The seven columns of the syndrome map are the seven nonzero vectors. Thus every invalid line table has exactly one line flip that kills its syndrome. This recovers the old Hamming one-error correction in the norm-identity language. It produces a valid table one line away, but does not establish which table was intended: two errors in an originally valid table can be miscorrected to another valid table.

Under §6's prime-selection dictionary, s(z) is the squareclass of `sigma(m(z))/30`, interpreted as parity subtraction even when that rational number is not integral. The same single line/prime toggle corrects both the norm-composition syndrome and the divisor-sum squareclass. Toggling uses inclusion/exclusion in the squarefree bank, not multiplying an already present prime again.

The quantitative conclusion strengthens the connection contract: the map preserves not only the valid set but also the exact number of failed norm rectangles and their mean-square energy, all determined by whether the shared three-bit syndrome is zero. It preserves neither individual signed coefficients nor the numeric value of sigma.
