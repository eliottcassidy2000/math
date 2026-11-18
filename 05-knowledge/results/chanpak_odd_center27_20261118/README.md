# Odd center 27: the next physical lift, two modular geometries, and a π threshold

**Main result.** Both counts 27·2^m−1 and 27·2^m+1 occur for every
integer m≥1. For m≥3 they have strong witnesses of the same order m+6,
with the same transitive flank, and two different marked cores. This is
the next centered step of the [slope-9 signed Collatz lift](../chanpak_machin_collatz_lift_20261118/README.md),
not arbitrary odd-center closure or a proof of Chan–Pak.

## Physical proof

For t≥1, put a transitive chain L of order t+2 before a five-vertex
tournament D: every L vertex beats every vertex of D. Add one apex v.
It beats exactly the first two vertices of L and the specified outmask
in D. All other L and D vertices beat v. Let h=H(D), κ=H(D+v) with the
indicated contacts, and

    R = sum over nonempty S⊆D of H(D−S) times
        (number of Hamilton paths of D[S] ending at a vertex beating v).

Use H(empty)=1. The following literal rows index D from 0 and give the
specified data:

| sign | outgoing rows of D | apex outmask on D | h | R | κ |
|---|---|---:|---:|---:|---:|
| − | (2,4,1,7,15) | 4 | 3 | 33 | 5 |
| + | (2,20,9,3,13) | 13 | 11 | 25 | 23 |

Set u=3·2^t, the binary tag `11` followed by t zeroes. By the
[odd-cycle packing formula](https://arxiv.org/abs/2307.05569), distinguish
packings with no cycle through v, with a v-cycle only on L, only on D,
or crossing the seam. Cycles away from v lie entirely in D. The
no-v and D-only terms sum to κ. The transitive L with apex has
H=u−1: for each of the two possible starts, ending the cycle at a
beaten-back vertex in L gives the geometric sums, with one adjacency
correction; equivalently F(3)=1 followed by t zero recurrences
F(2x)=F(x)+x+[x odd]. The L-only term is therefore h(u−2).
A crossing v-cycle begins at one of the first two L vertices and
then uses some nonempty S of D. At the two starts there are respectively
t+1,t spare L vertices; exactly half the subsets have the needed parity.
The weight 2 of an odd cycle cancels this half, giving
(2^(t+1)+2^t)R=uR. This counts residual packings on D−S through H(D−S).
Hence the whole host has

    H = (h+R)u + κ−2h = 36u ±1 = 27·2^(t+2) ±1.

Both hosts are strong: v reaches L and D via the first L vertex, and
each vertex of D reaches a marked in-contact. The tables are replayed
directly. The templates are changed, not a single-arc or native
metagraph edit. In each sign family, adding one L vertex immediately
before D is an induced extension, producing H⁻→2H⁻+1 and
H⁺→2H⁺−1. The small boundary pairs (53,55),(107,109) have
separate finite strong witnesses in the replay. In particular the
mixed arithmetic step from 71/73 lands on physical 107/109. By
contrast the ordinary positive sign takes 69→13, 71→107 and 73→55
at their respective odd shortcuts. All three are exact, but they
do not carry both ends of the same unit bracket to the same center.

There are two additional numerical landings with a structural cause.
Because 27=−1 mod 7 and 2 has order three, 7 divides the plus leg
exactly when m=0 mod 3, and never divides the minus leg. Thus the
plus ray includes 217=7·31 and the taxicab 1729=7·13·19, with its
same-order repair mate 1727. This is no theorem about infinitely
many prime quotients or about the next odd center 81.

The entire 7-adic filtration on the plus branch is explicit. Put
m=3j, so the number is 27·8^j+1. The condition v₇≥2 is exactly
j=4 mod 7, because 8^j=1+7j mod 49. More generally each condition
v₇≥d gives one residue j mod 7^(d−1). To lift it one digit use
8^(7^(d−1))=1+7^d c with c a unit mod 7 (by binomial induction).
Adding one step of 7^(d−1) to the exponent moves the
normalized error by a nonzero mod-7 increment; exactly one of the
seven digits cancels it. In particular every exact positive
7-adic valuation appears infinitely often on this physical ray.

## The two geometries at order 27

The prior [odd multiplication-table and stitched-sheet result](../collatz_spectrum_equality_and_stitched_sheets_20261005.md)
is inherited: the central entries are ±4^(-1) in the checkerboard
pattern, since their row and column labels are ±2^(-1). At modulus
27 these entries are 7 and 20. The unhalved affine sheets
3x−1 and 3x+1 take 7→20→7 modulo 27, reducing the rational
two-cycle {1/4,−1/4}. For a *specified exact* odd Collatz
valuation e, the modular target is 2^(-e)(3x+1) modulo 27.
The affine residue map on Z/27 is exactly three-to-one onto
the residue class 2^(-e) mod 3 and reads only x mod 9. The
unit 2 has order 18 mod 27. CRT shows that any such compatible
mod-27 residue can be combined with the separate single source
class imposing v₂(3x+1)=e modulo 2^(e+1). This does not read the
next valuation or imply integer convergence. For example, the
classical integer orbit from 27 has 111 steps to 1 and reaches 9232;
this finite trajectory is a different claim from the modular cycle.

In the tournament geometry write A for the skew ±1 off-diagonal
dominance matrix (zero diagonal), and J for the all-ones matrix.
The border

    [  1    1^T  ]
    [ −1    I+A  ]

is skew Hadamard of order 28 exactly when A1=0 and AA^T=27I−J,
by block multiplication. It is supplied by the
[Paley square-difference construction](https://doi.org/10.1016/j.ejc.2012.10.006)
on the additive group of F_27. The replay constructs the field by
x^3−x−1 and checks every integer dot product. The additive group
is (Z/3)^3, not Z/27. No ordinary circulant Hadamard matrix is
asserted by this bordering.

In fact a Z/27 circulant skew tournament cannot supply this same
border. Write A's difference sign f(s), f(−s)=−f(s), f(0)=0, and
let ζ be a primitive 27th root. Were the border Hadamard, all
nontrivial eigenvalues Σ f(s)ζ^(js) would be ±i√27, which lie
in Q(ζ^9) and are fixed by every automorphism ζ↦ζ^u with
u=1 mod 3. Thus Fourier inversion forces f(us)=f(s) for all u=1
mod 3. The signs on the nine units s=1 mod 3 must all be equal,
and those on s=2 mod 3 are opposite. For j=9 the nonunits cancel
by skewness and the resulting eigenvalue is ±9i√3, not ±3i√3.
Same cardinality, different Fourier constraints. This does not
address general circulant Hadamard matrices.

## Cubes, reciprocity, and exact approximation

The inherited [divisor balance theorem](../collatz_mod6_20260917_divisor_balance_family.md)
distinguishes p, p^3 and p^2qr by their exponent boxes; 27 is the
smallest nontrivial odd cube and has the chain 1,3,9,27. At the
cyclotomic layer relevant to Nagell–Ljunggren, if r≠3 is a prime
dividing x²+x+1 then x has order 3 mod r, hence r=1 mod 3. For
odd r≠3 this is exactly the splitting condition (−3/r)=1 by
quadratic reciprocity. The prime 3 occurs exactly to valuation
one when x=1 mod 3, since substituting x=1+3z gives
x²+x+1=3(1+3z+3z²). Thus a perfect cube of this form cannot
be 27. Seven, however, splits, and 18²+18+1=7³=343 is an exact
[Nagell–Ljunggren solution](https://tidsskrift.dk/math/article/view/15038).
This is a reciprocity *filter* plus an exact valuation, not a
generalization or proof of reciprocity; note also that the
repunit equation is not the quartic Ljunggren equation used to
explain Machin's 239.

For a precise example of an irrational threshold supplied by
tournaments, 355 and 113 are
[verified Hamilton counts](../chanpak_machin_blue_farey_20261118/README.md).
The ordinal products give the ratio of actual counts

    r_m = 355(27·2^m−1) / [113(27·2^m+1)].

It increases toward 355/113. Exact alternating rational bounds
in Machin's identity give r_19 < π < r_20. In other words the
actual dyadic count family finds a discrete crossing index, not just
two arbitrary close real numbers. This is not an optimal new
Farey bound nor evidence of exact absolute coverage for all integers.

Run `python3 verify.py` for the literal core data, small boundary
witnesses, direct counts, Hadamard dot products and π threshold.
