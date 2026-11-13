# The cited RH criterion and the precise Euler-product connection

**CITED:** Nicolas's equivalent formulation of RH.
**PROVED IN NOTE:** the finite residue-group comparison below.
No RH proof or numerical verification of its unbounded quantifier is claimed.

## Primary-source verification

Jean-Louis Nicolas, *Small values of the Euler function and the Riemann
hypothesis*, Acta Arithmetica155(3) (2012),311–321:
[primary HTML](https://arxiv.org/html/1202.0729v2),
[DOI](https://doi.org/10.4064/aa155-3-7).

Define

    c(n)=(n/phi(n)−e^gamma log log n) sqrt(log n),
    beta=2+gamma−log(4pi).

Theorem1.1(1.5) and Corollary1.1 state precisely

    RH iff, for EVERY n>=N_120569, c(n)<e^gamma(2+beta),
    N_120569=2·3·5···1591883.

Thus the user supplied the correct threshold, weak lower-bound inequality
on n, and strict upper bound on c(n). The coefficient is
e^gamma(4+gamma−log(4pi))=3.6444150964... . PDF text extraction can turn
the printed >= into >; the HTML preserves the TeX relation. The separate
limsup formulation is also equivalent, but is not a finite test replacing
the displayed universal condition.

The checker independently sieves through1591883, verifies that it is
prime, and counts exactly120569 primes. This confirms the finite threshold
identifier; it does not test the criterion for all integers above its
enormous primorial.

## Ambient units versus an affine power orbit

For any squarefree W, direct CRT counting gives

    |(Z/WZ)^*|/W=phi(W)/W=product_(p|W)(1−1/p).

The reciprocal Euler product W/phi(W) is the finite arithmetic quantity
entering the primorial extremal setting. It measures the density of
uniform residue classes that are units, without an exponent parameter.

For our A_k=99·2^k+23, let W be the product of the first25 primes.
Since k>=1, its2-coordinate is fixed and A_k is odd. At odd primes the
input x=2^k lies in the common cyclic orbit

    k ↦ (2^k mod p)_(p|W,p odd)

inside the product of multiplicative groups. A prime dividing99 or23
never rejects. For every other odd prime, the affine expression99x+23
excludes one nonzero input coordinate. If those group coordinates were
chosen independently and uniformly, the surviving proportion would be

    product_(p<100,p not in{2,3,11,23}) (p−2)/(p−1)
      =2851648031025/7696581394432.

The actual orbit has different marginal supports and shared exponent
phases, so its fraction is instead5392213/10201620. Both fractions are
exact, and their difference is not an error term in Nicolas's theorem.
The first counts an ambient product group; the second counts a particular
cyclic trajectory and affine forbidden coordinates.

The same loss appears already in the triples:30/phi(30)=15/4 and
42/phi(42)=7/2, while the235 and237 exponent-wheel fractions are3/4
and1. Restoring the multiplicative orders, discrete-log targets and the
shared exponent recovers the actual predicate.

This is the usable connection: prime-local factors supply finite groups;
the chosen orbit and its phase determine the sieve. Controlling the
all-integer extremal Euler-product error in the RH criterion does not,
by this comparison alone, control prime occurrence along99·2^k+23.
