---
id: THM-4497
title: "Coprime-atom weighted graphs preserve multiplicative relation kernels"
status: "PROVED + INDEPENDENTLY AUDITED; FINITE-EXACT Collatz controls"
source: "crossroads10-20260926"
depends_on: []
proofs:
  - 05-knowledge/results/crossroads10_20260926_arithmetic.md
scripts:
  - 04-computation/experiments/crossroads10_20260926_arithmetic.py
  - 04-computation/experiments/crossroads10_20260926_audit.py
---

# THM-4497: an exact arithmetic graph with a residual matrix

**PROVED + independently audited.** Reserved as an empty stub in pushed
checkpoint3130458ce. Root and geometry independently audited the kernel,
weighted-cycle, common-prime and sign-torsion arguments. The algebra is
self-contained and uses standard gcd-free bases; no algorithmic priority
or factoring-complexity claim is made. Collatz and G2 remain OPEN.

## All-family statement and proof

Let a_0,...,a_(N-1)>1 be integers. They admit a gcd-computable pairwise
coprime basis b_1,...,b_s>1 with a_i=prod_h b_h^E_(h,i), E_(h,i)>=0.
The atom matrix E has exactly the same rational and integral relation
kernel as the prime valuation matrix.

Indeed repeatedly replace overlapping distinct factors b,c by
gcd(b,c), b/gcd(b,c), c/gcd(b,c), discarding1 and duplicates. The product
of current distinct factors strictly decreases and all original values
remain products of current factors. Termination gives the claimed basis.
For p|b_h, its prime valuation row is v_p(b_h) times atom row h. Every
atom has at least one such prime and distinct atoms share none.

Turn one-support rows into grounded vertices and two-support rows
u c_i+v c_j=0 into edges transporting c_j=-(u/v)c_i. Retain parallel
edges and endpoint weights. A connected component has zero solution
if grounded or if a cycle has gain different from1. Otherwise choose
nonzero rational weights w_i so its solutions are c_i=t_C w_i.
Let b be the number of surviving components, including isolated vertices.
For each higher-support row h define

    R_(h,C)=sum_(i in C) E_(h,i) w_i.

Then the full relation kernel is exactly the image of ker R under this
component substitution, and

    rank_Q(E)=N-b+rank_Q(R).

The proof is ordinary elimination with an injective substitution from
the b free component parameters. Clearing denominators recovers integral
relations. An odd cycle always forces zero since its gain is negative;
an even cycle forces zero precisely when its exponent ratios are
unbalanced. Omitting R or the edge weights need not preserve rank.

## Factor-free sufficient subcertificate

Gcd stripping detects private primes and pair-exclusive primes without
factoring. In the resulting unweighted graph, each private vertex and
each component containing an odd cycle forces all its coefficients zero.
Remove forced coordinates and repeat: higher-support primes may become
private or pair-exclusive. This iterative test is sound but not complete;
the weighted graph with R is the exact completion.

## Paired values and Collatz scope

For pairs (-3m_i,3m_i+1), take one coprime basis of all2N positive
magnitudes and stack the atom rows of the two coordinates separately.
Full column rank is equivalent to multiplicative independence of the
pairs. A rational kernel vector clears to an integer vector and may be
doubled to eliminate the first coordinate's minus sign. Do not add the
two coordinate matrices or discard the common prime3.

The actual odd orbit of27 contains91=7*13,175=5^2*7,325=5^2*13.
Their exclusive-prime triangle proves independence although none has a
private prime and all three strict gcd-dominance tests fail at equality.
Seven other orbit nodes extend this to a ten-node example. An actual
consecutive ten-node prefix from199 contains the triangle299,253,143.

**FINITE-EXACT:** among9446 injective ten-odd-node prefixes from odd
sources3..20000, private primes certify6487, one graph pass8481,
iterative removal9061, exact first-slot rank9350, and exact paired rank9446.
The independent control universe of8835 small multisets agrees among
direct prime elimination, gcd-free atoms and graph compression.

The same unweighted C10 can have rank9 or10 after changing one endpoint
exponent. On one actual orbit, the odd values507 and13 are independent
but their first slots1521 and39 are dependent, since1521=39^2. These are
explicit failure boundaries for forgetting weights or the distinguished
factor3. The separate paired coordinates restore rank2 in that example.

No assertion covers every chronological prefix. Appending a node can
turn a pair-exclusive prime into a higher-support prime; certificates
cannot simply be repeated along an infinite orbit. Full rank alone also
does not imply descent. See the proof note for complete controls and
the carry-divisibility restriction on actual Collatz graph realizations.
