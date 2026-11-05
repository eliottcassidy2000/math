# Independent review: matching bandwidth and short content certificates

**PASS**, no corrections requested, 2026-11-05.

Frozen source in `session197-odd-cycle-path-polynomials/`:

- CONSEQUENCES.md: 5616a1d747edad5828f48fb36efb615cdb01d60eaf770bef9e00029d81b74e1d
- verify_consequences.py: 7bff94a30fbd7c1af8fb7eb4ab59886f97800e57b385ff0aeebc9d811cc992e6
- Normal/optimized output: ae77cc9b13ec2e1d3c861cea3da341c0867a094ca75b001dad5405c0b5b9f386

Both copied-source replays reproduce this output. Dependencies are the separately reviewed universal-context and companion verifiers.

## Exact bandwidth, including cancellation

Every exterior segment between successive clone runs starts in O and ends in I. Choosing an O-to-I edge in each such segment gives a matching because all segments have disjoint exterior vertex sets, including across different crossing cycles. Hence total marked-run count is at most nu, and every shift exponent lies in [-nu,nu].

A maximum matching supplies nu disjoint crossing-triangle skeletons, contributing X^(N-nu)delta^nu. This is a valid skeleton packing even though the quotient notation uses the same abstract root: its marked gaps receive disjoint actual clone sets. For each skeleton term with total run count nu, the extreme forward coefficient is positive. No lower-run term reaches that shift. Thus there is no hidden cancellation in p_nu.

Distinct shifts are linearly independent over Q[x]: applying a proposed relation to x^m for m=0,...,2nu gives the Vandermonde matrix on x-nu,...,x+nu, whose determinant is a nonzero constant. Therefore the nonzero extreme coefficient really determines the exact bandwidth, not merely one choice of representation.

The positive delta/Delta normal form and degree bounds follow from delta^2=Delta^2+4Delta without increasing the total shift budget. At nu=1, crossing terms are only delta or Delta, so the remaining coefficient P is exactly the deleted-root profile. The positive-m clone ray and its odd offset are correct; m=0 must retain the zero-base term, as stated.

## Content and congruences

The two normalized basis identities are exact integer-valued polynomial identities. In the even case the sum of adjacent binomial coefficients equals d/(k+1) times binom(d+k,2k+1), exactly the needed r normalization. Evaluation at d=1,...,s+1 gives an integer unit lower-triangular matrix on the gamma-weighted coordinates, so the finite gcd and full content agree. The size s+1 is ceil(n/2) in both parities, including n=1,2.

For the modular statements, q_k at nonzero height contains the factor d(d-1)(d+1), divisible by6 at every integer d. The corresponding r_k contains d^2(d-1)(d+1), divisible by12. Integer Laurent-shift coefficients at x=1 preserve these divisibilities. Thus the scalar context congruences follow, within the intact-module scope only.

## Independent controls beyond the source census

The source's rooted templates through order four only reach nu<=1. The new audit program `path-polynomial-audit/consequences_independent.py` therefore additionally tests the carousel-five, carousel-seven and Paley-seven roots, with nu=2,3,3 and respectively1,1,2 maximum cross matchings.

A brute-force matching enumeration is independent of the source augmenting-path routine. The full skeleton operator has precisely the predicted support, and its extreme coefficient is respectively

    x(x+1), x(x+1)(x+2), 2x(x+1)(x+2).

Fifteen full clone profiles (m=0,...,4 for these contexts) are independently reconstructed from actual directed odd cycles and disjoint cycle collections, not from skeleton shifts or path-cover DP, and agree with the operator. These controls exercise higher rank, deletion, and possible shift cancellations.

Finally all1,099 labelled tournament profiles through order five are independently reconstructed from literal cycle packings. Their gamma-weighted gcd, first-ceil(n/2)-evaluation gcd, and longer evaluation gcd agree. The supplied normal/-O replay separately passes its285 rooted matching controls and1,099 content profiles.

The balanced extreme-coefficient values above are additional controls, not an unstated extension of the frozen theorem. The general exact bandwidth and content claims follow from the proof, not the census.
