# Exact prefix formula and all-word dyadic gaps

Status: PROVED conditional only on the separately proved four-register physical construction. This is a family-specific obstruction, not a threshold counterexample or a universal finite-control obstruction.

Let w=b1...bn be any binary word and V(w) its ordinary binary value. For the physical family in the companion packet, put X=D+u and Y=D+v. Direct substitution in its exact register equations gives

0: (X,Y) -> (2X,13X/6-36),
1: (X,Y) -> (2X-36,13X/6-78).

Initially X=540 and Y=546. Therefore, with K0=15 and Kj=2K(j-1)-bj,

Xj=36Kj, Kj=15*2^j-V(b1...bj).

For n>=1 define

T(w)=7b1 + sum over j=2,...,n of bj*(K(j-2)-b(j-1)).

Then the remaining register and physical output have the exact formulas

u(w)=-23-6T(w),
D(w)=36Kn+23+6T(w),
s(w)=1,
B(w)=12Kn+7+2T(w).

Proof: digit 0 preserves u. Digit 1 sends u to u+X-Y-36. At the first digit this increment is -42. At any later digit j, the preceding output formulas for X,Y give X(j-1)-Y(j-1)-36=-6K(j-2)+6b(j-1). Summing these increments proves the displayed formula. The B formula follows from D=3B+2, since s=1.

Every summand in T is nonnegative: Kj>=14*2^j+1, and bj is zero or one. Also

0 <= T(w) <= 7+15*sum(j=2,...,n) 2^(j-2)
             =15*2^(n-1)-8.

This upper bound is valid for n=1 as well. Since

14*2^n+1 <= Kn <=15*2^n,

we obtain, uniformly over ALL words of length n>=1,

504*2^n+59 <= D(w) <=585*2^n-25,

168*2^n+19 <= B(w) <=195*2^n-9.

These simple bounds are not claimed sharp. They imply that the open interval

(585*2^n-25, 1008*2^n+59)

contains no D output from any word of any length, for every n>=1. The empty word gives D=563, below these intervals. Earlier lengths are bounded above by the left endpoint, and later lengths are bounded below by the right endpoint. Equivalently the B intervals

(195*2^n-9, 336*2^n+19)

are missing from the entire family. They contain arbitrarily long runs of odd integers.

Consequently no choice of a word for each target makes the physical B outputs of this UNCHANGED seed, these two operations, and this fixed closing context cover a tail of the odd integers. This is stronger than failure of a bounded additive correction. Changing seed, operations, available closing contexts, or combining other families is outside the claim. The broader eventual odd Hamiltonian-path spectrum question and F(8,137) remain open.

The prefix formula is also a positive exact synthesis rule: any chosen binary word gives a literal tournament with the displayed B and s=1. It makes the mixed-digit defect a weighted sum of earlier prefixes rather than an unexplained failure of commutation.

## Reproduce

Run `python verify.py`. This exact symbolic checker verifies the two quotient identities and the summand update algebra, and tests the closed prefix formula against the independently certified register dynamics for all 8191 words of lengths 0 through 12. The all-word proof is the algebra and induction above, not the finite tests. The physical realization depends on the companion multiport-control certificate.

## Reviewed physical dependency

The all-word physical realization and register equations are independently audited and published at [session204/multiport-control](https://github.com/eliottcassidy2000/math/tree/9e2699ee112e88966d540bdece7f68e2f0265f50/session204/multiport-control). The present independent audit is REVIEW.md.
