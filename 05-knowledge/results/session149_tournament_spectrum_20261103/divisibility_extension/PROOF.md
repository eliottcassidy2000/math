# A second bit of the even-order tournament count

For every even integer n>=8, the number A(n) of tournament isomorphism classes satisfies

A(n)/2^(n/2) = epsilon(n) phi(n)/2 (mod 4),

where epsilon(n)=+1 when n=2 mod 4, and epsilon(n)=-1 when n=0 mod 4. The quotient is an integer, by the earlier bound in ../divisibility/PROOF.md. In particular, v2(phi(n))=2 implies v2(A(n))=n/2+1, and v2(phi(n))>=3 implies v2(A(n))>=n/2+2. Together with the previously established first bit, if v2(phi(n))<=2, then v2(A(n))=n/2+v2(phi(n))-1.

## Proof

Use the notation E, z, d=v2(z), and c from the preceding proof. Every term in the Davis partition sum has valuation E-d; every odd integer is invertible modulo 4. We identify terms of valuation at most n/2+1.

If c>=6, the general inequality E-d >= (n-c)/2+binom(c,2)-v2(c!) and the bound v2(c!)<=c-1 give E-d-n/2 >= binom(c,2)-(c-1)-c/2 >=2. (Only even c arise.)

If c=4 and the four cycles are not all the same length, then d<=2: the largest 2-adic factorial valuation among nonconstant multiplicities of total four is attained by two pairs, and equals 2. Therefore E-d-n/2 >= 6-2-2=2. If all four cycles have common odd length i, then n=4i; n>=8 forces i>=3, and E=2(i-1)+6i=8i-2, d=v2(4!)=3, so E-d-n/2=6i-5>=2.

If c=2 with equal odd lengths i, then n=2i and E-d-n/2=i-2. For n>=8, i>=5, so this is >=3. If c=2 with distinct odd lengths i<j, then d=0 and E-n/2=gcd(i,j)-1. Since gcd(i,j) is odd, it equals zero for coprime pairs and at least two otherwise.

Consequently only the coprime pairs i<j with i+j=n survive modulo 4, with normalized contributions 1/(ij). Modulo 4 the inverse of an odd number is itself. If n=2 mod 4, i and j are congruent mod 4 so ij=1 mod4. If n=0 mod4, they have opposite residues 1 and 3, so ij=3 mod4. There are phi(n)/2 such unordered pairs: each i coprime to even n is odd and pairs with n-i, and n/2 is not coprime to n for n>=4. This proves the formula.

Source: Davis's formula as recorded at https://oeis.org/A000568, also deriving from Royle et al., Theorem 3.3, https://link.springer.com/article/10.1007/s10801-022-01197-0. This result is a deduction; no claim of priority is made. Note the lower bound is not an exact general formula when v2(phi(n))>=3: cancellation at still higher powers is not determined by this congruence.

## Reproduction

Run `python mod4_check.py` and `python -O mod4_check.py`. The script imports ../divisibility/partition_count.py and uses its exact rational Burnside sum for even n=8,10,...,80; errors are explicit exceptions. It verifies divisibility and the mod 4 formula. The common Burnside foundation means this is not independent evidence for the proof, but it exercises the finite partition formula beyond the source table's n=76.
