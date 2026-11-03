# Independent growth and sparse-prime audit

Status: PASS conditional on the inherited mixed-prefix tournament realization; it concerns only that named family.

The seed (Y,X)=(25,34) fails the cone 2X≥3Y: 68<75. Its two one-letter successors are q(25,34)=(411,688) and b(25,34)=(193,302), both in the cone. For arbitrary positive Y,X, the cone slack after q is 3Y+2X and after b is Y, so it is preserved. Inside the cone Y′_q=11Y+4X≥17Y and Y′_b=5Y+2X≥8Y. The first q meets 411≥17·193/8; the first b meets equality 193=8·193/8; the empty seed 25≥193/8. Hence all words satisfy Y≥(193/8)17^r8^s. The even-parity companions are 3Y+20 and 3Y+22, each at least (579/8)17^r8^s.

For any σ>δ, where 17^(−δ)+8^(−δ)=1, weighted counting over all words (even and odd parity, an upper bound) gives 2(8X/579)^σ/(1−17^(−σ)−8^(−σ)). An independent high-precision root is δ=0.28452464288487406704344… . Choosing σ<1 gives density zero among odd positive integers. Possible duplicate outputs only reduce the count.

For any fixed modulus M=7·13^b·19^c·31^d and reduced class j, the number of p≤x with 7p in the family is O(x^σ), while PNT for arithmetic progressions gives x/(φ(M)log x) primes in the class. Thus a proportion tending to one have 7p outside this family. This cannot exclude those 7p from the full tournament spectrum. The conclusion is asymptotic for each fixed M, not uniform as M grows with x; membership of the quotient in a residue class does not establish its primality.
