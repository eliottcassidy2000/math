# Complete quotient-ring reduction and finite classification certificate

25 October 2026. Proof assembled; the independent odd-ring reviewer has reproduced all 203 candidate results and checked the reduction proofs. See `session55-audit-odd-rings` for the independent certificate. This extends THEOREM.md. The only unbounded external mathematical dependency is the proved prime-field classification in `session55-finite-field-geometry/THEOREM.md` (including its cited Weil bound and audited bounded census). No extension-field result is substituted for quotient-ring arithmetic.

## 1. Universal eventual prime-power constructions

Let L be the ordinary Lo Shu integer magic array

    8 1 6
    3 5 7
    4 9 2.

An affine transform 1+tL remains magic, with line sum 3+15t.

For an odd prime p, every unit congruent to 1 modulo p is a square modulo p^k for all k≥1. Proof: start with root 1 modulo p. If r²=x modulo p^j, replacing r by r+s p^j changes its square by 2rs p^j modulo p^(j+1). The coefficient 2r is invertible modulo p, so a unique s modulo p corrects the next digit. Induction proves the claim.

Therefore 1+pL is a magic array of squares modulo p^k. Its nine entries are distinct whenever p^(k−1)>8, since their nonzero integer differences are p·d for 1≤|d|≤8. Thus:

- p=3: every k≥3 is non-Parker;
- p=5,7: every k≥3 is non-Parker;
- p≥11: every k≥2 is non-Parker.

For 2-adic powers, every unit congruent to 1 modulo 8 is a square modulo 2^k. For a self-contained induction start with root 1 modulo 8. For k≥3, replacing an odd root r modulo 2^k by r+2^(k−1) changes its square by 2^k modulo 2^(k+1), toggling precisely the next bit. Consequently 1+8L is a square-residue magic array modulo 2^k and is distinct for k≥7.

The remaining negative prime-power cases are finite and certified by complete local enumeration: 2,4,8,16,32,64 and 3,9,5,25,7,49. Here is an additional analytic obstruction for 25 and 49, independently supplied by the reviewer. Normalize any nonzero square center modulo p to 1 by square scaling. For p=5, the only opposite square pair summing to 2 is (1,1), so every nonzero-center array is constant. With center zero, the allowed offsets are 0,±1; requiring u,v,u+v,u−v all in this set forces one of u,v to vanish. Thus a nonconstant array modulo 5 has three zero cells. For p=7, center zero forces the all-zero array since −1 is nonsquare. Center 1 has allowed offsets 0,±1, and again one of u,v must vanish; if the array is nonconstant it has three zero cells. In a square-residue lift modulo p² every cell that was zero modulo p must still be zero modulo p². Three such cells prohibit distinctness. A constant unit reduction permits at most p distinct values modulo p², also fewer than nine. These alternatives exhaust the reductions and prove 25 and49 Parker. Finally Z/9Z has only four square residues, so is immediately Parker. For the other odd primes, the external prime-field theorem gives exactly the Parker primes

    3,5,7,11,13,17,19,23,31,43,47,67.

This determines every prime-power quotient ring. The proof uses neither arbitrary lifting of modular arrays nor a general divisor monotonicity assertion.

## 2. Every odd modulus

Modulo 3 every magic square-residue array is constant. To see this, its center c is 0 or 1 and every opposite pair sums to 2c. Among square residues {0,1}, a pair summing to 0 is necessarily (0,0), and a pair summing to 2 is necessarily (1,1). Thus all nine entries are c. In the CRT theorem its collision mask is the full relation, which cannot change an intersection. For gcd(m,3)=1,

    Z/(3m)Z is Parker iff Z/mZ is Parker.

If an odd n contains two coprime factors each at least 5, the orthogonal square-AP construction in THEOREM.md makes n non-Parker. This covers at least two distinct primes ≥5, and also n with 3-adic exponent at least 2 and another prime factor. The only remaining forms are 3^k, p^k, and 3p^k for p≥5. Applying Section 1 gives exactly these 28 odd Parker moduli:

    3,5,7,9,11,13,15,17,19,21,23,25,31,33,
    39,43,47,49,51,57,67,69,75,93,129,141,147,201.

The largest is 201.

## 3. Finite reduction for arbitrary moduli

Write n=2^a m with m odd. If a≥7, combine the non-Parker array modulo 2^a with a constant-zero array modulo m by CRT; distinctness is already supplied by the dyadic component. Likewise if m is non-Parker, its zero collision mask supplies distinctness regardless of the dyadic component. Hence any Parker n must have

    0≤a≤6,   m∈{1}∪{the 28 odd Parker moduli above}.

There are exactly 203 candidates, all at most 64·201=12864. `classify_candidates.py` computes the exact local antichains, applies the proved CRT criterion, and saves every result in `classification_candidates.json`. It finds 137 Parker moduli including the trivial modulus 1, equivalently 136 with n≥2. The largest is 3216=16·3·67. Every excluded candidate has a reconstructible CRT witness and explicit square roots; every included candidate has the complete minimal collision-mask obstruction. Together with the preceding reductions, this is a complete classification of Z/nZ, not an extrapolation from a bounded search.

The candidate checker has not searched all full moduli through 12864. That is unnecessary: all large-modulus results are certified by complete prime-power local tables and the exact CRT intersection theorem. The separate full-modulus census independently checks the composition mechanism at every n≤160.

## 4. Literature and claim boundary

Cain, *Gaussian Integers, Rings, Finite Fields, and the Magic Square of Squares*, arXiv:1908.03236v2, Section 9, states the odd maximum 201 as Conjecture 9.3 and the global maximum 3216 as Conjecture 9.4: https://arxiv.org/html/1908.03236v2#S9 . The established classification above proves both statements, subject to the stated prime-field theorem dependency and its separate audit. The identified maxima themselves are inherited conjectures, not newly discovered numerical predictions.

A targeted public-web search on 25 October 2026 found Cain's paper and the OEIS Parker-ring entry, but did not locate a later proof. This is a limited literature search, not proof of priority. No assertion about novelty should be made without further literature review. The integer magic-square-of-squares existence problem remains untouched: modular witnesses and CRT square roots do not supply nine integer square entries with equal unreduced line sums.
