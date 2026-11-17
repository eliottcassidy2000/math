# Separate the two-cover fringe, the Hamiltonian carry, and the odd phase

Write a one-apex tournament as a binary spine word b₀…b_(L−1): the spine is ordered forward, and the apex beats position i iff b_i=1. We build strong tournaments (first bit 1 and last bit 0) of **even** order (L odd). Write H for Hamiltonian paths and P for ordered two-path covers, with empty slots allowed. This result chooses actual tournaments with exact **prescribed residues**, not exact integer counts.

**Theorem.** At every power 2^a, every pair of residues with H odd and P≡0 mod gcd(8,2^a) is attained by such a word of length O(a²). Furthermore one can simultaneously prescribe H and P arbitrarily modulo any m coprime to 6, and prescribe H≡P modulo 3, with no change to the dyadic residues or even order. The guarded mixed construction has no useful absolute size bound of the kind that would prove spectrum coverage.

The dyadic restriction is necessary for even order ≥4: by the odd-pack formula each term in P has weight 2^(2c+s), where c is the number of odd cycles and s the number of remaining singleton vertices. Even order ≥4 gives divisibility by 8 termwise (if c=1, s≥1).

## 1. The endpoint fringe controls P₂

The one-apex odd-pack formulas are

    H = 1 + Σ (2^(d−1) + [d=1]),
    P = 2^(L+1) + Σ 2^(L−d)(3^(d−1)+1),

where the sums range over inversions b_i=1,b_j=0 with i<j and d=j−i. Indeed any odd directed cycle passes through the apex and chooses an even subset of the d−1 intervening spine vertices. No two such cycles are disjoint. The even-subset binomial sum gives the second read.

For large L, modulo 2^a the surviving inversions in P cross from one of the first a bits to one of the last a bits. Let c_j=1−b_(L−1−j) denote the reverse-complement final fringe. Then

    P ≡ Σ b_i c_j · 2^(i+j+1)(3^(L−2−i−j)+1)       (mod 2^a),

summing only fringe indices, with L long enough to keep these ends disjoint. An interior edit does not change this residue, and a length change by a multiple of the order of 3 mod 2^a preserves it if the fringes are kept.

Always take b₀=c₀=1. For L odd the term with i+j=2k or 2k+1 has exact 2-adic valuation 2k+3. This follows because 3^e+1 has valuation two for e odd and valuation one for e even. Program the radix-four digits of P/8, least significant first:

- At the first digit, (3^(L−2)+1)/4 is 1 if L≡3 mod4 and 3 if L≡1 mod4 (reduce 3^e mod16). The digit is that number plus b₁+c₁ mod4. With b₁+c₁ chosen to be at most one, all four digits are possible.
- For digit k≥1, the two bits b_(2k+1) and c_(2k+1) each increment it by 1 mod4; the corresponding exponent of 3 is even. In addition vary one of b_(2k),c_(2k) on the side for which the other side's index-one bit is zero. Its coefficient in this digit is 1 or 3 mod4. At least one such side exists by the first step. Interactions among these new high bits have index sum ≥2k+2 and hence only affect later digits. Coefficients from low opposite-side bits beyond index one likewise affect only later digits. The three available coefficients 1,1 and an odd number cover all four digit values by subset choice.

A truncated highest digit uses only its available low bits. This constructs every admissible P residue with fixed fringes and a chosen L mod4. When a≤3 there is no nontrivial P phase to select.

## 2. Interior carries adjust H without moving either fringe

Mod 2^a, the inversion formula for H only reads distances at most a. In a long all-zero guard, introducing a single isolated 1 contributes the sum over its next a zeros:

    Σ_(d=1)^a (2^(d−1)+[d=1]) = 2^a ≡0.

For two isolated-from-others 1s at internal distance j+1 (j≥1), superpose the singles but delete the first one's scheduled inversion to the second. The actual change is exactly −2^j modulo 2^a, provided the pair is surrounded by enough zeros. Existing 1s are too far away to detect these edits, and neither fringe changes. Make one guarded pair available for each j=1,…,a−1. Choosing the appropriate subset shifts H from its baseline to any prescribed odd residue. Separation at least a+1 and guard lengths O(a) per gadget suffice, so the dyadic word has length O(a²).

## 3. Insert the odd-phase word inside a protected zero interval

Use the finite word group from the [joint phase compiler](../chanpak_independent_primes_20261117/JOINT_PHASE.md). On a word ending in zero it has forward commands A=append 0 and B=append 10. It tracks normalized registers F,G,J,D,X,Y modulo m with H=F+X and P=D(Y+2G−2J). Vertical generators V₁,V₂ (commutators) and V₃ (the four-conjugate word) fix F,G,J,D and have translations in X,Y spanning the full plane mod m. In particular the determinant away from five is −5/72. At powers of five, V₃ changes Y+2X by 2, while V₁ fixes that sum and changes X by a unit. This is why simply omitting the five-prime part is unnecessary.

First construct the dyadic word with a long zero interval between its gadgets and its closing fringe. Split it in this interval. The prefix ends in at least a zeros and the fixed suffix begins in at least a zeros. We will insert literal positive commands at that point. In addition to the m-registers, track length modulo λ=lcm(4,ord_(2^a)(3)) and (F,X) mod3. The A and B actions remain bijective on this augmented finite state. Replace each forward B by append 10 followed by N zeros, where N≥a+1 is a multiple of the order of A on the augmented state. This has exactly the original modular action of B but places a long zero guard after every newly appended 1. Negative formal powers of either command are positive repetitions on its finite orbit. Consequently even inverse operations insert only isolated interior 1s, each which contributes zero to H mod2^a. No edit changes the two fringes. The vertical words have formal net length zero and therefore the implemented length change is a multiple of λ, preserving P mod2^a.

For mod3, A preserves X+F², and B sends F→F−1, X→X−F, increasing X+F² by one. Thus the commutators V₁,V₂ preserve H mod3. The V₃ word is the product of powers of conjugates W_j with exponents −6,19,−16,4. Their sum is 1; each increases X+F² by one and their total translation in F is zero (the translation polynomial vanishes at 2). So V₃ increments X and H by one mod3. Its exponent can be chosen mod3 independently of the choice mod m by CRT. At even order the odd-pack expansion gives P≡H (mod3), as required.

Finally append the original suffix. The fixed suffix acts identically on F,G,J,D whether or not we have inserted a vertical word; its changes to X and Y are additive functions of those registers and the literal suffix bits. Therefore initial vertical shifts in X,Y persist through the suffix. To meet the target P mod m, choose its Y shift to be the desired absolute P residual divided by the fixed final D modulo m. This programs both requested odd coordinates without breaking either dyadic invariant.

## Consequence for independent prime targets and limit

Feed B into the [unimodular induced pairing](../chanpak_independent_primes_20261117/UNIMODULAR.md) p=44P−13H, q=71P−21H. At every modulus of form 2^a·3·m, gcd(m,6)=1, this realizes arbitrary independent (p,q) phases subject to the visible even-source constraints: both odd, p≡q mod 2^min(a,3), and p≡−q mod3. It does not assert the two numbers are prime. Higher powers of 3 are not covered here.

This is still a local construction. In fact the same positive-P₂ reader on the one-apex family remains sparse in absolute counts by the [all-word span bound](../chanpak_eight_split_20261117/THIN_READERS.md). Thus highly flexible phase invariants can coexist with a density-one set of missed prime values. The next absolute repair must control size with the phase, or change the source family.

Replays: `python replay_dyadic.py` checks every joint dyadic phase through exponent seven and deeper samples. `python replay_mixed_phase.py` compiles literal guarded words and re-reads all three modulus factors using the independent unnormalized inversion recurrence. Both scripts also run in optimized Python.
