# Seven is the smallest prime for the three-family residue pattern

2026-10-19. **PROVED and FINITE-EXACT.** Optional addendum; the earlier frozen packet is unchanged.

For a prime p, seek three residues a,b,c such that 2a,2b,2c are nonsquares modulo p, while a+b,a+c,b+c are squares modulo p. Zero counts as a square. This is exactly the purely residue-based pattern that rules out all within-family edges and permits all three cross-family edge types in the simultaneous two-subdivision construction. The three residues must be distinct: equality of two would make their cross sum simultaneously square and nonsquare.

The smallest prime supporting this pattern is seven, and its unique unordered triple is {3,5,6}.

- Modulo two every residue is square, so no family is admissible.
- Modulo three the only residue a for which 2a is nonsquare is a=1, so three distinct families are impossible.
- Modulo five the admissible residues are exactly {1,4}, again too few.
- Modulo seven the admissible residues are exactly {3,5,6}. Their pair sums are 1,2,4 modulo seven, all squares. Since there are only three admissible residues, the unordered triple is unique.

The exact checker enumerates all residues and triples for p=2,3,5,7. This minimality concerns this particular three-family prime-field suppression pattern. It does not show that every subdivision construction needs modulus seven, that seven is minimal among arbitrary composite-modulus devices, or that fewer residue families cannot suffice. In fact the independent one-subdivision construction uses square scaling and a power of two instead.
