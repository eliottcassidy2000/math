# Addendum: exclusion of every distinct arithmetic-progression triple in the Lagrange A-part

This is a separate strengthening of the family-specific analysis in `proof.md`; it does not replace or modify the modulo-47 argument there. The source remains the Jean Lagrange family reproduced and attributed in Ajai Choudhry, arXiv:2508.07806v1, equation (1.1), https://arxiv.org/html/2508.07806v1#S1.E1 . The result here is independently derived, not attributed to the paper.

## Theorem

For every rational parameter pair (u,v), no three distinct values in the first quadruple A form an arithmetic progression. More precisely, for a nonzero pair (u,v), an equation

2 A_m = A_i + A_j

with three different indices holds if and only if the index triple is {0,1,3} and 2u+3v=0. In that event A_0=A_1=A_3. Here and throughout this addendum indices are zero-based, matching the machine certificate; source indices are one higher.

## Exhaustive factor proof

There are four unordered triples of indices and three choices of middle in each, hence exactly twelve possible AP equations. The appendix below factors all twelve integer polynomial balances 2A_m−A_i−A_j. Nine have no linear factor; the other three have the factor 2u+3v, precisely for triple {0,1,3}. Each balance is a nonzero integer constant times its listed factors. Every remaining factor is homogeneous of degree three or four.

For each residual factor F, the appendix provides a prime p and its complete evaluation vector

(F(0,1), F(1,1), …, F(p−1,1); F(1,0)) modulo p.

Every entry is nonzero. This is an exhaustive certificate that F has no projective zero over F_p: every projective point either has v≠0 and scales to [t:1], or is the point [1:0]. A homogeneous integer form with no such finite-field projective zero cannot have a nonzero rational zero. Indeed, clear denominators of a proposed rational pair, then divide the resulting integer pair by its gcd. Homogeneity preserves vanishing, and the primitive pair cannot have both coordinates divisible by p. Its reduction would therefore be a forbidden projective zero. This argument handles all rational denominators, including those divisible by p; no bounded search or integrality assumption on the original pair is involved.

It follows that only the explicitly retained linear factor can vanish rationally. Every rational point of 2u+3v=0 has (u,v)=(-3t,2t) with rational t. Direct substitution gives

A(-3t,2t)=(900,900,12996,900)t⁴,

where the notation on the left denotes the quadruple evaluated at that parameter pair. Thus its surviving index triple {0,1,3} consists of three equal values. At t=0 all entries are zero. Neither case is an AP of three distinct values, completing the proof.

## Consequence and exact boundary of the claim

This excludes an AP-K3,3 certificate obtained by choosing three original rows and three original columns of this square-sum K4,4 family, whenever that certificate requires the labels in each chosen part to be a three-term AP of distinct values. Whichever part is selected first, the selected triple from the original A part would violate the theorem. Exchanging the parts therefore does not evade it.

Arbitrary permutations are covered by exhaustion of the twelve equations. Translation of all labels in a part preserves each AP balance because its coefficients sum to zero. Multiplication by a nonzero rational square preserves both AP balance and distinctness. Consequently the obstruction persists under translation gauge and nonzero square homothety. Zero scaling collapses the configuration and is not a nondegenerate exception.

This is stronger information about the family's original labels than the whole-part parallelogram guard: it obstructs all distinct AP triples in A, not merely its four-label balance. It does NOT exclude every selection of nine square values from the ambient configuration, every possible regrouping into new additive coordinates, or every re-embedding into a larger ambient square-sum graph. Nor does it make a universal claim about other published families or MSOS existence.

## Twelve exact certificates


### Case 1: triple (0, 1, 2), middle 0

Exact factorization: 18*(160*u**4 - 160*u**3*v - 4492*u**2*v**2 - 10700*u*v**3 - 7335*v**4).

Residual F=160*u**4 - 160*u**3*v - 4492*u**2*v**2 - 10700*u*v**3 - 7335*v**4. Degree 4, prime p=13.

Affine residues, in order t=0,…,p−1: 10, 2, 12, 4, 12, 10, 3, 1, 6, 12, 5, 2, 12.

Infinity residue: 4.

### Case 2: triple (0, 1, 2), middle 1

Exact factorization: -9*(640*u**4 + 5408*u**3*v + 18068*u**2*v**2 + 27508*u*v**3 + 15705*v**4).

Residual F=640*u**4 + 5408*u**3*v + 18068*u**2*v**2 + 27508*u*v**3 + 15705*v**4. Degree 4, prime p=13.

Affine residues, in order t=0,…,p−1: 1, 2, 2, 5, 9, 6, 8, 8, 6, 9, 5, 2, 2.

Infinity residue: 3.

### Case 3: triple (0, 1, 2), middle 2

Exact factorization: 9*(320*u**4 + 5728*u**3*v + 27052*u**2*v**2 + 48908*u*v**3 + 30375*v**4).

Residual F=320*u**4 + 5728*u**3*v + 27052*u**2*v**2 + 48908*u*v**3 + 30375*v**4. Degree 4, prime p=13.

Affine residues, in order t=0,…,p−1: 7, 11, 4, 10, 11, 12, 2, 6, 7, 11, 8, 11, 4.

Infinity residue: 8.

### Case 4: triple (0, 1, 3), middle 0

Exact factorization: 2*(2*u + 3*v)*(1440*u**3 + 5856*u**2*v + 6314*u*v**2 + 1035*v**3).

Residual F=1440*u**3 + 5856*u**2*v + 6314*u*v**2 + 1035*v**3. Degree 3, prime p=7.

Affine residues, in order t=0,…,p−1: 6, 1, 6, 2, 5, 3, 5.

Infinity residue: 5.

### Case 5: triple (0, 1, 3), middle 1

Exact factorization: -(2*u + 3*v)*(1440*u**3 + 4704*u**2*v + 3626*u*v**2 + 1035*v**3).

Residual F=1440*u**3 + 4704*u**2*v + 3626*u*v**2 + 1035*v**3. Degree 3, prime p=7.

Affine residues, in order t=0,…,p−1: 6, 4, 4, 1, 4, 1, 1.

Infinity residue: 5.

### Case 6: triple (0, 1, 3), middle 3

Exact factorization: -(2*u + 3*v)*(1440*u**3 + 7008*u**2*v + 9002*u*v**2 + 1035*v**3).

Residual F=1440*u**3 + 7008*u**2*v + 9002*u*v**2 + 1035*v**3. Degree 3, prime p=7.

Affine residues, in order t=0,…,p−1: 6, 5, 1, 3, 6, 5, 2.

Infinity residue: 5.

### Case 7: triple (0, 2, 3), middle 0

Exact factorization: 2*(1440*u**4 - 672*u**3*v - 37484*u**2*v**2 - 93612*u*v**3 - 66015*v**4).

Residual F=1440*u**4 - 672*u**3*v - 37484*u**2*v**2 - 93612*u*v**3 - 66015*v**4. Degree 4, prime p=11.

Affine residues, in order t=0,…,p−1: 7, 7, 6, 6, 7, 7, 2, 8, 6, 8, 2.

Infinity residue: 10.

### Case 8: triple (0, 2, 3), middle 2

Exact factorization: 2880*u**4 + 53088*u**3*v + 249356*u**2*v**2 + 445548*u*v**3 + 273375*v**4.

Residual F=2880*u**4 + 53088*u**3*v + 249356*u**2*v**2 + 445548*u*v**3 + 273375*v**4. Degree 4, prime p=11.

Affine residues, in order t=0,…,p−1: 3, 4, 5, 1, 5, 4, 3, 3, 1, 1, 3.

Infinity residue: 9.

### Case 9: triple (0, 2, 3), middle 3

Exact factorization: -5760*u**4 - 51744*u**3*v - 174388*u**2*v**2 - 258324*u*v**3 - 141345*v**4.

Residual F=5760*u**4 + 51744*u**3*v + 174388*u**2*v**2 + 258324*u*v**3 + 141345*v**4. Degree 4, prime p=11.

Affine residues, in order t=0,…,p−1: 6, 7, 6, 2, 8, 7, 7, 8, 2, 6, 7.

Infinity residue: 7.

### Case 10: triple (1, 2, 3), middle 1

Exact factorization: -64*(45*u**4 + 498*u**3*v + 2023*u**2*v**2 + 3498*u*v**3 + 2160*v**4).

Residual F=45*u**4 + 498*u**3*v + 2023*u**2*v**2 + 3498*u*v**3 + 2160*v**4. Degree 4, prime p=13.

Affine residues, in order t=0,…,p−1: 2, 8, 8, 8, 2, 11, 5, 7, 2, 2, 7, 5, 11.

Infinity residue: 6.

### Case 11: triple (1, 2, 3), middle 2

Exact factorization: 128*(45*u**4 + 534*u**3*v + 2161*u**2*v**2 + 3624*u*v**3 + 2160*v**4).

Residual F=45*u**4 + 534*u**3*v + 2161*u**2*v**2 + 3624*u*v**3 + 2160*v**4. Degree 4, prime p=17.

Affine residues, in order t=0,…,p−1: 1, 7, 9, 3, 11, 13, 15, 15, 3, 12, 16, 15, 1, 9, 15, 4, 4.

Infinity residue: 11.

### Case 12: triple (1, 2, 3), middle 3

Exact factorization: -64*(45*u**4 + 570*u**3*v + 2299*u**2*v**2 + 3750*u*v**3 + 2160*v**4).

Residual F=45*u**4 + 570*u**3*v + 2299*u**2*v**2 + 3750*u*v**3 + 2160*v**4. Degree 4, prime p=13.

Affine residues, in order t=0,…,p−1: 2, 10, 8, 5, 11, 11, 4, 3, 9, 11, 12, 3, 2.

Infinity residue: 6.

## Reproduction

Run `python verify_ap.py` and `python -O verify_ap.py`. All checks raise explicit exceptions on failure and are independent of Python assertions. `ap-certificate.json` contains every exact factorization, all finite-field values, and the exceptional-line substitution. The original `verify.py`, `certificate.json`, and `proof.md` remain unchanged.
