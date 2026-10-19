# Independent audit: Berggren norm-one tree and absorber dictionary

Date: 2026-10-19. Result: **PASS**, subject to the minor leg-swap wording qualification below. This audits the mathematical claims directly; it makes no historical-priority claim.

## Proof checks

1. **All three matrices preserve the shell.** Exact integer matrix multiplication verifies Mᵀ diag(1,1,−1) M = diag(1,1,−1), and the displayed integer inverses are correct. Each row preserves oddness. For positive norm-one triples, c≥a,b, making every child positive and strictly increasing its third coordinate.
2. **Positive, strictly decreasing parent.** At any nonroot odd norm-one triple c≥3. The short uniform estimate

       9c² − 4(a+b)² ≥ 9c² − 8(a²+b²) = c²−8 > 0

   gives 3c>2(a+b), hence h>0. Also c<a+b, because (a+b)²−c²=2ab+1>0, so h<c. This avoids any boundary/interior case split.
3. **Exactly one inverse is positive.** r and s are odd, so neither vanishes. If a≥b then (2a+b)²−4c²=4ab−3b²+4>0, giving s>0; if b≥a, the symmetric calculation gives r>0. Thus the sign patterns are precisely (+,+), (+,−), (−,+), selecting B⁻¹, A⁻¹, C⁻¹ respectively. Descent terminates only at (1,1,1). It proves completeness, uniqueness, and absence of branch collisions.
4. **Boundary classification is exact.** a=1 forces b=c, and b=1 forces a=c. The branches are A^k(1,1,1)=(1,2k+1,2k+1) and C^k(1,1,1)=(2k+1,1,2k+1), with the root shared. No other address lies on either boundary. Once both legs exceed 1, their three children do too.
5. **Absorber criterion is exact for distinct positive x,y.** For x+y=t², the four sums are (at)²,t²,(bt)²,(ct)². When a,b>1, both hosts exceed max(x,y). Their difference is (a²−b²)t²−(x−y), which cannot vanish: if a=b it is y−x; otherwise |a²−b²|t²>|x−y|. If a=1 then u=y, and if b=1 then v=x, so either boundary is necessarily degenerate. Global backbone membership is an additional requirement, correctly excluded from the theorem.
6. **The PPT orientation matters.** The rooted (3,4,5) Berggren tree enumerates primitive triples with first leg odd and second leg even, not both possible leg orders. The draft states this correctly. The common-address correspondence is a rooted labeled-tree isomorphism, not a global linear map of quadratic shells.
7. **All diagonal nodes, not just some, are on the B ray.** At a=b>1, r=s=3a−2c>0, since c²=2a²−1<(3a/2)². Therefore the unique parent is again diagonal and is reached by B⁻¹. This proves every diagonal norm-one node has a pure B address, giving exactly the stated negative-Pell lane.

## Minor wording qualification

In section 3 of the draft, “Swapping legs corresponds to a different address in the second tree” is literally false at diagonal nodes. Use “For a≠b, swapping legs corresponds to a different address in the second tree.” Equivalently, leg exchange interchanges A and C throughout the address and fixes B. Pure B addresses are fixed.

No theorem statement needs correction.

## Seed identities

All are verified exactly:

- A³(3,4,5)=(9,40,41), and A³(1,1,1)=(1,7,7).
- B(9,40,41)=(171,140,221), and B(1,7,7)=(29,23,37).
- B(1,1,1)=(5,5,7), B²(1,1,1)=(29,29,41).
- The AP construction 2(c²−4ab),2c²,2(c²+4ab), applied to (9,40,41), gives (482,3362,6242).

Words in the verifier are chronological: AAAB means apply A three times, then B, i.e. the matrix product BA³.

## Reproducible finite checks

Run:

    python verify_berggren_absorber_bijection_20261019.py

The standard-library-only script uses exact integers and emits JSON. The saved receipt is `berggren_absorber_bijection_controls.json` and includes the script SHA-256.

Default run checks:

- Exhaustive independent enumeration through c=2,001: 5,395 odd norm-one triples, exactly matching the bounded generated tree.
- Exactly 2,001 boundary nodes and 3,394 nondegenerate coefficient triples.
- Independent enumeration of 319 parity-oriented primitive Pythagorean triples through the same height.
- Every one of 9,841 addresses through depth 8, with unique descent and the same address recovered on both shells.
- 370 distinct ordered positive square-sum edges with square root at most 10, for every nondegenerate coefficient triple: 1,255,780 absorber checks.
- Every diagonal node in the height bound is precisely on the central B ray.
- All requested seed identities and the relevant matrix/inverse identities.

These finite checks audit the proof and implementation; the infinite classification is established by the descent argument, not by enumeration.

## Extension: the complete positive integer norm-one forest

The extension is also **PASS**. Every ordered positive integer solution of a²+b²=c²+1 belongs to exactly one of three ternary trees, rooted at

    (1,1,1), (2,1,2), (1,2,2).

Their parity sectors are respectively odd/odd/odd, even/odd/even, and odd/even/even. These are the only possibilities modulo 4. Each matrix is the identity modulo 2, so the sectors are invariant.

The same parent proof works for every nonroot. The estimate h>0 holds whenever c≥3, without any parity assumption. The only points with c<3 are precisely the three roots. The possibility of a zero inverse coordinate is also isolated exactly: r=0 forces a=2k and

    k(3k−2b)=1,

so k=b=1 and (a,b,c)=(2,1,2). Symmetrically, s=0 occurs only at (1,2,2). Thus all nonroots have nonzero r,s, a unique positive inverse, and smaller positive c. The roots have h=−1,0,0 respectively and admit no positive inverse.

The complete degenerate boundary consists of four rays:

    A^k(1,1,1)=(1,2k+1,2k+1),
    C^k(1,1,1)=(2k+1,1,2k+1),
    C^k(2,1,2)=(2k+2,1,2k+2),
    A^k(1,2,2)=(1,2k+2,2k+2).

The first two share their root. There are no further boundary nodes. Therefore the full forest enumerates every ordered positive integer coefficient triple for the stated absorber form, with the same exact nondegeneracy criterion a,b>1. The whole forest is naturally in bijection with three labeled copies of the parity-oriented primitive-Pythagorean tree.

### Same-address absorber for the requested seed

The prior odd-tree boundary limitation is sector-specific:

    A³(2,1,2)=(8,31,32).

Thus the seed (9,40,41) has a nondegenerate coefficient triple at its unchanged AAA address in the even/odd/even sector. Its hosts are

    u=63x+64y, v=961x+960y,
    u+x=(8t)², y+v=(31t)², u+v=(32t)².

More generally,

    A^k(2,1,2)=(2k+2,2(k+1)²−1,2(k+1)²).

This is nondegenerate for k≥1, and the corresponding C ray in the other even sector is its leg swap. Every nonempty address gives a nondegenerate absorber in at least one sector: pure A or pure C addresses give exactly one; all other nonempty addresses give three. The empty address gives none.

### Extended exact receipt

The updated verifier and same-named JSON receipt check:

- All 10,987 positive norm-one triples with c≤2,001 against an independent exhaustive search.
- Sector counts 5,395, 2,796, and 2,796.
- All 4,001 degenerate boundary nodes and 6,986 nondegenerate coefficients.
- All 9,841 addresses through depth 8 in each of the three norm-one trees and the primitive-Pythagorean tree.
- 2,584,820 absorber checks over all sectors, including the earlier 1,255,780 odd-sector checks.
- Exact classification of r=0, s=0, and h≤0 in the bounded search.
- The same-address even-sector seed and both explicit opposite-boundary polynomial rays.
