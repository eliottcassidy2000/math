# Independent Session 27 lattice-controls audit

Date: 2026-10-22. Status: **PASS for the review's specified finite computations and cautions; no extremality certificate.** No repository operations, external downloads, or downloaded executable code were used. Source PDFs and all reviewed code/output are bound in `INPUT_SHA256SUMS.txt`. The three directly relevant PDFs were freshly extracted using pdftotext; all three resulting texts match the review's corresponding text byte for byte.

## 1. Exact rank-43 collision: PASS

The review's formulas agree with compression PDF equations (2.2)–(2.14), including the factors 4/3 and 8/3 and exceptional E6/E7/E8 values. Re-execution of the inspected symbolic controls verifies:

- R1 = 2A1 + 2A4 + A5 + 2A7 + A8 + D6; R2 = 2A2 + A3 + 2A6 + A9 + 3D5. Both have rank 43.
- (c1,c2,c3,c4) = (318,−4110,89608,−2433804).
- c5(R1)=93582144/5, c5(R2)=94375104/5, difference R2−R1 = 158592.
- Coxeter-rank masses d2,…,d10 = −2,4,3,−8,−5,12,1,−8,3 yield P(z)=z²(z−1)^5(z+1)(3z²+4z+2). Hence all power sums through degree 4 vanish, while the fifth power sum is 2160.
- Defining δj = (A per-rank value) − (D per-rank value), δ4(h)=−10h(h−2)(h−4). In particular 15δ4(8)=6δ4(10)=−28800, so the D correction to R2−R1 is zero. δ5(h)=(28/5)h(3h³−24h²+67h−58); the A-proxy fifth contribution is 72576 and D correction is 86016.

This explains the exhibited collision, not its global minimality or uniqueness. These are polynomial invariants of formal ADE multisets, unless their claimed identification with actual frame cumulants is separately established.

## 2. Exhaustive census and uniqueness: PASS

I inspected, compiled independently with warnings enabled, and reran `audit_ade.cpp`. Its full output matches the saved output byte for byte. Its recursion lists multisets of the 86 types A1..A43, D4..D43, E6/E7/E8, with nondecreasing type index; repeated copies are allowed by recursing at the same index. Sorting only by rank still supplies a fixed total indexed order among equal-rank types, so there is no omission/duplication from ties. Types above remaining rank are correctly skipped by the sorted-rank cutoff. Signatures are exact signed 64-bit integer (c1,c2,3c3,c4), not hashes. Conservative triangle bounds verified in `independent_controls.py` are safely below 2^63. Sorting exact signatures and counting adjacent duplicates therefore proves the reported multiplicities.

A second, independently implemented Euler-product dynamic program counts multisets without any signature recursion. All 43 individual rank counts agree with enumeration: ranks 1..42 contain 5,565,832 configurations; rank 43 contains 1,556,185; total 7,122,017. No signatures repeat within a rank below 43. Rank 43 has exactly one adjacent duplicate, which means exactly one doubleton fiber and no tripleton. The independently verified pair above has that signature, identifying the unique collision. This conclusion is specifically within fixed rank and the stated universe, not across unrestricted ranks or other root systems. The zero-rank empty multiset is intentionally absent from the totals.

## 3. Formal weight-48 theta series: PASS; old paper coefficient wrong

The exact truncated q-series control reruns successfully. The triangular basis E4^(12−3j) Δ^j has leading term q^j, so constant term 1 and q1..q4 = 0 uniquely determine coefficients:

1, −2880, 2324160, −465638400, 5611550400.

The q5 coefficient is 565866362880 and q6 is 2972108280960000. The old PDF p2 equation (2) actually prints 184139089776640 at q6, so the review's correction is accurate. The newer PDF's value agrees. Positivity was computed only for q5..q12; this is not a proof of all-coefficient positivity and no coefficient calculation proves realization by an extremal lattice.

## 4. Old96 matrix and minimum controls: PASS with important limits

Exact determinant and rational Schur/LDL elimination were rerun. The inert matrix data hash is:

7a6a3b8f761a9f02e1a687be8040aaf403c7ac75c33d92f010f2dd131f17bc24

The 96×96 integer matrix is symmetric, has even diagonal, determinant 1, and all 96 exact positive rational pivots. Thus it certifies a positive-definite even unimodular rank-96 lattice. Numerical eigenvalues play no role in that proof. First pivots are 28,21,56/3; last are 63/1742,22/567,1/22. Output matches the prior matrix output byte for byte.

A newly written independent C++ exhaustive sparse control checks every vector with entries 0,±1 and support at most 4, up to overall sign. Results by exact support:

- support 1: 96 candidates, minimum 20;
- support 2: 9120 candidates, minimum 20;
- support 3: 571520 candidates, minimum 24;
- support 4: 26575680 candidates, minimum 20.

These strengthen the reproducibility of the review's sparse result, but not its scope. They do not cover coordinates of magnitude at least 2, even at small support, or support exceeding 4. They cannot prove minimum 10, or rule it out.

Both saved heuristic witnesses were independently evaluated using unbounded Python integers: rounded-eigenvector witness norm 28 and LLL witness norm 20, each with maximum coordinate magnitude 1. The stored LLL output reports only 95 steps and its witness is a basis vector; it should not be advertised as a materially improved lower-minimum search or an exact reduced-basis certificate. Neither heuristic supplies a lower bound. A basis vector of norm 20 supplies an upper bound only. Combined with the standard even-unimodular rank-96 upper bound, the lattice's still-unresolved minimum lies in {2,4,6,8,10}. No complete shortest-vector enumeration, qfminim verification, or new96 matrix reconstruction was performed.

## 5. Orthogonal versus symplectic gluing: PASS, conditional split hypothesis

Old PDF p3 equation (10) explicitly describes a nondegenerate symmetric discriminant bilinear form, but p4 equation (12) uses the symplectic maximal-isotropic count. In odd characteristic a nondegenerate symmetric form cannot also be alternating. For a split orthogonal 2m-dimensional space over Fq, the total number of maximal totally isotropic m-spaces (both families) is 2∏(q^i+1), i=1..m−1. For m=24,q=7 this is exactly the review's corrected expression. The corresponding count in a fixed isotropic line's quotient gives inclusion probability 1/(7^23+1), denominator 27368747340080916344. Independently, double-counting lines per generator versus all isotropic lines yields the same ratio. Nonisotropic lines have probability zero. If the form were nonsplit, dimension-24 isotropic spaces would not exist; the stated count is explicitly conditional on split type, not an independent certificate of that type.

This correction alone would not destroy a genuinely established bound of 10^12 forbidden isotropic lines: the corrected union-bound estimate is still approximately 3.65×10^−8. The actual obstacle is that the paper has not supplied the required forbidden-set bound/certificate or justified its pure-tensor exclusion argument.

## 6. Schottky normalization: PASS, actual textual inconsistency

Compression PDF p9 section 5.2 explicitly defines J8 = Θ(E8²) − Θ(D16+) and immediately states Fh=(15mh/4)J8Θ(E8), with m30=896. For the h30 pair E8³ and D16+⊕E8, theta multiplicativity gives F30=Θ(E8)J8 with the corresponding orientation (or its negative if the pair is reversed), never 3360 times that expression. Since 15·896/4=3360, the normalization mismatch is real. To preserve the printed positive scalar formula for the stated orientation one would need J8(new)=J8(theta-difference)/3360, or instead replace 15mh/4 by mh/896. This identifies the normalization needed; it does not independently prove all five proportionality constants.

## Suggested presentation refinements, not changes to core results

1. Name the sparse test explicitly as a {0,±1}-coordinate test, not an arbitrary sparse-integer-vector test.
2. If mentioning LLL, say it found only an existing norm-20 basis vector and supplies no minimum certificate.
3. Keep split orthogonal counting conditional unless an exact Witt-index/gluing certificate is separately supplied.
4. Keep the formal ADE/frame distinction and the unique-collision/exhaustive-census distinction explicit.

No correction is required to the review's collision arithmetic, census totals, exact old96 Gram checks, theta q6 correction, orthogonal counting correction, or factor-3360 diagnosis.
