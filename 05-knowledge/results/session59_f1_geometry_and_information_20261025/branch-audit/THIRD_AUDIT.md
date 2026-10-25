# Session 59 independent audit: distinctness (third gate)

**Verdict: PASS**, 2026-10-25. No mathematical correction needed. This audit includes the final direct human proof, the doubled-image existence criterion, and the optional topology. No Git operations were performed.

## Frozen targets

- `distinctness-addendum.md`: SHA256 `407bd205a2f5bd5629dc6261db24a69f11dbccc3dace5fc7d3d0729111afd9e1`.
- `check_distinct_magic.py`: SHA256 `92e70e1536257dda3ad162abc7a3e1c7b020956f94a0b75f5bc51731df72fb0f`.
- `distinct-certificate.json`: SHA256 `ca67fd07c3463cdee9120444c8e9cf061aa44a08a3d40a7c0584703f85a1e2c4`.

All targets are in `/tmp/session59-f1-wildcards/`.

## Direct proof: all finite abelian groups

The direct proof is valid independently of any computation. Write B=2A, M=|B|. The map A[4] to B[2] given by doubling is surjective: if b=2a and 2b=0, then 4a=0; its kernel is A[2]. Consequently |B[2]|=h4/h2.

To make the four coefficient-grid offsets +/-x +/-z distinct, u=2x and v=2z must avoid u=0, v=0, u+v=0, and u-v=0. For nonzero u, exclude two values of v when u has order 2 and three otherwise. This gives M^2-4M+2+|B[2]| choices. Multiplying by the h2^2 lifts yields C=N^2-4h2 N+2h2^2+h2 h4.

The remaining collision equations are exactly z=+/-2x and x=+/-2z. On each, corner-goodness excludes A[4] union A[6]; their intersection is A[2], so each residual line has N-h4-h6+h2 good points. Pairing the two signs of a fixed orientation forces 4-torsion. Same-sign cross-orientation intersections force 3-torsion. These four intersections fail the corner conditions. The two opposite-sign cross-orientation intersections instead force 5x=0; every nonzero 5-torsion parameter is corner-good because 2,4,6 are units modulo 5. These two intersections have h5-1 points each, are disjoint on the corner-good locus, and no triple survives. Thus

`I = C - 4(N-h4-h6+h2) + 2(h5-1)`

is exactly the displayed D1 formula. Multiplication by N for unrestricted centers and N-9 for nonidentity centers is valid because translation preserves distinctness and exactly nine centers are forbidden for each successful pair.

## Independent finite certificates and consequences

- Independently generated all canonical nonzero grid differences, obtaining twelve sign classes. The four deleted primitive characters have kernels contained in their retained doubled counterparts for every abelian group. No invalid division by 2 occurs.
- A new generic-minor implementation certifies all 256 subsets of the reduced eight-character arrangement and recovers the stated Smith table, including the zero signed coefficient of type (1,3).
- `check_distinctness.py` directly checks offset injectivity and every center for 52 finite abelian groups: cyclic orders 1 through 40 plus twelve noncyclic direct products. It agrees with I, N I, and (N-9) I. Twelve additive-field specializations also agree. Normal and optimized outputs are identical.
- The h6=h2 h3 simplification is correct by coprime primary decomposition. Both cyclic parity formulas, period-60 quasipolynomial statement, C25/C5^2 comparison, F9 versus C9 contrast, and q=11 generic/additive comparison are correct.
- In odd characteristic, the eight projective collision directions specialize to four lines in characteristic 3 and six in characteristic 5; otherwise there are eight. In characteristic 2 opposite offsets coincide identically, so the additive injective locus is empty.

## Existence classification and attribution

The exact criterion `N >= 9 and |2A| >= 4` is correct. Necessity follows from the four coefficient-grid corners. Sufficiency exhausts finite exponents: an element of order at least 9 gives x=g,z=3g; prime exponents 3,5,7 supply two independent factors under the size assumption; exponent 4 requires two C4 factors; exponent 6 requires two C3 factors; and exponent 8 with N>=9 contains C8 x C2, with the stated valid witness x=(1,0),z=(2,1). Avoiding the identity changes N>=9 to N>=10.

The distinction between the image 2A and kernel A[2] is correctly maintained. The criterion correctly excludes arbitrarily large C4 x C2^k and C3 x C2^k, unlike an incorrect order/exponent-only shortcut.

Primary source checked directly: [Bowerman, Fleece, Insall, arXiv:2308.01858v1](https://arxiv.org/html/2308.01858v1), Definition 1.1 and Theorem 3.1, with Propositions 2.6, 2.7, 2.15, 2.17, 2.18. Distinctness is part of the source definition. Its Sylow-factor classification matches the stated doubled-image reformulation. The RECOVERED/REFORMULATED label and absence of priority claims are appropriate.

## Optional topology: independent 17-row exhaustive audit

`check_combined_topology.py` independently enumerates all 131072 subsets of the nine entry characters (1,i,j) together with the eight collision characters (0,a,b). Precomputed integer minors and gcds determine every subset's Smith factors. Weighting by 1 or by the product of Smith factors yields respectively:

- ordinary characteristic coefficients (degrees 3 through 0): `1,-17,79,-63`;
- arithmetic characteristic coefficients: `1,-21,164,-504`.

These independently establish both displayed factorizations. Direct coefficient transformations give additive Betti numbers `(1,17,79,63)` and toric Betti numbers `(1,24,209,690)`, and Euler characteristics 0 and -504. The four doubled characters each have two distinct connected components over C; the four extra toric divisors explain the degree-two coefficient change from -17 to -21.

The Poincare specializations were directly verified in [Liu, Tran, Yoshinaga, arXiv:1707.04551v3](https://arxiv.org/html/1707.04551v3), formulas (2.1), (2.2), and Corollary 5.6, during the preceding branch gate. The topology is therefore a valid cited application, not an unsupported inference from numerical resemblance. Full combined-census normal and optimized outputs agree byte-for-byte.

## Final execution and scope

The final author checker, including corner-good and residual-line controls, was replayed in isolated `distinct-replay-normal/` and `distinct-replay-optimized/` directories. Both modes pass, and their certificates agree with each other and the final published certificate. The author suite covers 4096 unreduced subsets, 256 reduced subsets, 75 directly enumerated groups, and 1060 invariant-factor models through order 512. Those model checks corroborate, rather than replace, the existence proof.

No blocker remains. All conclusions concern labelled group-magic arrays or the corresponding open base arrangements. Neither these formulas nor the existence classification produce an integer magic square of distinct squares, a rational section of the ordinary-square cover, a constant rational fiber size, or a broad F1-descent result.
