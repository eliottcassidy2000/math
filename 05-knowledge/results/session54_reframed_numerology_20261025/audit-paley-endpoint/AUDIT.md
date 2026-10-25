# Independent Paley endpoint audit

Verdict: PASS in the producer's explicit prime scope. No blocking mathematical issue. This is a source-bound elementary quadratic-order theorem plus finite-exact computations, not a modular coefficient law or literature novelty claim. The producer packet is frozen under inputs/session54-paley-endpoint-observer. No repository files were edited and no Git operations were performed.

The corrected 1729 reflection and MISTAKE-658 only refute the proposed uniform split-prime interpretation; they do not forbid other arithmetic carriers. The prior Session25 source already supplies (a,b)=(0,1),(1,8),(628,1101) at q=3,7,11. All seven source-copy hashes match both the producer source pins and the source origins inspected at audit time. Classical Rédei parity is used as an inherited theorem; the abbreviated alternative proof-route claims in its historical source are not independently endorsed.

## Mathematical audit

The square-affine group has order q(q-1)/2 and is regular on ordered arcs and on reverse arcs. The endpoint matrix therefore is E=aA+bA^T. Each path has distinct endpoints. A free action on ordered Hamilton words gives H=q r(a+b). The identified full prime-Paley automorphism group is inherited, not reproved here; the operator theorem needs only the square-affine subgroup. Because q and r are odd, Rédei's odd H gives odd s=a+b and odd nonzero delta=a-b.

The character identity S²=J-qI for S=A-A^T implies the rational decomposition into the constant line and augmentation subspace. E=(s(J-I)+delta S)/2 has eigenvalues rs and (-s±delta sqrt(-q))/2, with multiplicities 1,r,r. Its quadratic norm N=(s²+q delta²)/4 is integral and positive, and its quadratic discriminant -q delta² is negative. Thus Q[E]=Q×Q(sqrt(-q)), while on augmentation the generated integral order is Z[lambda]=Z+|delta|O_K, conductor |delta| and discriminant -q delta². This is an order in a varying imaginary quadratic field, not an elliptic-curve construction or a cross-prime modular coefficient identity.

Closed paths are exactly the reverse-arc endpoint channel, yielding C=r b. Translation acts on unrooted directed Hamilton cycles; a fixed cycle's successor permutation commutes with translation, hence is x↦x+d. Exactly r square steps give fixed cycles, so C≡r mod q and b≡1 mod q. No inference a≡1 or q|delta is valid: q19 gives a≡3 and delta≡2. Neither endpoint strict positivity nor spectral recovery of orientation is asserted: q3 has a=0; transposition preserves the spectrum and reverses the signed delta.

## Independent finite computations and integer safety

independent_count.cpp uses a separately written full-vertex-mask predecessor recurrence, with start 0 fixed. It does not call or reuse the producer recurrence. It agrees with every endpoint coordinate for q=3,7,11,19. At q19, a=3040229428, b=3817640437, H=1172695746915 and C=34358763933. Its peak DP allocation is 79,691,776 bytes, distinct from the producer's reduced-mask 37,748,736 bytes. Every partial DP count is at most 19!, hence below 2^64. This bound applies to counting, NOT to every spectral invariant: the q19 M and determinant exceed uint64 and require larger arithmetic, which the Python verifier uses. Exponential complexity and the q<=19 certificate limit remain explicit; q23 was not run.

Both producer replay modes match the pinned receipt. The independent normal and Python -O receipts agree byte-for-byte. Mutation, missing-file and extra-file rejection tests pass in both modes. Mathematical hostiles separately reject a full integral direct product, false q19 congruence, scalar substitution for the quadratic eigenvalue, and all-positive endpoint claims.

## Exact full integral refinement (separate proof)

Let t=rs, Q(X)=X²+sX+N, and M=Q(t)=q(qs²+delta²)/4. In the rational product identify E with (t,lambda). The image of Z[E] is exactly

{(u,A+B lambda) in Z × Z[lambda] : u≡A+B t (mod M)}.

Proof: divide any integer polynomial P(X) by the monic Q(X). Its remainder A+BX is integral and P(t)=A+Bt+Q(t)R(t), so the congruence is necessary. Conversely, if u-A-Bt=kM, the polynomial A+BX+kQ(X) maps to the desired pair. The quotient map from Z×Z[lambda] to Z/MZ sends (u,A+B lambda) to u-A-Bt; it is an additive surjection with precisely this kernel. Therefore the index is M. Equivalently, the basis 1,E,E² has coordinate determinant M. This is a fiber product of rings, with Z[lambda]→Z/MZ given by lambda↦t, well-defined because Q(t)=M.

The moduli for q=3,7,11,19 are 3,1078,91045845,4247362459861422810691. Since M>1, the full generated integral algebra is NOT the full product Z×Z[lambda]. In particular, (1,0) is excluded already at q3. The rational direct product and augmentation order remain completely correct. This refinement is an audit-side corollary, not a modification of the frozen producer.

## Closed replay

Run python3 audit.py and python3 -O audit.py, comparing output with normal.json and optimized.json. Run python3 hostiles.py and python3 -O hostiles.py, comparing with hostiles.normal.json and hostiles.optimized.json. Run python3 verify_manifest.py to verify exact closed inventory and all SHA-256 pins. These scripts require only the packet, Python, and a C++17 compiler. A separate clean-directory replay receipt records copying only this packet to a temporary directory and reproducing all four outputs.
