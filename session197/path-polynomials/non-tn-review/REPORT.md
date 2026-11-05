# Independent review: an actual context with a negative TN minor

**PASS, no correction.** This binds the frozen `session197-odd-cycle-path-polynomials/` files:

- PHYSICAL_NON_TN_CONTEXT.md:57ca1e3a4e7f2047217d5541383b39e2d078b8935a9f3472970368d7e1b2aea9.
- verify_physical_non_tn_context.py:1de40feb975b542690397bb40f99f9f909c2015c61ef69d9f8c672d2198fc205.
- non_tn_context_normal.stdout and optimized counterpart:1a5e3c14745d2a6162721654160b0736e51119e047c875cbd5fd995fe6eb64d4.

## Proof and scope

The literal context C3[T,Tr3,Tr2] is a tournament substitution with five exterior vertices. Its unique whole-profile linear operator exists by the named reviewed cover-linearity theorem. Consequently the actual input profiles x,x^3,x^5 determine the first three q columns exactly, since q0=x, q1=x^3-x and q2=x^5-5x^3+4x. They are not being used to infer the entire operator from three samples.

The monic triangular change to the r basis gives exactly the displayed columns. Their first three rows have determinant-468310. Total nonnegativity requires every finite minor to be nonnegative, so this single exact negative minor suffices, despite every column coefficient being nonnegative.

Appending a universal sink is multiplication by x. The identity Xr_k=q_(k+1)+(k+1)^2 q_k gives the exact3-by-3 lower bidiagonal left factor with diagonal1,4,9. Higher-index columns of this left factor cannot enter its first three output rows. Thus the padded leading minor is precisely36 times the original determinant, namely-16859160. The padded context is physical and parity-preserving.

I checked the named PERIODIC_STIELTJES_WORDS.md dependency at SHA82b941c9ce53b41ca73ee0224b240fe59c61df701df5fde6f9d29283f0fad166. Its Section4 gives TN column matrices for each cross-parity X operation and for D in either parity, and hence for every parity-returning D/X word. The padded context therefore cannot equal such a word as a full-profile operator. This conclusion does not depend on a choice of initial vacuum. It does not imply a negative Hankel minor for an output sequence, a scalar-count obstruction, or minimality of the physical counterexample. The source states these distinctions correctly.

## Independent physical verification

The standalone `check.py` imports no source helper. It generates the actual cyclic-block graphs directly. For each graph it enumerates directed odd cycles as literal words anchored at their least vertex, then enumerates disjoint cycle packings by the least remaining vertex. This computes the entire physical B polynomial independently of the source's path-cover computation and cyclic-skeleton operator.

It checks eleven actual profiles: every clone length0 through6, and the separately constructed padded graphs at lengths0,1,3,5. It verifies both exact singleton row lists, the three physical monomial columns, both parity-basis coordinate changes, every displayed coefficient, both3-by-3 determinants and all coefficients' nonnegativity. The zero-clone graph is handled as the actual ordered join Tr3->Tr2.

Normal and Python -O runs produce identical receipts:

- check.py:be6bb6788d30077bd37ae0f87d383634886c1d1a126d890111251e7fc1c1ec79.
- independent-normal.json and independent-optimized.json:b9c9a5ff7527a595cd6c4b2ab2b1edf9782ec5caec6822ee4a34f8a8592df662.

The source checker and its two named Python helpers were copied into an audit-only directory before replay. Both source normal and optimized outputs reproduce the frozen receipt exactly. Original files and dependency hashes remain unchanged.

Reproduce the independent audit with `python3 check.py` and `python3 -O check.py`. The negative minors are exact integers, not floating-point numerical tests or a finite-census extrapolation.
