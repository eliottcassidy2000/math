# Independent audit: planar affine cancellation and the positive-word boundary

2026-11-10. **PASS after notation clarifications.**

Reviewed: `../planar_actions/PROOF.md`. This audit does not import or call the producer checker. The source proof was read directly and compared with THM-5085's inherited base-x physical transfer.

## Verdict and scope

The affine commutator, Hamiltonian formulas, forbidden group-orbit outputs, positive-region obstruction, and positive-word normal form are sound. They do not supply a tournament realization of inverses or prove Chan–Pak completeness. The permanent-gap interpretation relies explicitly on THM-1370 and THM-4940; those all-order exclusion proofs were not independently re-audited here.

Two clarity repairs were requested and applied by the producer:

1. The opening names C=BS while declaring rightmost-first matrix products. Since the intended chronological steps are B then S, write C=S∘B, or label BS explicitly as a chronological word rather than a matrix product.
2. The binary code for C is01 in ordinary high-to-low positional notation: C alone has n=2,m=1. The printed10 needs an explicit low-to-high convention and reversed block reading. Prefer A=0,C=01, with chronological concatenation giving the padded high-to-low binary expansion of m. The normal-form formula and uniqueness argument themselves are correct.

## Independent proof review

On Delta=0, THM-5085 gives B:(h,R,0)->(h+xR/2,2R-d,-d), then S yields(h+xR/2,4R-2d,0). S alone gives(h,2R,0). Deleting the root appends only sinks, so d stays fixed. For x>0 every packing polynomial has nonnegative coefficients, and the singleton rooted path contributes d to R; after a sink R>=2d>0. This justifies the region for arbitrary positive real x, not only bases2 and3.

Set k=x/2. Directly composing C A^-2 gives D(h,R)=(h+kR/4,R-2d). Directly composing A D A^-1 D^-1 gives K(h,R)=(h-kR/8-kd/4,R-2d). The elementary two-step nilpotent affine subgroup law then yields [D,K](h,R)=(h-3kd/4,R). Inverses reverse that sign. The Hamiltonian equations R'=-2d integrate exactly to D,K using the printed linear terms; their Poisson bracket has derivative-3kd/4. A,C themselves have determinants2,4, so they are correctly not called symplectic.

At d=1, J^-8 sends(1,2) to(1+3x,2), while J^-24 sends(1+x,6) to(1+10x,6). The two seeds are physical after a sink. At x=2 these yield7,21; at x=3 they yield10,31. They are positive integral terminal pairs on the selected sheet but not thereby physical states. Every nonempty positive word strictly increases R in the stated invariant region, already excluding those unchanged-R targets from that positive orbit.

For the normal form, chronological append A maps(n,m) to(n+1,2m), and append C maps it to(n+2,4m+1). Thus high-to-low padded bits concatenate0 or01. No adjacent1 bits and highest position<=n-2 are exactly the parsing condition: pair each1 with the zero immediately preceding it; unmatched zeros are A. The slope2^n and intercept-2dm distinguish all maps when d!=0. For initial R=2d, the intervals in equation(5) are disjoint even for successive n, as the next lower endpoint exceeds the previous upper endpoint. Hence that seed's final R uniquely determines its positive word. This is not freeness of arbitrary S/B operators or all physical tournament histories.

## Reproducible independent controls

Run from this directory:

    python3 verify.py
    python3 -O verify.py
    python3 verify.py --inject-failure
    python3 -O verify.py --inject-failure

Both normal runs passed with byte-identical JSON: 8,343 explicit checks, 1,596 positive words of total weight<=14, and48 reconstructed physical instances. Both injected-failure runs exited nonzero using explicit RuntimeError, independently of Python assertions.

- Rational commutator controls cover x in{1/3,2,3,7/2}, d in{1,3/2}, multiple h and R, including outside the positive region. D,K,J and inverse-J signs are checked from fully expanded chronological A/C words.
- Word normal forms, completeness of admissible bit masks, seed injectivity, and interval inequalities are checked through weight14.
- Physical instances start from a singleton or directed triangle and append an initial sink, then all A/C words of weight<=4, in bases2 and3. The audit constructs actual tournament adjacency, enumerates directed odd cycles, computes disjoint-cycle packing polynomials by subset recursion, and independently enumerates rooted directed paths and their complement packing counts. No transfer recurrence is used to compute the actual states.

Finite controls support, rather than replace, the uniform algebraic and combinatorial proofs. The checker does not assert that all physical-state inequalities hold at forbidden output points. No external literature novelty is claimed.

## Padded cone-hostile supplement

The revised proof adds A^3 padding to the singleton seed and A^2 padding to the triangle seed. Since J is horizontal it commutes with A; the outputs become(1+3x,16) and(1+10x,24). Independent exact checks confirm both outputs and h<=xR/2 for x=2,3. Monotonicity is an obstruction relative to the corresponding padded seed; the original-seed exclusion instead follows from inherited global gaps or the unique R decoder. A final wording request makes that seed explicit.
