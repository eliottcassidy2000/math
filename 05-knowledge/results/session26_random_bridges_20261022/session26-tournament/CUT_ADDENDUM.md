# One-apex counts as directed cuts and Ferrers hook weights

Session26. PROVED below; independent audit requested. No literature-priority claim. Existing source/correction search was interrupted by a shared executor outage, so this packet does not assert that the extrema are new to the repository or literature.

Let w=b1...bm encode a tournament with one marked apex and a transitive deletion, as in PROOF.md. Define c1=2 and cd=2^(d-1) for d>=2. Then

    H(w)-1 = sum_(j<i) bj(1-bi)c_(i-j).

This follows from Session25's exact path formula: the boundary contribution b1-bm plus the number of adjacent01 transitions equals the number of adjacent10 transitions. Thus the usual inversion weight2^(distance-1) receives one additional unit precisely at adjacent10. The expression is exactly a directed cut capacity on the forward complete DAG, with the positions carrying bit1 on the source side. Strong words have b1=1,bm=0, so it is an s-t cut.

## Minimum and equality

For m>=3 the mandatory pair1,m contributes2^(m-2). Each interior position j contributes c_(j-1) if bj=0, or c_(m-j) if bj=1, through its pair with the opposite endpoint. All remaining interior inversions have positive weight. Hence

    H-1 >= 2^(m-2) + sum_(j=2)^(m-1) min(c_(j-1),c_(m-j)).

Since the weight sequence2,2,4,8,... is nondecreasing, the pointwise minima are attained by initial interior zeros followed by interior ones; this also kills every interior inversion. For m>=5 the minimizing words are exactly

    1 0^a 1^b 0,  a+b=m-2, |a-b|<=1,
    H_min=(2^a+1)(2^b+1).

Equality forces no interior inversion and then the balanced threshold, proving necessity as well as sufficiency. The small plateau cases must be separate: m2 gives10,H3; m3 gives100 or110,H5; m4 gives1000,1010,1110,H9. The remaining m4 word1100 hasH11. A zero exponent must not be substituted into the product formula.

## Maximum and equality

Swapping adjacent01 to10 strictly raises the cut by at least2. Every earlier1's edge to the local0 gets longer; every later0's edge from the local1 gets longer; other pairs are unchanged, and the new adjacent inversion contributes2. Consequently, for fixed positive counts a ones,b zeros, the unique maximizing word is1^a0^b. Its value is

    H=2+(2^a-1)(2^b-1)=2^m-2^a-2^b+3.

Over all strong words of lengthm>=2, the maximum is attained exactly when a+b=m and |a-b|<=1. The global maximum is2^m-2^floor(m/2)-2^ceil(m/2)+3.

## Exact Ferrers map

Rows are indexed by ones left-to-right, columns by zeros right-to-left. The cell at a pair exists precisely when the one precedes the zero. Row lengths are nonincreasing. Its arm counts intervening zeros and its leg intervening ones; therefore its hook length equals the positional distance. Thus

    H=1+sum_cells 2^(hook-1)+number_of_corners.

Strong words correspond exactly to diagrams in their a-by-b rectangle having a full first row and a full first column. The smallest such shape is the hook(b,1^(a-1)); the largest is the rectangle. Adding a Young-lattice cell increases the positive hook statistic, exactly reflecting the adjacent swap. This is not the standard hook-length tableau-count formula. Shape[3] (word1000) and shape[2,1] (word1010) both giveH9 although their hook multisets differ. ScalarH therefore loses marked shape data; no distinct-unrooted-isomorphism assertion is made about this pair.

## An all-order absent value and infinitely many gaps

Remove leading zeros and trailing ones from any word. These are source/sink vertices and do not changeH. The remainder is strong or transitive. Therefore all feedback-vertex-at-most-one counts equal1 or a strong word's count. The small lists above and the minimum15 fromm5 show that13 never occurs anywhere in this family.

Let A(X) count distinct H<=X in this whole family. For k>=4 put R=floor((k+1)/2). Then

    A(2^k) <= 7*2^(k-4)+2^(R-1).

Indeed any strong core has lengthm<=k+1, since H>=2^(m-2)+1. There are2^k-1 strong words of lengths2..k+1. Reverse-complement preservesH; its fixed words occur only at even lengthm=2r and number2^(r-1), totaling2^R-1. Thus the number of such orbits, plusH1, is2^(k-1)+2^(R-1). At top lengthm=k+1 the2^(k-3) words starting11 and ending00 are a reverse-complement-invariant excluded set: four forced cross pairs contribute9*2^(m-4)=9*2^k/8>2^k. They remove at least2^(k-4) orbits, proving the bound.

Consequently at least2^(k-4)-2^(R-1) odd integers through2^k are missing. There are infinitely many holes, despite Session25's complete modular coverage even within a subfamily. For arbitrary X, monotonicity of the missing-count function between dyadics gives at leastX/32-O(sqrt(X)) missing odd integers. This is a lower-density bound, not existence or evaluation of a natural-density limit.

The cut and Ferrers maps preserve actualH but not individual Hamilton paths, general insertion/substitution ports, or merged-metagraph wiring. None of these results settles the full tournament H-spectrum.

## Controls completed during the executor outage

Pure-JavaScript BigInt checks completed: all8191 words of lengths0..12 satisfy the cut identity; all20481 adjacent01 swaps increase capacity by at least2; all511 words of lengths0..8 agree with independently implemented vertex-level subset DP. Every Ferrers cell's hook equals its positional distance throughout the8191-word universe, and hook-plus-corner counts agree. Minima, maxima and equality words agree throughm12. Parent independently checked strong words throughm16. The combined standalone verifier is now saved as verify_cut.js, with matching cut-output.json; run node verify_cut.js. The later product-law review in ASYMPTOTIC_REVIEW.md records a stronger asymptotic missing-value bound.
