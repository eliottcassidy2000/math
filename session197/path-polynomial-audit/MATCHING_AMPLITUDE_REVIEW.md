# Independent review: leading matching amplitude

**PASS**, no corrections requested, 2026-11-05.

Frozen source in `session197-odd-cycle-path-polynomials/`:

- MATCHING_AMPLITUDE.md: 443097af71c4d1c364c4bc8ad5d417453df6aea71d734ccb5700b74b93bf48ee
- verify_matching_amplitude.py: e062734acf04e628a5fee3ed15d57123a22bcb6f7a2f1058edd1972f0566773c
- Normal/optimized output: 343c5aee64fb83ae9d6c4173b59f7c10a251ce8d68dc992449aab07ba7a47c81

Both copied-source replays reproduce the frozen output. No source was edited.

## All-order proof

Each of the nu marked clone runs needs an exterior segment with distinct O and I endpoints. For a crossing skeleton with t runs, its exterior size is therefore at least2t. The deficit between total exterior size N and its bin-profile exponent is at least nu. Equality forces no pure exterior cycle, exactly one run and two exterior vertices in every crossing cycle, and every remaining vertex unused. These are precisely maximum matchings, proving degree N-nu and leading coefficient mu.

For a fixed maximum matching, its directed O-to-I pairs are distinct labelled objects. A permutation partitions these objects into cyclically ordered cycles. Joining each pair to the next through a marked clone gap is always legal. Conversely cutting at the marked gaps recovers both the matching and permutation. Thus neither rotation nor unordered-cycle multiplicities have been omitted. The rising-factorial cycle enumerator gives exactly the claimed coefficientwise subfamily. Positivity of the extreme shift layer, already proved separately, prevents cancellation of this lower bound.

When N=2nu, the nu exterior segments exhaust all vertices and each has exactly two vertices. There is no other kind of extreme-shift skeleton, giving balanced equality. The case nu=0 includes the single empty matching and is consistent with the deleted-root profile.

At x=1, the unique largest absolute exponential base is nu+1. Its coefficient is p_nu(1), so it is the leading Stirling run coefficient a_(nu+1). This identification does not require assuming that all negative-base coefficients vanish separately. The ratio limit already identifies the unique largest term. The factorial lower bound follows.

## Independent finite evidence

The new `path-polynomial-audit/matching_independent.py` checks all383 source control templates. It enumerates matchings by subsets of cross edges, separately from the source's left-vertex backtracking. It then computes the highest run coefficient from literal directed path words on every exterior subset, followed by disjoint-subset convolution of initial, internal and final fragments. This is independent of the source's ordinary Hamilton subset DP and Stirling inversion.

All383 literal top amplitudes agree with the sum of the extreme shift coefficients and satisfy mu*nu!. The leading degree/coefficient and balanced equalities also agree. Earlier independent controls reconstructed the full physical operator columns at matching ranks two and three, including the rising-factorial polynomials for the carousel-five, carousel-seven and Paley-seven roots.

The theorem concerns the fixed rooted template's highest clone-response layer. It correctly leaves lower-mode cancellation and inverse target construction unresolved.
