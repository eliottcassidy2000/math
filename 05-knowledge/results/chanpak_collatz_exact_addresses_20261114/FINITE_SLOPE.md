# Terminal-zero cut transfers do not support a finite menu of exact slopes

**PROVED, scoped.** A finite set of exact affine source/reader templates at fixed spine count can produce only logarithmically many different counts, if the moves are transitive suffixes transporting templates *as lines*. This rules out a particular all-size digit-ready architecture, not the full one-apex family, sampling a line at isolated points, or genuinely changing cores.

## Set-up and inheritance

Use the [four-coordinate cut transfer](../chanpak_selective_target_20261114/CUT_STATE.md) and the [fixed-cycle congruence theorem](../chanpak_cross_problem_theorems_20261114/FIXED_CYCLE_CONGRUENCES.md). At a terminal-zero seam the signed bias is zero. Fix a positive spine count h. The two variable coordinates are S (start weight) and N (full path count). A family might be described by exact readiness relation N=rS+c, in which individual sources realize some discrete parameter values along the line. For this argument a transition must transport that relation identically as an affine formula, not merely succeed at one chosen source state.

Appending a bit 1 doubles S and adds h, makes the bias h, and leaves N unchanged; appending 0 doubles S, adds S plus the former bias to N, and resets the bias. In particular, for a word w of length L starting at bias zero,

    S' = 2^L S + h s(w)
    N' = N + t(w) S + h p(w),

where s,p depend only on the word and
```
t(w) = sum over zero positions i of 2^(i−1), with positions numbered from 1.
```

The formula for t is also immediate by differentiating with respect to the initial S: the derivative of S just before position i is 2^(i−1), and only a 0 uses it in the N update. Thus at a specified length L, t determines the *entire physical word*; its binary digits record the zero positions. The same conclusion holds with any conditions restricting the permitted words.

## The finite-template theorem

Take a finite seed list with S>0 and bias zero, a finite slope set R, and any collection of transitive terminal-zero blocks. Require each block to transport an exact affine relation of slope in R into an exact affine relation of slope in R. Intercepts may range arbitrarily. For a concatenation of total length L, telescoping the slope laws gives

    r_out = (r_in + t(w))/2^L.

For a specified pair of start and end slopes, this forces t(w)=2^L r_out−r_in. At most |R|² complete words of any given length can be used, regardless of the way the word was split into a succession of moves or how many intercepts were considered. For a physical nonempty word ending 0, the last update adds at least 2^(L−1)S_initial to the count because S has at least doubled at every step and the bias at that step is nonnegative. Other count increments are nonnegative too. Thus H≤X gives L=O(log X), uniformly over the finite positive-S seed list. The number of exact outputs ≤X is O(|seeds|·|R|² log X), including the empty choices independently.

No rationality of r is needed. As a special case, any nonempty words preserving the *same* slope have t/(2^L−1)=r. The repeating binary fractions 0.overline(d), d the high-to-low digits of t padded to L, coincide. Uniqueness of the infinite nonterminating periodic binary expansion implies that all such words are powers of one primitive word (the all-zero and all-one streams may be handled directly; all-zero t is not allowed for terminal zero). This is about slopes only; the affine intercept imposes additional restrictions.

## What changes the conclusion

This is a finite-profile no-go for a particular exact-invariance proposal, distinct from the [57-generator ray result](../session236_chanpak_offsets_orbits_20261109/last_letter/PROOF.md). The obstruction is a reverse-address uniqueness, akin to a Collatz itinerary which fixes the branch when both ends and the total clock are given. It places no restriction on all changing-core seams, a nonlinear invariant, an unbounded slope deck, or target-tested isolated line intersections. In particular, changing a finite set of offsets cannot by itself repair this strategy: intercepts never entered the proof. Any successful target descent has to bypass the one-address-per-length bottleneck.
