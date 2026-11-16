# Why these carries do not yet fill a larger low digit

Two quantitative bounds distinguish a physical rollover from
a complete bank of remainders.

**Unmarked spacer.** At fixed cut widths and high product P, an
admissible left factor d determines both masks: its partner is
P/d, then b=2^k−d and a=P/d−2^r. Thus at most τ(P) exact low
values occur at that product. The divisor bound τ(P)=P^{o(1)}
is enough to exclude a complete bank of K/2 odd digits when
K=2^t≥P^ε for any fixed ε>0. For completeness, for each η>0
choose a prime threshold beyond which 2≤p^η; for those primes
v+1≤2^v≤p^(ηv). For the finitely many small primes,
(v+1)/p^(ηv) is bounded uniformly in v, giving τ(P)≤c_η P^η.
On the variable-radix branch the high product is of order D²,
and K=D/2, so a full low bank at one fixed product is impossible
just by increasing those factor choices.

**With spacer contacts in the three-sheet branch.** Use the
[variable radix notation](VARIABLE_RADIX.md), and write
H=K m²+C for the unmarked high state. The downshifted outer
state reads the same entire palette as that high state by
converse, so there are only two different palettes to count.
Within the integer interval [H−K+1,H], the high palette can
only use symmetric spacer masks: a nonsymmetric mask has
F−R nonzero, and the estimate in VARIABLE_RADIX.md bounds
its displacement below in magnitude by m−K>K. In the
upshifted palette, nonsymmetric masks have |F−R|≥2.
If positive, the displacement from H−2 is at least
2m−3K+4>K+4; if negative, its magnitude is at least 2m.
So here too only symmetric masks are eligible. Each palette
has at most 2^ceil(t/2) such masks. Together they therefore
hit at most 2·2^ceil(t/2) counts in that interval of length
K, not its K/2 odd targets for t≥6.

This does not rule out a growing interval from other cuts,
other outer factors, or multiple apices. It locates the
next bottleneck: find a new source of independent low digits,
not a larger nominal radix for the same masks.
