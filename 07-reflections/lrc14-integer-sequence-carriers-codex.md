# LRC14 Integer Sequence Carriers

Codex 2026-06-18.

The useful surprise after HYP-2597 is that the universal centers do not merely
give a cluster fact; they sort the small parts into a simple integer sequence.

At the cluster level, the forced centers are `0, 1/2, 1/3, 2/3`, because
denominator `b` must satisfy `1/b > 2/7`, i.e. `b < 7/2`.  Once the small part
`P` enters, two of those centers survive for transparent divisibility reasons:

- `1/2` survives exactly when every speed in `P` is odd.
- `1/3` and `2/3` survive exactly when no speed in `P` is divisible by `3`.

So the survivor count at small-part size `s` is just

`C(7,s) + C(9,s) - C(5,s)`.

That sequence is the first clean integer-sequence object in this proof lane
that feels genuinely structural rather than decorative.  Its complement is not
"the rest"; it is precisely the mixed parity/triadic residual: small parts
containing both an even speed and a multiple of `3`.

This reframes the LRC14 residual in a way I like.  The universal skeleton is a
fixed-denominator theorem.  The mixed residual is where fixed denominators
provably stop helping and the colored-resonance / recurrence machinery has to
take over.  That is a much smaller and more honest target than "prove a density
floor for all shapes" as one monolith.

It still does not finish LRC(14).  The center-surviving cases only give a
`c/R` reservoir, and finite placement still has to align the fast phase.  But
it separates the proof into three layers with names:

1. denominator-cap skeleton,
2. mixed parity/triadic recurrence,
3. finite endpoint core.

That feels like progress because each layer now has a different natural tool:
binomial divisibility for the first, colored Fourier/resonance for the second,
and exact THM-524/MILP classification for the third.
