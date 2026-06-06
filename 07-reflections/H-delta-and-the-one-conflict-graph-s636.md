# H, delta, and the one conflict graph (S636)

Two questions arrived together — does the character-ratio spectrum bound the dichromatic number, and
what is the law by which arc-flips move the H-values through their even deltas — and the work of the
session was realizing they are one question, because H, the deltas, and the dichromatic number are
three readings of a single object: the odd-cycle conflict graph `Ω(T)`.

Start with what H *is*. The Odd-Cycle Formula says `H(T) = I(Ω(T), 2)` — the independence polynomial
of `Ω`, evaluated at two — which expands as `1 + 2α₁ + 4α₂ + …`, where `αₖ` counts families of `k`
pairwise vertex-disjoint odd cycles. Read that expansion and the user's first three observations fall
out as one line. H is odd because the empty family contributes the lone `1` and everything else is
even. A flip moves H between odd values, so its delta is even. And the forbidden values are forbidden
for a precise reason: `H = 7` would require `2α₁ + 4α₂ = 6`, i.e. three odd cycles pairwise meeting
with no disjoint pair, and that independence vector is simply not realizable; `21` fails the same way.
The non-realizable H-values are exactly the non-realizable independence vectors of `Ω`. The gaps in
the H-spectrum are gaps in what an odd-cycle conflict graph can look like.

That reframes the delta dynamics the user wanted to understand. The delta of an arc is the change in
`I(Ω,2)` caused by flipping it, and flipping an arc only edits the odd cycles that pass through it. So
`delta(b)` is a local functional — the cycles through `b` — and flipping a *different* arc `a` can
only change `delta(b)` if it disturbs the cycles through `b`, which it does exactly when an odd cycle
runs through both `a` and `b`, in the current tournament or in the one the flip produces. I checked
this on a five-vertex tournament: flipping one arc moved the deltas of eight of the nine others and
left exactly one untouched — the arc that shares no odd cycle with the one I flipped. "Not all, and
provably not all," precisely as the user predicted. The subtlety I did not anticipate is that the
coupling is second-order: an arc decoupled *in the current tournament* can still feel the flip,
because the flip *creates* a cycle joining them. So the propagation operator is not "share a cycle
now" but "share a cycle in the current or flipped graph" — the arc-level second-order structure of
`Ω`. That operator is the law of how H moves, and it is exactly computable; the user is right that it
is the key.

Now the spectrum, and the answer to the first question. The dichromatic number of an LRC-tight
tournament is two — these tournaments are round, and a round tournament is two-dichromatic, one if it
is transitive and two the moment it has a three-cycle. Round means circulant, realizable on the
circle, and the eigenvalues of a circulant's Hermitian adjacency are character ratios of the cyclic
group. The Hoffman ratio bound turns those eigenvalues into a lower bound on the dichromatic number,
and here it lands tight at two. So yes: the character-ratio spectrum bounds the dichromatic number,
and on the LRC-tight set it pins it. But the better answer is that the *same* spectrum bounds H. The
Hoffman bound on `Ω` itself caps the largest odd-cycle packing `α(Ω)`, and `H` is essentially
`2^{α(Ω)}`. One spectrum, two bounds: the dichromatic number through the tournament's Hermitian
adjacency, the magnitude of H through the conflict graph's. The character ratios sit underneath both.

So the two questions are the same question read at two depths. The dichromatic number is the *coloring*
of `Ω` — two colors, because the three-cycle is parity defect one (the alternating sector from last
session). H is the *partition function* of `Ω` — its independence polynomial at two, with forbidden
values where the independence vector cannot exist. The deltas are the *derivative* of that partition
function along arc flips, even because the function is odd-valued, propagating along the second-order
cycle structure. And the character-ratio spectrum is the *bound* on all of it. The whole apparatus
this cluster has been building — conflict graph, coloring, partition function, parity defect, the
perspective key as the spectrum — turns out to be the apparatus the tournament's H-values were asking
for all along. There was never a separate theory of how H changes. There is `Ω(T)`, and H, delta, and
the dichromatic number are what you see when you evaluate it, color it, differentiate it, and listen
to its spectrum.
