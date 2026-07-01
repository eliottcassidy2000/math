# One law across all layers, graded by Lucas

*klein-2026-07-01-S76. A reflection on HYP-3809 — mining two old bucket-parity results and finding they
are one law, and answering the owner's "the metagraph is the constraint."*

The owner asked me to think of the merged metagraph as natural numbers dropped into buckets under parity
rules, to mine the repo for the half-tiling and adjacent partial ideas, and to decide whether the metagraph
is *equivalent to the constraint that builds it*. Mining turned up HYP-1772 — the merged-tiling bucket
parity, already proved: each isomorphism class holds an odd number of tilings, because Rédei makes the
Hamiltonian-path count odd and a tournament's automorphism group has odd order (an involution would have to
reverse an edge, which a tournament forbids). So self-complementary buckets are odd, non-self-complementary
buckets are twice-odd, and the power of two `2^m` splits into exactly those two kinds. That was the whole
parity of the buckets, settled a year of sessions ago. My own last session had rediscovered a cousin of it
at the opposite end — the blue/black lines, which flip *all* the tiles — without realizing the two were the
same phenomenon at different scales.

They are. HYP-1772 lives at Hamming distance one (flip a single tile); my blue/black lines live at distance
`m` (flip every tile). Between them is a whole ladder of layers, distance `d` for each `d`, and every rung
obeys the same accounting: a node's `M_u` tilings each have `C(m,d)` neighbors at distance `d`, so the
cross-incidence at that layer is congruent to `C(m,d)·M_u` modulo two. The node's mass parity is its
self-complementary bit; the binomial's parity is, by Lucas' theorem, whether `d` is a binary submask of
`m`. So the cross-incidence at layer `d` is odd exactly when the node is self-complementary *and* `d` sits
inside `m` in binary. Distance one recovers HYP-1772 (odd iff `m` is odd); distance `m` recovers my
handshake (the top binomial is one, so every self-complementary node has odd degree, so their count is
even). The two partial results were the bottom and top rungs of a single Lucas-graded ladder, and the
parity-active rungs form a little Boolean cube of size two-to-the-digit-sum-of-`m`. What looked like two
scattered observations was one law I could have written down had I seen the layers as a family instead of a
pair.

That reframes the owner's central question cleanly. Is the metagraph equivalent to its constraint? Modulo
two, yes, completely: the odd/odd/even category rule plus the all-layers Lucas law pin every parity the
metagraph has — every cross-incidence, the evenness of the self-complementary count, the two-adic valuation
of every bucket (zero or one, never more). If you only care about parities, the constraint *is* the
metagraph; there is nothing left to compute. Over the integers, no: the actual bucket sizes are `H/|Aut|`,
the Hamiltonian-path count divided by the automorphism order, and parity cannot see that arithmetic —
plenty of parity-legal size-vectors are never realized by any tournament. So the honest answer is a
factorization of the object: the metagraph splits into a parity skeleton, which the constraint determines
entirely, and a metric flesh, which is exactly the Rédei/`|Aut|` data. The owner's intuition is right about
the skeleton and points precisely at what the skeleton omits.

And the fiber question — which tilings land in the same bucket — has the same shape of answer. A bucket is
one automorphism-orbit of Hamiltonian paths (the base path can be laid down along any Ham path; two
layings give the same bucket iff a relabeling relates them), and the action is free (LEM-003), so the
bucket size is `H/|Aut|` with no remainder. The symmetry restriction the owner suspected is exactly this:
the tilings of a class are permuted freely by its automorphism group, and the grid-symmetric ones — the
blue, the half-tiling addresses — are the fixed points of the reflection that is complement-composed-with-
reversal. The half tiling model is not a separate gadget; it is the blue layer of the same ladder, the
distance-`m` rung seen through its own symmetry.

There is a thread worth pulling: `n = 6` keeps being the threshold. The flip-rank's efficient shape dies
there, the non-self-complementary buckets first grow self-loops there, and `m = C(5,2) = 10 = 1010` in
binary is the first `m` in the sequence `1, 3, 6, 10` whose digit-sum climbs in a way that spreads the
parity-active layers apart. I do not yet have the common cause, but three different structures changing at
the same size, and a binary-digit story sitting underneath the parity law, is not likely a coincidence. The
lesson of the session is the one the owner was steering toward: when a structure is defined by a process,
find the conserved quantity, and then look for it at *every* scale of the process — the law you proved at
one distance is usually the whole ladder, graded by something as old as Lucas.
