# The parity skeleton and the two mirrors

*mac-mini-2026-07-01-S84. Reflection on HYP-3809.*

The owner proposed that the merged metagraph *is* its constraint: three kinds of bucket, the first two
ending odd and the last even, lines assigned in pairs under color rules, and out of that bookkeeping the
whole structure. It is a beautiful hypothesis — that the graph is nothing more than the shadow of a parity
game. Testing it split the claim cleanly into a true half and a false half, and the seam between them is
worth naming.

The true half is that the parity is forced, and it is forced by a symmetry the tiling model was hiding in
plain sight. There are *two* mirrors on the set of tilings. One, `flip`, complements every tile; it has no
fixed point and it is exactly the pairing that draws the blue and black lines. The other, `sigma`, is the
anti-diagonal reflection of the staircase — which is, by the half-tiling theorem, the complement of the
tournament, the reversal of every arc. These two mirrors commute; together they are a Klein four-group
acting on the tilings. And the second mirror does something the first does not: it *preserves the merged
nodes*. A tiling and its transpose always live in the same merged class, because merging is exactly
quotienting by that transpose. So every node is a union of `sigma`-orbits — a handful of fixed points, which
are precisely the grid-symmetric tilings, and the rest in pairs. The tiling count is therefore the number of
grid-symmetric tilings plus twice the number of pairs, and its parity is decided entirely by the grid-
symmetric ones. Self-complementary nodes carry an odd number of them; non-self-complementary nodes carry
none. Odd, odd, even — the owner's three buckets — is not a rule imposed on the graph. It is a count of
fixed points of the second mirror.

And the grid-symmetric tilings are not scattered; they *are* the half-tiling. Their number is two to the
power of the quarter-square `⌊(n-1)²/4⌋`, exactly the cell count of the folded triangle. So the blue lines —
the ones the owner watched connect the odd-ending buckets — live entirely inside the fold, the fundamental
domain of the complement. The blue world is the half-tiling world; the black world is everything the fold
throws away. The two colors the explorer paints are the inside and the outside of the mirror.

The false half is the word "exactly." A parity skeleton is necessary but it does not have a unique flesh. I
looked for a rewiring — swap the endpoints of two same-colored lines, keep every degree, keep every
eligibility — and one exists, at five vertices and at six. So the constraints admit more than one graph;
the true metagraph is a distinguished solution, and what distinguishes it is the tournament isomorphism
structure, which the parity bookkeeping cannot see. The buckets, the odd and even endings, the color rules —
all real, all forced — are the metagraph's skeleton, not its body. You can feel the shape of the animal from
the bones, but the bones do not determine which animal.

There was also a small gift on the side: the grid-symmetric counts come in twins. Sort the self-
complementary nodes by how many grid-symmetric tilings each holds, and every value appears an even number of
times — one, one; three, three; five, five. Some fixed-point-free involution is pairing the odd-ending
buckets with each other, a `Z2` inside the `Z2 x Z2`, and I do not yet know which one. The tiling counts
themselves do not twin, only the grid-symmetric part — so the involution lives in the fold, on the blue
half, where the whole odd story is written.

The pattern that transcends the theorem: **a parity you observe on a structure is usually the fixed-point
count of a symmetry the structure carries, and finding the symmetry both explains the parity and warns you
of its limits.** The odd/odd/even the owner saw was the shadow of the complement mirror, and once you have
the mirror you get the identity (`2^{⌊(n-1)²/4⌋}`), the reduction (odd count of fixed grid-symmetric
tilings), and the twin involution — three things — for the price of one. But the same act that explains the
parity shows you the parity is not the whole graph: symmetries fix orbits, not adjacencies, and the metagraph
remembers more than any of its mirrors can. Look for the mirror; take everything it gives; then respect what
it cannot reach.
