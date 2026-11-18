# Four rigid isomorphism classes at infinitely many counts

The [rigid-fibre construction](RIGID_CHIRAL_COUNT_FIBRES.md)
gives two rigid nonisomorphic, nonconverse tournaments L,R at
each level j≥1, of equal order 8j+14 and equal count N. It
proves L non-self-converse. In fact R is non-self-converse
too, so L,L*,R,R* are **four distinct ordinary isomorphism
classes** at that order and count. This does not assert an
edge in the merged metagraph.

In R the unique indegree-one vertex is original vertex 4,
and the unique outdegree-one vertex is the final added root.
Deleting original vertex 4 from the five-vertex R source
leaves rooted profile (1,0,8,1). Write its even appended word
as W. By the retained-root transfer, deleting that vertex
from the full R has count A=1+8W−K(W).

Let 2^s be W's highest bit and P=W−2^s. Deleting the final
root gives the literal preceding tournament, with count
d=9+25P−5K(P)−(P mod2). A converse isomorphism would
have to match the two unique degree-one vertices and preserve
their deletion counts, requiring A=d.

Here W=u+16·16^j V, where
u=(7·16^j−2)/5 and V=(299·16^j+36−10j)/10.
For every j≥1, 28·16^j≤V<32·16^j, by direct inequalities
(19·16^j+36−10j≥0). Consequently the top three bits
of V, and of W, are 111, and P≥3·2^s/4.
The elementary checksum bound 0≤K(z)≤z then gives
d≥8+20P and A≤1+8(P+2^s), so
d−A≥7+12P−8·2^s≥7+2^s>0. This contradiction proves
R is chiral at every level. The existing proof already
distinguishes L from R and R*, and L from L*; conversing
these comparisons completes the four-class claim.

Replay: `python3 four_replay.py`, with `fixtures.py` and
`unit_gain_replay.py` from this packet. The replay independently
checks the four-vertex deletion profile and exact inequalities
through j=30. The displayed argument covers all j.
