# Merged Line Parity Even/Odd

The answer is both yes and no, and that is the useful part.

As an addendum to HYP-2250's boundary-vector audit, if "black/blue" means the
corrected complement-line convention inside the
fixed-base tiling cube, then the black portion is always an even graph in the
strong Eulerian sense, and the active blue portion is an odd-degree graph.  The
S675b audit checks this exactly through `n=7`.

If "black/blue" means the older simple merged-metagraph edge coloring by
SC/non-SC endpoint type, the conjecture is false.  The black subgraph already
has odd-degree witnesses at `n=4`, and the blue active subgraph stops being
all-odd at `n=5`.

That distinction is a good sign, not a bookkeeping annoyance.  It says parity
lives one level up from the simple graph, in the line lift where each tiling
knows its complement partner.  Projecting to `G_n/Z_2` keeps adjacency but
forgets the number of line endpoints and how self-loops contribute two
half-edges.  The missing coordinate is exactly the kind of address HYP-2245
warned us to retain.

The proof mechanism is small.  In the fixed-base-path cube `Q_m`, complement
pairs tilings as `{x,C(x)}`.  Grid-reflection-fixed tilings form the blue
support.  Over an SC class, non-fixed tilings pair under reflection and the
fixed-base tiling count is odd, so the blue endpoint count is odd.  Over the
black side, SC nodes have `odd - odd`, and NS merged pairs have `odd + odd`;
both are even.  Therefore black is Eulerian and blue is the odd-boundary
coset.

This is the clean relation to even and odd numbers.  An even graph is a
boundary-zero object over `F_2`; the black line lift is in that kernel.  The
blue line lift has boundary one on the SC support, so it is an affine coset of
the same cycle space.  Odd is not "bad even"; odd is even plus a named boundary.

The even-graph equinumerosity thread should use this form.  Bare count equality
between fixed-path tilings and labelled even graphs is true but weak.  The
useful statement is that a specific line-lift inside the tournament quotient
already lands in the even-graph cycle-space language.  The next question is
whether the black Eulerian carrier has stable cycle basis, component rank, or
transport invariants that survive further quotienting.

For LRC14, the transfer pattern is tempting: find a line lift where the black
side is an Eulerian certificate and the blue odd boundary marks the exact floor
atoms.  If the visible quotient loses that parity, attach the owner/carry
address that makes the boundary visible again.
