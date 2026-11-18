# Read the marked register without squaring it: an exact chiral target

Use X -> TT_s -> b with just one right singleton b, and reverse
r -> b for the chosen left port r. This gives a one-sided version
of the marked mirror; the right slot factor is 1, not q_r^start.
When s>=1 its exact count is

    h + 2^s q_r^start,

and at s=0 it is h+q_r^start-e_r. The proof is the same two-segment
decomposition: if b->r is used, b ends the first segment and r
starts the second. If it is not used, the baseline loss occurs
only at s=0. With strong X the output is strong: r reaches the
other left vertices, all of them reach b directly or via the
center, and b returns to r.

For the first order-six fixture A from README.md, the ports r=1
and r=2 have q_r^start=18 and 17. At s=1 the two order-eight
outputs have exactly 53 and 51 Hamilton paths. They differ only
on the two arcs from the last vertex b to those two ports.
In particular the marked blue probes compute an exact input
needed for an output that could never itself be blue: 53 is an
all-order self-converse hole in the preceding packet's complete
range below 75. This is a concrete use of the blue construction
as a reader, not as a demand that the output preserve symmetry.
No claim is made that the two outputs are native metagraph
neighbors or that an arbitrary target admits the required port.

For independent reproduction, outgoing row masks for the two
outputs are respectively (250,68,193,198,206,222,128,2) and
(250,196,65,198,206,222,128,4). The direct subset DP in verify.py
returns 53 and 51 on these rows.
