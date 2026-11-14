# A cycle fan recovers the ordered SCC sizes of its source

**All orders.** Let (R) be a nonempty tournament on (m) vertices,
with ordered strongly connected component sizes (s_1, …, s_k).
Form the cyclic wrap (T=C_3[b,a,R→x]): (b→a),
(a→R,x), (R,x→b), and (R→x).
Every Hamiltonian cycle of (T) is
(b→a→r_1→⋯→r_m→x→b) for exactly one
Hamiltonian path of (R). Reverse the cycle, then restore one cut arc
to its original orientation. Call the result the cut mate.
This is the closed path-line operation of the
[cycle fan](../chanpak_h9_cycle_metagraph_20261114/CYCLE_NECKLACES.md).

**The boundary is finer than its size.** In each fan, the nonstrong mates
have exactly the following ordered SCC-size compositions, with one
incidence for each row or eligible index:

| Original arc restored | Ordered SCC sizes |
| --- | --- |
| (x→b) | (m+2,1) |
| (a→r_1) | (1,m+2) |
| (r_m→x) | (m+1,1,1) for (m≥2), and (1,1,1,1) for (m=1) |
| (r_i→r_{i+1}), (i=s_1+⋯+s_j<m), (i=1) | (1,1,m+1) |
| (r_i→r_{i+1}), same condition, (i≥2) | (i+1,m-i+2) |

Every other cut mate is strong. From one fan with ordinary orientations,
discard the three portal rows and read (i) from the surviving prefixes;
the consecutive differences of (0,i_1,…,i_{k-1},m) recover the
source composition. This does not require keeping the labels of (R)
or knowing which arc was cut. It does require the SCC sizes in order,
not an independently converse-merged signature for each mate.

## Proof

In the full reversal, (b→x) is the only exit of (b),
(x→r_m) is the only exit of (x), and (r_1→a) is the
only entry of (a). Restoring the corresponding arcs leaves (b) a
sink, (x) a sink followed by (b) as the next sink upon deletion,
or (a) a source, respectively. After removing the indicated singleton
vertices, the remaining vertices lie on a directed cycle, except when
the remaining block is ({r_1,a}) at (m=1); there (r_1→a).
The reversed (a→b) can always be bypassed along
(a→x→r_m→b).

For an internal arc of the reversed source path, let (P=(r_1,…,r_i)),
(S=(r_{i+1},…,r_m)). The full reversal adds only one arc
from (S) to (P). If there was an original arc (s→p)
across that cut, the reversed arc can be bypassed: from (r_{i+1})
use (b,x,r_m) and the reversed suffix to reach (s); after crossing
to (p), use the reversed prefix to reach (r_1), then (a) and
finally (r_i) if needed. Thus removal, hence also restoration,
leaves a strong mate. If there was no original crossing, restoring the
arc gives the ordinal split (P+a) → (S+b+x). The right block has
cycle (b→x→r_m→⋯→r_{i+1}→b);
the left has cycle (a→r_i→⋯→r_1→a) for (i≥2),
and otherwise (r_1→a). A Hamilton path follows the SCCs in ordinal
order, so the no-crossing case occurs precisely at the component boundaries.

## Counts and scope

Write (h=H(R)), Q=P_R(2) for the ordered two-cover count allowing
empty slots, and (A=|Aut R|). A cyclic triangle
with two singleton blocks and a third block (B) has (3H(B))
single-visit paths and (P_B(2)-2H(B)) return paths, so the wrap
has (H(T)=h+2Q): the singleton (x) doubles the two-cover count.
The two outdegree-one vertices (b,x) are distinguishable (only
(x) points to the other one); this fixes (a) as well, so
(Aut T=Aut R), acting freely on cycles
because it fixes (a). Per ordinary node the three incidence counts
are therefore

- closed: (m+3)h/A;
- closed with nonstrong mate: (k+2)h/A;
- open: [2Q-(m+2)h]/A.

The cyclic-wrap fan can read order and SCC structure that the scalar
(h+2Q) alone is not designed to record. For example, ordinally
permuting the SCCs of (R) preserves (h) and (Q) by their
multiplicativity, yet changing the ordered size composition changes
this fan signature. In particular (R=C_3→TT_1) and
(R'=TT_1→C_3) give two wraps with (H=51), (h=3), and (Q=24),
but the extra boundary signatures are (4,3) and (1,1,5).
This is not a reconstruction of all orientations or a coverage proof.

The general boundary-size theorem was suggested by the background
research pass; this note adds the ordered decomposition and inverse.
Run `python cut_replay.py` for the 207-fan direct check through source
order four; no finite check is required by the all-order proof.
