# The last one-apex internal wiring: independent verification

22 October 2026. This concerns only the induced one-apex stratum, not bridges in the full metagraph. No novelty claim is made.

## Structural proof at order eight

Use positions 1,...,8 along the pinned directed Hamiltonian path. Put T'=kappa(T). Every three-consecutive-position window is a cyclic triangle in exactly one of T,T'. Consequently the union of any deletion-transitive apex of T and any such apex of T' hits every three-window. A pair hitting every window must be {3,6}: its first element is at most 3, its second at least 6, and their gap at most 3. No singleton suffices. Thus both endpoints are nontransitive, their apices are unique, and these occur at 3 and 6 in opposite roles.

First fix T's apex as 3 and T's companion's apex as 6. On B={1,2,4,5,7,8}, both tournaments are transitive and their agreement arcs are precisely 1->2, 4->5, 7->8. Each agreeing pair must be adjacent in either transitive order: if u->x->v intervened between an agreeing u->v, flipping the two cross-pair arcs would make a cyclic triangle in the companion. Hence the three pairs form consecutive blocks, with internally forward orders, and the block orders reverse between the endpoints.

The path constraints 5->6->7 force block (45) before (78) in T. The constraints 2->3->4 in T' force block (12) before (45) in T', equivalently block (45) before (12) in T. Thus T's B-order is either

- (45),(12),(78), with three cuts for inserting 6 in T-3 and one cut for inserting 3 in T'-6;
- (45),(78),(12), with one cut for 6 and three cuts for 3.

Every such choice determines all arcs except the free arc between 3 and 6. These constructions satisfy both transitivity conditions and every path arc. They give (3+3)*2=12 distinct tilings; exchanging the endpoint roles gives another 12. This proves completeness of the 24-element internal tiling set.

Neither endpoint can be self-converse. Any anti-isomorphism must fix its unique apex and reverse the transitive order on the other seven vertices, so also fixes that order's middle vertex. But an anti-isomorphism cannot fix two distinct vertices, since it would reverse their mutual arc. All internal edges are therefore black.

## Actual edge-orbits, not merely endpoint pairs

Encode a tiling by the non-path arcs (i,j), j>=i+2, listed lexicographically, bit 1 meaning i->j and the first pair being the least significant bit. Then kappa XORs with 2^21-1; reflection sends the bit at (i,j) to that at (9-j,9-i). The fundamental edges are the orbits of these two commuting involutions. All six order-eight internal orbits have size four:

- {2247,32767,2064384,2094904}
- {2511,8191,2088960,2094640}
- {2503,24575,2072576,2094648}
- {4095,6607,2090544,2093056}
- {6343,28671,2068480,2090808}
- {6599,20479,2076672,2090552}

The endpoint word is the apex's out-neighborhood indicator along the transitive deletion order, minimized over deletion apices and converse (reverse-complement). The complete nontrivial induced components are:

    0001010 --- 1000010 --- 0010110
    1000110 --- 0010010 --- 1001010
    0011010 === 1010010

Here === denotes two parallel black edges, specifically the distinct orbits with minima 6343 and 6599. There are 47 additional isolated nodes, giving 55 nodes in all. There are no loops. The four path edges are bridges only in this induced graph; neither double edge is a bridge there. Nothing here classifies full-metagraph bridges.

## Complete small-order catalogue

n | ordinary classes | merged nodes | selected tilings | blue incidences | internal tilings | internal blue edges | internal black edges | loops | isolated nodes
3 | 2 | 2 | 2 | 2 | 2 | 1 | 0 | 0 | 0
4 | 4 | 3 | 8 | 4 | 8 | 2 | 1 | 1 | 0
5 | 9 | 7 | 43 | 9 | 30 | 3 | 6 | 1 | 1
6 | 21 | 12 | 209 | 9 | 52 | 0 | 13 | 2 | 3
7 | 48 | 29 | 926 | 34 | 48 | 0 | 12 | 0 | 15
8 | 106 | 55 | 3898 | 18 | 24 | 0 | 6 | 0 | 47

The JSON contains every edge orbit, its literal tilings and merged endpoints for n=3,...,8, together with all selected nodes. Loops are counted as edges, once in the loops column. At n=8 there are also exactly 18 blue boundary edges and 1928 black boundary edges, obtained by subtracting the 12 internal black incidences from the total (3898-18)/2=1940 black incidences.

## Reproduction and completeness of the independent computation

Run `python catalogue.py`. It enumerates all 2^(n-1) apex words on a transitive base, enumerates their actual directed Hamiltonian paths, relabels each path to the fixed path, and deduplicates the resulting bitmasks. Every one-apex tournament has at least one such representation, and every fixed-path presentation arises from one of these Hamiltonian paths; thus this enumerates exactly the selected tiling set, without enumerating all tournaments. Ordinary canonical labels minimize over all deletion-transitive apices. Merged labels additionally minimize over reverse-complement. An assertion checks consistency when multiple rooted representations yield the same literal tiling. Edge orbits are taken under the actual kappa and reflection maps, preserving multiplicities. The enumeration does not rely on the structural argument or the root's twelve-configuration implementation.
