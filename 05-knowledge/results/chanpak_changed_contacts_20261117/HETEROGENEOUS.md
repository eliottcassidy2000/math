# Heterogeneous two-vertex exact restoration

The homogeneous hypothesis in the first proof is essential to its uniqueness statement. Here are two other exact induced restorations H' = 7(H−8), with nonidentical contact sets on the two newest vertices. The second is a prime-capable ray.

In either case delete the two final labels of the stated child, obtaining the source. Then replace the designated old slot I in both by TT_k, keeping all old arcs. Since each slot has one outgoing neighbor outside itself, it can be visited at most twice in a Hamiltonian path. Let h and w count the valid directed quotient skeletons with one and two occurrences of I (no two occurrences adjacent), respectively. Enumerate by ordering all outside vertices and inserting I into the requisite number of their gaps, then checking every arc. The count after inflation is h+w(2^k−2), since the two-run case gives 2^k−2 ways to distribute the ordered slot material between nonempty runs. Matching w' = 7w and h' = 7h−56 therefore proves the desired induced restoration for every k≥1.

## Fixture 1: H'/7 = 15·2^k−9

Child rows (252,445,504,32,136,208,154,8,249) from source (124,61,120,32,8,80,26). Slot I=3, outgoing to vertex 5 only. Source (h,w)=(29,15), child (147,105). Both quotients and all their inflations are strong. This quotient ray has permanent factor 3.

## Fixture 2: H'/7 = 50·2^k−17

Child rows (350,300,376,112,98,449,2,351,88) from source (94,44,120,112,98,65,2). Slot I=6, outgoing to vertex 1 only. Source (h,w)=(91,50), child (581,350). Both quotient tournaments are strong. The new vertices 7 and 8 have old outgoing masks 95 and 88, respectively, rather than both beating only the slot. Therefore this fixture is outside the homogeneous D→A→B classification. At k=1,3, the literal quotient is 83,383 respectively; no prime infinitude is asserted. For every odd k it is coprime to 3, and it is always 7-free (the powers of 2 modulo 7 avoid 3).

The one-/two-cover counting also gives H'=7H−56H(B) when TT_k is replaced by an arbitrary common tournament B. This is a structural identity about independent old-block substitution, not about scalar division of arbitrary sources. As in the original homogeneous witness, neutrality H(B)=1 is what makes the intercept absolute.

Replay with replay_heterogeneous.py for literal quotient skeletons and subset-DP checks on inflated tournaments.

## One physical read, two different repairs

Write h=H(B) and P=P₂(B) with empty slots allowed. The homogeneous fixture and the two above give the following quotients H'/7 for the **same choice of arbitrary B**:

| Geometry | H'/7 |
|---|---:|
| Homogeneous | 5P−3h |
| First heterogeneous | 15P−9h |
| Second heterogeneous | 50P−17h |

Thus the first heterogeneous geometry realizes literal **strong** tripling of the homogeneous output for every B, without ordinal-appending an isolated directed triangle. The second is 10 times the homogeneous quotient plus 13h. Both identities are comparisons of physical constructions, not claims that knowing the quotient lets one reconstruct B.

In the prime-capable heterogeneous family, a prescribed quotient p asks for a physical pair satisfying `P=(p+17h)/50`. It is strictly a target equation, not just a congruence. For a module of order n≥2, P is divisible by 4 (pair each subset with its complement, using Rédei oddness), so `h≡−153p (mod 200)`. Also P≥(n+1)h by splitting Hamiltonian paths and counting the two empty-slot cases, so h≤p/(50n+33). These are selectors and pruning conditions; existence of such a pair at every prime is not proved.
