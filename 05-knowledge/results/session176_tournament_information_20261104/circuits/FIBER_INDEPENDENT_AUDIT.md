# Independent audit of exponential response fibers

**PASS.** The reviewer checked the updated exponential-fiber proof and the all-order positivity argument in PROOF.md. No correction is required.

For positivity, a Hamilton path on the other vertices admits the free-orientation pair immediately before its first vertex not dominating both members, or at the end if no such vertex exists. The entering and exiting arcs exist by the stated choice. Deleting only the free mutual edge gives a valid endpoint cover. Tournament parity then upgrades positive total cover count to a positive integer capacity.

For the fiber, ordered triangles cannot be traversed backwards between blocks. If both endpoints of a two-path cover lie in the same block, earlier blocks lie entirely in the first path and later ones entirely in the second, giving internal capacity 3^(k-1). For forward endpoints in different blocks, intervening blocks make the cover impossible; adjacent blocks give 3^(k-2). Reverse endpoints admit three choices at each boundary triangle, twelve at each intervening triangle, and three at each exterior triangle. Hence K(v,u)=3^k*4^t and the displayed pair capacities follow. Internal weighted outdegree has one equal-capacity edge at every vertex, while all cross-block orientations are unchanged. This proves equal entire one-insertion responses.

The distinction between labels and isomorphism is essential and is now explicit: all 2^k bases are isomorphic as unlabelled tournaments. The result is exponential loss of labelled orientation information, not an exponential number of nonisomorphic response twins.

## Independent computation and a continuation witness

`fiber_independent_audit.py` imports no producer code. It constructs K by literal vertex permutations and every split, permitting a failed edge only at that split. It checks all six triangle-orientation choices at k=1,2. These are finite controls on the proof, not its justification.

It also checks all128 sequential two-vertex extensions of the opposite labelled triangles, using literal Hamilton-permutation enumeration. This supplies an explicit continuation counterexample to treating the one-insertion response quotient as an autonomous insertion state:

- Start with T on labels0,1,2 and arcs0→1→2→0, and its reverse T'. Their entire one-insertion responses agree.
- Insert x with outgoing bits to0,1,2 equal to (0,0,1).
- Insert y with outgoing bits to0,1,2,x equal to (0,1,0,0).
- The final Hamilton counts are respectively13 and11.

Thus identical complete one-step response does not suffice for two-step insertion prediction. This is an exact finite witness with all vertices and arcs specified, not an assertion about arbitrary continuation depth. All normal/optimized audit receipts agree.
