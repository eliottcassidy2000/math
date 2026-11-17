# A Bézout repair giving consecutive odd counts with identical marked interface

This combines the exact independent ordinal switches with a newly constructed switch at path count 97. A motif suggested by Catalan's adjacent powers is useful here as an exact determinant condition rather than an assertion about congruence search: if two blocks have ordinary counts h1,h2 and signed independent changes δ1,δ2 in their K channel, the joint correction is h2δ1+h1δ2. Matching it to a small integer is a condition on two actual physical blocks.

An old 8-vertex block A has outgoing rows (158,164,56,242,98,1,167,52). Tag choices 144 and 6 have common interface (H,P₂,W,Z,E,S)=(389,5600,1168,1376,102,88) and respective K=217,225. The newly matched module of count 97 changes its arcs as well as its tags: choices
    (254,240,162,54,4,144,60,88), tag 61;
    (142,136,18,20,131,31,191,44), tag 79.
These have identical (H,P₂,W,Z,E,S)=(97,2560,1867,1715,69,97), with K=−1093,−1095. Here W/Z are weighted starts/ends; E/S are unweighted endpoint counts, and K=H−U, with U the two-contact ordered cover read.

Ordinally concatenate A and the new module. The transfer fixes J and W, and scales each signed K move by the other ordinary count: common mixed terms are unchanged exactly because the weighted starts and ends and both endpoint counts agree. Thus the four variants have
    J=37733, W=3716343,
    K=−2454875 +776 i−778 j, i,j in {0,1}.
The two interior switched states (i,j)=(0,0),(1,1) differ by two, due to 389=4*97+1 (so 8*97−2*389=−2). The independent switches may be made in either order.

Attach the one-exit q and B as usual: q beats precisely the tagged vertices, B beats q and loses to every other vertex. For arbitrary nonempty B, the two corresponding output counts are W P₂(B)+K H(B), so the common move changes the absolute count by precisely −2 H(B). In particular B=TT_k produces adjacent odd Hamiltonian counts at every k. Direct reachability on all four quotient graphs proves strength; it survives all substitutions at B. At k=11, the consecutive odd counts are the actual twin primes 7608615587 and 7608615589. Trial division certifies both. This neither asserts infinitely many prime twins nor that every target belongs to this bank.

Replay: python replay_bezout.py certifies the two full switch interfaces, the four ordinal states and their strengths, and both numerical primes.
