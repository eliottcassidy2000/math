# Changing the reader after a tournament carry

The radix-four reset reopens **physically** on the infinite
m=11·2^v branch. In [the main construction](BANK.md) the same
boundary rewrite gives 4m²+1→4m²+3 at every v≥2. It reverses three
old arcs, adds two vertices, and retains every old spacer arc.
[Five more specific reversals](LOCAL_MOVES.md) then give 4m²+5
without adding a vertex. A [second square rail](SECOND_RAIL.md),
m=19·2^v, also reopens to offsets +3 and +5.

[The same construction](BANK.md) gives a verified bank of
**38 consecutive odds, 4m²−5 through 4m²+69**, at every height
of this branch. After one fixed boundary setup, the entire bank
uses the same growing tail and a 13-vertex mutable boundary.
Each +2 can be achieved by bounded rewiring with no added vertices.
[Two anchors rule out](NEED_A_NEW_READER.md) a universal repair
based on moving the cut within the transitive one-apex family.

[The Collatz side](COLLATZ.md) explains why the two low digits
want different signed readers, and supplies an exact decreasing
odd-count selector avoiding 7 and 21. The exceptional 7 landings
have explicit alternate smaller parents on real tournament rays.
What would complete a proof is the physical inverse of the signed
branches and of the bypass, not improvement of the convergence of
ordinary Collatz. The 38-count bank does not yet give those
inverses at an arbitrary prescribed target.

## Reproduce

Run verify.py, check_local_moves.py and check_collatz.py
with Python 3. Compile no_one_apex.cpp for the reduced-spine checks.
root_atlas.cpp takes the complete order-8 CSV from the
[fixed-count census](../actual_hamilton_multiplicities_20261021/session24-multiplicities/)
to replay reader discovery; verifying the bank does not require it.
