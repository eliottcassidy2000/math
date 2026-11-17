# One arc suffices even for the continuable, locally prescribed prime-pair seed

This strengthens the main [fourteen-step construction](README.md). Its two-contact relay is unnecessary if the edited arc is **not incident to the marked root** and the deletions need only have equal counts rather than identical arcs.

Take rows (8,17,35,54,101,195,143,31), mark vertex 3 and reverse the edge 4↔7. The marked states in the order (H,E,O,g) are

    (315,198,235,37),   (329,198,235,37).

Thus the state that determines the continuation is already equal, and the states are normalized. Apply N,S,S,S,N in chronological order. Independent subset DP gives

    (7343,6595,6632,37),  (7357,6595,6632,37).

The two tournaments are strong, have prime quotients 1049 and 1051, and still differ on exactly one arc. Every subsequent common S/N word preserves the count gap and that one-arc relation. This argument does not contradict the earlier impossibility for a *single root recontact on the same deletion*: the edited arc is in the deletion this time, and the two deletions only have equal g.

For pure N continuations, the reads after t additional steps are 13190·2^t−5847 and 13190·2^t−5833. At t≡0 mod 42 both reads are congruent to their seeds modulo 147, ensuring v₇=1 and quotients coprime to 3. The integers need not be prime.

More generally the first read for a continuation ending N, with address z, is

    H(z) = −5884 + 13190(z+1) + 37 F(z).

The main joint phase/scale proof applies with constants −5884,13190,37 in place of −482,1506,25: 37 is odd and gcd(13190,37)=1, all that the phase proof requires. Apply it at full size 7x with modulus 7M, imposing H≡7c mod 7M for any parity-compatible quotient residue c. Thus the **one-arc** pair may be chosen strict v₇=1, avoids any fixed finite simultaneous prime sieve, and has its lower quotient within O_M(√x) of x while giving the requested exact residue. This does not show the error is zero or that a pair has prime quotients at all scales.

The constraint that both marked *deletions* be identical is what made the two-contact lower bound apply in the principal construction. Equal continuation does not require identical deletion, because marked path-cover counts lose information about the underlying arcs. This distinguishes a root recontact from a fixed-root change inside a count fibre.

Replay: `python single_arc_relay.py` (also with `-O`), using the literal DP and recurrence checker of `replay.py`.
