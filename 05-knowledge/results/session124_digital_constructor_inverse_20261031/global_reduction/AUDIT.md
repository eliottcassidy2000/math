# Independent audit: reduction to the non-3-divisible lane

**PASS.** The bridge is an exact conditional equivalence, not a proof of the remaining universal existence statement.

Assume the established exclusions7 and21. The assertions

1. every positive odd N other than7,21 is attained by a tournament;
2. every positive odd m with3 not dividing m and m≠7 is attained by a tournament

are equivalent once an actual H=63 tournament is supplied.

The forward implication is immediate. For the reverse implication, write the unique factorization N=3^k m with3 not dividing m. If m≠7, an attaining tournament for m, followed by k ordinal C3 factors, has count N. Here k=0 requires no factor and m=1 is permitted. If m=7 and k≥2, use the explicit H=63 tournament followed by k−2 ordinal C3 factors. The remaining cases k=0,1 are exactly7,21 and are excluded, rather than assumed attainable.

Ordinal-sum multiplicativity requires no strongness hypothesis: if every vertex of A beats every vertex of B, a Hamiltonian path cannot return from B to A, so it is uniquely a path of A followed by a path of B. This gives H(A⊕B)=H(A)H(B). A directed C3 has exactly three Hamiltonian paths. In particular the proof does **not** factor an H=63 witness through a hypothetical H=7 witness.

Universal positive answers from the specified restricted non-3 inverse would imply assertion2, because every accepted word reconstructs an actual tournament. They are therefore **sufficient** for the unrestricted conjecture. They are not shown necessary: an unrestricted witness may lie outside that seed/constructor language. No assertion is made that these inverse answers are all positive, or that resulting tournaments are strong.

## Independent H=63 witness verification

The supplied seed[1,0,0] and singleton word[1,0,0,1,1,0,0] give an eight-vertex transitive spine and one external vertex. With the spine first and the external vertex last, the independently reconstructed outgoing masks are

`[510,252,504,496,224,192,384,256,50]`.

`independent_audit.py` imports no constructor code. It checks all9! vertex permutations and separately runs adjacency/subset DP; both counts are63. Normal and optimized outputs match. The receipts additionally check the factorization branches for5,000 small odd integers as controls, not as the proof of the universal reduction.

Only files inside this new addendum directory were written. No sealed digital-package file, million-formatting patch, existing code or existing receipt was changed.
