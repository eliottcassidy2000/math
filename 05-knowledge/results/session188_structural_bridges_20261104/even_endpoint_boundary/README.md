# An even prime quotient with balanced endpoint parities

This is an exact obstruction to a proposed endpoint-only extension of the odd-pair quotient theorem. It is not a full-response collision and is not asserted to have minimum possible order.

The ten-name tournament with outgoing bit rows

    (598,508,344,865,40,69,656,317,113,438)

is prime: its only nontrivial module is the whole vertex set. Let P_uv be the number modulo two of all simple directed paths from u to v, including the singleton path when u=v. Its row sums and column sums are both the bit vector 296, supported on physical names {3,5,8}. The same statement holds for Hamilton paths only. Thus the all-simple-path start-minus-end parity vector vanishes at even order.

Nevertheless P is not symmetric: P_01=1 and P_10=0. This also refutes the proposed implication that (P+P^t)1=0 forces P=P^t for tournament path matrices. The familiar total-path check is consistent: the sum of all entries of P is one modulo two because every induced tournament has odd Hamilton count, while the sum of its ten diagonal entries is zero.

`verify.py` computes both full parity matrices separately for each starting name by literal subset path extension. It checks all pairs are oriented, all nontrivial subsets for modularity, the endpoint vectors, and the asymmetric entry. `normal.json` and `optimized.json` agree and retain both complete matrices.

No equality of complete two-port functions, integer endpoint counts or modulo-eight mixed traces is claimed. The proved odd-order pair-child theorem is unaffected. The failed extrapolation concerns a proposed route for even quotient order only.
