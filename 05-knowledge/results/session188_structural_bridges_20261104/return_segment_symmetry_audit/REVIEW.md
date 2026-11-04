# Independent review: return-segment symmetry

PASS. The frozen all-order proof and numerical controls were reviewed without changing the source packet.

The return decomposition is unique, and the permutation action preserves every physical label, every exterior join and both actual port incidences. Each nonempty i-return needs at least two distinct non-i symbols: one symbol would force a loop or a tournament two-cycle. Therefore, with at most one other repeated child, equal return segments are impossible and the symmetric-group action is free. The factorial sharp construction exhausts every other label in its returns, leaving no extra prefix/suffix choices.

For a power-of-two number R of returns, the multinomial orbit valuation is sum s2(m_j)-1. Odd orbit size forces a single repeated word class. Such a return contains no port because each port has multiplicity one and R is at least two. Contraction to the unique remaining i and expansion back are inverse, retaining arbitrary prefix/suffix joins. The first out-to-in transition gives the claimed heavy-label triangle. The condition is necessary only, and the proof does not extend the Hamilton return action across arbitrary cover cuts. For a non-power-of-two return count, mixed classes can have odd orbit size, as 3!/(2!1!)=3 already shows.

The independent standard-library checker imports no producer code. Its forward used-inventory dynamic program reproduces every complete parity and factorial record, while an explicit word enumeration independently reconstructs every retained return-permutation orbit. It separately enumerates every contracted return word and compares each exact summand, not only the final residue. It also checks 256 additional four-child cases outside the producer parity domain, both sharp constructions, and the necessary-not-sufficient control.

All 13,248 parity cases, 768 factorial cases, 144 explicit word sets and 1,062 explicit words pass. The original producer was replayed in both normal and optimized Python and compared as complete output bytes. Both independent audit modes also agree byte for byte. Source pins and the inherited coarse-modulo-eight proof pin are retained in the receipts. This is independent mathematical and exact-computation review, not a formal proof-assistant certificate or a response-separation theorem.

Reproduce with python audit.py and python -O audit.py from a checkout containing the sibling source and inherited proof. Exact audit allowlist: REVIEW.md, audit.py, normal.json, optimized.json, SHA256SUMS.
