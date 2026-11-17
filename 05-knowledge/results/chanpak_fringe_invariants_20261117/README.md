# Separate the information a construction needs

This session produced two different invariants for exact tournament constructions:

- [Fringe and carry compiler](FRINGES.md): the two-cover count modulo powers of two depends only on the fringes and a length phase. Independent interior changes program the Hamiltonian count. The odd-modulus word can then be inserted without disturbing either. For the earlier mixed read (p,7q), all compatible joint phases are now constructively reachable modulo 2^a·3·m for gcd(m,6)=1. This is not an absolute-target theorem.
- [A growing one-arc unit repair](UNIT_SEAM.md): paired binary addresses give exact +2 windows and, without a carry, preserve every longer odd-cycle count and all low bits of the two-cover read. A second unbounded switch moves H by 12 while preserving the full P₂ exactly.
- [Complete one-port substitution invariant](ONE_PORT.md): four marked path-cover sequences are necessary and sufficient for substituting one tournament core for another when the exterior sees all but one distinguished vertex uniformly. Explicit finite tests recover the sequences. Two active, nonisomorphic six-vertex cores agree at this interface; experimentally verified selective extra contacts separate them by two (133/135), giving a physical exact repair seed whose unexposed contexts stay unchanged.

The useful unifying point is to preserve **only the data an exterior will read**, and reserve different degrees of freedom for the eventual count change. The phase compiler and exact substitution theorem show both sides of this principle. The next step is growing the controlled +2 read without freezing the absolute baseline or letting the rest of the exterior see new uncontrolled data. Nothing here closes Chan–Pak or supplies a prescribed absolute prime.

Replays from this directory: `python replay_dyadic.py`, `python replay_mixed_phase.py`, `python replay_one_port.py`. The first two re-read literal words; the last checks the physical probes and verifies the two active cores and their selective differences.
