# Independent audit: 7, 13 and 19

Outcome: PASS for the finite base and p-adic extension, conditional on inherited gadget legality. Independently written forward BFS in `audit_7_13_19.c` reaches 5,806,080 states of 5,978,882 ambient states, and explicitly tests all states: there are no primitive misses and no nonprimitive hits. The expected count is 2(7²−1)(13²−1)(19²−1)=5,806,080, including both B-parity layers.

`tangent_7_13_19.py` independently recomputes matrices, orders, and tangent vectors. For L=B^144 and M=Q L Q^(-1), both loops are identity modulo each prime; B orders mod 7,13,19 are 16,12,18. The two tangent vectors at z=(25,34), with columns in (L,M) order, and their determinants are:

  7: (3,4), (1,3), det 5;
  13: (8,1), (8,10), det 7;
  19: (6,15), (16,5), det 18.

The conjugation action order and even parity are as in `PREFIX_AUDIT.md`. For target prime p at current precision p^r and other primes q at q^s_q, use E=p^(r−1)∏q≠p q^(s_q−1). Since both loops are I mod every q, raising to E preserves each q^s_q. Mod p^(r+1), the next digit shifts by ∏q≠p q^(s_q−1) times the same two basis vectors; the scalar is invertible mod p. Thus arbitrary next digits can be chosen while preserving previously fixed coordinates and parity. Pullback by a base word transfers the seed-fibre result to any primitive target in either parity. Odd modulus still makes the integral X₂≡g mod 2 condition compatible with all target residues.

The proof claims congruence coverage and inherited strong-tournament realization, not exact integer coverage or primality. The source packet and independent audit code and output are in this directory.
