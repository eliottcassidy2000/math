# An exact FRACTRAN operator and an ear-image filter

Let F_0(x)=x² and, for a series F=∑_{n≥1}c_n x^n, let C_d F(x)=∑_{k≥1} c_{dk}x^k. Given an **ordered** FRACTRAN program (a_i/b_i)_{1≤i≤m} in lowest terms, define, for J⊆{1,…,i−1}, d(i,J)=lcm(b_i,{b_j:j∈J}). Put

T F(x) = ∑_i ∑_{J⊆{1,…,i−1}} (−1)^|J| (C_{d(i,J)}F)(x^{a_i d(i,J)/b_i}).

Then F_{t+1}=T F_t. Proof: inclusion–exclusion selects exactly integers divisible by b_i and by none of b_1,…,b_{i−1}. Their image is n a_i/b_i. An input for which no fraction applies contributes nothing. This works coefficientwise, and preserves the priority rule, unlike a naive sum over all applicable fractions. The state-and-time formal series G(y,x)=∑_{t≥0}y^t F_t satisfies G=x²+y T G. Its pure-power observation is Q(y,z)=∑_{p≥1}z^p [x^{2^p}]G(y,x). The seed contributes y^0 z and is not prime output. This is an exact recurrence for the machine, not a new prime formula or a claimed equivalence with ear insertion.

For Conway's standard PRIMEGAME, the positive-time pure-power exponents are the primes (the established theorem). If one forgets time and retains support, the resulting indicator is P(z)=∑_{p prime}z^p. The printed standard program has coordinate primes 2,3,5,7,11,13,17,19,23,29; 27=3³ is composite and cannot replace the coordinate 29. Both standard and Guy variants and their actual first 15 output exponents and arrival times are recorded in `math/05-knowledge/results/procgen_atlas_20260924.out`, L1. In particular the tenth output is 29, not 27.

For the specified order-five ear carrier, E(z)=2z^15+4z^17+2z^19+6z^23+2z^25+2z^27+2z^29+4z^33+4z^37+2z^41. The coefficientwise (Hadamard) product gives the exact prime-filtered ear enumerator

E⊙P = 4z^17+2z^19+6z^23+2z^29+4z^37+2z^41.

Its total coefficient 20 counts labeled nonconstant cuts whose path count is prime; the remaining 10 have composite path count. This is a canonical filter once both data sources and the equal-integer comparison are chosen, not a bijection of machine states with tournaments. For example, the four cuts with path count 17 and the one PRIMEGAME visit with exponent 17 already preclude a multiplicity-preserving bijection at that value. Four equal counts can be related to the output only by adding four arbitrary labels to the machine event. No method converting a FRACTRAN transition into ear addition, or vice versa, has been exhibited.

The ear polynomial and carrier are verified by `../ear/test.py` and `../ear/README.md`. PRIMEGAME details and exact output times: `math/05-knowledge/results/procgen_atlas_20260924.out` lines near the PRIMEGAME section. The program's coordinate list and Conway source are identified in [`HISTORY.md`](../history/HISTORY.md) and `05-knowledge/results/session148_square_sum_20261103/prime_alignment/README.md`; primary historical source: https://gwern.net/doc/cs/computable/1987-conway.pdf.
