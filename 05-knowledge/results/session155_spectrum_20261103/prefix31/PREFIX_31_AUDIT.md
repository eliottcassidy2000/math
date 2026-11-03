# Independent audit of 31 kernel fibre and lift

Status: PASS for algebraic statements, conditional on inherited legality and source-word realization. Reproduction: `python audit_31.py`, recorded in `audit_31.out`.

Set P=7·13·19=1729, Q=[[11,4],[18,7]], B=[[5,2],[8,3]], T=19152, A=Q^T and C=B A B⁻¹. Computation independently verifies A=I mod P; consequently C=I mod P. In the parity-augmented finite group, A and C each have B-parity zero. At modulus 31, A=(27,5;7,22), C=(12,15;13,6), and forward BFS from (25,34) mod 31 obtains every one of the 960 nonzero vectors. Since A and C are invertible modulo 31, their finite positive-word monoid is a group; the full seed orbit makes it transitive on all nonzero vectors.

To produce any four-prime base state, first use the established 1729-orbit theorem for the required old-base state and parity, applying positive word W to the seed. Its image modulo 31 remains nonzero because both generators are invertible. Choose a word in A,C moving this image to the target mod 31; these matrices fix the old-base state and parity. If target modulus N is larger, B⁻¹ is a positive word with the correct odd parity: choose the even order S of (B,1) in GL₂(Z/N)×Z/2 and use B^(S−1). Temporal order for C on columns is B⁻¹, then A, then B. Analogously Q⁻¹=Q^(ord(Q)−1) has parity 0.

For all higher prime powers, use L=B^576 and M=Q L Q⁻¹. Both are identity modulo each of 7,13,19,31. The displacement-vector determinants at z₀ are respectively 3,8,3,13, all nonzero. At target prime p and current exponents r_q≥1, E=∏_q q^(r_q−1) raises the p-adic precision by one via the two independent tangent vectors, while preserving every other q^r_q. Pull back the final target to the seed fibre by an inverse of a parity-correct base word, fill the seed fibre with even loops, and reapply that word. All inverses are justified in the augmented finite group; positive words can be chosen separately for each modulus.

This verifies residue coverage, not equality for exact target integers or prime generation.
