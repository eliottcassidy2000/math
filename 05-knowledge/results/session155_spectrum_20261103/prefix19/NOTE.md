# Mixed-prefix orbit at all 7-, 13-, and 19-powers

**Result.** For every a,b,c≥1, put N=7^a 13^b 19^c. In each parity layer for the number of b prefixes, the pairs (Y,X) realized modulo N by the actual source words are exactly those primitive at each of 7,13,19. In each layer there are

    (7^(2a)−7^(2a−2)) (13^(2b)−13^(2b−2)) (19^(2c)−19^(2c−2))

such pairs. The same statement for any exponents a,b,c≥0 follows by projecting from positive exponents. This is a count of residue states, not words, tournaments, or exact Hamiltonian-path counts.

In the even layer, the inherited literal strong-tournament companions have Hamiltonian-path counts H=3(Y+6)+2=3Y+20 and H=3(Y+6)+4=3Y+22, respectively. Each companion separately hits every residue modulo N infinitely often, and thus every odd residue modulo 2N infinitely often. Choose X=1, and then Y by the appropriate affine congruence; 3 is a unit modulo N. The source lattice restriction X≡g mod 2 causes no conflict: N is odd and each residue is realized by actual prefix words, whose integer states automatically satisfy the lattice condition.

A stronger valuation corollary: let t7,t13,t19≥0, D=7^t7 13^t13 19^t19, and let j be any unit modulo N. Each companion has infinitely many counts H with exact p-adic valuations t_p at all three primes and H/D≡j mod N. Prescribe H≡D j mod D N. For example, any desired 13- and 19-power residue of H/7^t can be achieved while requiring v7(H)=t, without requiring restrictions on the quotient at 13 or 19. These are congruence results, not primality or equality to a chosen integer.

## Conventions and finite base

Actions are on columns z=(Y,X)^t, with starting state z0=(25,34)^t:

    Q = [[11,4],[18,7]],    B = [[5,2],[8,3]].

Q and B denote the q2 and triangle prefix actions, and e is the parity of B actions. In the inherited positive-companion construction g=6 when e=0. The coordinate action and tournament interpretation are inherited from `golden_inverse_defects_20261103/TRANSPORT.md` at f6d5304849e50b29df6e368239e92fffa92e57e3 and the S154 mixed-prefix packet at fb8bc7775a0678ecf6ef4169faff452220573572. Matrix multiplication is right-to-left on columns; when we speak of an action sequence, we mean its chronological application order, not a literal prefix spelling.

The supplied `orbit.c` exhausts the forward orbit modulo m=1729 from (e,z)=(0,z0), independently exhausts the reverse orbit, and separately scans all 2m² ambient states for each direction. Both are exactly the primitive predicate, with 5,806,080 states total, i.e. 48*168*360 per parity. Since det Q=5 and det B=−1, neither action can destroy primitivity. This is a finite exhaustive certificate for the only base computation used below, not extrapolation from a sample. Reverse reachability agrees with forward reachability because the generators are invertible elements of the finite group GL₂(Z/m)×Z/2 and positive powers realize their inverses.

## Two loops and p-adic lifting

Let L=B^144 and M=Q B^144 Q^(-1), with inverse in the appropriate modular ring. B has orders 16,12,18 modulo 7,13,19, respectively, so L=M=I modulo 1729; each has even B-parity. Exact arithmetic gives

    p       (L−I)z0/p mod p       (M−I)z0/p mod p       determinant
    7             (3,4)                 (1,3)                 5
    13            (8,1)                 (8,10)                7
    19            (6,15)                (16,5)                18.

The two displayed columns form a basis of F_p² for each p. All quotients are computed modulo p², where the numerator is divisible by p. The supplied `lift.py` independently computes the table.

For the final modulus N, let T be the order of Q in GL₂(Z/N). Then Q^(-1)=Q^(T−1), so M is a positive action sequence: apply Q^(T−1), then B^144, then Q. L and M and any of their nonnegative powers are therefore realizable by source-prefix actions, and they preserve parity. T depends on N and need not be computed to prove existence.

Suppose v=z0 modulo 1729 and currently matches a chosen target modulo p^r for each p in P={7,13,19}, where the depths are r_p≥1. To advance the p-digit while preserving all the other established residues, use the exponent

    E_p = p^(r_p−1) * product_{q in P, q != p} q^(r_q−1).

For A=L or M, write A=I+p C_p. The binomial theorem gives

    A^E_p v = v + p^r_p (product_{q != p} q^(r_q−1)) C_p z0
              mod p^(r_p+1),
    A^E_p v = v mod q^r_q for every q != p.

The scalar product over q≠p is a p-adic unit. Because the two C_p z0 vectors form a basis modulo p, appropriate powers of L^E_p and M^E_p select any desired next p-digit. Their first-order changes add modulo p^(r_p+1); take coefficients 0,...,p−1, so negative powers are unnecessary. Inductively this reaches every target in the fibre of z0 modulo 1729, preserving even parity.

For any primitive target w modulo N and either parity e, the finite base certificate gives a positive word W of parity e taking z0 to w modulo 1729. Its invertibility implies W^(-1)w lies in the seed fibre modulo 1729. Reach that pullback using the even loops, then apply the same word W. This proves the entire orbit and the cardinality formula by CRT.

For infinitude, once a legal word reaches a desired companion residue, repeat Q^T arbitrarily many times. It is the identity modulo N and does not change b parity. On actual positive integer states, Q strictly increases Y, hence the companion counts are unbounded and distinct. This step uses the inherited legality of repeated q2 prefixing; modular algebra alone does not supply tournaments.

## Reproduction

    cc -O3 -std=c11 -Wall -Wextra -o orbit orbit.c && ./orbit
    python3 lift.py

`orbit.c` checks both forward and reverse base orbits in the entire ambient space. `lift.py` verifies all loop displacements and runs 120 deterministic, constructive seed-fibre lift tests with depth triples (2,2,2), (3,2,2), (2,3,2), (2,2,3). The tests corroborate the induction but do not replace it. No claims are made about other prime factors, exact spectral coverage, independence of the two companions, or primality.
