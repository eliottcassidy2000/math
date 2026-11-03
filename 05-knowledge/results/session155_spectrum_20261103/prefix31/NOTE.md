# Mixed-prefix orbits at powers of 7, 13, 19 and 31

**Theorem.** For a,b,c,d≥1 and N=7^a 13^b 19^c 31^d, the mixed q2/triangle prefixes of the literal source seed realize, in each parity layer of the triangle-prefix count, precisely every pair (Y,X) modulo N that is nonzero modulo each of 7,13,19,31. Thus the number of residue states in each parity layer is

    product over p^r in {7^a,13^b,19^c,31^d} of (p^(2r)−p^(2r−2)).

Projections give the same conclusion if some exponents are zero. This counts residue pairs, not words, tournaments or exact integer counts.

In the even layer g=6, h=Y+6, and each of the two inherited legal strong-tournament companion families, H=3h+2=3Y+20 and H=3h+4=3Y+22, separately takes every residue mod N infinitely often. Choose X=1 and solve the corresponding affine equation for Y, as 3 is a unit modulo N. Since H is odd, every odd class modulo 2N is attained infinitely often. The original integer condition X≡g mod 2 is obeyed by the constructed source words and does not obstruct residues modulo odd N.

For any t_p≥0 and D=product p^t_p, and any unit j modulo N, each companion family has infinitely many H with exact v_p(H)=t_p for each p in {7,13,19,31} and H/D≡j mod N: request H≡D j mod D N. This does not assert primality or exact prescribed-integer coverage.

## Actions and finite base

Use columns z=(Y,X)^t, the literal seed z0=(25,34)^t, and the chronological prefix actions

    Q = [[11,4],[18,7]],     B = [[5,2],[8,3]].

Keep parity e of the number of B actions. The source interpretation, legal words, g and companions are inherited from `golden_inverse_defects_20261103/TRANSPORT.md` at f6d5304849e50b29df6e368239e92fffa92e57e3 and the S154 mixed-prefix packet at fb8bc7775a0678ecf6ef4169faff452220573572. Matrix products act right-to-left. Q and B are invertible at all four primes, so primitivity is necessary.

Set m=7*13*19=1729. `orbit.c` computes the forward and reverse orbits from (0,z0), scans the entire ambient space and checks they comprise all primitive pairs in either parity; each direction visits exactly 5,806,080 states. Since the augmented generators are invertible in the finite group GL₂(Z/m)×Z/2, inverses are positive powers; the finite forward certificate itself suffices.

We add prime 31 without enumerating all (1729*31)² pairs. Set

    T=19152,    A=Q^T,    C=B A B^(-1).

Exact arithmetic shows A=I modulo 1729. Therefore both A,C fix every old-base vector and have even B parity. Modulo 31 they are

    A=[[27,5],[7,22]],    C=[[12,15],[13,6]].

The independent finite BFS `kernel.py` starts at z0 mod31 and applies these two matrices; its orbit is exactly all 960 nonzero vectors. Since these are invertible maps on a finite set, their positive semigroup is a group. One full nonzero orbit implies transitivity from any nonzero starting vector.

Given any primitive desired vector modulo 1729*31 and either parity, first use the old-base certificate to choose a positive prefix-action word W giving its old-base residue and parity. Its 31-component is nonzero. Then append a positive word in A,C taking that 31-component to the desired one. It fixes the old-base vector and parity, proving the full four-prime base orbit. No assumption of independence of unrelated phase indices is used.

All inverse symbols have positive word realizations at every finite modulus N: Q^(-1)=Q^(ord(Q)−1). For B, take S the order of (B,1) in GL₂(Z/N)×Z/2; S is even and B^(-1) is represented by B^(S−1), which has odd B parity. Hence C=B A B^(-1) is represented by an even-parity positive sequence. Temporal order for C is B^(-1), then A, then B.

## Lifting all four exponents

Let L=B^576 and M=Q B^576 Q^(-1). B has respective orders 16,12,18,64 modulo 7,13,19,31, whose least common multiple is 576. Thus L,M=I at each prime and both have even B parity. The exact first-order columns at z0 are

    p       (L−I)z0/p mod p       (M−I)z0/p mod p       determinant
    7             (5,2)                 (4,5)                 3
    13            (6,4)                 (6,1)                 8
    19            (5,3)                 (7,1)                 3
    31            (0,23)                (21,0)                13.

All determinants are units. The inverses in M again have positive realizations: if R is the order of Q mod N, first apply Q^(R−1), then B^576, then Q.

For a current vector v≡z0 modulo all four primes, suppose its residue is fixed to the target modulo p^r_p for each p. For one chosen p, put

    E_p=p^(r_p−1) product_{q != p} q^(r_q−1).

Writing U=I+p C_p for U=L or M and expanding binomially gives

    U^E_p v = v + p^r_p (product_{q != p} q^(r_q−1)) C_p z0
              mod p^(r_p+1),
    U^E_p v = v mod q^r_q  for q != p.

The scalar product is a p-adic unit and the two columns form a basis, so nonnegative powers of L^E_p and M^E_p select any next p-digit while fixing all other prime-power residues. Iterate to realize all lifts of z0. To reach an arbitrary primitive target w at given parity, select a positive base word W giving w modulo 7*13*19*31, pull w back by the invertible W, build that seed-fibre lift with the even loops, and apply W. CRT gives the cardinality above.

To obtain infinitely many counts for any attained residue, repeat Q^R for R=ord(Q mod N). This preserves parity and residues, and strictly increases the actual positive Y coordinate; inherited source legality gives distinct actual tournament counts.

## Reproduction and scope

    cc -O3 -std=c11 -Wall -Wextra -o orbit orbit.c && ./orbit
    python3 kernel.py
    python3 lift.py

The first program exhaustively checks the 1729 base in both directions. The second checks Q^19152=I mod1729, the matrices A,C and their entire 960-point fibre orbit. The third checks the four tangent bases and runs 150 deterministic constructive seed-fibre lift tests at depth vectors (2,2,2,2), (3,2,2,2), (2,3,2,2), (2,2,3,2), (2,2,2,3). Tests corroborate, rather than replace, the proof. No claims are made for other primes, independent selection of the two companions, exact spectral coverage, or primality.
