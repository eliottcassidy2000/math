# The exact 15/17 bridge: equal multiplier-orbit counts, unequal symmetry

Status: elementary PROVED corollary of CITED inputs, with exact checks
and same-agent review. No independent audit, natural isomorphism, or new
classical circulant classification is claimed.

For odd n>1 let S range over subsets of Z/n\{0} containing exactly one
member of each pair{x,-x}. There are 2^((n-1)/2) such connection sets.
The unit group U(n) acts by multiplication. These are **multiplier orbits
of presentations**, with the cyclic coordinate retained. At composite n
we do not silently identify them with full tournament isomorphism classes.

Akin's [Theorem 3.13](https://arxiv.org/html/1806.11241v4)
proves that this unit action is free iff n is a squarefree product of
Fermat primes. Its fixed-subset construction (Theorem 3.14) also gives the
odd-stabilizer obstruction outside that family. The repo already imports
these results in THM-4644; they are not new claims here.

## F1. A finite ladder of exact equinumerosities

Put F_m=2^(2^m)+1. Suppose every F_0,...,F_m is prime, with m>=1.
By the classical Fermat product identity,

    n_-=F_m-2=product_(j<m) F_j,
    n_+=F_m.

Both unit actions are free. Their totients are

    phi(n_-)=2^(2^m-1), phi(n_+)=2^(2^m).

Hence the number of multiplier orbits is the SAME at both orders:

    2^((n_--1)/2)/phi(n_-)
      =2^((n_+-1)/2)/phi(n_+)
      =2^(2^(2^m-1)-2^m).                             (F1)

One extra negation pair doubles the connection-set count; the totient
doubles as well. This exact cancellation explains the equality.

| m | orders | multiplier orbits at each order |
|---:|---:|---:|
|1|3,5|1|
|2|15,17|16|
|3|255,257|2^120|
|4|65535,65537|2^32752|

This initial-segment ladder stops unconditionally: F5 is composite, so the
hypothesis fails for every m>=5 regardless of undiscovered later Fermat
primes. This uses the same initial-segment principle already proved in
THM-4801, with a different observable. The following carry argument goes
further: it classifies every collision of this orbit count within the
free-action family, including any undiscovered Fermat primes.

## F1b. Complete collision classification by binary carries

**Theorem.** Among distinct odd integers greater than 1 for which the full
unit action on tournament connection sets is free, equal numbers of
multiplier orbits occur exactly at the four pairs in the table above.
No assumption about the completeness of the known Fermat primes is used.

For an auxiliary integer s>=0, define

    P(s)=product_(i: bit_i(s)=1) F_i,
    E(s)=(P(s)-1)/2-s.

These auxiliary formulas use all Fermat NUMBERS, whether prime or not.
Their binary representation is the Pascal/Sierpinski row: expanding the
product gives P(s)=sum_(t a binary submask of s) 2^t, with no carries.
Indeed in characteristic 2, (1+X)^s is the product of 1+X^(2^i) over
the set bits i of s. This supplies the exact row map behind the earlier
Sierpinski/initial-segment connection in THM-4801.
If n is a squarefree product of Fermat primes, it has a unique index s
with n=P(s), phi(n)=2^s, and orbit count 2^E(s).

Suppose s has exactly m trailing binary ones. Its factors below bit m
multiply to F_m-2, while its factors above bit m have some product A.
Incrementing s replaces precisely those lower factors by F_m. Therefore

    P(s)=A(F_m-2), P(s+1)=A F_m,
    E(s+1)-E(s)=A-1 >= 0.                            (F2)

Equality holds exactly when A=1, that is s=2^m-1. For s>=1 these zero
increments are isolated. Hence E(s)=E(t) for 1<=s<t iff
(s,t)=(2^m-1,2^m), m>=1. This implies
(P(s),P(t))=(F_m-2,F_m), and both belong to the free-action family exactly
when every F_0,...,F_m is prime. The known factorization of F_5 rules out
all m>=5, even if later individual Fermat primes exist.

The exclusion of n=1 matters: E(0)=E(1)=E(2)=0 gives an auxiliary triple
before that exclusion. This proof classifies collisions of the multiplier
orbit count, not full composite-order tournament isomorphism classes.
It is an elementary corollary of Akin's action criterion and the product
identity; no external novelty is claimed.

## F2. Why equal counts do not supply the desired symmetry map

U(15)=C2 x C4, whereas U(17)=C16. Even after adding a one-bit factor,
U(15) x C2 is not C16 (their exponents are 4 and 16). No surjective group
homomorphism C16->U(15) exists because a quotient of a cyclic group is
cyclic. Thus a one-bit refinement that is required to lift the FULL unit
actions by a quotient homomorphism is impossible. An arbitrary bijection
between the two 16-element orbit sets exists; a port-compatible or canonical one has
not been supplied.

There is a second difference. The circulant on Z15 with

    S={1,3,4,6,7,10,13}

is the lexicographic tournament C3[C5], with blocks given by residues mod 3.
Here C3 and C5 denote the regular cyclic tournaments, not bare cycle graphs.
Independent cyclic rotations of the three five-vertex blocks, together
with cyclic block rotation, give a subgroup C5^3 semidirect C3 of order 375.
Its multiplier stabilizer is nevertheless trivial. Reversing the roles,
C5[C3] is circulant on Z15 with

    S={1,2,5,6,7,11,12},

and has a subgroup C3^5 semidirect C5 of order 1215. Exact full automorphism
orders are not needed or newly claimed. At prime 17, the inherited
THM-871/4644 theorem gives precisely C17 for every circulant tournament.

The distinction is therefore between a coordinate multiplier action and
all automorphisms of the underlying tournament. This is the same kind of
mark-forgetting issue that the A000568 attachment census must retain.
