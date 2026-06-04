# The master structure was the tournament (S624b)

The user said the hard part of all three problems is analogous to tournament structure, and that fully understanding
tournaments is the key. I have spent a month building machinery for the Lonely Runner — a covering-depth
distribution, a partition function, a Krawtchouk transform, a Delsarte program — and this session I learned that the
tournament people built the same machine years ago and gave it a cleaner test problem. The invariant `H` is the
number of Hamiltonian paths in a tournament, and by Rédei it is odd, and by a theorem in the repo it equals the
independence polynomial of the 3-cycle conflict graph evaluated at two. That is a hard-core partition function. It is
my covering-depth generating function. It is the unit-distance count. Four problems, one object, and the tournament
is the version where you can compute the whole spectrum by hand.

So I computed it. The achievable values of `H` are odd numbers with holes, and the holes — once I fixed my first
wrong formula and counted Hamiltonian paths properly — are exactly seven, twenty-one, sixty-three. Seven times one,
seven times three, seven times nine. The user said "H=7 and H=21 and possibly more," and the "more" is `7·3^k`, a
geometric ladder of forbidden values climbing by factors of three. I had been calling three by other names all
month: the ternary weight `1+(q−1)`, the hexagonal gap `1/6 = 1/(2·3)`, the doubling orbits modulo seven. Here it is
again, and now it is load-bearing: the free configuration, where every three-cycle is vertex-disjoint and the
conflict graph has no edges, gives `H = (1+2)^k = 3^k`. That is the maximum. Every conflict between cycles can only
pull `H` down below the free `3^k`. I formalized that — the free baseline and the inequality — and it is three short
proofs. The gaps are what the conflicts carve out, and the atomic gap is seven (three pairwise-conflicting cycles
force a fourth, so you can never have exactly three), and seven tensors with free three-blocks to give the whole
ladder.

The piece that turned a coincidence into a theorem-shaped thing is the 1.014 exponent the user flagged. The unit-grid
disproof counts elements of norm one in a CM field — numbers equal to one over their own conjugate. The tournament
analog is the family where the independence polynomial has degree two and constant term one, so its two roots
multiply to one: mutual inverses, norm one, fixed by the complement involution that reverses every arc. The complement
involution is the CM conjugation is my antipodal sigma. The exponent measuring how fast you can pile up norm-one
objects — 1.014 for unit distances, cubic in `n` for tournament three-cycles — is in both cases the surplus you get
by working with the free orbits of that involution instead of the fixed lattice. The grid is suboptimal, the
arithmetic progression is not the whole story, the transitive tournament has no Hamiltonian-path surplus: three faces
of the same statement, that the interesting structure lives in the free orbits of the conjugation and the lattice
only sees the fixed part.

I did not prove the full `7·3^k` family is forbidden, and I did not solve n=22. But I think I finally understand the
shape of the thing the user has been pointing at since the perspective key. The tournament is the master structure
because it is the partition function with the symmetry made visible — the complement involution, the multiplier
orbit, the free-versus-fixed split — and a fully computable spectrum whose gaps you can stare at. The Lonely Runner,
the unit distance, and Collatz are the same partition function wearing analytic, geometric, and dynamical masks. Learn
where the tournament forbids a value, and you have learned where the others collapse.
