# The three prime sets as fixed-point regimes, not only order factors

Status: PROVED finite-group lemma and exact application using the cited Monster/Baby Monster orders and centralizer2.B. No new finite-simple-group theorem is claimed.

## General lemma

Let a finite group G act transitively on Ω, with point stabilizer H. Fix a prime p dividing |G| and write a=v_p(|G|), b=v_p(|H|). A subgroup K fixes a point of Ω if and only if K is contained in a conjugate of H.

Consequently:

1. If b=a, every Sylow p-subgroup of G fixes at least one point. Choose a Sylow p-subgroup inside H, and conjugate it to any chosen Sylow subgroup.
2. If 0<b<a, some subgroup of order p fixes a point, but no Sylow p-subgroup fixes any point. The first assertion follows from Cauchy's theorem in H; the second follows from orders. This does not say every element of order p has a fixed point.
3. If b=0, no subgroup of order p fixes any point.

More precisely, b is the largest exponent k for which some subgroup of order p^k fixes a point. A Sylow p-subgroup of H attains b; containment excludes larger k. Subgroups of all intermediate orders exist inside a finite p-group. The deficit a−b equals v_p([G:H]).

## Monster action on2A involutions

Take G=M, Ω the2A involution class, H=C_M(t)=2.B. The actual index is |M|/(2|B|). The orders used are

|M|=2^46 3^20 5^9 7^6 11² 13³ 17·19·23·29·31·41·47·59·71,
|B|=2^41 3^13 5^6 7² 11·13·17·19·23·31·47.

Therefore the user's partition has a precise action-theoretic interpretation:

- {17,19,23,31,47}: b=a=1. Every Sylow p-subgroup centralizes some2A involution.
- {2,3,5,7,11,13}: 0<b<a. Some order-p subgroup centralizes a2A involution, but no full Sylow p-subgroup does.
- {29,41,59,71}: b=0. No order-p subgroup centralizes any2A involution.

The centralizer deficits on the six small primes are (4,7,3,4,1,2), not (5,7,3,4,1,2); the difference at2 is exactly the double cover in2.B. Dividing by the order of B instead produces the user's integer after stripping the four missing primes. A numerical order ratio must not be called the index of a subgroup B without an embedding; here2.B is the actual subgroup.

Thus the three sets encode three distinct fixed-point regimes of a named group action. This is a stronger interpretation than saying that the same primes occur in two factorizations. It does not give a tournament action, identify7 with21, or preserve Hamiltonian path counts.

## Relation to the modular-curve calculation

For odd good-level primes, the companion note gives a separate cited geometric criterion using supersingular rationality on X_0^+(2). The four primes absent from B have a conjugate pair of coarse2-isogeny points even though their underlying supersingular j-vertices are rational. This supplies a concrete example of information retained by an edge/level structure and lost by forgetting it.

The fixed-point lemma above distinguishes all three prime sets; level-two rationality by itself distinguishes the absent set from the two present sets, not their differing p-adic stabilizer depths. These are compatible views with different resolution, not identical invariants.

Primary source trail for group orders, the2.B centralizer, and modular criteria is in `monster/MONSTER_BABY_PRIME_PARTITION.md`; the elementary action proof above is complete on its own once those inputs are supplied.
