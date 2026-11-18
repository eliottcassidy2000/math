# A complete projective switching chart, and its exact target limitation

**PROVED all-order:** all cut switches of a fixed bounded tuple can be physically realized inside one Paley host by moving a pole. **FINITE-EXACT:** two different five-vertex switching histograms. **OPEN:** selecting the switching class from an arbitrary Hamiltonian target.

## Inherited geometry

THM-4602, `01-canon/theorems/THM-4602-tournaments-as-chirotopes-pfaffian-criterion-and-legendre-frieze-strips.md`, already identifies quadratic-character projective configurations with Paley switching classes and representative rescaling with vertex switching. We do not claim the Möbius switching identity as new. The following adds the explicit bounded-host extension argument from `PALEY_TARGET_REFORMULATION.md` to realize **every** cut, and identifies the remaining discrete target obstruction.

Let q=3 mod 4 be an odd prime power, chi its quadratic character, and x_i distinct finite coordinates. For a fractional-linear map f(x)=(ax+b)/(cx+d), of nonzero determinant D=ad-bc, assume cx_i+d is nonzero for every i. Then

    chi(f(x_j)-f(x_i))
      = chi(D) chi(x_j-x_i) chi(cx_i+d) chi(cx_j+d).

This follows directly from f(y)-f(x)=D(y-x)/((cx+d)(cy+d)); a nonzero quadratic character equals its inverse. Thus a square determinant realizes cut switching, and a nonsquare determinant additionally takes the converse. This is generally not a tournament automorphism and does not preserve H.

## Every switch is physically available

Take the bound q>m^2 4^m from the previous note. For any desired signs s_i, the proved mixed-pattern estimate gives a pole a outside {x_i} with chi(x_i-a)=s_i. (Replacing chi(a-x_i) by its negative adjusts all prescribed signs and uses chi(-1)=-1.) Choose

    f_a(x)=-1/(x-a).

Its matrix is [[0,-1],[1,-a]], determinant 1, independently of the chosen signs. No pole is among the tuple, every image is finite, and injectivity preserves distinctness. Its induced tournament is exactly the cut switch s_i s_j of the original tournament. Global sign reversal gives the same cut, so there are 2^(m-1) possible labelled outputs. Conversely every square-determinant Möbius image avoiding infinity belongs to this switching class. Thus these outputs are exactly the full switching class, not an arbitrary family of tournaments.

This is an explicit finite construction: search poles until the desired sign pattern occurs. It chooses a physical move from prescribed cut data. It does **not** choose the cut from a prescribed Hamiltonian total.

## The complete invariant that cannot be changed by a pole

Write b_ij=+1 for i->j and -1 otherwise. For each distinct triple, the signed cyclic product b_ij b_jk b_ki is invariant under cut switching. These triple products completely determine the labelled switching class. Indeed, if b and b' have the same products, fix vertex 1, set s_1=1 and s_i=b'_1i/b_1i. Equality on triples (1,i,j) gives b'_ij/b_ij=s_i s_j. Therefore a pole move changes the reading within a class, but cannot change this triple-product datum. A new source is required when that invariant excludes the requested count.

## A conservation law inside each chart

For every fixed n-vertex tournament T,

    sum_{cuts s modulo global negation} H(T^s) = n!.

For each ordered permutation pi, orienting its n-1 consecutive edges forwards imposes equations s_pi(j) s_pi(j+1)=b_pi(j),pi(j+1). Because these edges form a spanning tree, the signs are determined uniquely up to simultaneous negation. Thus exactly one cut makes that permutation a Hamiltonian path. Sum over permutations. This is the repository's THM-1467 (`oriented-spanning-tree-switching-sum`), coupled to THM-474 (`tilings-are-switching-classes`), not a new conservation law or assertion of independence. It gives a mean n!/2^(n-1) but no surjectivity of the integer-valued reading.

The finite failure is already sharp at n=5. All 1,024 labelled tournaments give one of two switching histograms (count: multiplicity among 16 cuts):

    A: 1:5, 9:5, 11:5, 15:1
    B: 3:3, 5:6, 9:3, 13:3, 15:1.

Example A is the transitive tournament with outgoing bit rows [30,28,24,16,0]. Example B has rows [26,28,25,16,0]. Thus the complete projective chart of A misses the physically realized totals 3,5,13. The chart of B misses 1 and 11. Both have weighted sum 120. These finite obstructions refute the proposal that a single fixed-source Möbius chart is automatically a universal target selector, even when every cut is physically available. They do not obstruct larger orders or changing source classes.

## Connections with the requested older geometry

The older `planar_quintic_s5_petersen_resolvent_20261019.md` realizes the Petersen graph on the ten pairs of five labels, with disjointness adjacency, and recovers K5 from that incidence structure. Its A5 reduction is a quadratic-discriminant Galois reduction, not a physical A5 action on a tournament. Here the analogous useful discipline is retaining the character multiplier: projective coordinates alone forget the cut signs needed to identify a particular tournament. The numerical occurrence of five is not itself a target construction.

Pythagorean parametrization and the present pole maps both replace a quadratic restriction by rational coordinates, but the character signs are additional data that ordinary rational parametrization does not select. No implication from primitive triples to target coverage is claimed. The constructive deliverable is the pole search for an already specified cut, together with the triple-product invariant explaining exactly which freedom that search cannot provide.

## Reproduction

Run `python projective_replay.py` in this folder. It uses F_(3^9), q=19683>4^2 4^4, embeds the transitive four-vertex tournament, finds all eight cuts via poles, and directly checks every transformed edge and H by subset DP. It then checks all 1,024 labelled five-vertex tournaments, all their cuts, the conservation law and the two histograms above. `projective_replay.json` includes field polynomial, literal coordinates, poles, matrices, and H values. Finite tests do not replace the all-order proofs.
