# The regular-module bound and its sharp three-cycle structure

**Status: PROVED.** Elementary argument, reviewed within this research session. No prime-degree group classification, odd-order classification theorem, or novelty claim is used.

## 1. The regular-module bound

Let T be a regular tournament on d vertices. Let A be a proper nontrivial homogeneous module: 1<m=|A|<d, and every vertex outside A either beats all of A or loses to all of A. Then

    d ≥ 3m.

At equality the whole tournament is a cyclic blow-up of three regular m-vertex tournaments, with A as one of the three modules.

**Proof.** Put r=(d−1)/2. All vertices of A have the same external outdegree, so regularity of T forces the induced tournament T[A] to be regular. Thus m is odd and its internal outdegree is (m−1)/2.

Let B be the d−m outside vertices. Let L consist of vertices that beat A, and U those beaten by A. The outdegree of a vertex of L within T[B] is r−m; that of a vertex of U is r. The sum of these outdegrees is the number of arcs of T[B], hence

    |U|r + |L|(r−m) = (d−m)(d−m−1)/2,
    |U|=|L|=(d−m)/2=:h.

Since d and m are odd and m<d, h≥1. The h vertices of L induce a tournament with h(h−1)/2 arcs. Their total outdegree in T[B] is h(r−m), so

    h(r−m) ≥ h(h−1)/2.

Cancel h and substitute r and h. This gives d≥3m.

If d=3m, the inequality is an equality, so no vertex of L sends an arc to U. Consequently U→L, while by definition L→A and A→U. All three blocks have m vertices and all cross arcs are uniform. Each vertex beats exactly m vertices outside its block, so the three induced tournaments are regular, with internal outdegree (m−1)/2. Thus A,L,U are homogeneous modules arranged in a directed triangle. QED.

Conversely, the cyclic blow-up of any three regular m-vertex tournaments is regular and has those three homogeneous modules, each of size one third of the whole. This proves sharpness of the module bound for every odd m≥3.

## 2. The single-cycle support corollary

Suppose an automorphism of a regular d-vertex tournament is a single cycle of odd length m, fixing every other vertex, where 1<m<d. Its support A is homogeneous: every fixed outside vertex has a constant relation to the cycle orbit. Therefore d≥3m by §1. At equality the same directed-triangle decomposition holds, and the induced tournament on A is cyclic.

The excluded case m=d is essential: a regular cyclic tournament has a full cycle on its entire vertex set.

In particular, an odd transitive permutation group of degree d cannot contain a single m-cycle when m<d<3m. To see this without any classification theorem, an odd group cannot interchange the endpoints of an unordered pair: such an element would have even order. Orient each pair orbit consistently to obtain an invariant tournament. Transitivity makes its outdegrees equal, so it is regular, and the support bound applies.

## 3. Sharp transitive realizations

Take three identical regular cyclic tournaments of odd order m and arrange their blocks in a directed triangle. Rotation of one block is a single m-cycle automorphism. Independent rotations of the three blocks and a cyclic block permutation generate the transitive group C_m wr C3, of odd order 3m³.

Thus the single-cycle threshold is sharp for every odd m≥3, not only for primes. At m=5 its first possible degree is fifteen. At degree thirteen it excludes a single 5-cycle in an odd transitive action, providing the source-side Sylow input for [the exact double-triple boundary](SHARP_N13_BOUNDARY.md).

This is an elementary regular-tournament module observation, not a claim of a new theorem in the literature. The prime-five specialization is used for Sylow bookkeeping; primality is unnecessary in the geometric bound.
