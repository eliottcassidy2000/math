# Independent affine-reflection parity addendum

Status: PASS, explicitly restricted to affine-reflection-symmetric Hamiltonian paths in prime-order Paley tournaments. Reproduce `python session35-singer-audit/check_reflection.py`; exact results are in `reflection_receipt.json`.

## Theorem and proof

Let p≡3 mod4 be prime and let a Hamiltonian path in the Paley tournament on Fp satisfy p_(p−1−i)=b−p_i. Then the number of cyclic triangles containing no edge of the path is even.

Translate the unique reflection center b/2 to zero, so the reflection is r(x)=−x. It reverses every tournament arc, since −1 is a nonsquare. It preserves the set of cyclic triples and preserves the set of unordered path edges. Thus it acts as an involution on the admissible cyclic triples, and their number modulo two equals the number of fixed admissible triples.

A fixed three-element set is exactly {0,a,−a}. Choose the sign of a so that a is a square. The arcs −a→0→a hold, so the triple is cyclic exactly when −2a is a square, equivalently χ(−2)=1. For p≡3 mod4, the supplementary quadratic-character formula gives χ(−2)=1 iff p≡3 mod8, and −1 iff p≡7 mod8.

If p≡7 mod8, there are no fixed cyclic triples and the admissible count is even.

If p≡3 mod8, all (p−1)/2 fixed triples are cyclic. Exactly one contains path edges: the central consecutive triple {−a,0,a}. Indeed zero occupies the central path position, its two neighbors are opposite, and those are its only path edges. An edge between opposite nonzero vertices cannot occur in the path: their reflected positions have sum p−1, which is even, whereas adjacent indices have odd sum. Therefore exactly (p−3)/2 fixed cyclic triples are admissible. This number is even when p≡3 mod8. The involution proves the theorem.

This does not assert that every anti-automorphism or every blue path in every Paley tournament is induced by an affine reflection. Such a claim would need a separate automorphism theorem; none is used here.

## Independent checks and failed converse

The checker enumerates every reflected vertex ordering with center zero, filters actual directed Hamiltonian paths, then directly enumerates admissible cyclic triples. Other centers are translates, so this gives exhaustive reflection-symmetric controls at each stated order:

- p=3: one center-zero path, admissible degree 0, fixed admissible count 0.
- p=7: nine center-zero paths, degrees 2 (three paths) and 4 (six paths), fixed count 0.
- p=11: 185 center-zero paths, degrees 26,28,30,32,34, fixed admissible count exactly 4 throughout.

The converse already fails at p=11. The supplied path

(6,4,5,3,7,0,1,10,8,9,2)

is independently verified to be directed Hamiltonian and black under the full reversed-adjacency test, yet has 32 admissible cyclic-triangle reversals. Thus even degree cannot characterize blue paths beyond the special P7 result.
