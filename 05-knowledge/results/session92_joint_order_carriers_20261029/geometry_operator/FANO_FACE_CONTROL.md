# Seven Fano forbidden triples give 64 surviving sectors

**Session92, 2026-10-29. PROVED application plus FINITE-EXACT control.**
Companion to `HOLONOMY_FACE_CARRIER.md`; the original frozen files are unchanged.

Take seven named modules indexed by the nonzero vectors of F2^3. Let the
minimal forbidden active sets be its seven lines `{a,b,a+b}`. Apply exactly
the signed boundary construction of the companion theorem. The face counts,
indexed by cardinality0 through4, are

```
1, 7, 21, 28, 7.
```

Consequently the leakage Gram on3^7=2187 feasible states has

```
nullity = 64,       rank = 2123.
```

This is not an identification with the tournament folded-path64-state model.
It is a separate concrete application of the phase/face theorem.

## Counting proof

There is one empty face, seven singleton faces, and all21 pairs. Of the35
triples exactly seven are lines, leaving28. A nonline triple `{a,b,c}` is a
basis of F2^3. To extend it to a line-free four-set, its fourth point cannot
be any of `a+b,a+c,b+c`, so must be `a+b+c`. This extension is unique. Every
four-set contains four triples, giving28/4=7 four-faces.

For this four-face the omitted three points are
`{a+b,a+c,b+c}`, a Fano line. Conversely its complement is line-free. No
five-set is line-free: starting from a four-face, every omitted point is a
sum of two of its points. This proves the entire face vector without relying
on a matrix computation.

The seven four-faces also equal the supports of the seven nonzero linear
forms on F2^3 evaluated at its seven nonzero points. They therefore are
precisely the nonzero supports of the binary simplex[7,3,4] code, dual to the
classical Hamming[7,4,3] code. Here this statement is checked directly from the
linear-form definition, not inferred just from equal cardinalities. It gives
an actual set-level bridge to the repository's Session85 Fano/Hamming work.
The present64-face family is **not itself a linear code**: it contains all
seven singleton vectors but omits every line vector. Indeed a singleton and
a disjoint two-subset of a line are both faces and their XOR is not. Its size
64 does not imply six independent or commuting bit toggles.

## Exact scope and reproduction

`check_fano_face.py` performs gain propagation on all128 occupied-module
sectors, independently testing balancedness against line avoidance. It never
forms the dense2187-by2187 Gram. Each cube vertex is a named port-choice
assignment; the propagation ratio is the negative of the actual boundary
coefficient. A conflict along a square or cycle kills that sector's kernel.
The checker also reconstructs the four-faces from nonzero linear forms and
checks a literal XOR failure.

Both normal and `-O` execution exit0 with **264 checks**, face vector
`(1,7,21,28,7)`, nullity64 and rank2123. The rank conclusion uses the proved
real-positive Gram theorem; gain enumeration alone is not a finite-field
Gram-rank claim. There is no asserted tournament/code isomorphism, no
six-qubit interpretation, and no enumeration speedup claim.
