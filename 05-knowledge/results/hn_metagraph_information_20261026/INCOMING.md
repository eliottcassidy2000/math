# Final incoming work: a positive metric escape and prime-family counting

2026-10-26. Commits075c9f922c,15825436d2 and62430ddb33 arrived during final
sync and were rebased before publication. Concurrent results below are
recovered, not claimed as new. See [the integration checker](audits/check_incoming.py)
and its [exact receipt](audits/incoming.json).

## A different pairing escapes the antipodal-frame theorem

The concurrent [planar packet, U1--U3](../rank_geometry_metagraph_20261026/unit_distance/RESULT.md)
independently excludes the original two-loop candidate by a stronger rank
argument without antipodality: double flips force three nonzero pair-difference
vectors to be mutually orthogonal. Our MT1--MT2 instead give the canonical
section for arbitrary networks of antipodal planar frames.

Its positive alternative matters. On a unit regular hexagon with vertices
z0,...,z5, use pairs `{z0,z3}`, `{z1,z5}`, `{z2,z4}`. Reflection fixes the first
pair pointwise and swaps both others, acting as even translation110 with no
fixed endpoint frame. Squared pair lengths are4,3,3. Two pairs are not
antipodal, so this lies outside MT1. All15 metric distances and all8 frame
actions are independently replayed here.

This translation is also one of the four deterministic maps passing our
unprescribed two-frame S6 gate. Both necessary gates leave it alive; neither
makes it a forced colouring relation. The unit-wheel table is still
`w(beta)=[proper on C6](6-|image beta|)`, with14910 nonzero boundary assignments
and24600 six-colourings. Its zero-weight rainbow chart does not cover all
colourings. This gives a concrete next carrier with non-antipodal pairing;
the gadget forcing its complete relation remains missing.

## Prime families and clone families are different decompositions

The incoming [family census](../a000568_intransitivity_families_20261026/FAMILIES.md)
uses prime modular-decomposition quotients, automorphism cycle indices, and
ordered strong-component assembly. Our bank uses pair-free clone cores.
Order9 has136655 prime tournaments but147095 pair-free cores.

Our independent standard-library evaluator reads the23 incoming cycle-index
types through9, checks each group mass, and evaluates
`T=SEQ(S), S=X+P(T_+)` coefficient by coefficient. With the incoming order10
prime count7906161 it recovers A(10)=9733056 and S(10)=9355949. Every coefficient
through9 agrees with our complete saved bank. This verifies the species
evaluation **conditional on the supplied census inputs**; it does not rerun
nauty or supply an order10 representative bank. The incoming packet retains
its own census provenance and bounded status.

## Unsigned cyclic triples do not replace signed curl

The incoming intransitivity object C3(T) is the hypergraph of cyclic triples.
The exact forgetful map is `C3(T)={triples with |K_ijk|=3}`. It discards signs
and the signed transitive-triple values, so C3(T)=C3(T^op).

Even the **aligned** pair `(C3(T),R_T)` therefore loses global converse, unlike
`(K_T,R_T)`. A four-vertex cyclic triangle dominating a singleton and its full
converse are a nonisomorphic witness, with score multisets(0,2,2,2) and
(1,1,1,3). The checker verifies identical cyclic-triple sets and every named
flip response. These information quotients need different repairs.

The incoming THM-895 correction (unit-flow potential d/n) is retained too.
No scalar, geometry or finite-field coincidence is exported as a colouring
obstruction. All links preserve the incoming packets' own audit scopes.
