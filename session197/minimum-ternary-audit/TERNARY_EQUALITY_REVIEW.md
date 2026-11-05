# Independent audit of amplitude-60 equality rigidity

Result: PASS, with no change to the frozen source packet.

Inputs and SHA-256:

- `session197-prime-atoms/TERNARY_EQUALITY.md`: `09d40d8ae91a183cfa66452184d9ee110562b544295d2506936745ea0624fc2f`
- `session197-prime-atoms/verify_ternary_equality.py`: `e5f67960d210c90b4953ce4122f29f7913d714937e63db78805fa24e8a5f456f`
- `session197-prime-atoms/verify_ternary_equality.cpp`: `b1c0d43c08cdeb5cd0fb8936f7778a6670769746c99b1385e04a0323bc347b6d`
- `session197-prime-atoms/receipt-ternary-equality.json`: `0a5cb8b13df5e64e3c53bee85d481f901f46b1b5010a7ce98a19d613c59eaec2`

The all-order reduction correctly reuses the independently reviewed minimum-amplitude proof. At amplitude 60 the exact strong core has H179, outside SCC factors are one, and the exact Busch comparison still bounds core order by ten. All candidate roots must have matching number two, so the two exact clone equations suffice, including for an eventual rather than initial identity.

The cap-179 generator uses the same safe complete strong-deletion construction and exact canonicalization as the reviewed minimum packet. I checked that it includes all rooted candidates of exact count 179, not only those previously found at the lower cap. My complete normal and Python-optimized replays match the supplied receipt byte-for-byte: 1,595 strong graphs, 152 candidate roots, four survivors of the first equation and exactly two survivors of both equations. As an additional inventory check, filtering the separate complete uncapped order-seven bank gives exactly the same 352 capped strong representatives.

The root-preserving converse permutation, isomorphism from the original hybrid clock, and full automorphism checks are correct. Both representatives are canonical and different, which establishes their nonisomorphism. Exactly one tested root in each graph survives, and the matching theorem excludes any untested higher-rank alternative root.

The arbitrary-order count formulas follow rigorously: the SCC chain consists of a transitive prefix, one of the two rigid nontrivial cores, and a transitive suffix. Different prefix lengths cannot be isomorphic because SCC order is invariant. Each resulting source is rigid and has a unique clock root. Thus there are 2(n-7) unlabelled and 2(n-7)n! labelled sources at every n>=8, and none below eight. This classifies only the amplitude-60 transitive-clone stratum, as stated.

All frozen hashes remained unchanged. No correction or unresolved blocker was found.
