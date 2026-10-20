# Independent audit: hand-local-configuration supplement

Date: 2026-10-20. Verdict: **APPROVED**.

Audited source: `../session15-deficit-centers/HAND_LOCAL_CONFIGURATION_SUPPLEMENT.md`.
SHA-256: `acdde10d8f7d1eb76a3fd9bc45fd7b556f91b40f5a4981b68c2dee08b28984f9`.

The hand proof correctly replaces the finite six-vertex local lemma in the self-converse ambient application. It does not claim that an arbitrary strong non-self-converse ambient tournament must acquire an exit from the hand argument: the argument first gives an exit or entry, and uses the ambient anti-automorphism to turn an entry into an exit. That conversion is valid even if the image is a different induced subset.

Detailed checks:

1. Moon's theorem supplies a strong five-vertex core in every strong ambient tournament of order at least six, and supplies a four-cycle within that core. The strong four-vertex tournament is the stated S4.
2. If the fifth vertex has other than two outneighbors on S4, it gives the required entry or exit immediately. The six balanced extensions exhaust all remaining possibilities.
3. All four displayed cycle words consist of actual directed cycles, and the listed terminal is beaten by at least three of their vertices. The two remaining balanced cases have exactly the claimed structures: the regular C5, and Q5 with a=2, b=1, and directed triangle on 0,3,4.
4. Every four-subset of regular C5 is strong: after deleting one vertex, its internal scores are all 1 or 2. Any outside vertex has a majority of at least three in one direction on the five-core, and a four-subset containing that majority gives an entry or exit.
5. In Q5, any pair x,y of triangle vertices makes {a,b,x,y} strong: orient the pair x→y and use the cycle a→x→y→b→a. If all three such sets are balanced against the outside vertex, the equations p+t_x+t_y=2 force all three indicators equal. The two remaining cases are precisely those stated, and both displayed cycles and all three exits are correct.
6. An anti-automorphism reverses the induced cycle and reverses all three majority arcs, converting any entry into an exit. The subsequent connectivity/prefix argument remains unchanged.

An independent optional control script, `audit_hand_supplement.py`, verifies the explicit witnesses and exhaustively checks all labelled strong five-vertex cores. Of 544 strong inputs, 480 already have an entry or exit; the other 64 are exactly 24 labels of C5 and 40 labels of Q5. All 32 ways to adjoin a vertex to either core have an entry or exit. For each core, 31 extensions have an exit, 31 have an entry, and 30 have both. The two extreme extensions illustrate why the ambient self-converse conversion is genuinely needed in this proof.

These computations are hostile controls, not a replacement for the short exhaustive hand case argument. Exact counts are in `hand_supplement_checks.json`.

The accepted complete result D(z)=z^4 can therefore use this hand lemma instead of the local six-vertex enumeration, while still relying on the independently checked complete residual periodic-presentation census.
