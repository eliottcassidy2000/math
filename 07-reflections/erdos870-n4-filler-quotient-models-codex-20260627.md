# Erdos-870 And n=4 Filler Quotient Models

**Session:** codex-2026-06-27
**Hypothesis:** HYP-3145
**Artifacts:** `04-computation/lrc14_erdos870_n4_filler_models_codex_20260627.py`, `05-knowledge/results/lrc14_erdos870_n4_filler_models_codex_20260627.out`

The two n=4 tables separate two kinds of modeling.

The fixed-path tiling table is a representative calculus.  It is useful
because it exposes the three off-path toggles, but the quotient to
isomorphism classes is not congruent.  The `S` class is the tell: it contains
`c`, `ab`, `ac`, `bc`, and `abc`, so multiplying from an `S` representative
can land in different classes.  That is exactly the controlled-forgetting
warning the LRC14 work keeps rediscovering: a visually clean table may still
hide a destroyed coordinate.

The partial-score table is different.  Four fixed arcs force the score
scaffold `(0,1,1,2)` and leave two arcs `x=(0,1)` and `y=(2,3)`.  The four
states `E,x,y,xy` hit `T,+,-,S` exactly and close as a Klein square.  This is
closer to a proof interface than the tiling cube because the quotient is
legal before the table is used.

The Erdos-870 connection is the proof architecture: a small core plus finite
fillers, with deletion/nonminimality controlled at the interface.  For LRC14,
this suggests a sharper next move.  Do not ask the fixed-path tiling cube to
be the algebra.  Use it as an atlas and ambiguity detector, then choose
finite fillers that force a residue/score scaffold and leave a smaller signed
core for HYP-3129/HYP-3132/HYP-3136.  That is the bridge from the user's table
to the current edge-witness frontier.

The practical packet field is `filler_core_interface`: it records the fixed
filler scaffold, the remaining core variables, whether the quotient is a true
congruence, and which deleted coordinate becomes named debt if it is not.
This field should sit next to HYP-3134's global-consistency class and
HYP-3135's middle resolvent/SPEC payload.
