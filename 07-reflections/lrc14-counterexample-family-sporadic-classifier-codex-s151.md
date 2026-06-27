# LRC14 Counterexample Family/Sporadic Classifier - Codex S151

S151 turns the current LRC14 counterexample hunt into a gate plus taxonomy:
a true row must be primitive, survive the q-clock exclusion, sit at `qdiv>14`,
and be exactly `covered` under the S146 threshold classifier.  Everything else
is not a counterexample; it is either a boundary atom or has positive Haar-open
witness mass.

The bounded audit covered AP, one AP swap through `add<=420`, two swaps through
`add<=60`, and three swaps through `add<=30`.  Across `68368` exact
`qdiv>=14` rows, the true candidate bucket was empty.  The two closed/tight
rows are AP and GW.  The rest split into positive q14 migration, covering-core
positive rows, unit-petal/GW strips, K33/state-lift packets, single-14-tail
combs, and unlabelled covering repairs.

The important conceptual move is that "sporadic" is no longer a vibe.  In this
interface it means `qdiv>14`, exact status `covered`, and no named family label.
The observed sporadic-like rows are only a reservoir: `16253` unlabelled
covering repairs appear, but all of them are positive-open in this audit.

Connections back through the repo:

- HYP-2955 supplies the qdiv/Haar migration bank and the fact that no hard
  covered row appeared there.
- HYP-2960 explains why AP/GW are the only q14 boundary atoms and why K33,
  petal, and Jacobsthal labels must travel as a packet.
- HYP-2947/HYP-2937 explain the low-frontier split into AP/GW, unit-petal, and
  K33/state-lift obligations.
- HYP-2953/HYP-2954 provide the source-spectrum/functor view that makes this
  taxonomy quotient-preserving rather than just a scalar census.

Assumption challenge: the tournament vertices were not runners.  I considered
runners, residues, divisor clocks, cover skeletons, Haar fronts, boundary
owners, C27/K33 packet labels, and proof families.  The chosen tournament uses
classifier gates/families as vertices because that is what preserves the LRC
predicate.  It destroys row identity after the divisor skeleton and packet
labels are recorded, so it cannot replace the exact interval computation.

Next useful theorem shape: prove every unlabelled distributed cover skeleton
has positive strict Haar mass by constructing a witness interval from its cover
skeleton.  That would collapse the current sporadic reservoir and leave only
the named family theorem plus K33/state-lift impossibility.
