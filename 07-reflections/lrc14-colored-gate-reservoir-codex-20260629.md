# LRC14 Colored Gate-Reservoir Reflection

The useful mental shift is that "coloring" is a compression operation, not a
decoration.  The question is what theorem-predicate survives the compression.

The old phase-color reservoir had colors `b mod 14` at the CRT grid.  HYP-3471
pulls that color down to the local survivor gate: each gate has endpoint
labels like `B1:7|E:84`, so the color word is typed residue
`B1:7|E:0`, plus sidecars saying which branch survives, which bad block is
adjacent, and how the minimal covers change across the gate.

The HYP-3462 AP84 corridor-splice certificate and HYP-3470 AP84 color-grid
bridge are siblings, not replacements.  HYP-3462 keeps the structural AP
corridor carrier, HYP-3470 keeps exact CRT placement for the canonical
`q=14V` AP tail, and HYP-3471 keeps the typed gate word after the obstruction
has become local.  In information terms, HYP-3462 preserves the AP splice,
HYP-3470 preserves grid witnesses, and HYP-3471 preserves graph-composable
boundary payload.

That is a satisfying little unification.  The AP84 transient packet is exactly
the four-color picture:

```text
B1:7|E:0
E:0|B1:5
B0:5|E:0
E:0|B0:7
```

But the bank refuses to be compressed to that packet.  Only `67/130` dead rows
contain the AP typed palette.  The good news is stronger and cleaner: every
dead row still has some rank-`<=2` E/branch gate.  The color changes; the
gate type does not disappear.

This also explains the explicit failure of "proper coloring" as a proof
metaphor.  There are many same-branch gates:

```text
B0|B0: 741
B1|B1: 741
B0|B1 / B1|B0: 182 total
```

So a naive three-color graph story is false.  What survives is a typed
compression ladder:

```text
endpoint kind only      -> 8 gate colors
numeric residue mod 14  -> 147 gate words
typed residue mod 14    -> 360 gate words
structural sidecar      -> 161 structural words
full colored gate word  -> 1727 words
```

That ladder has information-theoretic content: the row predicate stays alive at
the E/branch-gate level, but random rows need typed residues and delta sidecars
to prevent illegal identifications.  The AP packet is a low-entropy terminal
case, while `random_covering_031` is the high-entropy gluing stress.

The proof route I would try next is a current/cut proof of the strengthened
finite lemma:

```text
dead component => E/branch rank-2 gate.
```

In graph language, a dead component is a two-branch blocker island.  A survivor
gate is a boundary edge where the island leaks to an even wall.  If a full
two-colour cover claims there is no such E/branch leak, it should force a
closed current with impossible boundary divergence, or a Menger cut whose
terminal set is exactly the branch-coloured blocker support.  HYP-3471 supplies
the color alphabet for that cut: typed residues name the terminals, while
branch mask, adjacency, and deltas keep the compression honest.

This is not LRC14 closed, but it is a better target than "use colorings" in the
abstract.  Prove the E/branch gate theorem, use HYP-3462 for the closed AP84
corridor splice, and use HYP-3470 only where actual AP84 CRT witnesses are
needed; the color-reservoir story then has a theorem-shaped place to stand.
