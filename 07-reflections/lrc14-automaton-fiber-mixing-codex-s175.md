# LRC14 Automaton Fiber-Mixing - S175

The prompt's finite-automaton lane has now become more precise.  Moser,
fibbinary, Hurwitz doubling states, Fermat-Catalan power guards, and
Ostrowski-Hadamard gap warnings are not proof scalars.  They are packet side
channels.  S175 asks what happens when we actually group LRC14 rows by those
side channels and compute exact safe measure.

The answer is sharp:

```text
residue + Moser/fibbinary terminal state mixes boundary and open rows.
```

In the exact single-swap atlas through tail `180`, AP and Goddyn-Wong are the
only zero-open rows.  Yet AP shares its `residue_mfc_pairs` and
`residue_terminal_pairs` fibers with open rows such as `12->26` and `12->96`.
Goddyn-Wong shares its one-dipole fiber with open tails `12->38`, `12->52`,
and `12->150`.

That makes the missing coordinate visible.  It is not another bit in the
fibbinary DFA.  It is magnitude: Farey scale, tail height, endpoint ownership,
and safe-component geometry inside a fixed residue-language fiber.  I called
the lost coordinate the `magnitude_cocycle` because this is exactly the
HYP-2990/HYP-3002 rule in a new costume: if a quotient forgets a coordinate, it
must prove fiber constancy, reconstruct the coordinate, annihilate it by a
dual certificate, descend to a family theorem, or route it to named residual
debt.

The outside ideas land cleanly now:

```text
2-adic Littlewood / Hurwitz:
  keep native transition state, but do not expect terminal DFA state to decide
  LRC boundary-vs-open behavior.

Moser-de Bruijn:
  a base-4 no-carry label, good for subfibers, unsafe as a tightness scalar.

Fibbinary / Zeckendorf:
  a better dyadic carry normal form, still not enough without magnitude.

Ostrowski-Hadamard:
  a warning that lacunary supports need boundary labels, not smoothing away.

Fermat-Catalan / perfect powers:
  a finite-exception/no-lift ledger; S175's `perfect_power_word` also mixes.

Sticks/potatoes:
  approximate visibility is a guardrail; exact safe components are the proof
  object.
```

The most useful theorem target from this pass:

```text
Within each residue-automaton fiber, every nonzero magnitude cocycle either
opens a strict safe component, descends to a proved single-swap formula, or
emits K33/F7/THM-572 residual debt.
```

That is close to the repo's current proof grammar and much less misty than
"finite automata might help."  They help by telling us exactly which coordinate
they are not allowed to forget.
