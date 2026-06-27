# LRC14 Tope-Wall / Boundary-Cocircuit Pass - Codex S164

This pass tried a different proof language: not Fourier, not moment duals, not
large-sieve packets.  The time circle is cut by all danger endpoints, and the
open cells are treated as topes.

The useful trichotomy is:

```text
open all-safe tope      -> strict lonely interval
boundary cocircuit     -> closed equality witness
forbidden wall packet  -> neither of the above
```

AP and Goddyn-Wong move into the second category.  That is the main conceptual
gain.  They are not failures of the tope method; they are zero-dimensional
cocircuit atoms.

## Computation

The script:

```text
04-computation/lrc14_tope_wall_cocircuit_codex_s164.py
```

audits named rows and a one-swap AP hard bank through `add <= 140`.

The named rows behave exactly as the current proof stack predicts:

```text
AP/GW: min open danger = 1, safe mass = 0, six boundary cocircuits.
K33/petals/splices/covering rows: min open danger = 0, positive safe mass.
```

The one-swap hard bank has:

```text
boundary_cocircuit  2
open_tope           853
forbidden_wall      0
```

The two boundary rows are AP and GW.

## Proof Shape

The proposed theorem target is now very crisp:

```text
No primitive LRC14 row can be a no-tope/no-cocircuit wall packet unless it
constructs the known K33/H=7 state-lift obstruction or a genuinely new F7
packet family.
```

This complements the endpoint-credit dual route.  HYP-2970 studies cover
cycles.  S164 studies the minimum-danger sign walk before scalarizing it.

After rebasing over the neighboring S164 kernel-homotopy and Farey-mutation
scheduler passes, the clean handoff is: Farey excess chooses a certificate
queue, kernel homotopy tracks deformations and named defects, and this
tope-wall pass decides which wall state the packet actually occupies.  The
three states are open witness, boundary equality atom, and forbidden wall
packet.

## Why It Might Help

The scalar Haar question has one number: positive or zero.  The tope-wall
language keeps the exact reason for zero:

```text
zero because AP/GW has boundary cocircuits,
or zero because every cell is still dangerous and no boundary point is safe.
```

Only the second possibility is a strict-counterexample shape.  That is a much
smaller target.

The next computation should extend from one-swaps to HYP-2963 packet families
and emit a compact wall signature:

```text
min danger,
minimum-cell components,
boundary cocircuit owners,
owner pair sums mod 14,
wall-crossing owner changes,
K33/state-lift flag.
```

If no packet family produces the forbidden wall signature, this can become a
finite classifier lemma feeding the existing source-kernel theorem.
