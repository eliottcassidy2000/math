# LRC14 Colored Gate Partition-Lattice Reflection

codex-2026-06-29

## What Changed

HYP-3471 made the colored gate-reservoir concrete:

```text
dead_components(row) > 0
  => rank <= 2 E/branch survivor gate
```

on the current `135`-row bank, with `130/130` dead rows satisfying the
strengthened E/branch gate condition.  HYP-3474 asks the next question: after
we have a colored gate, what is a quotient allowed to forget?

The answer is pleasantly severe.  The theorem-failure predicate is identically
false on the bank, so every quotient preserves the implication as a sampled
fact.  But most quotients do not preserve the route information needed to turn
the sampled fact into a proof.  The right target is a partition-lattice
statement:

```text
quotient Q is legal for target P
  iff every fiber of Q is P-pure.
```

This reframes the proof problem as a labelled packet theorem rather than a
search for one blessed scalar.

## Main Exact Lessons

The singleton axis audit gives the first guardrail:

```text
K endpoint-kind set:
  8 fibers, mixes route flags across 132 rows

S structural sidecar set:
  107 fibers, still mixes route flags across 32 rows

N numeric mod-14 set / T typed mod-14 set / F full colored-gate set:
  123 fibers, route-pure on the current bank
  but each still mixes 12 rows for the dead minimum structural gate family

C count profile:
  122 fibers, route-pure on the current bank
  but mixes 15 rows for the dead minimum structural gate family

M minimum E/branch structural gate:
  17 fibers, mathematically meaningful but route-mixed across 124 rows
```

The strongest new finite lead is:

```text
(C,M), (N,M), (T,M), and (F,M)
```

are pure for both the dead minimum structural family and the dead minimum full
family, with `125` fibers and max fiber size `7`.  The minimum structural gate
is not enough by itself, but it becomes a legal obstruction-family coordinate
once paired with a high-cardinality row shadow.

## The Count-Profile Trap

The most dangerous finding is also the most interesting one: the low-rank gate
count profile `C` is route-pure on the finite bank and wins the tournament
score tie.

That does not make count a proof quotient.  The duplicate-fiber report shows
why it behaves so well here.  The large duplicate classes are AP-tail copies:

```text
covering_AP_with_84,
ap_omit_12_tail_84x01,
...,
ap_omit_12_tail_84x12
```

with the same route flags and minimum structural gate.  The only random
duplicate in the displayed count fibers is the nondead pair
`random_covering_044/random_covering_053`.

So `C` is a lead for a hidden gate-count recursion or AP-tail normal form, but
it is not yet an admissible proof quotient.  To use it, we would need a
reconstruction theorem:

```text
count profile -> row shadow / min structural gate / route flags
```

or a dual certificate showing the lost interval geometry and owner-current
coordinates cannot affect the conclusion.

## Proof Reframe

The current proof frontier should be phrased as:

```text
dead cover component
  -> E/branch survivor gate
  -> min structural gate M
  -> row shadow C/N/T/F
  -> route flags and named debts
```

This looks like a Myhill-Nerode theorem for proof states.  Two rows may be
identified only when every future proof continuation sees the same route
predicate.  That ties back cleanly to the tournament-tiling program: a
Hamiltonian path or fixed color order is only legal when the forgotten edge
orientation belongs to a pure fiber for the chosen target.

## Connections To Older Threads

Irreducibility: evaluation of a polynomial at one integer is a quotient from
coefficient/factor geometry to an integer.  Singh-style factor-count witnesses
are legal only for factor-count targets; they are not automatically legal for
coefficient-lift irreducibility.  HYP-3474 is the same warning in LRC clothing.

Unital designs and C27 labels: a design quotient forgets point identity but
keeps block incidence.  HYP-3474 says the LRC quotient must say which block
incidences, gate colors, and sidecars survive.

Faulhaber moments: odd moments suffice only because the balance equation
annihilates the even moment directions.  Here a quotient may forget an axis
only after the target predicate annihilates that direction or after a sidecar
restores it.

Colored discrepancy / Haar zipper: scalar discrepancy is too blunt unless the
color fiber and boundary cocycle are retained.  The AP84 bit `A` is the clean
example: it perfectly preserves AP-palette presence, but it mixes every other
route target.

Tournament metagraphs: the vertices in HYP-3474 are not runners or arcs.  They
are quotient carriers.  The tournament orders quotient candidates by how much
proof-predicate purity they preserve, with scalar compression as a secondary
bonus.

## Next Session Hooks

1. Prove the strengthened HYP-3471 finite lemma structurally:

```text
dead_components(row)>0 => rank<=2 E/branch survivor gate.
```

2. Try to prove why AP-tail duplicate fibers are count-profile pure.  If this
is just the HYP-3456/HYP-3458 mod-`35` recursion in another costume, it could
become a compact AP normal-form lemma.

3. Generalize the partition-lattice script to new row banks.  If `C` stops
being route-pure outside the current bank, keep it as telemetry.  If it
continues to be route-pure, look for a count-profile reconstruction theorem.

4. Build the explicit proof-state automaton:

```text
state = (route flags, min structural gate, AP packet, gluing debt)
transition = gate deletion / AP splice / owner-current discharge
```

The desired theorem is that every dead row reaches an accepting E/branch state
before any illegal quotienting step occurs.
