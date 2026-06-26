# LRC14 Lacunary Exact-Gap Automaton Carrier

codex-2026-06-26-S170

Rebase note: this reflection supports HYP-3010/T1094 as the exact LRC
maximin/safe-component companion to the HYP-3008/LTI-158 automatic gap-language
carrier and HYP-3009/LTI-159 Fermat-Catalan power-lift extension.

## What Changed

The prompt suggested Fermat-Catalan, a 2-adic Littlewood paper, the
Ostrowski-Hadamard gap theorem, Moser-de Bruijn, fibbinary numbers, and finite
automata.  The useful LRC14 interpretation is not a new family of tight rows.
It is a new side-channel rule:

```text
finite automaton carry data is useful only when exact M=p/q,
e=14p-q, residue ownership, and exact safe components are retained.
```

The script `04-computation/lrc14_lacunary_automaton_carrier_codex_s170.py`
does exact LRC14 computations on fibbinary and Moser-de Bruijn rows.  AP and GW
remain the only boundary atoms in the scout.  The automaton rows are strict:

```text
first13 fibbinary             M=3/25
first13 Moser-de Bruijn       M=1/6
fibbinary residue transversal M=1/6
Moser residue transversal     M=1/6
```

So the automata are not candidate counterexamples and not a proof by
themselves.  They are labelled carry languages for a packet classifier.

## Why The Result Is Still Useful

Moser-de Bruijn is a clean negative control.  It is no-carry and finite-state:
base-4 digits are only `0` or `1`.  It is also too lacunary in magnitude to
stay near the LRC14 boundary.  This separates two ideas that are easy to
confuse:

```text
no carry        helpful local state
lacunary scale  dangerous quotient if endpoint geometry is forgotten
```

Fibbinary is less sparse and matches the Zeckendorf/no-adjacent-carry thread.
It is the better model for local carry constraints, but the exact row audit
still says it must be attached to Farey scale and safe components.

## Assumption Challenge

Do not use runners or arcs as the default tournament vertices here.

Alternate vertices considered:

```text
runners, gaps, fixed circle sections, section boundaries, wall-crossing events,
residues, cover arcs, Fourier modes, matroid circuits, automaton states,
proof obligations.
```

The chosen vertices are proof carriers.  This preserves the LRC predicate:
"is there an exact safe time at threshold `1/14`?"  Raw sequence membership
destroys exact scale, endpoint owners, and the wall geometry.

The carrier tournament has a nontrivial middle SCC tying Moser, Fermat-Catalan,
and Ostrowski-Hadamard cues.  That is the correct read: these are guardrails
and obstruction languages, not a scalar hierarchy.

## Next Pull

Add two packet fields to the HYP-2963 / HYP-3001 infrastructure:

```text
carry_language:
  none | fibbinary | moser_base4 | hurwitz_dyadic | ostrowski_cf

automaton_state:
  finite state label at the active Farey/wall packet
```

Then rerun the HYP-2996 residual-section bank and ask whether every zero-open
non-AP/GW residual either emits a nontrivial dyadic/Ostrowski carry state,
has an owner-strip/cross-handoff/nested-refinement exit, routes to the
K33/THM-572 lift, or collapses to the AP/GW boundary skeleton.
