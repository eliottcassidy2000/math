# Codex Synthesis: AP/GW As Base And First-Derived Tournament Atoms

- Created: 2026-06-24T08:40:20Z
- Agent: codex-synthesis
- Post: 20260624-084020Z-derived-tournament-ap-gw-conditions

## Session Meat

The AP/Goddyn-Wong common core now has a sharper research shape:

```text
AP = base endpoint identity.
GW = first admissible nilpotent derivative of that identity.
```

The Tao-derived-multiplicative-functions analogy is not decorative.  It gives
the right type signature for an AP perturbation.  A first perturbation should
obey a Leibniz rule across coprime clock factors, and a second perturbation
should either vanish or produce a new obstruction ledger.  For `14 = 2 * 7`,
the only first perturbation that survives the Jacobsthal gate is `12 -> 24`.

The Collatz/Farey analogy also earned its keep.  A nontrivial Collatz cycle is
blocked by rational approximation data for `log_2(3)`.  An LRC14 endpoint atom
is blocked by rational approximation data around `1/14`; the familiar
`3/41` resonance is the first place to audit.

The new tournament technique is the part to operationalize.  Four tournaments
were defined:

```text
T_press  = runner endpoint-pressure tournament.
T_gap    = residue multiplicity/gap tournament.
T_gate   = acceleration-site Jacobsthal tournament.
T_square = node-squared local tournament, each runner carrying a 13-node
           internal tournament.
```

The summit target is a finite set `A_14` of achievable isomorphism signatures:

```text
AP class,
GW class,
plus unit/reflection/renaming symmetries.
```

If every primitive endpoint candidate maps into `A_14`, the LRC14 tight locus
collapses to AP/GW after the existing quotient/divisibility ledgers are paid.

## Random Repo Niche

The strongest new "fun" condition is `T_square`.  It takes the user's
node-squaring prompt literally:

```text
each tournament node becomes a tournament of the original size.
```

The first computation found a useful warning: if `T_square` keeps only inner
score sequences, then AP, GW, and all single folds look transitive and
indistinguishable.  The full relative profile

```text
Counter_u(distance(xu), distance(yu), distance((y-x)u))
```

is the keeper.  AP has a `6 | 1 | 6` local-profile histogram.  Single even
folds, including GW, move into a first-fold profile class.  The Jacobsthal gate
then selects `12 -> 24` from that class.  A candidate with two independent
folds should show a larger local-profile defect and should be routed to
support-6/Kuratowski debt rather than AP/GW.

## Concrete Next Computation

Build `scripts/lrc14_tournament_atoms.py` with:

```text
input: residue multiplicity vector m_0,...,m_13 plus divisibility ledger
output:
  pressure profile blocks,
  residue-gap tournament fingerprint,
  acceleration-gate source list,
  node-squared inner-class histogram,
canonical tournament hash/isomorphism class.
```

Prototype now exists as:

```text
04-computation/lrc14_derived_tournament_atoms_codex_s145.py
```

Start with the finite family:

```text
AP,
all single accelerations AP \ {v} union {2v},
all double accelerations,
the AP-residue false look-alike {1,...,11,13,26},
known support-6 wall examples.
```

Expected first result:

```text
only v = 12 has the admissible T_gate source fingerprint;
the AP-residue false look-alike fails covering_ok;
coarse pressure shells and inner score sequences are too weak;
node-squared relative-profile hashes preserve the single-fold defect class;
double accelerations still need the next support-6/Farey audit.
```

That would make the tournament method genuinely useful instead of just poetic:
it would explain which finite isomorphism classes are achievable and which
pieces of continuous arithmetic each class forgot.

## Connections

Connect this post to:

- `07-reflections/anatomy-of-a-tight-runner-set.md`: original AP/GW eight-layer
  anatomy.
- HYP-2917/HYP-2918/HYP-2919: covering quotient, doubling operad,
  Jacobsthal/Farey constraints.
- HYP-2946: Kuratowski/perfect-carrier mediant guardrails.
- HYP-2947: measurable rank recombination and support-6 address debt.
- `20260624-083010Z-haar-baire-any-angle-lrc14`: Borel/Baire/Haar endpoint
  owner route.
