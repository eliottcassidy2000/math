# LRC14 Gap Automaton Carrier - S173 Reflection

This pass tries to keep the user's finite-automaton and lacunary prompts from
collapsing into another scalar sequence shadow.

The useful object is not "fibbinary" or "Moser" alone.  It is a packet field:
which automaton language is being used, which state word the packet occupies,
which transition is native, and which LRC predicate survives after forgetting
the rest of the row.

## What Changed

S172 added the observer-frame ledger: boosts and Poincare language are safe
only after the observer velocity, tube metric, lattice, and sign debt are
retained.  S173 adds the same kind of ledger for gap languages:

```text
fibbinary no-adjacent language      -> native doubling x -> 2x
Moser-de Bruijn base-4 0/1 language -> native scaling x -> 4x
Hadamard lacunary support           -> natural-boundary warning
2-adic Littlewood/Hurwitz doubling  -> continued-fraction transition state
Fermat-Catalan exponent budget      -> finite-exception ledger
visibility core                     -> largest labelled source kernel
```

The finite audit makes one point hard to ignore: Moser-de Bruijn support embeds
inside fibbinary support through `N=4096`, and structurally this is because
base-4 digits `0,1` become binary pairs `00,01`.  But Moser has zero
violations under `x -> 4x` and `63` violations under `x -> 2x` through the same
bound.  Therefore "dyadic" is not specific enough.

## Concrete Output

The S173 script records:

```text
fibbinary_count=378
moser_count=65
intersection_count=65
fibbinary_mixed_residue_classes_mod14=14
moser_mixed_residue_classes_mod14=14
fibbinary_double_closure_violations=0
moser_times4_closure_violations=0
moser_double_closure_violations=63
```

On the unit-excess lane `q=14p-1`, there are `21` fibbinary/Moser hits with
`p<=384`.  These are not counterexample candidates.  They are route labels
that can be attached to exact `M/qdiv` packets.

The tournament signal is better than expected.  Conservative rank-priority
gauges are all transitive, but fieldwise majority produces one nontrivial SCC:

```text
{fibbinary_no_adjacent_language,
 moser_base4_digit_language,
 ostrowski_hadamard_lacunary_boundary}
```

That SCC is the missing warning: finite-state dyadic closure, base-4 sparse
support, and genuine lacunary natural-boundary force are three different
proof currencies.  A proof route can rank them only after choosing a retained
predicate.

## Assumption Challenge

Vertices considered:

```text
runners
speed gaps
automaton states
base-4 digit masks
Zeckendorf carry states
doubling transitions
valuation exponent triples
visibility cores
proof carriers
```

The chosen tournament vertices are proof/gap carriers.  This preserves the
question "which proof predicate survives the quotient?" and destroys direct
runner geometry.  That destruction is acceptable only because the proposed
packet fields reattach exact `M/qdiv`, endpoint, route, and certificate data.

## Proof Direction

The theorem shape I would now try is:

```text
Finite-State Residual Dichotomy.
Given a primitive non-AP/GW LRC14 packet after the HYP-2963 labels, augment it
with the S173 automaton/gap fields.  If every active sequence-shadow quotient
is finite-state and transition-compatible, then the packet must route to an
existing q, Fejer, Ramanujan, Haar, moment, twist, or THM-572 state-lift
certificate.  If not, the first failed transition defines the named F7
nonregular/lacunary residual.
```

This is not yet a proof, but it is a sharper fork.  It says F7 should not be a
bag of unknowns; it should be a failure of finite-state packet closure, a
natural-boundary obstruction, or a state-lift debt with exact labels.

## Next Pulls

1. Add the S173 automaton fields to the HYP-2963 packet bank or a sidecar
   manifest.
2. Build product automata for `q=14p-1` under fibbinary, Moser, Zeckendorf, and
   exact-period Ramanujan labels.
3. For each HYP-2963 hard non-AP/GW packet, compute the induced carrier
   tournament class and compare it to S173's `n=4,5,6` canonical class words.
4. Decide whether the fibbinary/Moser/Hadamard SCC is killed by an existing
   Fejer/Ramanujan/Haar certificate or must become a named F7 sub-sector.
