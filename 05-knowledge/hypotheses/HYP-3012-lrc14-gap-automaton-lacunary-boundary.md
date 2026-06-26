---
id: HYP-3012
title: LRC14 gap automaton carrier and lacunary-boundary guardrail
status: SYNTHESIS / necessary-condition and proof-interface carrier; not a proof
source: codex-2026-06-25-S173
script: 04-computation/lrc14_gap_automaton_carrier_codex_s173.py
result: 05-knowledge/results/lrc14_gap_automaton_carrier_codex_s173.out
related:
  - HYP-3011
  - HYP-3010
  - HYP-3009
  - HYP-3008
  - HYP-3007
  - HYP-3006
  - HYP-3004
  - HYP-3003
  - HYP-3002
  - HYP-2998
  - HYP-2997
  - HYP-2983
  - HYP-2982
  - HYP-2963
  - THM-572
  - OPEN-Q-108
---

# HYP-3012: LRC14 Gap Automaton Carrier

## Claim

This is the S173 extension of the incoming HYP-3008 automatic-gap scout, the
HYP-3009 Fermat-Catalan automatic-gap / power-lift ledger, the HYP-3010
exact-gap/safe-component carrier, and the HYP-3011 automatic lacunary
safe-component filter.  It keeps the
Moser/fibbinary automaton base but adds lacunary-boundary, exponent-budget,
visibility-core, and induced tournament isomorphism-class ledgers.

The prompt's gap-language family gives a new necessary condition for any
sequence-shadow proof of LRC14:

```text
If a proof quotient uses fibbinary, Moser-de Bruijn, lacunary, 2-adic,
Fermat-Catalan, Skolem-Mahler-Lech, or visibility-core data, then the quotient
must retain the automaton state / gap support / valuation budget that makes the
data meaningful, or prove that the LRC14 predicate is constant on the forgotten
fiber.
```

The finite-state lesson is sharp.  Moser-de Bruijn support embeds inside the
fibbinary language because base-4 digits `0,1` become binary pairs `00,01`.
But the native transition is different:

```text
fibbinary:  x -> 2x closes by appending binary 0
Moser:      x -> 4x closes by appending base-4 0
Moser:      x -> 2x fails generically
```

So a dyadic or 2-adic LRC route must declare which transition it is using.  A
raw statement such as "sparse sequence membership" is too lossy.

## Computation

Script:

```text
04-computation/lrc14_gap_automaton_carrier_codex_s173.py
```

Output:

```text
05-knowledge/results/lrc14_gap_automaton_carrier_codex_s173.out
```

Finite audit through `N=4096`:

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

The residue readout is the guardrail: both languages meet all residue classes
modulo `14`, and every residue class is mixed with nonmembers.  Therefore an
LRC14 proof cannot use residue plus language membership as a certificate unless
the exact packet data are still present.

On the LRC14 unit-excess chain `q=14p-1`, the script finds `21` hits with
`p<=384` in the fibbinary/Moser languages.  Examples:

```text
p=3   q=41    fibbinary
p=5   q=69    fibbinary+moser
p=23  q=321   fibbinary+moser
p=79  q=1105  fibbinary+moser
p=93  q=1301  fibbinary+moser
```

These are not counterexamples.  They are prompts for packet fields: the same
unit-excess Farey lane can intersect several finite-state languages, so the
proof object must remember exact `M/qdiv`, endpoint ownership, automaton state,
and doubling/base-4 transition state.

## Tournament Analysis

Vertices are proof/gap carriers, not runners:

```text
labelled_lrc_gap_automaton_packet
two_adic_littlewood_hurwitz_doubling
fibbinary_no_adjacent_language
moser_base4_digit_language
ostrowski_hadamard_lacunary_boundary
fermat_catalan_exponent_budget
skolem_mahler_lech_zero_set_guard
large_stick_potato_visibility_core
raw_sequence_membership_scalar
```

Pairwise observable:

```text
finite_state
gap_support
doubling_transfer
valuation_budget
lrc_packet_transfer
boundary_guardrail
quotient_safety
proof_maturity
```

The rank-priority gauges are all transitive.  The proof-priority path is:

```text
labelled_lrc_gap_automaton_packet
> two_adic_littlewood_hurwitz_doubling
> moser_base4_digit_language
> fibbinary_no_adjacent_language
> large_stick_potato_visibility_core
> ostrowski_hadamard_lacunary_boundary
> skolem_mahler_lech_zero_set_guard
> fermat_catalan_exponent_budget
> raw_sequence_membership_scalar
```

The fieldwise-majority gauge is not transitive:

```text
score_hist={0:1,1:1,2:1,3:1,5:3,7:1,8:1}
directed_3cycles=1
scc_sizes=[3,1,1,1,1,1,1]
nontrivial_scc={fibbinary, Moser, Ostrowski-Hadamard lacunary boundary}
hamiltonian_path_count=3
```

This is the important tournament signal.  Fibbinary wins by finite-state
dyadic closure; Moser wins by base-4 sparse support; Ostrowski-Hadamard wins by
true lacunary boundary force.  None dominates the other two under all retained
dimensions.  A proof that collapses them to a scalar "gap sequence" is hiding a
3-cycle.

The induced tournament isomorphism-class census on this carrier surface gives:

```text
n=4: distinct_classes_generated=3
n=5: distinct_classes_generated=4
n=6: distinct_classes_generated=4
```

Representative canonical words include:

```text
n=4: 000000, 000010, 001001
n=5: 0000000000, 0000000010, 0000001001, 0001000011
n=6: 000000000000000, 000000000000010, 000000000001001, 000000001000011
```

This is not a universal enumeration of tournament classes.  It is an exact
achievable-subset census for the present proof-carrier surface, and future
agents can compare new carriers against these words.

## Source Synthesis

- `arXiv:2506.04110` contributes the 2-adic Littlewood/Hurwitz idea that
  doubling of continued-fraction data is itself a retained transition state.
- Ostrowski-Hadamard contributes the warning that truly lacunary supports can
  force natural boundaries; for LRC14 this is boundary-residual language, not a
  scalar proof.
- Fermat-Catalan contributes the finite-exception pattern
  `1/p + 1/q + 1/r < 1`; the S173 audit sees `150` hyperbolic triples in
  `2<=p<=q<=r<=10`, with first examples
  `(3,3,4)`, `(2,4,5)`, `(2,3,7)`.
- DOI `10.5555/1109557.1109610` resolves to a geometry paper on large sticks
  and potatoes in polygons, so it is used here only as a visibility-core
  analogy: the largest inscribed/visible structure is a labelled core, not a
  number-theoretic theorem.

## Necessary Conditions

1. **Automaton-fiber condition.**  Any sequence/gap quotient must retain
   `automaton_language_id` and `automaton_state_word`, or prove LRC14 predicate
   constancy on automaton-language fibers.

2. **Transition-state condition.**  A doubling route must say whether it uses
   `x -> 2x`, `x -> 4x`, or a recentered 2-adic/Hurwitz transition.

3. **Lacunary-boundary condition.**  A Hadamard-gap claim is a boundary
   obstruction until it is attached to a finite labelled LRC source-packet
   family.

4. **Exponent-budget condition.**  Fermat-Catalan style exponent budgets are
   finite-exception ledgers over packet fibers.  They cannot replace exact
   `M/qdiv`, endpoint, route, or certificate labels.

5. **Visibility-core condition.**  "Largest core" analogies may guide a
   labelled source-kernel search, but they must carry the geometry they used:
   visibility graph, core boundary, and packet route.

6. **Exact safe-component condition.**  HYP-3010's exact maximin and
   safe-component labels are the negative-control side channel.  A sequence
   language can be a strict positive row while AP/GW remain zero-open boundary
   atoms, so the automaton state must be read beside exact `M`, safe mass,
   endpoint ownership, and component geometry.

7. **Induced class condition.**  If a quotient replaces packet data by a
   tournament isomorphism class, it must name the vertex set, pairwise
   observable, tie rule, and induced class word, then prove that the LRC14
   predicate is constant on that class or identify the first class collision.

## Packet Schema Additions

Add these fields to HYP-2963-style records or to a sidecar manifest:

```text
automaton_language_id
automaton_state_word
gap_support_ratio_label
hadamard_boundary_warning
doubling_transition_state
base4_digit_mask
zeckendorf_carry_state
valuation_exponent_budget
finite_exception_budget
visibility_core_label
safe_component_label
induced_tournament_class_word
```

## Candidate Theorem Target

A primitive non-AP/GW LRC14 residual whose exact packet data are finite-state
under all declared carry/doubling languages must either:

```text
emit a known q / Fejer / Ramanujan / Haar / state-lift certificate,
or expose a named nonregular or natural-boundary residual sector.
```

This reframes F7.  It is not "whatever remains"; it is a specific failure of
finite-state labelled packet closure, or a specific nonregular/lacunary
boundary object.
