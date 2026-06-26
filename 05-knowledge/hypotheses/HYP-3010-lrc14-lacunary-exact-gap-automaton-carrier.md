---
id: HYP-3010
title: LRC14 lacunary exact-gap automaton carrier
status: EVIDENCE / proof-interface carrier; not a proof
source: codex-2026-06-26-S170
tangent: T1094
script: 04-computation/lrc14_lacunary_automaton_carrier_codex_s170.py
result: 05-knowledge/results/lrc14_lacunary_automaton_carrier_codex_s170.out
related:
  - HYP-2963
  - HYP-2981
  - HYP-2996
  - HYP-2998
  - HYP-2999
  - HYP-3000
  - HYP-3001
  - HYP-3008
  - HYP-3009
  - LTI-089
  - LTI-149
  - LTI-150
  - LTI-158
  - LTI-159
  - OPEN-Q-108
---

# HYP-3010: LRC14 Lacunary Exact-Gap Automaton Carrier

## Claim

Fibbinary and Moser-de Bruijn languages are useful for LRC14 only as labelled
carry side channels.  They should not be scalarized into "sequence rows" or
"lacunary rows" without retaining exact LRC data.

The admissible carrier is:

```text
exact M=p/q and e=14p-q
residue ownership mod 14
strict safe-component geometry at threshold 1/14
finite automaton state / carry language
```

The prompt's Fermat-Catalan and Ostrowski-Hadamard cues are guardrails.  Power
collisions should be treated as scarce primitive events with side conditions,
not as dense representation fibers.  Hadamard-lacunary supports warn that
forgetting packet labels can create natural-boundary behavior rather than a
smooth analytic continuation route.

## Computation

Script:

```text
04-computation/lrc14_lacunary_automaton_carrier_codex_s170.py
```

Stored output:

```text
05-knowledge/results/lrc14_lacunary_automaton_carrier_codex_s170.out
```

The script builds two small automata:

```text
fibbinary: no adjacent 1s in binary
Moser-de Bruijn: base-4 digits in {0,1}
```

Then it computes exact LRC14 gaps and exact safe components for AP, GW, the
first 13 terms of each language, and the first mod-14 residue transversal in
each language.

Key readout:

```text
AP boundary                         M=1/14   safe_mass=0
GW boundary 12->24                  M=1/14   safe_mass=0
first13 fibbinary                   M=3/25   safe_mass=66077/399840
first13 Moser                       M=1/6    safe_mass=4264747/40348854
fibbinary residue transversal       M=1/6    safe_mass=8317919/55327860
Moser residue transversal           M=1/6    safe_mass=4264747/40348854
```

The exact result is a negative-control result: these automata restore useful
carry and residue labels, but the tested automaton rows are strict positive
rows, not new tight atoms.  Moser is especially useful as a warning because it
is finite-state and no-carry yet too lacunary in magnitude to preserve the
LRC14 boundary.

## Tournament Analysis

Vertices are proof carriers, not runners:

```text
farey_exact_M_scheduler
large_stick_safe_component
two_adic_hurwitz_transducer
fibbinary_no_adjacent_carry
moser_de_bruijn_base4_sparse
fermat_catalan_power_collision
ostrowski_hadamard_gap_guard
raw_sequence_scalar
```

Pairwise observable:

```text
(exact scale retention,
 LRC predicate preservation,
 finite-state checkability,
 carry-collision control,
 quotient-loss warning)
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,4:2,5:2,7:1}
directed_3cycles=2
scc_sizes=[1,1,1,1,4]
hamiltonian_path_count=5
```

One high-score Hamiltonian path:

```text
farey_exact_M_scheduler
> large_stick_safe_component
> two_adic_hurwitz_transducer
> fibbinary_no_adjacent_carry
> moser_de_bruijn_base4_sparse
> fermat_catalan_power_collision
> ostrowski_hadamard_gap_guard
> raw_sequence_scalar
```

The middle SCC is the important warning: Moser, Fermat-Catalan, and
Ostrowski-Hadamard cues all carry useful obstruction information, but none is
safe as an LRC quotient unless the exact scale and endpoint geometry are
reattached.

## Source Links

The 2-adic Littlewood paper states the bounded-partial-quotient problem for
all dyadic multiples and uses Hurwitz multiplication by 2; it also constructs
long chains of equivalent quadratic irrationals.  This is the finite-transducer
analogy for LRC14's dyadic carry lane:

```text
https://arxiv.org/abs/2506.04110
```

The Ostrowski-Hadamard gap theorem is used only as a lacunarity guardrail: a
power series with exponents separated by a uniform ratio greater than 1 has no
regular boundary point at the unit circle:

```text
https://en.wikipedia.org/wiki/Ostrowski%E2%80%93Hadamard_gap_theorem
```

The ACM DOI resolves to the SODA paper on computing large sticks/potatoes in
polygons.  For this LRC pass it contributes the exact "largest safe component"
analogy, not a Fermat-Catalan theorem:

```text
https://doi.org/10.1145/1109557.1109610
```

## Proof Target

Add two fields to the HYP-2963 / HYP-3001 packet records:

```text
carry_language:
  none | fibbinary | moser_base4 | hurwitz_dyadic | ostrowski_cf

automaton_state:
  finite state label at the active Farey/wall packet
```

This is the exact-gap companion to HYP-3008/LTI-158 and HYP-3009/LTI-159:
those notes verify the automatic-language and Fermat-Catalan/power-lift packet
fields, while this note attaches those languages to exact LRC maximin gaps and
strict safe-component geometry.

Then test the finite HYP-2996 residual-section packet bank:

```text
Every zero-open non-AP/GW residual either
  emits a nontrivial dyadic/Ostrowski carry state,
  has an owner-strip/cross-handoff/nested-refinement Haar exit,
  routes to K33/THM-572 state lift,
  or collapses to the AP/GW boundary skeleton.
```

This would turn the automaton idea into a packet-preserving zipper tooth rather
than a new scalar sequence shadow.
