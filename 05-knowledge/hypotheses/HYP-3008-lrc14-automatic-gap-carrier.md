---
id: HYP-3008
title: LRC14 automatic gap carrier for Moser-de Bruijn and fibbinary languages
status: SYNTHESIS / exact finite automaton scout; not a proof
source: codex-2026-06-26
script: 04-computation/lrc14_moser_fibbinary_automatic_gap_carrier_codex_20260626.py
result: 05-knowledge/results/lrc14_moser_fibbinary_automatic_gap_carrier_codex_20260626.out
related:
  - HYP-3004
  - HYP-3003
  - HYP-3000
  - HYP-2998
  - HYP-2990
  - HYP-2963
  - LTI-150
  - LTI-153
  - LTI-154
---

# HYP-3008: LRC14 Automatic Gap Carrier

## Claim

Moser-de Bruijn and fibbinary numbers give a useful finite-automaton carrier for
the LRC14 `2`-adic proof route only when the automaton state is retained.

The useful split is:

```text
Moser-de Bruijn:
  binary 1s only in even positions
  base-4 digits only 0 or 1
  stable under n -> 4n
  not stable under n -> 2n unless the even/odd bit-position phase is remembered

Fibbinary:
  binary expansion has no adjacent 1s
  stable under n -> 2n
  exposes carry debt under n -> n+1
```

Thus Moser-de Bruijn is an exact automatic subset of fibbinary, but it is a
`4`-adic / even-position gap carrier rather than a theorem-safe `2`-adic scalar
on its own.  Fibbinary is the safer `2`-adic normal-form language, matching the
Zeckendorf/path-rank lane of HYP-3000 and the summand/multiplicand guardrail of
HYP-3003.

## External Source Discipline

The cited arXiv paper `2506.04110` is about the `2`-adic Littlewood conjecture:
it improves Badziahin's lower-bound obstruction from `C >= 8` to `C >= 15` and
uses Hurwitz multiplication-by-`2` machinery.  This motivates the `n -> 2n`
closure test here; it does not prove an LRC14 statement.

The Ostrowski-Hadamard gap theorem concerns lacunary power-series exponents
with ratios bounded below by some `lambda > 1`.  In this note the theorem is
used only as a warning: the primitive Moser atoms `4^j` are lacunary, while the
sorted value set is a finite automaton language with dense local clusters.

The ACM DOI supplied in the prompt resolves to "Finding large sticks and
potatoes in polygons", not a Fermat-Catalan paper.  I treated it as a
largest-safe-component / inner-approximation analogy only.  It contributes no
arithmetic closure fact.

The Fermat-Catalan input is therefore handled conservatively as a valuation
gate: any p-adic or linear-recurrence zero-set argument must emit its valuation
or eventual-periodic coordinate explicitly before it can interact with these
automatic languages.

## Computation

Script:

```text
04-computation/lrc14_moser_fibbinary_automatic_gap_carrier_codex_20260626.py
```

Result:

```text
05-knowledge/results/lrc14_moser_fibbinary_automatic_gap_carrier_codex_20260626.out
```

The script builds three exact DFAs:

```text
moser_binary_even_lsb      states: even, odd, dead
fibbinary_no_adjacent_lsb  states: last0, last1, dead
moser_base4_digits_0_1     states: ok, dead
```

Checks through `2^14-1`:

```text
Moser subset of fibbinary: True
binary Moser DFA agrees with base-4 Moser DFA: True

bits=14:
  Moser count      = 128 = 2^7
  fibbinary count  = 987 = F_16
  intersection     = 128
```

Closure under the operations most relevant to the current LRC14 mode atlas:

```text
moser     n -> 2n:   1/128 preserved
moser     n -> 4n:  64/64 preserved
moser     n -> n+1: 64/128 preserved

fibbinary n -> 2n: 610/610 preserved
fibbinary n -> 4n: 377/377 preserved
fibbinary n -> n+1: 610/987 preserved
```

Residue profiles mod `14` are not uniform enough to replace packet labels:

```text
Moser residues:     min=9,  max=10
fibbinary residues: min=52, max=88
```

The raw characteristic prefixes also fail a small eventual-period stress:

```text
no period <= 64 after cutoff <= 512, through 2^14-1
```

This is not a theorem about non-periodicity.  It is a quotient warning: raw
automatic value-set membership should not be silently identified with a
Skolem-Mahler-Lech zero set, whose index set is eventually periodic for linear
recurrences.

## Tournament Analysis

Vertices are proof-language carriers, not runners:

```text
product_phase_automaton
fibbinary_2adic_shift_dfa
fermat_catalan_valuation_gate
sml_eventual_periodic_zero_gate
moser_base4_gap_dfa
ostrowski_hadamard_atom_gap
stick_potato_safe_geometry
raw_sequence_scalar_shadow
```

Pairwise observable:

```text
packet_label_retention
2adic_shift_control
carry_boundary_memory
gap_lacunarity
sml_eventual_period_guard
finite_automaton_exactness
lrc_transfer_safety
```

Tie Hamiltonian path is the listed order.  The exact fingerprint is:

```text
score_hist=[(0,1),(1,1),(2,1),(3,1),(4,1),(5,1),(6,1),(7,1)]
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1,1]
hamiltonian_path_count=1
score_order:
  product_phase_automaton
  > fibbinary_2adic_shift_dfa
  > fermat_catalan_valuation_gate
  > moser_base4_gap_dfa
  > sml_eventual_periodic_zero_gate
  > stick_potato_safe_geometry
  > ostrowski_hadamard_atom_gap
  > raw_sequence_scalar_shadow
```

## Assumption Challenge

Candidate vertex sets considered:

```text
runners, gaps, fixed circle sections, section boundaries, wall crossings,
residues, cover arcs, Fourier modes, matroid circuits, proof obligations,
automatic-language states, and valuation gates
```

I chose automatic-language / valuation proof carriers because the predicate
under test is not loneliness itself.  It is whether a quotient is allowed to
forget a `2`-adic phase, a carry boundary, or an eventual-periodic recurrence
coordinate.

Preserved LRC predicate:

```text
quotient safety for sequence-shadow packet labels under 2-adic shifts and carry moves
```

Destroyed information:

```text
exact endpoint owners, safe interval lengths, Farey binding scale, and
state-lift labels unless attached from the HYP-2963 packet bank
```

Challenged assumption:

```text
"Moser-de Bruijn is sparse/lacunary, so it should automatically be a good
2-adic LRC carrier."
```

The scout refutes that scalar version.  Moser's base-4 gap structure is real,
but the `2`-adic LRC route sees the missing odd-position phase immediately:
`n -> 2n` destroys Moser membership for every nonzero audited Moser value.

## Next Proof Target

Add an `automatic_gap_carrier` field to HYP-2963-style packet records:

```text
none
fibbinary_2adic_normal_form
moser_even_phase_gap
product_phase_automaton
sml_eventual_periodic_gate
valuation_gate_only
```

Then test whether any hard non-AP/GW packet has an automatic-language label
that predicts its existing certificate route:

```text
Fejer / Toeplitz
Ramanujan exact-period
endpoint-owner bridge
C27 / unital
K33 / state-lift
F7 residual
```

If the automatic label mixes routes inside one fiber, it becomes a guardrail
against scalar sequence-shadow arguments rather than a proof shortcut.
