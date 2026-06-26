---
id: HYP-3011
title: LRC14 automatic lacunary safe-component filter
status: SYNTHESIS / exact finite scout and proof-interface guardrail; not a proof
source: codex-2026-06-26-S171
script: 04-computation/lrc14_automatic_lacunary_packet_filter_codex_s171.py
result: 05-knowledge/results/lrc14_automatic_lacunary_packet_filter_codex_s171.out
related:
  - HYP-3008
  - HYP-3002
  - HYP-3000
  - HYP-2999
  - HYP-2997
  - HYP-2996
  - HYP-2995
  - HYP-2990
  - HYP-2963
  - HYP-1902
  - THM-572
  - LTI-160
  - LTT-062
---

# HYP-3011: LRC14 Automatic Lacunary Safe-Component Filter

## Claim

Finite automata and lacunary sequence languages are useful for the LRC14 proof
only as labelled packet filters.  They should not be used as scalar sequence
analogies.  Complementing HYP-3008's membership and closure audit for
Moser-de Bruijn and fibbinary languages, this note asks whether the first
natural automatic rows are actual LRC14 boundary threats.  The packet record
should retain:

```text
automaton state,
residue mod 14,
exact-period denominator label,
gap-block / lacunary block label,
safe-component certificate or named residual.
```

The transfer from the user's prompts is:

- The 2-adic Littlewood paper's multiplication-by-2 algorithm is a local
  finite-window machine over continued-fraction digits.  For LRC14 this is a
  model for the even-fold / dyadic packet branch, not a replacement for real
  gap geometry.
- The Ostrowski-Hadamard gap theorem says genuine Hadamard gaps create a rigid
  analytic boundary phenomenon.  In LRC language, large gaps should route to
  scale-separated interval peelers or exact-period packet certificates.
- Fibbinary numbers and the Moser-de Bruijn sequence are regular / automatic
  bit languages.  Their value for LRC14 is that the DFA state is a compact
  carry/support label.
- Fermat-Catalan contributes the proof shape: discharge infinite structured
  families, then audit a finite named exception ledger.

## Finite Scout

The script computes exact threshold-open safe components at LRC14 threshold
`1/14` for four 13-runner rows:

```text
AP = {1,...,13}
GW = {1,...,11,13,24}
first 13 positive fibbinary numbers
first 13 positive Moser-de Bruijn numbers
```

The exact output shows:

```text
AP: positive_safe_components=0, safe_measure=0
GW_12_to_24: positive_safe_components=0, safe_measure=0
fibbinary_first13: positive_safe_components=38, safe_measure=66077/399840
moser_de_bruijn_first13: positive_safe_components=64, safe_measure=4264747/40348854
```

Thus the two automatic/lacunary rows are not boundary threats in this finite
stress test.  AP and GW remain the zero-open denominator-14 boundary atoms.

## Automaton Guardrail

The fibbinary DFA reads binary words least-significant bit first and forbids
`11`.  The Moser-de Bruijn DFA reads least-significant bit first and forbids a
`1` in odd bit position.  The latter is a sublanguage of the former, but the
language inclusion is not a proof certificate: it says only that a packet has a
sparse carry/support label.

The proposed HYP-2963 packet fields are:

```text
binary_language_state,
moser_even_bit_support,
fibbinary_no11_state,
gap_block_profile,
hadamard_block_ratio,
automatic_filter_exit,
```

and each field must be paired with exact `M=p/q`, endpoint owners,
exact-period labels, and a certificate route.

## Tournament Analysis

Vertices are proof carriers and automatic packet filters, not runners:

```text
labelled_packet_dfa_filter
hurwitz_2adic_window_automaton
fibbinary_no11_carry_filter
moser_even_bit_support_filter
hadamard_gap_block_certificate
fermat_catalan_exception_ledger
raw_growth_or_sequence_name
```

Pairwise observable: retention vector over

```text
predicate, finite_state, exact_period, gap_block, exception_ledger,
certificate, anti_scalar
```

Switch/gauge: orient `A -> B` if `A` has more strictly larger retention
coordinates than `B`; ties use the listed Hamiltonian path.  The scout reports
a transitive carrier tournament with one Hamiltonian path and raw sequence
name/growth last.

## Proof Target

Add a DFA-labelled packet field to the HYP-2963 classifier.  A primitive LRC14
packet may use automatic data only if the automaton state either:

```text
gives a strict-safe component certificate,
routes to AP/GW boundary equality,
descends to a scale-separated family,
is annihilated by Fejer/Ramanujan/endpoint data,
or is emitted as named F7/THM-572 residual debt.
```

This makes automaticity part of the controlled-kernel zipper rather than
another scalar shadow.
