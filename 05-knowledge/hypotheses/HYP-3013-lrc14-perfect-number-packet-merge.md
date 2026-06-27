---
id: HYP-3013
title: LRC14 perfect-number packet merge and divisor-lattice guardrail
status: SYNTHESIS / packet-schema guardrail and proof-interface carrier; not a proof
source: codex-2026-06-26-S174
script: 04-computation/lrc14_perfect_number_packet_merge_codex_s174.py
result: 05-knowledge/results/lrc14_perfect_number_packet_merge_codex_s174.out
related:
  - HYP-2221
  - HYP-2220
  - HYP-2941
  - HYP-2945
  - HYP-2946
  - HYP-2999
  - HYP-3003
  - HYP-3008
  - HYP-3009
  - HYP-3011
  - HYP-3012
  - HYP-2963
  - THM-572
  - OPEN-Q-108
---

# HYP-3013: LRC14 Perfect-Number Packet Merge

## Claim

This pass merges the older perfect-number / aliquot fixed-point carrier
(`HYP-2221`, `HYP-2941`, `HYP-2945`) into the current LRC14 automatic-gap
packet stack (`HYP-3008` through `HYP-3012`).

The safe theorem-facing statement is:

```text
Even perfect numbers are the exact n=2 unit-excess control chain.
The LRC14 lane q=14a-1 is only a prime-q deficient shadow, and composite q
rows prove that divisor factorization is load-bearing.
```

So the merge is not "perfect products solve LRC14."  It is a packet-schema
guardrail: retain exact `M=p/q`, the unit-excess address, prime/composite `q`,
factorization, abundancy defect, product/Kpq incidence, and automaton state
before using perfect-number analogies.

## Computation

Script:

```text
04-computation/lrc14_perfect_number_packet_merge_codex_s174.py
```

Output:

```text
05-knowledge/results/lrc14_perfect_number_packet_merge_codex_s174.out
```

The even-perfect controls through Mersenne exponent `17` are the expected
Euclid-Euler rows:

```text
N=6,28,496,8128,33550336,8589869056
sigma(N)/N=2
defect=0
```

For the LRC14 power-of-two shadows `a=2^k`, `q=14a-1`, `0<=k<=12`, the script
finds:

```text
lrc14_shadow_rows=13
prime_q14_shadow_rows=3
composite_q14_shadow_rows=10
prime_q14_defect_formula_ok=1
```

The prime `q=14a-1` rows occur at:

```text
a=1,   q=13,   defect=12/13
a=16,  q=223,  defect=12/223
a=256, q=3583, defect=12/3583
```

The exact formula is the HYP-2945 bridge:

```text
sigma(aq)/(aq) = n(2a-1)/(na-1) = 2 + (2-n)/(na-1)
```

At `n=14`, prime `q=14a-1` gives deficiency `12/(14a-1)`.  At `n=2`,
prime `q=2a-1` gives defect `0`, the perfect-number control.

The composite `q14` rows are the warning.  In the bounded scan every composite
`q14` row is abundant rather than deficient, for example:

```text
a=2, q=27=3^3, defect=-2/9
a=4, q=55=5*11, defect=-16/55
a=8, q=111=3*37, defect=-21/37
```

Thus the primality / divisor-lattice label is not a cosmetic annotation.  A
quotient that keeps only `p*q` or only the power-of-two address erases the
sign of the abundancy defect.

## Packet Fields

Add these fields to any HYP-2963-style packet or sidecar that imports the
perfect-number route:

```text
unit_excess_apex
perfect_control_status
abundancy_defect
divisor_lattice_factorization
prime_q_flag
product_incidence_rank
automaton_transition_state
```

Use `perfect_control_status` only for the `n=2` calibration lane.  For LRC14,
the theorem-facing data are the exact defect and its factorization witness.

## Tournament Analysis

Vertices are proof carriers and side channels, not runners, arcs, or raw
sequence entries:

```text
labelled_lrc_packet_sheaf
divisor_lattice_abundancy_packet
exact_farey_unit_excess
lrc14_n14_deficient_shadow
perfect_n2_fixed_control
k33_product_incidence
automatic_power_of_two_state
fermat_catalan_power_budget
raw_product_scalar
```

Pairwise observable:

```text
lrc_predicate
exact_scale
divisor_lattice
graph_minor
automaton_state
power_exception
anti_scalar
proof_maturity
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1,1,1]
hamiltonian_path_count=1
```

The resulting order puts the labelled packet sheaf first and raw product
scalar last.  The only safe scalar use of the perfect-number lane is as a
negative control after divisor, exact Farey, Kpq, and automaton labels have
already been attached.

## Assumption Challenge

Alternate vertices considered: runners, gaps, fixed circle sections, section
boundaries, wall-crossing events, residues, cover arcs, Fourier modes, matroid
circuits, divisor atoms, and proof obligations.

The chosen quotient preserves the LRC predicate only through labelled packet
fields.  It deliberately destroys raw runner identity and refuses to collapse
divisor or automaton information to a product scalar.

## Next Pulls

1. Add the seven packet fields above to the HYP-2963 sidecar schema for any
   unit-excess / product-incidence route.
2. Audit hard non-AP/GW packets on `q=14p-1` by `(prime_q_flag,
   abundancy_defect, automatic_transition_state, Kpq_route)`.
3. Compare composite `q14` abundant rows against the HYP-3009 perfect-power
   no-lift ledger; `q=27=3^3` is the first shared warning row.
4. Keep HYP-2946's Kuratowski split attached: perfect products give edge-load
   controls, but K33 state-lift debt begins at the graph-minor packet, not at
   raw edge count.
