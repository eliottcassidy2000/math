---
id: HYP-3032
title: LRC14 analytic sieve-clock bridge into the side-channel repair ladder
status: EVIDENCE / finite analytic-clock quotient audit and zipper-theorem target; not a proof
source: codex-2026-06-26-S196
tangent: T1113
script: 04-computation/lrc14_analytic_sieve_clock_bridge_codex_s196.py
result: 05-knowledge/results/lrc14_analytic_sieve_clock_bridge_codex_s196.out
related:
  - HYP-3027
  - HYP-3026
  - HYP-3024
  - HYP-3023
  - HYP-3020
  - HYP-2985
  - HYP-2984
  - HYP-2983
  - HYP-2982
  - HYP-2979
  - HYP-2978
  - HYP-2963
  - HYP-2997
  - HYP-2995
  - HYP-2992
  - THM-572
  - OPEN-Q-108
---

# HYP-3032: LRC14 Analytic Sieve-Clock Bridge

## Claim

The analytic-number-theory packet from HYP-2982/HYP-2983 should be attached to
the HYP-3027 repair ladder as a family of clocks, not as scalar proof data.
In particular:

```text
sum mu(n)/n        = cancellation / tail guard
sum mu(n)^2/phi(n)= inverse primitive-unit capacity
large sieve        = minor-arc family majorant
exponential sums   = resonance checksums
smoothing choices  = kernel-homotopy defect ledger
Kaczynski labels   = boundary approach-class data
```

The central guardrail is:

```text
mu^2/phi is a capacity meter with a blindness certificate.
```

It sees squarefree primitive-unit capacity, but it kills prime powers and
repeated-prime packets.  Therefore it is useful only if the squarefree
blindness report, exact denominator, smoothing policy, exponential-sum
checksum, Kaczynski boundary approach, and LRC packet labels remain attached or
are discharged by an existing dual/zipper theorem.

## Computation

Script:

```text
04-computation/lrc14_analytic_sieve_clock_bridge_codex_s196.py
```

Stored output:

```text
05-knowledge/results/lrc14_analytic_sieve_clock_bridge_codex_s196.out
```

The audit uses the HYP-3026 named carrier bank and adds the two concrete
residual mixed pairs from HYP-3027:

```text
drop(10,13)->add(20,26)  vs  drop(8,12)->add(16,24)
drop(12,13)->add(26,36)  vs  single swap 12->72
```

It computes exact `M=p/q`, safe status, arithmetic factorization, `mu(q)`,
`phi(q)`, `mu^2(q)/phi(q)`, prefix sums `sum mu(n)/n`,
`sum mu(n)^2/phi(n)`, `sum phi(n)`, `theta(q)`, and prime reciprocal sums.

## Key Readout

The squarefree-blind rows are:

```text
loose_12_to_26
GW_fiber_tail_12_to_38
GW_fiber_tail_12_to_52
GW_fiber_tail_12_to_150
petal_10_to_20
petal_13_to_26
P10_plus_GW
P10_plus_K33
fibbinary_first13
```

This is exactly the warning from HYP-2982 in theorem-facing form: the C27
`q=27=3^3` petal branch and the `q=25` fibbinary control vanish under
`mu^2/phi`.  They cannot be managed by a squarefree capacity scalar unless the
Ramanujan/Fejer/p-adic side channels are reattached.

The HYP-3027 residual pair survives the analytic quotient stress:

```text
residual_petal_drop10_13_add20_26:
  route=BOUNDARY-PETAL-SPORADIC, M=2/23, q=23

residual_cover_drop8_12_add16_24:
  route=COVERING-MOMENT, M=2/23, q=23
```

Both are prime-squarefree `q=23`, both have `mu^2/phi=1/22`, and even the
non-route analytic packet

```text
(prime_squarefree_unit, q=23, open, bar_count=6,
 boundary_count=18, zero_sum_pairs=6, first_chart_den=14)
```

still mixes their theorem route.  Thus the analytic bridge does not bypass
HYP-3027; it identifies the exact place where the next packet label must fire.

The second HYP-3027 pair separates arithmetically:

```text
drop(12,13)->add(26,36): M=3/37, q=37, K33 state-lift
single swap 12->72:      M=6/77, q=77, covering moment
```

So exact denominator/factorization can already distinguish that mixed
topology pair, while the first residual pair needs a more geometric or
packet-family coordinate.

## Quotient Stress

The S196 audit reports:

```text
mu2_phi_blind_flag:
  fibers=2, mixed_route=2, mixed_status=1

mu2_phi_value_only:
  fibers=9, mixed_route=2

sieve_clock_only:
  fibers=6, mixed_route=4

sieve_clock_plus_exact_denominator:
  fibers=11, mixed_route=2

nonroute_clock_packet:
  fibers=12, mixed_route=1
```

The only non-route analytic-packet mix left in the small audit is the
HYP-3027 petal/covering residual at `q=23`.  Declared smoothing/exit labels
split it, but those labels are route handoffs, not independent proof quotients.

## Tournament Analysis

Vertices are analytic proof clocks and repair packets, not runners or arcs.
The pairwise observable is:

```text
status retention
route retention
denominator retention
squarefree blindness
smoothing retention
exponential-sum declaration
boundary approach
finite checkability
noncircularity
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1,1,1,1]
hamiltonian_path_count=1
```

High-retention path:

```text
labelled_repair_ladder_packet
> analytic_sieve_clock_bridge
> kaczynski_boundary_approach
> smoothing_explicit_formula_packet
> exponential_sum_checksum
> circle_method_major_minor_split
> large_sieve_minor_arc_gate
> mu2_phi_inverse_unit_clock
> mobius_mu_over_n_tail
> raw_prime_count
```

## Proof Target

Analytic repair-clock lemma:

```text
Inside a fixed automatic/residue/fusion fiber, the first nonzero analytic
clock among
  mu/n tail,
  mu^2/phi capacity,
  large-sieve minor-arc budget,
  exponential-sum checksum,
  smoothing defect,
  Kaczynski approach class
either opens a strict component, descends to AP/GW/C27/K33/covering,
is dual-annihilated by Fejer/Ramanujan/Haar, or emits F7/THM-572 residual debt.
```

This turns prime sums, Mobius sums, quadratic-sieve/large-sieve language,
circle-method major/minor arcs, smoothing, saddle/explicit-formula packets, and
Kaczynski boundary ambiguity into typed repair teeth inside the existing LRC14
packet theorem target.

## Assumption Challenge

Alternate vertices considered: runners, primes, denominators, residues,
Fourier modes, smoothing kernels, exponential phases, boundary approach
classes, local obstructions, and proof obligations.

Chosen vertices are analytic proof clocks and repair packets because those are
the objects whose quotient safety is being tested.  They preserve
open/boundary status and theorem-route purity only when exact denominator,
squarefree blindness, smoothing choice, exponential checksum, Kaczynski
approach class, and packet labels remain attached or are discharged.

The challenged assumption is that `sum mu(n)/n`, `sum mu(n)^2/phi(n)`, or
prime-sum data can itself be an LRC14 proof quotient.  S196 says no: these are
high-value clocks, but they become theorem data only inside the HYP-3027
first-nonzero repair ladder.
