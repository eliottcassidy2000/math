---
id: HYP-2983
title: LRC14 analytic sieve/Kaczynski exponential-sum packet synthesis
status: SYNTHESIS / analytic proof-template integration, not a proof
source: codex-2026-06-24-S162
script: 04-computation/lrc14_analytic_sieve_kaczynski_s162.py
result: 05-knowledge/results/lrc14_analytic_sieve_kaczynski_s162.out
related:
  - HYP-2982
  - HYP-2980
  - HYP-2979
  - HYP-2978
  - HYP-2977
  - HYP-2976
  - HYP-2975
  - HYP-2974
  - HYP-2973
  - HYP-2972
  - HYP-2971
  - HYP-2970
  - HYP-2969
  - HYP-2968
  - HYP-2966
  - HYP-2965
  - HYP-2964
  - HYP-2963
  - HYP-2956
  - HYP-2946
  - HYP-2938
  - HYP-2900
  - HYP-2899
  - HYP-2679
  - HYP-2227
  - HYP-1963
  - HYP-1962
  - HYP-1953
  - THM-548
  - THM-523
  - THM-566
  - THM-571
  - THM-572
  - OPEN-Q-108
---

# HYP-2983: LRC14 Analytic Sieve/Kaczynski Exponential-Sum Packet Synthesis

This hypothesis imports the prompt's analytic-number-theory cluster into the
current LRC14 source-kernel program:

```text
prime sums / Goldbach templates
Mobius and squarefree-totient sums
large-sieve and circle-method decompositions
upper-bound sieve majorants
exponential-sum control
smoothing, saddle, and explicit-formula kernels
Kaczynski/Bagemihl boundary ambiguity
```

It is complementary to HYP-2982.  HYP-2982 computes the larger analytic packet
weight atlas and separates primitive denominator capacity from inverse-unit
large-sieve/Selberg weights.  HYP-2983 focuses on the Kaczynski boundary and
exponential-sum proof-template role those weights should play after the LRC
packet labels are retained.

The import is deliberately not a claim that LRC14 has an Euler product or that
prime sums are the right variables.  The claim is that the successful ternary
Goldbach/Helfgott-style proof architecture has the right *shape* after being
translated into labelled LRC packets:

```text
arithmetic main packet
  -> major/minor decomposition
  -> explicit exponential-sum bound
  -> finite exceptional verification
  -> smoothing/boundary defect ledger.
```

For LRC14, the arithmetic main packet is not a singular series over primes.  It
is the retained packet

```text
qdiv + exact M/Farey + Ramanujan/totient exact-period data
+ Haar/Baire open-vs-boundary status
+ endpoint-owner and C27/K33/state-lift labels.
```

The analogue of the minor-arc estimate is not "generic cancellation in one
scalar sum."  It is an explicit exponential-sum/resonance lemma saying that
off-resonance packets have a certified safe interval, while resonant packets
localize into the finite AP/Goddyn-Wong, petal, K33, covering, or state-lift
atlas.

## Computed S162 Artifact

Script:

```text
04-computation/lrc14_analytic_sieve_kaczynski_s162.py
```

Stored output:

```text
05-knowledge/results/lrc14_analytic_sieve_kaczynski_s162.out
```

The script scans the local repo for eight motifs:

```text
prime_sum_goldbach
mobius_mu_sums
totient_phi_packets
large_sieve_circle
exponential_sum_core
smoothing_saddle_formula
kaczynski_boundary
lrc_source_kernel
```

The top counts are:

```text
lrc_source_kernel     2920
prime_sum_goldbach    2364
exponential_sum_core  1094
totient_phi_packets    471
mobius_mu_sums         387
smoothing_saddle_formula 367
large_sieve_circle      86
kaczynski_boundary      24
```

The strongest co-occurrences are between LRC source-kernel language, prime /
Goldbach templates, exponential sums, and totient/Mobius packets.  That is
exactly the desired signal: the repo already contains pieces of the analytic
proof architecture, but they are distributed across quotient guardrails,
Fourier/Toeplitz duals, Ramanujan packets, source spectra, and boundary-value
notes.

## Explicit Weight Ledger

S162 records four bounded sums:

```text
A(N) = sum_{n<=N} mu(n)/n
B(N) = sum_{n<=N} mu(n)^2 / phi(n)
C(N) = sum_{n<=N} mu(n)^2 / n
D(N) = (1/N) sum_{n<=N} phi(n)/n
```

At `N=20000` the computed values are:

```text
A(N) =  0.00140031
B(N) = 11.23619615
C(N) =  7.06451593
D(N) =  0.60793954
```

Interpretation:

```text
A(N): cancellation diagnostic, not positive density.
B(N): squarefree primitive-packet sieve capacity.
C(N): squarefree mass ledger.
D(N): coprime-density floor, close to 6/pi^2.
```

The rule for LRC use is strict.  `sum mu(n)/n` may certify cancellation or
kernel removal, but it cannot by itself prove a positive lonely interval.
`sum mu(n)^2/phi(n)` is more valuable for this project because it measures
primitive squarefree exact-period capacity.  It should be attached to
Ramanujan/Farey endpoint-owner packets before any scalar collapse.

## Kaczynski Boundary Translation

THM-548 and HYP-2679 already place the true-wide route in a Fatou /
Bagemihl-Kaczynski boundary-value frame.  In S162 language:

```text
ordinary approach arcs       -> decorrelated/Fatou boundary limits
ambiguous approach arcs      -> resonance exceptions
resonance exceptions         -> AP/GW/K33/petal finite packets
boundary defect under kernel -> smoothing-choice ledger
```

This gives a necessary condition for any AP/GW-like tight atom:

```text
It must be a Kaczynski ambiguous boundary atom for the LRC phase functions.
All nonambiguous approaches must either yield a strict safe interval or route
to a named state-lift / finite-packet obstruction.
```

AP and Goddyn-Wong have the expected symptoms:

```text
zero strict-open Haar mass,
q=14 boundary equality,
endpoint-owner zero-credit pairs r+s == 0 mod 14,
no strict primitive exact-period witness through the tested q<=42 range,
neutrality under the current moment/twist/Fourier dual tests.
```

The conjectural classification upgrade is that every other packet loses at
least one of these symptoms: a positive bridge, a twist witness, a Toeplitz
negative direction, a danger-count dual, a lift packet, or a K33/state-lift
debt appears.

## Analytic Proof-Template Translation

The suggested proof skeleton is:

```text
1. Choose a smoothing family K_eta for the danger indicator.
2. Form a smoothed deficit F_{S,eta}(t) from C_S(t)-1.
3. Expand F_{S,eta} into finite Fourier/Ramanujan/exact-period modes.
4. Split modes by q/Farey major arcs and off-resonance minor arcs.
5. Use large-sieve / upper-bound-sieve estimates only as labelled majorants.
6. Use a Kaczynski boundary ledger to account for smoothing-choice defects.
7. Show every resonant residue packet is AP/GW equality or finite K33/petal/
   covering/state-lift debt.
```

The "explicit explicit formula" obligation should name every term:

```text
smoothed main packet
+ Ramanujan/totient exact-period contribution
+ large-sieve minor-arc remainder
+ endpoint-owner boundary defect
+ Kaczynski ambiguous-approach residual
+ finite exceptional-packet ledger.
```

A proof that hides one of these terms repeats the old scalar quotient failure.

## Tournament Analysis

S162 uses proof modules as tournament vertices, not runners, primes, or arcs.

Vertices:

```text
lrc_labelled_source_kernel
exponential_sum_engine
kaczynski_boundary_resonance
mobius_totient_packet_ledger
circle_large_sieve_decomposition
smoothing_saddle_explicit_formula
upper_bound_quadratic_sieve
raw_scalar_prime_sum
```

Pairwise observable:

```text
(LRC predicate retention,
 exact arithmetic packet retention,
 local/global decomposition strength,
 exponential-sum control,
 boundary/resonance exception handling,
 smoothing/adaptivity,
 auditability)
```

Switch:

```text
A -> B if A wins more coordinates than B.
Ties follow the declared Hamiltonian path:
lrc_labelled_source_kernel > exponential_sum_engine
> kaczynski_boundary_resonance > mobius_totient_packet_ledger
> circle_large_sieve_decomposition > smoothing_saddle_explicit_formula
> upper_bound_quadratic_sieve > raw_scalar_prime_sum.
```

Fingerprint:

```text
score_hist: {0:1, 1:1, 3:3, 5:1, 6:1, 7:1}
directed_3cycles: 1
SCC_sizes: [3,1,1,1,1,1]
nontrivial SCC:
  mobius_totient_packet_ledger
  smoothing_saddle_explicit_formula
  circle_large_sieve_decomposition
Hamiltonian_path_count: 3
```

The nontrivial SCC is the main warning.  Mobius/totient packet weights, the
circle/large-sieve split, and smoothing/saddle explicit formulas should be
treated as a coupled analytic middle layer.  Ordering them linearly is a
convenience, not a theorem.

## Assumption Challenge

The vertices in the S162 tournament are proof modules.  Other useful vertex
sets remain open and should be tested:

```text
Fourier/Ramanujan modes
major-arc denominators
endpoint-owner pairs
boundary approach arcs
smoothing kernels
wall-crossing events
finite packet obligations
proof-failure modes
```

For each quotient, the preserved predicate must be stated.  The live predicate
is either:

```text
strict LRC14 failure,
positive strict safe interval,
AP/GW boundary equality,
or named state-lift/finite-packet debt.
```

The quotient is not admissible unless it preserves that predicate or records a
recoverable defect.

## New Necessary Conditions

Any strict LRC14 counterexample packet surviving the modern stack should
satisfy all of the following:

1. It is invisible to raw denominator witnesses and lies in the `qdiv>14`
   branch.
2. Its Mobius-cancellation ledger cannot be converted into a positive density
   without endpoint-owner residual terms.
3. Its `mu^2/phi` primitive-packet capacity is fully saturated or exactly
   neutralized by a Ramanujan/Farey boundary packet.
4. Its off-resonance Fourier/Ramanujan modes are large-sieve exceptional, not
   generic.
5. Every upper-bound sieve majorant that proves positivity must name the packet
   label it majorizes; scalar majorization is inadmissible.
6. Changing the smoothing function either preserves the same family label or
   produces a recorded boundary defect.
7. Its Kaczynski approach set is ambiguous: no ordinary approach arc alone
   discharges it by decorrelated boundary values.
8. It avoids endpoint potentials, moment duals, twist ladders, Toeplitz PSD
   violations, Fejer/spectral shadows, lift packets, and taut-bridge positive
   curvature.
9. It avoids AP/GW equality while also avoiding K33/H=7 state-lift debt.

This list is intentionally overdetermined.  The hoped-for theorem is that the
conditions are mutually inconsistent.

## Candidate Breakthrough Statement

The analytic-sieve version of the source-kernel theorem is:

```text
For every primitive LRC14 source packet, either an off-resonance smoothed
exponential-sum estimate gives a strict lonely interval, or the packet belongs
to the finite resonant atlas: AP/GW boundary equality, unit-petal/two-block
discharge, K33/state-lift debt, covering/lift packet with positive bridge, or a
classified fixed-margin sporadic.
```

The missing proof is the explicit exponential-sum/resonance estimate with
adaptive smoothing and a labelled Kaczynski boundary defect term.
