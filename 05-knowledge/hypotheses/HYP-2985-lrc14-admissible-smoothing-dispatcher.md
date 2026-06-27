---
id: HYP-2985
title: LRC14 admissible smoothing dispatcher and packet-exit lemma
status: PROOF-INTERFACE / routing theorem target, not a proof
source: codex-2026-06-24-S164
script: 04-computation/lrc14_smoothing_dispatcher_codex_20260624.py
result: 05-knowledge/results/lrc14_smoothing_dispatcher_codex_20260624.out
related:
  - HYP-2984
  - HYP-2983
  - HYP-2982
  - HYP-2981
  - HYP-2979
  - HYP-2978
  - HYP-2974
  - HYP-2973
  - HYP-2972
  - HYP-2963
  - HYP-2908
  - HYP-2901
  - HYP-2679
  - THM-548
  - THM-572
  - OPEN-Q-108
---

# HYP-2985: LRC14 Admissible Smoothing Dispatcher

This pass merges HYP-2982's analytic weight atlas and HYP-2983's
Kaczynski/exponential-sum proof template into a single dispatcher.  It was
rebased after HYP-2984 was claimed for the sibling kernel-homotopy lane; the
two are complementary.  HYP-2984 asks when a kernel deformation preserves a
certificate or emits a boundary defect, while HYP-2985 asks which packet family
is allowed to use which kernel or smoothing policy.

Working claim:

```text
Every LRC14 smoothing or analytic quotient must declare the packet family it
is allowed to act on, the labels it retains, the labels it forgets, and the
handoff for any forgotten label.
```

The theorem target is an admissible-smoothing lemma on HYP-2963 packet fibers:

```text
Every primitive packet is one of:

1. AP/GW boundary equality;
2. labelled Fejer interval certificate;
3. Ramanujan exact-period / prime-power handoff;
4. labelled Selberg or large-sieve late-denominator minor-arc packet;
5. off-resonance Kaczynski/Abel decay packet;
6. Freiman-reduced resonant finite atlas row;
7. HYP-2908/THM-572 state-lift obligation.
```

The point is not to add one more scalar estimate.  The point is to decide which
clock is being used and which information that clock is allowed to forget.

## Computed Dispatcher

Script:

```text
04-computation/lrc14_smoothing_dispatcher_codex_20260624.py
```

Stored output:

```text
05-knowledge/results/lrc14_smoothing_dispatcher_codex_20260624.out
```

The dispatcher treats smoothing policies as proof carriers:

```text
raw_scalar_density
mertens_mu_over_n_tail
selberg_large_sieve_squarefree
ramanujan_exact_period_projector
fejer_toeplitz_interval
kaczynski_abel_boundary
hybrid_labelled_packet_sheaf
```

and packet families as proof obligations:

```text
AP_GW_boundary_atom
K33_near_12_to_36
P10_plus_GW_splice
petal_q27_unit_strip
covering_q41_or_q63_front
late_prime_power_denominator_wall
truewide_off_resonance_far_packet
truewide_resonant_freiman_packet
```

The useful local-policy readout is:

```text
AP/GW boundary              -> Kaczynski/Abel boundary labels
K33 near 12->36             -> Fejer/Toeplitz interval certificate
P10+GW and q=27 petals      -> Fejer plus Ramanujan prime-power side channel
covering q41/q63 front      -> Fejer/Toeplitz interval certificate
late prime-power wall       -> Kaczynski/Ramanujan; Selberg only as precondition
true-wide off resonance     -> Kaczynski/Abel signed decay
true-wide resonance         -> Kaczynski/Freiman finite reduction or state lift
```

This sharpens the HYP-2982 warning.  The squarefree weight `mu^2/phi` is a
Selberg/large-sieve normalizer, not a final LRC14 certificate.  It sees `q=14`
and prime `q=41`, but erases prime-power or repeated-prime packets such as
`25,27,36,63,84,98,168,280,4312`.  Those packets need exact-period,
endpoint-owner, Fejer, Kaczynski, or state-lift side channels before any
scalar quotient is admissible.

## Four Clocks That Matter

1. Endpoint-owner clock.
   This is the taut-wall/bridge clock.  It decides AP/GW boundary equality and
   guards against forgetting which endpoint pair owns a closed witness.

2. Exact-period denominator clock.
   This is the Ramanujan/Farey clock.  It catches `q=27`, `q=41`, `q=4312`,
   and other packets that squarefree weights can erase.

3. Smoothing/certificate clock.
   This is the Fejer/Toeplitz clock.  It turns positive-open packets into
   interval-enclosed dual certificates.

4. Far-approach boundary clock.
   This is the Kaczynski/Abel clock.  It separates true-wide off-resonance
   decorrelation from resonant Freiman finite models.

The clocks that can usually be demoted are raw prime counts, raw `sum mu(n)`,
and raw `sum mu(n)/n`.  They are useful diagnostics or tails, but not
lonely-interval certificates without packet labels.

## Proof Use

HYP-2985 gives a way to use the analytic-number-theory analogy without
overclaiming it:

```text
major arcs    -> exact-period / endpoint packet families
minor arcs    -> labelled Selberg/large-sieve family estimates
smoothing     -> Fejer or Kaczynski kernel chosen by packet type
exceptions    -> boundary approach classes, not averaged-away errors
finite check  -> AP/GW, petal, K33, covering, or state-lift atlas
```

The promising next lemma is not "large sieve proves LRC14."  It is:

```text
After HYP-2963 packet classification, each packet has an admissible smoothing
policy.  If all admissible policies fail, the failure data itself constructs a
named HYP-2908/THM-572 state-lift obligation.
```

That is a sharper interface than another global bound, because it says exactly
where Fejer, Ramanujan, Selberg, and Kaczynski are allowed to act.

## Tournament Analysis

Vertices are smoothing/quotient policies, not runners.

Pairwise observable:

```text
retained LRC predicate payload,
exact-period data,
boundary approach,
interval formalizability,
state-lift handoff.
```

Stored fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1}
directed_3cycles=0
hamiltonian_paths=1
Hamiltonian path:
hybrid_labelled_packet_sheaf >
kaczynski_abel_boundary >
fejer_toeplitz_interval >
ramanujan_exact_period_projector >
selberg_large_sieve_squarefree >
mertens_mu_over_n_tail >
raw_scalar_density
```

## Assumption Challenge

Alternate vertex sets considered: runners, denominators, endpoint-owner pairs,
exact-period modes, smoothing policies, far-approach classes, packet families,
and proof obligations.  The chosen quotient uses smoothing policies as
vertices because the question is no longer "which runner wins" but "which
proof clock is allowed to certify this packet."  It preserves the LRC predicate
only when packet labels are retained.  It destroys raw runner identity, so any
final theorem must attach endpoint-owner, exact-period, interval, boundary, and
state-lift labels before applying a scalar estimate.
