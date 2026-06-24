---
id: HYP-2978
title: LRC14 Ramanujan-divisor quotient guardrails
status: INQUIRY / quotient-admissibility proof lane with finite collision audit; not a proof
source: codex-2026-06-24-S161
related:
  - HYP-2979
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
  - HYP-2963
  - HYP-2956
  - HYP-2946
  - HYP-2938
  - HYP-2887
  - HYP-2886
  - HYP-2264
  - THM-406
  - THM-572
  - OPEN-Q-108
artifacts:
  - 04-computation/lrc14_ramanujan_divisor_quotient_guardrails_codex_s161.py
  - 05-knowledge/results/lrc14_ramanujan_divisor_quotient_guardrails_codex_s161.out
---

# HYP-2978: LRC14 Ramanujan-Divisor Quotient Guardrails

This hypothesis reserves the divisor/Ramanujan quotient-admissibility lane
requested on 2026-06-24.  The guiding principle is:

```text
A quotient is admissible for an LRC14 proof only if it preserves the predicate
needed by the next implication, or records an explicit certificate explaining
what information was intentionally forgotten.
```

Equivalently, each quotient must declare:

```text
preserved predicate
forgotten labels
compensating transform or side-channel
defect certificate when a forgotten label is load-bearing
```

The intended external seed is the divisor-function neighborhood:
`sigma_k(n)`, Dirichlet convolution, multiplicativity, Ramanujan sums
`c_q(n)` as primitive-root power sums, and the bridge from divisibility data to
Fourier/cyclotomic packets.  The intended internal seed is the repeated repo
lesson that scalar quotients are useful only after labelled fibers are retained:
irreducible cores, unital designs, Faulhaber moment positivity, Pollock degree
jumps, unit-distance norm layers, tiling/solid analogies, and the current LRC14
dual stack.

The immediate LRC14 hook is HYP-2974's divisor-curried Fourier coefficient:

```text
hat F_S(k) = sum_{v in S, v|k} sin(pi*(k/v)/7)/(pi*(k/v)).
```

Mode `k` sees both the divisor fiber `v|k` and the quotient `k/v`.
Ramanujan sums `c_q(n)`, the sums of nth powers of primitive q-th roots, are
therefore candidate exact-period unit characters: they retain primitive
residue-period data after averaging over units instead of collapsing to a bare
divisor count.

HYP-2979 is the first concrete child route: exact-period Ramanujan projectors
for q-ladders, endpoint sums, and primitive unit phase packets.

Finite audit added by S161:

```text
qdiv_only route-mixing collisions                 1
open_state_only route-mixing collisions           1
mod14_residue_multiset route-mixing collisions    1
ramanujan_14_profile route-mixing collisions      1
unit_counts_14_27_41 route-mixing collisions      2
divisor_lcm_scalars route-mixing collisions       1
guarded_packet_signature route-mixing collisions  0
```

The named-row table shows the main warning sharply.  AP, GW, the residue liar
`12->26`, near/K33 `12->36`, petals, and `P10+K33` all share the same coarse
lcm divisor scalars; several also share the same `c_14` profile.  Those
channels mix AP/GW boundary, q-witness, K33, petal, and covering proof routes.
The over-labelled packet signature avoids route mixing precisely because it
keeps q/Farey, Haar-open status, endpoint/state labels, and the relevant dual
certificate exits attached.

Tournament Analysis over quotient channels is transitive:

```text
labelled_lrc_packet_sheaf
  > exact_period_packet
  > ramanujan_primitive_shell
  > gcd_strata
  > totient_jordan_unit_capacity
  > squarefree_psi_support
  > raw_divisor_counts
```

with `score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1}` and no directed 3-cycles.
The readout is not that Ramanujan shells are complete.  It is that they are the
first reasonable arithmetic side-channel above scalar divisor counts, while
endpoint owners and state-lift labels remain load-bearing.

Remaining work:

1. Read the divisor-function page and one-hop related pages around multiplicative
   functions, Ramanujan sums, Dirichlet convolution, abundant/perfect numbers,
   and divisor summatory problems.
2. Extend the collision audit from named rows to packet families in the
   HYP-2963 bank.
3. Turn the following admissibility criterion into a theorem-facing lemma:

```text
Any quotient used to rule out an LRC14 counterexample must retain enough
divisor/cyclotomic phase data to distinguish AP/GW boundary atoms, positive
Toeplitz/Ramanujan exits, and K33/state-lift debts.
```

Expected falsifier shape: two rows with the same scalar divisor signature but
different LRC route.  The S161 audit already exhibits such pairs for multiple
scalar channels, so those channels are demoted to features.  Any final quotient
must include Ramanujan packet labels or exact endpoint-owner/Farey data.

Broader theorem target:

```text
Every post-THM-571 Moon-core packet either
  has a quotient whose preserved data forces a known dual certificate,
  or exhibits a Ramanujan/divisor defect showing exactly which label
  the quotient illegally forgot.
```
