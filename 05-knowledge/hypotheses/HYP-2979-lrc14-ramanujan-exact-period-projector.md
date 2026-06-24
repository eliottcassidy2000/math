---
id: HYP-2979
title: LRC14 Ramanujan exact-period projector packets
status: PROOF-INTERFACE / exact-period primitive-phase audit; not a proof
source: codex-2026-06-24-ramanujan-projector
related:
  - HYP-2977
  - HYP-2978
  - HYP-2975
  - HYP-2974
  - HYP-2973
  - HYP-2972
  - HYP-2970
  - HYP-2969
  - HYP-2963
  - HYP-2901
  - HYP-2908
  - THM-572
  - OPEN-Q-108
artifacts:
  - 04-computation/lrc14_ramanujan_exact_period_projector_codex_20260624.py
  - 05-knowledge/results/lrc14_ramanujan_exact_period_projector_codex_20260624.out
  - 07-reflections/lrc14-ramanujan-exact-period-projectors-codex-20260624.md
---

# HYP-2979: LRC14 Ramanujan Exact-Period Projector Packets

This file reserves the exact-period Ramanujan-sum route requested in the
2026-06-24 Ramanujan prompt. It is deliberately separate from HYP-2978's
Ramanujan-divisor quotient guardrail: HYP-2978 says what a quotient is allowed
to forget, while this hypothesis proposes one retained exact-period quotient.
The proposed object is not a proof yet. It is a packet interface between:

- HYP-2972's rational twist ladder,
- HYP-2973's danger-count moments,
- HYP-2974's Fourier-Toeplitz dual,
- HYP-2975's endpoint taut-current layer, and
- HYP-2977's spectral-shadow route.

For a denominator `q`,

```text
c_q(n) = sum_{1 <= a <= q, gcd(a,q)=1} exp(2*pi*i*a*n/q)
       = sum_{d | gcd(q,n)} d*mu(q/d)
       = mu(q/gcd(q,n))*phi(q)/phi(q/gcd(q,n)).
```

The working idea is to replace scalar "q is blocked" data by primitive-unit
exact-period packets. For LRC14, the layer `q=14` is especially sharp:

```text
gcd(14,n)=14 -> c_14(n)= 6
gcd(14,n)= 7 -> c_14(n)=-6
gcd(14,n)= 2 -> c_14(n)=-1
gcd(14,n)= 1 -> c_14(n)= 1
```

Thus `c_14(r+s)` is a primitive-unit trace of the AP/Goddyn-Wong zero-credit
taut relation `r+s == 0 mod 14`, while `c_14(r-r')` records exact-period
residue coincidence with parity/seven-adic labels still attached.

The theorem-facing hope is:

```text
Any LRC14 residual that is invisible to raw qdiv, Haar/Baire boundary support,
danger-count moments, Toeplitz PSD, and spectral shadows must be visible in a
Ramanujan exact-period packet; otherwise it is an AP/GW equality atom or it
routes to the HYP-2908/THM-572 K33/H=7 state lift.
```

The computation still needs to decide whether the useful packet is a primitive
phase witness profile, a shifted Carmichael autocorrelation profile, a
Ramanujan-subspace projection of `N_S(t)`, or a supercharacter table over the
unit-group quotient.

## Computed Evidence

The first HYP-2979 audit adds an exact integer projector on functions over
`Z/qZ`:

```text
E_q(f) = sum_{a,b mod q} f(a) f(b)c_q(a-b).
```

For the danger count `N_S(a/q)` and the weak-safe indicator, this is the
primitive-frequency shell energy.  It keeps exact-period variation of the
phase function, rather than only the speed multiset.

Named rows:

```text
row                  qdiv first weak q first strict q safe_mu
AP                     14           14          None        0
GW 12->24              14           14          None        0
residue liar 12->26    12           12            12 426/35035
near 12->36            14           14            41   1/1260
petal 10->20           14           14            27    1/980
petal 13->26           14           14            27    1/182
P10+GW                 14           14            27    1/980
P10+K33                14           14            27   4/2205
covering 12->84        15           41            41 563/105105
covering 12->168       15           41            41  263/30030
covering 6->98         15           25            25 1543/294294
```

For q=14, AP/GW and many q14-front rows share the same primitive boundary
packet `(((0,2),6),)`: each primitive phase is weak-safe with two boundary
speeds.  Covering rows instead show `(((1,2),6),)`: every primitive q=14 phase
has one danger speed and two boundary speeds.  This is a useful local guard,
but it still mixes proof routes.

Stress pass:

```text
rows audited          21906
no weak q<=42             0
no strict q<=42           2
no strict examples        AP, GW 12->24
```

The bank consists of named rows plus AP one-swaps through replacement `180` and
two-swaps through added values `36`.  The striking readout is that the weak
primitive-phase ladder reproduces HYP-2972's bounded-rational witness pressure,
while the stricter primitive packet separates AP/GW as the only rows in this
stress bank with no strict primitive witness by `q<=42`.

Tournament Analysis over proof carriers is transitive:

```text
labelled_packet_sheaf > toeplitz_fejer_dual > spectral_shadow
  > ramanujan_danger_energy > primitive_phase_packet
  > carmichael_autocorrelation > c14_endpoint_sum_trace
  > raw_qdiv > raw_runner_residues
```

Fingerprint: `score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}`,
`directed_3_cycles=0`, `hamiltonian_paths=1`.

## Proof Target After The Audit

The useful theorem is not "Ramanujan profile determines the route".  HYP-2978
and the q=14 packet collisions refute that.  The useful theorem is a packet
handoff:

```text
After current Moon-core reductions, every LRC14 residual either
  (1) has a primitive weak or strict q-phase packet certifying the row;
  (2) has a negative Toeplitz/Fejer or positive spectral-shadow dual;
  (3) has an AP/GW c_14 endpoint-sum boundary packet;
  (4) has a Ramanujan danger-energy defect that forces a labelled handoff; or
  (5) carries the K33/HYP-2908/THM-572 state-lift debt.
```

Next task: extend the strict primitive-witness audit from the AP-neighborhood
stress bank to the full HYP-2963 bank, then interval-enclose the finite list of
late primitive packets (`q=25,27,34,40,41`) so that they become exact
certificates rather than just computed witnesses.

Assumption challenge: the tournament vertices should not be fixed too early.
Candidate vertices are speeds, residues, primitive denominator phases,
denominators, Ramanujan modes, endpoint owner pairs, taut-current atoms,
Ramanujan-subspace components, and proof obligations. This route will start
with proof obligations / exact-period modes because they preserve the LRC
predicate "primitive weak witness, strict open mass, or boundary-only equality"
and destroy raw runner ownership unless endpoint/C27/K33 labels are reattached.
