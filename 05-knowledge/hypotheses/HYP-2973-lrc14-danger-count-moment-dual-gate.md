---
id: HYP-2973
title: LRC14 danger-count moment-dual gate
status: PROOF-TARGET / exact count-distribution dual certificate; not a proof
source: codex-2026-06-24-S158
script: 04-computation/lrc14_danger_count_moment_dual_codex_s158.py
result: 05-knowledge/results/lrc14_danger_count_moment_dual_codex_s158.out
related:
  - HYP-2968
  - HYP-2967
  - HYP-2966
  - HYP-2965
  - HYP-2964
  - HYP-2963
  - HYP-2962
  - HYP-2961
  - HYP-2956
  - HYP-2949
  - HYP-2948
  - HYP-2908
  - THM-572
  - OPEN-Q-108
---

# HYP-2973: LRC14 Danger-Count Moment-Dual Gate

## Claim

Let

```text
N_S(t) = #{s in S : ||s t|| < 1/14}.
```

Then `safe_mu(S)=mu(N_S=0)`.  If `S` were a strict LRC14 counterexample,
then `N_S(t)>=1` for every `t`.  Therefore any polynomial `P` on
`{0,...,13}` satisfying

```text
P(0)=1
P(n)<=0 for n=1,...,13
```

gives an exact dual certificate:

```text
safe_mu(S) >= E[P(N_S)].
```

Using the factorial basis

```text
P(n)=1 + sum_k y_k binom(n,k),
```

the lower bound depends only on the danger-count moments
`E[binom(N_S,k)]`.  A positive value of `E[P(N_S)]` proves a positive safe set
without locating the safe interval, without choosing a denominator-14 apex,
and without using the `u=14t` lift packet.

## Evidence

Script:

```text
04-computation/lrc14_danger_count_moment_dual_codex_s158.py
```

Stored output:

```text
05-knowledge/results/lrc14_danger_count_moment_dual_codex_s158.out
```

Default run:

```text
python3 04-computation/lrc14_danger_count_moment_dual_codex_s158.py \
  --max-degree 13 --one-swap-limit 80
```

The audit computes the exact distribution of `N_S` by sweeping rational danger
arc endpoints, then enumerates exact polynomial-dual LP vertices.  The feasible
dual vertex counts by degree are:

```text
{1: 1, 2: 1, 3: 11, 4: 10, 5: 45, 6: 36,
 7: 84, 8: 56, 9: 70, 10: 35, 11: 21, 12: 6, 13: 1}
```

Named-row readout:

```text
AP             safe_mu=0        first positive dual degree: none
GW 12->24      safe_mu=0        first positive dual degree: none
12->36         safe_mu=1/1260   first positive dual degree: 9
10->20         safe_mu=1/980    first positive dual degree: 9
13->26         safe_mu=1/182    first positive dual degree: 8
P10+GW         safe_mu=1/980    first positive dual degree: 9
P10+K33        safe_mu=4/2205   first positive dual degree: 9
12->84         safe_mu=563/105105 first positive dual degree: 8
12->168        safe_mu=263/30030 first positive dual degree: 7
6->14          safe_mu=7/858    first positive dual degree: 7
6->28          safe_mu=7/858    first positive dual degree: 7
```

One-swap AP bank through `add<=80`:

```text
audited rows: 871
zero safe-measure rows: 1
positive safe-measure rows: 870
first positive dual-degree histogram:
  {5: 109, 6: 345, 7: 400, 8: 14, 9: 2, none: 1}
rows certified by degree <=9 dual: 870
```

The lone zero row is the AP/Goddyn-Wong equality family.  The worst degree-9
lower bound among positive rows is `10753/28108080` at `drop(10)->20`, whose
exact safe mass is `1/980`.

## Negative Signal

This is not a cheap second-moment proof.  The named hard rows remain invisible
through degree `6`.  The count-only dual begins separating positive rows from
AP/GW at degrees `7`, `8`, and `9`, and degree `13` recovers exact safe mass by
interpolation.  Thus any global proof from this route must either:

1.  prove a structured degree-8/9 dual family on fixed-margin labelled packets,
2.  show AP/GW are the only packets defeating every low-degree count dual, or
3.  route the remaining count-indistinguishable packets to C27/K33/state-lift
    structure.

## Tournament Analysis

Vertices are proof carriers, not runners:

```text
full danger-count distribution
degree-13 interpolating dual
degree-9 moment dual
endpoint/lift packet
gK8/Delsarte sector moment
Bonferroni pair moment
raw Haar safe measure
raw runner set
```

Pairwise observable:

```text
cover predicate retention, exact count distribution, low-degree certificate
strength, AP/GW equality visibility, finite-atlas fit, anti-scalarization
```

The conservative carrier tournament is transitive in the audit, with
Hamiltonian path:

```text
full danger-count distribution
> degree-13 interpolating dual
> degree-9 moment dual
> endpoint/lift packet
> gK8/Delsarte sector moment
> Bonferroni pair moment
> raw Haar safe measure
> raw runner set
```

## Proof Target

Prove a **danger-count moment gate**:

```text
Every primitive LRC14 packet outside AP/GW has a positive
degree <= 9 count-dual expectation, unless it carries a labelled
C27/K33/HYP-2908/THM-572 state-lift obstruction.
```

This would give a proof branch that is genuinely different from HYP-2968:
HYP-2968 finds explicit lift intervals; HYP-2973 forgets the geometry and asks
whether the integer multiplicity distribution already forbids full cover.
