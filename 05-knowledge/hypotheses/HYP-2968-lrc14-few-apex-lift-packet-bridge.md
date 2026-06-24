---
id: HYP-2968
title: LRC14 few-apex lift-packet bridge
status: PROOF-TARGET / exact bounded lift-packet evidence for the covering boundary-moment branch; not a proof
source: codex-2026-06-24-S152
script: 04-computation/lrc14_few_apex_lift_packet_probe_codex_s152.py
result: 05-knowledge/results/lrc14_few_apex_lift_packet_probe_codex_s152.out
related:
  - HYP-2967
  - HYP-2966
  - HYP-2965
  - HYP-2964
  - HYP-2963
  - HYP-2962
  - HYP-2961
  - HYP-2956
  - HYP-2955
  - HYP-2954
  - HYP-2953
  - HYP-2929
  - HYP-2911
  - HYP-2909
  - HYP-2908
  - THM-571
  - THM-570
  - THM-568
  - THM-523
  - THM-572
  - OPEN-Q-108
---

# HYP-2968: LRC14 Few-Apex Lift-Packet Bridge

## Claim

The LRC14 covering residual should be split into two labelled packets:

```text
S = R union 14Q
```

1.  `|Q|>=7`: the apex-majority branch is already closed by THM-571,
    modulo the accepted below-frontier LRC input used there.
2.  `1<=|Q|<=6`: the active covering branch is a **few-apex lift packet**.

The few-apex theorem target is:

```text
primitive qdiv(S)>14, S=R union 14Q, 1<=|Q|<=6
  -> some labelled lift packet has positive regular-open Haar mass
     or the residual state-lifts to the K33 / THM-572 endpoint.
```

This corrects the older shorthand in HYP-2961 that still listed
`|Q|>=7` as live.  The actual live covering side is the `|Q|<=6`
scale-separated / finite-core package.

## Exact Lift Packet

Set `u=14t`.  For each `u`, the original time circle has fourteen lifts:

```text
t = (u+k)/14,  k=0,...,13.
```

For a multiple `14m in 14Q`,

```text
||14m t|| = ||m u||.
```

For a residual speed `r in R`, the danger set in lift `k` is exactly:

```text
||r(u+k)/14|| < 1/14
  <=> exists n in Z with |r(u+k)-14n| < 1
  <=> u in ((14n-1)/r-k, (14n+1)/r-k) cap [0,1].
```

Thus each row emits fourteen finite rational Borel/Baire interval fronts.
A positive complement interval in any lift is a strict LRC14 witness.

The quotient preserves:

```text
qdiv>14,
the Q/R split,
|14Z cap S|,
exact rational lift-safe mass,
whether a strict Haar/Baire witness exists.
```

It destroys fine residual-interval ownership and the exact maximizing Farey
node, so the tightest rows are checked against direct Haar intervals and, when
needed, exact `M`.

## S152 Computation

Script:

```text
04-computation/lrc14_few_apex_lift_packet_probe_codex_s152.py
```

Stored output:

```text
05-knowledge/results/lrc14_few_apex_lift_packet_probe_codex_s152.out
```

The bank is structured rather than exhaustive: AP-replacement rows with
`1<=|14Z cap S|<=6`, compact and stretched multiplier profiles, and `qdiv>14`
forced by restoring every divisor witness `d=2..14`.

Default stored run:

```text
audited qdiv>14 rows          8190
zero strict lift mass         0
no positive lift              0
smallest safe_t mass          563/105105
```

By `k14=|14Z cap S|`:

```text
k14  mode       rows   min safe_t        min open lifts
1    compact      13   563/105105        2
1    stretched    13   7/858             2
2    compact      78   4999/315315       4
2    stretched    78   9577/420420       4
3    compact     286   151/4620          6
3    stretched   286   46363/1261260     6
4    compact     715   126767/2522520    6
4    stretched   715   681/12740         6
5    compact    1287   178999/2522520    8
5    stretched  1287   180017/2522520    8
6    compact    1716   26366/315315      8
6    stretched  1716   21173/252252      8
```

The tightest direct-Haar checks all matched the lift-packet mass.  Optional
exact `M` fallback on the four tightest rows gave:

```text
drop(12)->14*6       M=7/89
drop(6)->14*1        M=2/23
drop(6)->14*2        M=2/23
drop(12)->14*12      M=14/173
```

All exceed `1/14`.

## Tournament Analysis

Vertices are proof carriers, not runners:

```text
qdiv scalar witness
THM-571 apex-majority descent
few-apex lift packet
exact M/Farey fallback
boundary-moment bridge
K33/state-lift endpoint
raw apex residue tournament
raw runner set
```

Pair observable:

```text
LRC predicate retention,
Q/R split retention,
lift interval exactness,
boundary-label retention,
state-lift compatibility,
anti-scalar guard.
```

The stored tournament is transitive:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
c3=0
scc=(1,1,1,1,1,1,1,1)
hp=1
```

Hamiltonian path:

```text
qdiv scalar witness
> THM-571 apex-majority descent
> few-apex lift packet
> exact M/Farey fallback
> boundary-moment bridge
> K33/state-lift endpoint
> raw apex residue tournament
> raw runner set
```

## Proof Use

HYP-2968 turns the covering bucket F5/F6 from HYP-2956/HYP-2963 into an
exact interval theorem target:

```text
covering zero-open non-migration
  -> no positive lift interval in all fourteen packets.
```

The structured bank found none.  The next proof should try to promote the lift
count to a theorem:

```text
If qdiv(S)>14 and 1<=|14Z cap S|<=6, then the Q-safe set in u cannot be
covered by the fourteen residual-lift danger ledgers unless a retained
nonunit packet state-lifts to the HYP-2908 / THM-572 obstruction.
```

This is the boundary-moment bridge in a labelled-packet form: the proof should
not compare raw runners, raw residues, or raw scalar safe mass until the Q/R
lift label and possible K33 endpoint have been retained.
