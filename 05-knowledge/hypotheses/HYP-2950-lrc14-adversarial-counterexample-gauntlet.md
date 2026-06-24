---
id: HYP-2950
title: LRC14 adversarial counterexample gauntlet
status: PROOF-DISPROOF INTERFACE / bounded adversarial gauntlet; not yet a proof
source: codex-2026-06-24-S148
related:
  - HYP-2952
  - HYP-2951
  - HYP-2949
  - HYP-2948
  - HYP-2947
  - HYP-2946
  - HYP-2944
  - HYP-2940
  - HYP-2937
  - HYP-2908
  - THM-572
  - OPEN-Q-108
results:
  - 04-computation/lrc14_adversarial_counterexample_gauntlet_codex_s148.py
  - 05-knowledge/results/lrc14_adversarial_counterexample_gauntlet_codex_s148.out
---

# HYP-2950: LRC14 Adversarial Counterexample Gauntlet

S148 searches for LRC14 counterexample-shaped rows while retaining the current
packet language:

```text
exact M/Farey branch
-> regular-open Haar/Baire safe mass
-> C27 owner/carry shell code
-> q=3 unital / affine-depth packet
-> K33/state-lift flag
-> PH-style rank
```

The test should deliberately include adversarial families:

```text
AP/GW/petal/K33 rows,
S138 two-swap frontier rows,
divisor-loaded lcm tails,
GW floor-odd iso impostors,
random covering-repair rows,
and packet-preserving mutations.
```

The proof target is a gauntlet statement:

```text
Every counterexample candidate either
  (a) has M > 1/14 by exact witness,
  (b) is endpoint-only AP/GW tight,
  (c) is a unit-petal discharge,
  (d) or preserves a K33/state-lift address for HYP-2908/THM-572.
```

The disproof target is equally explicit:

```text
find any primitive 13-speed row with M < 1/14,
or find a threshold-safe strict-Haar-zero row that is not AP/GW,
or find a low packet whose labels evade the S145 rank split.
```

## Computation

The script
`04-computation/lrc14_adversarial_counterexample_gauntlet_codex_s148.py`
stores output at
`05-knowledge/results/lrc14_adversarial_counterexample_gauntlet_codex_s148.out`.

Default run:

```text
single_limit=360
two_swap_limit=42
alias_depth=6
lcm_tail_max=8
```

## Main Audit

Named adversaries and shell-alias impostors:

```text
audited = 42
named counterexamples = 0
named low unknown packets = 0
```

AP single-swap scan through replacement `v<=360`:

```text
raw rows            = 4512
direct q-safe       = 2175
hard q>=14          = 2337
below threshold     = 0
tight hard rows     = 2
M<=3/41             = 3
M<=2/27             = 5
```

AP two-swap scan through added values `<=42`:

```text
raw rows            = 32046
direct q-safe       = 22128
hard q>=14          = 9918
below threshold     = 0
tight hard rows     = 2
M<=3/41             = 3
M<=2/27             = 7
```

The only hard tight rows in both scans are AP and Goddyn-Wong.  The best
single-swap and two-swap hard rows are both AP with `M=1/14`; no row with
`M<1/14` appears.

## Low Frontier Stability

The two-swap scan extends the S145 ceiling from added values `<=40` to `<=42`.
The low frontier is unchanged:

```text
M<=3/41: AP, GW, near/K33 12->36
M<=2/27: those three, P10, P13, P10+GW, P10+K33
```

Thus the HYP-2947 measurable rank split remains stable in this larger bank:

```text
rank 0: AP, GW, P10, P13, P10+GW
rank 1: K33, P10+K33
unknown low packets: 0
```

## Disproof Pressure

The gauntlet also tests false proof quotients.  Shell aliases repeat the same
coarse C27 transfer labels but become safely loose:

```text
GW-shell aliases:   best M=11/137, worst M=1/12
K33-shell aliases:  best M=12/149, worst M=1/12
P10-shell aliases:  M=1/10
P13-shell aliases:  M=1/13
```

The floor-odd Goddyn-Wong tournament impostor

```text
{1,2,3,4,5,6,7,8,9,10,11,13,360}
```

is also safely loose:

```text
M=6/73
strict_safe_mu=6667/630630
```

Readout: shell labels and tournament shadows are useful only as retained
addresses.  They do not characterize tightness or counterexamples.

## Tournament Analysis

Tournament vertices are proof/disproof routes:

```text
exact M counterexample search
q-threshold direct witnesses
regular-open Haar/Baire mass
S145 packet rank split
K33 state-lift endpoint
C27 shell labels alone
raw tournament/impostor iso
raw scalar count
```

Pair observable:

```text
exactness,
branch retention,
measurable witness retention,
packet fit,
endpoint closure,
scalar-decoy resistance.
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
c3=0
singleton SCCs
Hamiltonian paths=1
```

Hamiltonian path:

```text
exact M counterexample search
> q-threshold direct witnesses
> regular-open Haar/Baire mass
> S145 packet rank split
> K33 state-lift endpoint
> C27 shell labels alone
> raw tournament/impostor iso
> raw scalar count
```

## Counterexample Constraints

No disproof row was found in this gauntlet.  Any future counterexample must
evade all of the following:

1. all `q<=13` direct witnesses, so it must have `q_threshold>=14`;
2. AP/GW endpoint-only tight classification in the tested single/two-swap banks;
3. S145 rank split for every tested `M<=2/27` low packet;
4. HYP-2948/HYP-2949 regular-open safe-mass separation;
5. shell-alias and tournament-impostor decoys, which are loose under exact `M`.

The remaining theorem is:

```text
Every primitive LRC14 counterexample reduces to one of these gauntlet families,
or to a packet-preserving mutation that still enters HYP-2947.
```

## Rebase Connection: Boundary And Derived-Tournament Filters

Concurrent HYP-2951 and HYP-2952 add two compatible filters to this gauntlet.

HYP-2951 gives a Haar-Baire boundary-owner finite check:

```text
one-swap added<=160: only AP/GW are boundary-only
two-swap added in 14..40: no boundary-only rows
smallest positive masses: 1/1260, 1/980, 4/2205
```

HYP-2952 gives a derived-boundary tournament necessary-condition stack:

```text
AP/GW-kind boundary candidates should pass the transitive apex-pressure class
and first-derived AP profile before C27/Farey labels route them.
```

The concurrent `comms/POKE-COORDINATION.md` update promotes that HYP-2952
filter to the current cluster checkpoint.  In that reading, this HYP-2950
gauntlet is the adversarial backstop: try to falsify the apex-pressure /
Jacobsthal-repair front filter by exact `M`, Haar/Baire mass, shell aliases,
and two-swap packet scans before trusting it as a proof interface.

Together with S148, the emerging proof ladder is:

```text
exact M gauntlet
-> Haar/Baire boundary-owner rigidity
-> derived AP/GW tournament filter
-> HYP-2947 measurable rank packet
-> HYP-2908/THM-572 state-lift endpoint.
```
