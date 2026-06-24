---
id: HYP-2981
title: LRC14 Robbins-Robin bridge guardrails and interval Fejer packet certificates
status: PROOF-INTERFACE / named-row interval-certificate scaffold; not a proof
source: codex-2026-06-24-S162
related:
  - HYP-2979
  - HYP-2978
  - HYP-2977
  - HYP-2976
  - HYP-2974
  - HYP-2973
  - HYP-2965
  - HYP-2963
  - HYP-2956
  - HYP-2908
  - THM-572
  - OPEN-Q-108
artifacts:
  - 04-computation/lrc14_packet_fejer_interval_scaffold_codex_s162.py
  - 05-knowledge/results/lrc14_packet_fejer_interval_scaffold_codex_s162.out
  - 07-reflections/lrc14-robin-robbins-fejer-interval-scaffold-codex-s162.md
---

# HYP-2981: LRC14 Robbins-Robin Bridge Guardrails And Interval Fejer Packet Certificates

This hypothesis records the 2026-06-24 follow-up route combining:

- Robbins' graph theorem: an undirected graph has a strongly connected
  orientation exactly when it is connected and has no bridge.
- Robin's number-theory theorem from the divisor-function page: the divisor
  inequality for `sigma(n)` past `5040` is equivalent to RH.
- Neville Robbins-style divisor/partition side readings: divisor functions
  appear as fixed-extreme partition fibers and cyclotomic partition products.
- HYP-2974's current primary proof obligation: replace floating Fejer
  quadratic-form signs by rigorous interval certificates anchored to labelled
  packet `P(S)` fibers.

The spelling collision is useful but should not be blurred: the graph theorem
is Robbins, while the divisor-function RH criterion is Robin.  The common
mathematical lesson is the same controlled-kernel rule from HYP-2978:

```text
connection = map + declared kernel + certificate for every forgotten bridge.
```

In graph Robbins, a bridge is exactly the obstruction to orienting all edges
while preserving two-way reachability.  In Robin's divisor criterion, a possible
RH violation must pass through the extremal divisor-density kernel
`sigma(n)/(n log log n)`, with superabundant/colossally abundant numbers as the
load-bearing fiber rather than arbitrary integers.  In LRC14, the analogue of a
bridge is a packet fiber where a coarse quotient identifies rows with different
proof routes.

## Interval Certificate Target

For a row `S`, HYP-2974 uses

```text
F_S(t)=sum_{v in S} 1_{||v t||<1/14}-1.
```

A strict counterexample would have `F_S >= 0` almost everywhere, so every
Toeplitz moment section is PSD.  The Fejer full-bank audit found, for every
positive HYP-2963 packet-bank row, an explicit rational center `x` and degree
`d` with floating

```text
Q_d(S,x)=c_0+2*sum_{1<=k<=d}(1-k/(d+1))*c_k*cos(2*pi*k*x) < 0,
```

where

```text
c_0=6/7,
c_k=sum_{v in S, v|k} sin(pi*(k/v)/7)/(pi*(k/v)).
```

HYP-2981 makes the remaining obligation exact:

```text
For each labelled packet fiber P(S), attach a certificate tuple
  (family, route, safe component [a,b], center x, degree d, interval Q)
with upper(Q)<0.
```

The certificate is admissible only if it states what it forgets:

```text
preserves:     cover => F_S>=0 => Toeplitz PSD necessary condition
forgets:       raw runner names and most endpoint-owner data
reattaches:    labelled packet family, exact safe component, center, degree,
               divisor-curried coefficient fiber, and residual route
residuals:     AP/GW equality atoms; otherwise refine interval, center, or
               packet label until upper(Q)<0 or route to HYP-2908/THM-572.
```

Thus the Fejer certificate is a Robbins-style orientation of the packet proof
graph: every non-bridge packet receives a directed dual exit; AP and GW are the
only allowed unoriented bridge/equality atoms.

## Tournament-Adjacent Use

Tournament vertices should be proof obligations, not runners:

```text
interval Fejer certificate,
labelled packet fiber P(S),
Toeplitz PSD cone,
Ramanujan exact-period projector,
danger-count moment dual,
endpoint taut bridge,
K33/state-lift debt,
raw divisor or runner quotient.
```

Pairwise observable:

```text
A -> B if A preserves the cover/noncover predicate while discharging more
route-mixing fibers that B leaves as bridges.
```

Tie Hamiltonian path:

```text
interval Fejer certificate
> labelled packet fiber
> Toeplitz PSD cone
> Ramanujan exact-period packet
> moment/endpoint duals
> state-lift debt
> raw quotient.
```

The immediate computation is
`04-computation/lrc14_packet_fejer_interval_scaffold_codex_s162.py`, a named-row
interval certificate scaffold for the hardest floating Fejer rows: `K33
12->36` at degree `159` and `P10+GW` at degree `280`, plus smaller-margin
named rows.  Its stored output certifies interval upper endpoints `<0` for all
five selected rows and attaches each row to its labelled packet fiber.  This is
not yet the production certificate emitter: the current rational `pi` enclosure
and Taylor interval engine must be replaced by a formally sourced backend, and
the final pass should group the HYP-2963 bank by packet family rather than
treating `21911` rows individually.
