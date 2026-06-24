# LRC14 certificate handoff atlas: the proof is a zipper now

S164 packages the current LRC14 proof stack as a certificate handoff atlas:

```text
script: 04-computation/lrc14_certificate_handoff_atlas_codex_s164.py
output: 05-knowledge/results/lrc14_certificate_handoff_atlas_codex_s164.out
hypothesis: HYP-2986
```

The Tournament Analysis vertices are not runners.  They are proof carriers:

```text
labelled source packet
exact interval Fejer certificate
Ramanujan exact-period projector
endpoint bridge graph
twist-ladder blocker
danger-count moment dual
analytic-sieve/Kaczynski packet
tournament state lift
raw scalar shadow
```

The pairwise observable is retained proof payload:

```text
LRC predicate, exact scale, phase/period, topology, endpoint owners,
packet family, dual certificate, formal checkability, residual routing.
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,4:3,6:1,7:1,8:1}
directed_3cycles=1
SCC_sizes=[1,1,3,1,1,1,1]
Hamiltonian_path_count=3
```

The nontrivial SCC is the point: Ramanujan exact-period packets, endpoint
bridges, and twist ladders can each be the right carrier depending on which
predicate is preserved.  A handoff has to say what it keeps and what it is
allowed to forget.

## Zipper theorem target

If the following six arrows are proved, LRC14 is done in this framework:

```text
O1 source-kernel exclusion:
   every primitive row enters the atlas.
O2 formal interval backend:
   prototype Fejer/trig intervals become checkable balls.
O3 family compression:
   selected Fejer/twist certificates lift to packet templates.
O4 admissible smoothing:
   analytic-sieve/Kaczynski quotients retain approach labels.
O5 state-lift construction:
   zero-open non-AP/GW residuals build HYP-2908/THM-572.
O6 F7 definition:
   any Johnson-harmonic residual sector is named with its preserved predicate.
```

Then every primitive LRC14 row either has a strict witness/dual certificate,
is the AP/Goddyn-Wong equality atom, or constructs the forbidden tournament
state lift.

## Fixed-margin transfer

The useful import from `arXiv:2606.22636` is structural: fixed margins are not
statistics to average away; they are fibers.  The proof pattern is:

```text
preserve margin fibers
reduce to a low-row heat-bath core
split count sectors from Johnson harmonic sectors
```

LRC14 translation:

```text
preserve packet fibers
reduce arbitrary primitive rows to a finite labelled packet core
split ordinary count/moment sectors from a named F7 harmonic residual
```

So F7 should not be an anonymous "everything else" bucket.  It has to be a
labelled harmonic residual sector whose predicate is strong enough either to
empty it or to build the THM-572 state lift.

Best next attacks: O3 family compression for the Fejer certificates, O5
state-lift construction for zero-open non-AP/GW residuals, and O6 as a
Johnson-sector predicate.
