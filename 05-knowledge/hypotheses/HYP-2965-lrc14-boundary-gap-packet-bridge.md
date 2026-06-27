---
id: HYP-2965
title: LRC14 boundary-gap packet bridge
status: PROOF-BRIDGE / exact local certificate for the covering branch; not a proof of LRC14
source: codex-2026-06-24-S152
related:
  - HYP-2964
  - HYP-2963
  - HYP-2962
  - HYP-2961
  - HYP-2956
  - HYP-2955
  - HYP-2954
  - HYP-2953
  - HYP-2951
  - HYP-2950
  - HYP-2949
  - HYP-2947
  - HYP-2908
  - THM-523
  - THM-566
  - THM-572
  - OPEN-Q-108
results:
  - 04-computation/lrc14_boundary_gap_packet_bridge_codex_s152.py
  - 05-knowledge/results/lrc14_boundary_gap_packet_bridge_codex_s152.out
---

# HYP-2965: LRC14 Boundary-Gap Packet Bridge

## Claim

The F6 covering kernel from HYP-2956/HYP-2962/HYP-2963 should be expressed as
an exact endpoint-collision theorem.

For a finite row `S`, the strict safe set

```text
O(S) = { t : ||s*t|| > 1/14 for every s in S }
```

is the complement of a finite union of rational danger arcs.  Every positive
component `(lo, hi)` of `O(S)` is an exact **boundary bridge** between two
adjacent danger-arc endpoints, with rational length `hi-lo`.  Therefore:

```text
positive boundary bridge  =>  M(S) > 1/14
strict counterexample     =>  qdiv(S)>14 and every boundary bridge is pinched
```

The covering branch is not primarily a row search.  It is an endpoint-collision
problem with labels: q-covering clocks, endpoint owners, source-spectrum
identity, C27/unital packets, K33/state-lift flags, and the gK8/L_y moment
image.

Rebase note: HYP-2964 now names the moon-core proof skeleton after
apex-majority elimination.  HYP-2965 is the local F6/covering-face complement:
inside that moon core, a zero-open covering residual must be a zero
boundary-bridge packet.

## Computation

Script:

```text
04-computation/lrc14_boundary_gap_packet_bridge_codex_s152.py
```

Stored output:

```text
05-knowledge/results/lrc14_boundary_gap_packet_bridge_codex_s152.out
```

Default bank matches the HYP-2963 bounded classifier scale:

```text
single_limit=180
two_swap_limit=36
alias_depth=4
lcm_tail_max=5
```

Readout:

```text
qdiv>14 covering rows audited      = 1187
positive strict-open rows          = 1187
zero-open covering rows            = 0
rows with zero net endpoint current= 1187
smallest safe_mu                   = 1543/294294 at single swap 6->98
smallest max bridge                = 1/728 at single swap 6->98
```

Thus the bounded covering bank has no F6 row.  More importantly, every audited
covering row has zero net first endpoint current.  The boundary-moment bridge
cannot be a raw divergence proof.  It has to be a localized transition-packet
or second-moment proof, consistent with the gK8/L_y direction.

## Boundary-Gap Packet

For each strict safe component `(lo, hi)`, attach:

```text
length            = hi-lo
left_end_owners   = speeds whose danger arc ends at lo
right_start_owners= speeds whose danger arc starts at hi
transition label  = left owner packet -> right owner packet
```

Example from the lcm-tail covering row `12->84`:

```text
safe_mu = 563/105105
components = 8

[29/70, 163/392]       len=3/1960   5:g1 -> 0:g14
[491/1176, 41/98]      len=1/1176   0:g14 -> 7:g7
[33/392, 13/154]       len=1/4312   0:g14 -> 11:g1
[15/182, 97/1176]      len=1/15288  13:g1 -> 0:g14
```

The components occur in paired transitions, so a first-current sum cancels.
The retained object must be the whole transition multiset with lengths.

## Theorem Target

The useful bridge theorem is:

```text
Boundary-Gap Packet Theorem.

Let S be a primitive 13-speed row with qdiv(S)>14.
If every strict boundary bridge is pinched, then the localized endpoint-owner
transition packet has either:

  (A) positive gK8/L_y moment image, producing a strict source interval; or
  (B) a nonunit K33/H=7 state-lift label, contradicted by THM-572 once lifted.
```

This would close the F6 covering zero-open non-migration kernel in HYP-2956.

## Relation To HYP-2963

HYP-2963 classifies the bounded bank into q-witness, AP/GW boundary,
unit-petal, K33/state-lift, and covering-moment routes, with no unknown packet.
HYP-2965 zooms into the covering-moment route and replaces "positive
Haar-open" with exact rational bridge certificates.

In HYP-2963 terms:

```text
COVERING-MOMENT positive row
  -> one or more positive boundary-gap packets

SOURCE-SPECTRUM-UNKNOWN strict candidate
  -> qdiv>14 plus no positive boundary-gap packet
```

Thus the unknown bucket can be searched for as a zero bridge-packet row.

## Assumption Challenge

Candidate vertex sets considered:

```text
runners,
danger arcs,
safe components,
component endpoints,
endpoint owner events,
divisor clocks,
C27 shells,
K33 incidence,
exact-period residues,
Fourier/moment modes,
fixed-margin fibers,
proof obligations.
```

Chosen vertices:

```text
proof coordinates of a covering boundary-gap packet.
```

Preserved LRC predicate:

```text
qdiv>14 candidate status plus exact strict safe intervals;
any positive bridge length is a rigorous witness that M(S)>1/14.
```

Destroyed information:

```text
raw row identity once the endpoint owners and transition lengths have been
attached.
```

Challenged assumption:

```text
the covering branch is not a raw runner tournament or first-current problem.
The exact audit shows all covering rows have zero net endpoint current, so the
bridge must use localized transitions or second moments.
```

## Tournament Analysis

Vertices:

```text
q_cover_gate
exact_boundary_gap
endpoint_owner_current
source_spectrum_pullback
C27_unital_packet
K33_state_lift_packet
gK8_Ly_moment_image
fixed_margin_swap_fiber
raw_runner_row
```

Pair observable:

```text
which coordinate preserves the covering LRC predicate and prevents a positive
boundary gap from being scalarized away.
```

Switch/gauge:

```text
lexicographic retention of q-covering, exact gap, endpoint owner, source
identity, state-lift visibility, moment image, and anti-scalar guard.
```

Hamiltonian path:

```text
source_spectrum_pullback
> exact_boundary_gap
> endpoint_owner_current
> gK8_Ly_moment_image
> K33_state_lift_packet
> C27_unital_packet
> fixed_margin_swap_fiber
> q_cover_gate
> raw_runner_row
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}
directed_3cycles=0
SCCs=nine singleton SCCs
Hamiltonian paths=1
```

## Status

This is proof progress, not a proof of LRC14.  It promotes F6 from an informal
"zero-open covering kernel" into a finite endpoint-collision packet theorem.
The next hard step is to prove that a zero bridge packet has positive localized
gK8/L_y image or constructs the K33/H=7 state-lift.
