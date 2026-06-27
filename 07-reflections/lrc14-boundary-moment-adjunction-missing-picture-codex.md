# LRC14 Boundary-Moment Adjunction: Missing Picture Attempt

This pass searched the recent LRC14 packet/tournament work against older
gK8, exact-period, and boundary-currency artifacts.  The missing picture I see
is not another invariant.  It is a bridge:

```text
adaptive exact-period boundary packets
    -> signed boundary/depth quotient
    -> gK8 / L_y moment dual
```

If this bridge can be made into an honest monotone chain map, LRC14 has a
plausible final proof shape.

## Three Regimes

Let

```text
qdiv(S) = min{d >= 2 : no speed in S is divisible by d}.
```

The q-witness splits primitive 13-row candidates:

```text
qdiv < 14:  immediate loose witness, M(S) >= 1/qdiv > 1/14.
qdiv = 14:  endpoint boundary core, where AP/GW live.
qdiv > 14:  covering residual, the real LRC14 danger zone.
```

The recent AP/GW census work is a classification of the middle branch.  It
does not by itself close the covering residual.  The older gK8/L_y work is a
way to kill the covering residual, but it lives after scalarization into
miss-count moments.  The exact-period packet atlas sits between them.

## What The History Says

The useful artifacts line up as follows.

```text
HYP-2951:
  In the AP one-swap neighborhood up to added value 160, only AP and GW are
  strict-Haar-zero boundary-only rows.  Two-swap AP rows are all positive-open.

HYP-2952:
  Apex-pressure tournaments on U14 realize only 4 of 56 six-vertex classes in
  AP single-residue mutations.  AP, GW, K33-near, petals, and splices are all
  transitive, so this is a front filter, not a classifier.

HYP-2950:
  The adversarial gauntlet found no counterexample in named rows, shell aliases,
  lcm tails, AP single-swaps through 360, or AP two-swaps through added value
  42.  It also refuted shell-label-only and raw tournament-impostor shortcuts:
  those aliases are loose under exact M.

HYP-2886:
  Fixed finite denominator bases are false closures.  The correct discrete
  object is an adaptive exact-period unit-packet atlas, with mod-7, chi_7,
  affine-pair, and CRT-defect labels retained.

HYP-2704:
  The seven-sector survival currency is a boundary/depth quotient.  The
  death-chain recurrence turns the high p6 debt positive after two decorrelated
  far hits.

HYP-2840 and the gK8/L_y route:
  The wide branch should be closed by far-element decorrelation and
  consecutive extremality in the moment dual, not by a divergent spectral
  envelope.

HYP-2666:
  The endpoint tax is ordered: shell gate first, then p1/missed-sector tax.
  Scalarizing before paying the packet gate creates false constants.
```

These are not competing stories.  They are the same boundary currency at
different resolutions.

## Proposed Bridge

For a row `S` and denominator `D`, let `P_D(S)` be the exact-period unit-packet
hypergraph:

```text
vertices:  a in (Z/DZ)^*
runner s forbids a when 14*min(sa mod D, D-sa mod D) < D.
```

An LRC14 counterexample covers every exact-period packet in every relevant
scaled denominator chart.  A minimal such cover should have a boundary class:

```text
beta_D(S) = minimal signed boundary of the packet cover.
```

The missing map is a signed quotient

```text
B_D: beta_D(S) -> seven-sector missed-depth currency
```

with output in the same space as HYP-2704/gK8:

```text
p0,...,p6
C = p1+p2+p3+p4 - 4p6
L_y = 10q_0 + q_3 + 10q_6
```

Here `B_D` should forget exact denominator labels but retain enough boundary
orientation to preserve the LRC predicate:

```text
positive regular-open witness
versus
boundary-only packet
versus
covered residual.
```

The old exact-period atlas says where `beta_D` comes from.  HYP-2704 says the
right scalar quotient is a boundary/depth character, not raw count.  gK8 says
the moment image should be positive away from the bounded core.

## The Kernel Claim

The bridge should have a tiny kernel on zero-open packet boundaries:

```text
ker(B_D) among zero-open packet boundaries
  = AP/GW source orbit
    plus explicitly named K33/state-lift debts.
```

This is the connection between the new AP/GW tournament census and the old
wide proof.

If `qdiv(S)=14`, the boundary-support part of the kernel is allowed and gives
the tight AP/GW census.  If `qdiv(S)>14`, the AP/GW kernel is forbidden by
divisibility: the row has a multiple of 14 and is in the covering residual,
not the denominator-14 boundary orbit.  Then a minimal counterexample would
have to sit in the K33/state-lift kernel or in a genuinely new zero-open
kernel.  That is exactly where HYP-2908/THM-572 and the HYP-2950 gauntlet
should be aimed.

## LRC14 Proof Skeleton

The theorem target becomes:

```text
Let S be a primitive 13-speed row.

1. If qdiv(S) < 14, the q-witness gives M(S) > 1/14.

2. If qdiv(S) = 14, the boundary-source-core theorem says:
     S is AP/GW, or S has positive regular-open witness, or S has petal/K33
     escape.  Hence M(S) >= 1/14.

3. If qdiv(S) > 14, choose the adaptive exact-period packet cover of S and
   its minimal boundary beta.  Apply B_D.  Either:
     a. B_D(beta) has positive gK8/L_y slack, giving a witness;
     b. beta lies in AP/GW kernel, impossible because qdiv>14;
     c. beta lies in named K33/state-lift kernel, discharged by HYP-2908/THM-572.
```

Everything difficult is in step 3, but it is now one missing bridge rather than
many unrelated analogies.

## Four Lemmas To Prove Or Falsify

1. **Minimal Packet Boundary Lemma.**

```text
Every primitive qdiv>14 counterexample has an adaptive exact-period chart D
with a nonzero minimal cover boundary beta_D(S).
```

2. **Boundary-Moment Chain Map.**

```text
B_D(beta_D(S)) equals the HYP-2704 boundary/depth quotient plus a signed
resonant deviation controlled by retained mod-7/affine/CRT labels.
```

3. **Kernel Classification.**

```text
If B_D(beta)=0 and the strict-open witness is empty, then beta is AP/GW-owned,
or it carries a named K33/state-lift packet, or the gauntlet has found a new
kernel class.
```

4. **Far Positive Image.**

```text
For qdiv>14 packets outside the named kernel, B_D(beta) has positive gK8/L_y
slack after the shell gates and exact-period labels are paid.
```

This last lemma is where HYP-2840, THM-534, HYP-2704, and HYP-2666 should
merge.

## Tournament Analysis

The useful tournament vertices are not runners.  I considered:

```text
runners,
unit witnesses U14,
exact-period packets a/D,
packet-cover boundary faces,
missed-depth states 0..6,
C27 owner transfers,
K33/state-lift obligations,
proof lenses.
```

Chosen layer for this pass:

```text
vertices = proof lenses in the boundary-to-moment bridge.
```

Pairwise observable:

```text
which LRC predicate is preserved before scalarization:
  q-witness,
  boundary-only support,
  exact-period packet ownership,
  mod-7/affine/CRT defect,
  missed-depth boundary character,
  gK8/L_y positivity,
  state-lift obstruction.
```

Switch/gauge:

```text
A -> B when A preserves strictly more boundary ownership before the first
scalar comparison, tie-broken by the proof order below.
```

Tie Hamiltonian path:

```text
qdiv gate
> exact-period packet atlas
> boundary-owner/C27-unital packet
> boundary-moment map B_D
> missed-depth death-chain quotient
> gK8/L_y moment dual
> K33/THM-572 state lift
> raw scalar M or residue data
```

Fingerprint:

```text
transitive proof-lens tournament,
singleton SCCs,
unique Hamiltonian path.
```

This quotient preserves the proof predicate "where can a counterexample still
hide?"  It destroys row identity and many denominator labels, so it is only
valid after exact-period packet ownership and C27/unital/K33 labels have been
attached.

## What To Compute Next

Upgrade the HYP-2950 gauntlet with a boundary-to-moment ledger for each named
row:

```text
AP
GW
near/K33 12->36
petal 10->20
petal 13->26
S138 splices
divisor-loaded lcm tails
random qdiv>14 covering repairs
```

For each row, report:

```text
qdiv,
strict-open Haar mass,
boundary-owner skeleton,
exact-period first surviving D and unit packets,
mod-7/chi_7/affine/CRT labels,
missed-depth vector image,
gK8/L_y slack,
kernel label: AP/GW, petal, K33, wide-positive, unknown.
```

The search target is no longer "find a pretty invariant."  It is:

```text
find a qdiv>14 row whose exact-period boundary maps to zero gK8/L_y slack
without being AP/GW-owned or K33-owned.
```

That is a much smaller and more falsifiable missing picture.
