---
id: HYP-2947
title: LRC14 measurable rank recombination packet
status: PROOF-INTERFACE / finite low-frontier recombination; not yet a proof
source: codex-2026-06-24-S145
related:
  - HYP-2949
  - HYP-2948
  - HYP-2946
  - HYP-2945
  - HYP-2944
  - HYP-2942
  - HYP-2940
  - HYP-2937
  - HYP-2908
  - HYP-2248
  - THM-572
  - OPEN-Q-108
results:
  - 04-computation/lrc14_measurable_rank_recombination_codex_s145.py
  - 05-knowledge/results/lrc14_measurable_rank_recombination_codex_s145.out
---

# HYP-2947: LRC14 Measurable Rank Recombination Packet

S145 combines the recent POKE carriers into one finite proof interface:

```text
exact M/Farey branch
-> Haar/Borel witness carrier
-> C27 shell transfer
-> q=3 unital local chart
-> affine depth/order signature
-> K33/Kuratowski packet
-> HYP-2248 measurable address tax
-> PH-style bad-child rank
-> HYP-2908/THM-572 state-lift endpoint
```

The target is not a scalar proof.  It is a recombination lemma for the known
low-frontier rows:

```text
Every low-frontier packet either
  (a) is tight AP/GW,
  (b) discharges through a unit C27 petal/two-swap strip,
  (c) or carries the nonunit K33/depth-14 state-lift address.
```

The S145 computation confirms this split for the exact AP two-swap frontier
with added values `<=40`, after the rigorous `q_threshold>=14` filter:

```text
q>=14 exact rows audited = 8674
M<=3/41 rows            = 3
M<=2/27 rows            = 7
unknown low atoms       = 0
```

## Packet Rank

Define the proof-interface rank by retained nonunit obligations:

```text
rank(packet) = uncharged nonunit/K33 obligations after exact M/Farey,
               C27 owner/carry labels, q=3 unital chart choice,
               and product-measurable witness addresses are attached.
```

For the S145 low frontier this gives:

| Row | M | Components | Route | Rank |
|---|---:|---|---|---:|
| AP | `1/14` | none | tight AP floor | 0 |
| GW `12->24` | `1/14` | `GW` | tight GW floor | 0 |
| near-miss `12->36` | `3/41` | `K33` | K33 state-lift obligation | 1 |
| `drop(10,12)->add(20,24)` | `2/27` | `P10,GW` | petal plus tight-GW discharge strip | 0 |
| `drop(10,12)->add(20,36)` | `2/27` | `P10,K33` | petal plus K33 state-lift strip | 1 |
| `10->20` | `2/27` | `P10` | unit-petal discharge | 0 |
| `13->26` | `2/27` | `P13` | unit-petal discharge | 0 |

Thus the exact local rank histogram is:

```text
rank 0: 5 rows
rank 1: 2 rows
```

The only state-lift obligations are:

```text
near-miss 12->36
drop(10,12)->add(20,36)
```

These are precisely the rows that keep a nonunit K33/Kuratowski address.  If
those rows functorially produce the complete tournament conflict packet required
by HYP-2908, then THM-572 closes the endpoint.

## Depth-14 Signature

HYP-2944's affine-depth grammar becomes a useful order-sensitive invariant.
The linked q=3 unital chain

```text
GW -> K33 -> P10
```

has component depths

```text
[3, 4, 1]
```

and suffix-depth profile

```text
[8, 5, 1]
```

with total depth `14`.  Among all permutations of the same component multiset,
this is the unique order with depth sum `14`:

```text
K33,GW,P10  -> 13
GW,K33,P10  -> 14
K33,P10,GW  -> 15
GW,P10,K33  -> 17
P10,K33,GW  -> 18
P10,GW,K33  -> 19
```

Readout: depth `14` is not a scalar coincidence.  It appears only when the
tight AP/GW block is followed by the nonunit K33 block and then by the unit
petal discharge address.

## Measurability Guardrail

The quotient audit is:

```text
global phase             safe         compact group translation on T=R/Z
C27 shell labels         keep label   finite labelled residue packet
q=3 unital chart         keep chart   branch-local pair-completion chart
Kpq/K33 minor flag       keep flag    graph-minor predicate
owner-private deletion   mandatory    HYP-2248 address tax
raw row scalar M only    unsafe       forgets witness address
```

The preserved LRC predicate is:

```text
positive Haar measure of GOOD cap G_P after exact packet labels survive.
```

Only after a measurable invariant action is proved may runner identity or raw
residue labels be discarded.  This is the Borel/Baire/Haar import from the
current POKE thread.

## Rebase Connection: Baire-Haar Boundary Fronts

Concurrent HYP-2948/HYP-2949 adds the regular-open boundary-front companion to
this rank packet.  Exact strict-safety calibration at threshold `1/14` gives:

```text
AP                 safe_mu = 0        endpoint-only residual
GW 12->24          safe_mu = 0        endpoint-only residual
near/K33 12->36    safe_mu = 1/1260   positive open witness
petal 10->20       safe_mu = 1/980    positive open witness
petal 13->26       safe_mu = 1/182    positive open witness
```

This is compatible with S145's rank split: rank-0 tight packets AP/GW are
boundary-only, unit petals are Haar-positive discharge packets, and the K33
near-miss is Haar-positive but still carries a nonunit state-lift address.

Use HYP-2948/HYP-2949's Haar-Baire front only after the S145 packet address
has been retained.  Otherwise the positive open mass has lost the owner/carry
data needed by HYP-2908/THM-572.

## Tournament Analysis

Tournament vertices are proof interfaces, not runners:

```text
exact M/Farey branch
Haar/Borel witness carrier
C27 owner/carry shell code
q=3 unital local chart
affine order/depth code
K33/Kuratowski state-lift flag
PH bad-child rank
raw scalar count
```

Pairwise observable:

```text
branch retention,
address retention,
measurability,
state-lift fit,
extension-rank control,
anti-scalar guard.
```

The tournament is transitive:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
c3=0
singleton SCCs
Hamiltonian paths=1
```

Hamiltonian path:

```text
exact M/Farey branch
> Haar/Borel witness carrier
> C27 owner/carry shell code
> q=3 unital local chart
> affine order/depth code
> K33/Kuratowski state-lift flag
> PH bad-child rank
> raw scalar count
```

## External Analogies

Iravanian's real-factor recombination paper (`arXiv:2410.15880`) is a useful
analogy only at the address level: local factors become useful when the
recombination subset is retained as data.  Euler/Bell's translation of Euler's
divisor-sum recurrence (`arXiv:math/0411587`) gives the companion warning: a
sparse signed support becomes a theorem after a product/log-derivative address,
not after scalar counting.

For LRC14, the analogue is:

```text
local low-frontier atoms
-> retained C27/unital/K33/affine/Haar addresses
-> finite rank recombination
-> state-lift endpoint
```

## Remaining Hard Theorem

S145 proves no global LRC14 theorem by itself.  It proves the exact local packet
classification for the known two-swap ceiling.  The remaining theorem is:

```text
Every primitive LRC14 counterexample, after AP-tail and q-threshold reductions,
enters the measurable packet language above.
```

Equivalently, every surviving bad atom must either:

1. reduce to AP/GW tightness,
2. enter a unit C27 petal/two-swap discharge,
3. or carry the nonunit K33/depth-14 state-lift address.

The honest risk remains that low-frontier labels recur in safely loose rows.
Therefore exact `M`/Farey data and measurable owner/carry addresses must stay
attached throughout the proof.
