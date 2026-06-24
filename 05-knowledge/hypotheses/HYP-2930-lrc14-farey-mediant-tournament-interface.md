---
id: HYP-2930
title: LRC14 Farey-mediant tournament proof interface
status: PROOF-INTERFACE / local theorem target; not a proof of LRC14
source: codex-2026-06-23-S128
related:
  - HYP-2932
  - HYP-2925
  - HYP-2926
  - HYP-2928
  - HYP-2929
  - HYP-2924
  - HYP-2923
  - HYP-2920
  - HYP-2921
  - HYP-2917
  - HYP-2621
  - THM-523
  - THM-568
---

# HYP-2930: Farey-mediant tournament interface

This hypothesis makes the prompt "the mediant of a Farey sequence; a
tournament LRC14 proof" exact.

For any exact value

```text
M(S) = p/q,    threshold = 1/14,
```

define the **Farey excess**

```text
e(S) = 14p - q.
```

Then `M(S) > 1/14` iff `e(S) > 0`.  The special case `e(S)=1` is exactly the
Farey-neighbor condition against `1/14`:

```text
det [[1,p],[14,q]] = q - 14p = -1.
```

Equivalently,

```text
p/q = p/(14p-1).
```

For `p>=2`, this is the iterated Stern-Brocot mediant child

```text
p/(14p-1) = mediant(1/14, (p-1)/(14(p-1)-1)).
```

So the LRC14 near-miss value is not numerology:

```text
2/27 = mediant(1/14, 1/13)
3/41 = mediant(1/14, 2/27)
3/41 - 1/14 = 1/(14*41).
```

## Computation

The script
`04-computation/lrc14_farey_mediant_tournament_s128.py` stores its output at
`05-knowledge/results/lrc14_farey_mediant_tournament_s128.out`.

It uses three explicit Tournament Analysis carriers.

1. **Apex winding tournament.**
   Vertices are the thirteen speeds.  Pairwise observable is
   `frac((s_i-s_j)/14)`.  The switch is the half-circle cutoff
   `(0,1/2)`.  Ties at `0` or `1/2` are oriented by the increasing-speed
   Hamiltonian path.  This preserves the denominator-14 winding class and
   destroys magnitude.
2. **Optimum Farey winding tournament.**
   Vertices and switch are the same, but the scale is the exact optimizer
   `t=M(S)`.  This preserves the first exact escape class above the floor and
   sees the Farey-mediant child.
3. **Row-level proof-priority tournament.**
   Vertices are candidate rows.  Pairwise observable is exact escape height
   `M(S)-1/14`; smaller escape beats larger escape, with the listed row order
   as the tie Hamiltonian path.  This preserves proof priority and destroys
   runner geometry.

The row-priority tournament is transitive in the audited set:

```text
score_hist = ((0,1),(1,1),...,(14,1))
c3 = 0
scc = fifteen singletons
hp = 1
```

This is useful as a proof ledger, not as a geometric invariant.

## Exact local readout

The exact optimum/Farey classes in the audited local set are:

```text
F0: AP
F1: GW 12->24
F2: residue-liar 12->26
F3: near-miss 12->36, 12m family m=4,5,6,7,8
F4: petal 8->16
F5: petal 10->20
F6: petal 7->14
F7: petal 9->18
F8: petal 11->22, petal 13->26
```

The tight classes are exactly `F0` and `F1` in this local audit.  The
unit-excess Farey children are already separated from them:

```text
near-miss 12->36    M=3/41, e=1, F3 = med(1/14,2/27)
petal 10->20        M=2/27, e=1, F5 = med(1/14,1/13)
petal 13->26        M=2/27, e=1, F8 = med(1/14,1/13)
12m family m=5      M=1/13, e=1, F3 = right Farey parent
```

The apex classes are still magnitude-blind.  For example AP, the residue-liar
`12->26`, and `12m family m=8` share the same apex class `A0`, even though
their exact `M` values differ.  Thus a raw apex tournament cannot prove the
census.

## Proof target

The plausible tournament proof is multi-scale:

```text
candidate row
  -> apex winding class at 1/14
  -> first Farey-mediant / optimum winding class
  -> q-threshold and off-apex data
  -> forbidden tight class unless AP or GW.
```

The local theorem target is:

```text
Every non-AP/GW q=14 survivor either has nonunit excess e>1,
or enters a unit-excess Farey child class outside the tight optimum classes.
```

Either alternative certifies `M(S)>1/14`.  In that form, the Farey mediant is a
proof certificate: `e=1` says the row escapes at the first possible
Stern-Brocot child above the floor, while the optimum winding tournament says
which exact isomorphism class realizes that escape.

S131/HYP-2932 refines the unit-excess alternative by the complete-bipartite
product ledger:

```text
p=1: K_{1,13}   coarse q-threshold parent,
p=2: K_{2,27}   planar two-block strip,
p>=3:           K_{3,3} product-minor wall.
```

Thus the unit-excess proof branch should split into a two-block/petal rigidity
lemma for `p=2` and a finite three-owner obstruction packet for `p>=3`.

## Guardrail

This does not prove LRC14.  It proves the right local dictionary for the
mediant route.  It is compatible with the incoming tournament-spectrum verdict
HYP-2928: the useful object is not one fixed tournament class, but the
Farey-indexed spectrum together with its binding scale.  It is also compatible
with HYP-2925/HYP-2926: fixed non-periodic tournament discriminators can be
fooled in unbounded families, so the mediant datum must be retained as a scale
label rather than collapsed to one class.

The missing theorem is a realizability/rigidity statement: after the global
reductions, every remaining non-AP/GW candidate must be forced into one of the
non-tight Farey classes recorded here, or into an even larger nonunit-excess
class.
