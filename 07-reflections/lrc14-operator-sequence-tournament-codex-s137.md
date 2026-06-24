# LRC14 Operator Sequence Tournament (S137)

This session pulled in the new S131-S136 work and treated it as signal rather
than restarting the prompt.  The useful object is no longer "which of the four
Farey mutations is best?"  S130 already answered the first-order version:
`q` is the binding denominator, `p+q` is the least destructive additive ledger,
`p*q` is a product/coimage side channel, and powers are magnitude stress tests.

S137 adds a sequence and tournament warning around that answer.

## Seven Topics Kept Together

The computation explicitly keeps seven carriers in one frame:

```text
1. Farey unit-excess sequence p/(14p-1)
2. additive n+2 recursion versus product n*2/coimage growth
3. C=27 shell transfer
4. K_{p,q}/K33 incidence
5. octahedral L(K4) support-six current
6. Clebsch/halved-cube folded covariance
7. Paley-Zygmund second-moment existence
```

The point is not breadth for its own sake.  These are the objects that current
OPEN-Q-108 notes keep naming as possible carriers for the remaining
`|M14|<=6` finite-core atom.

## Sequence Readout

On the unit-excess chain:

```text
q        = 14p - 1       linear, Delta=14
p+q      = 15p - 1       linear, Delta=15
p*q      = 14p^2 - p     quadratic, Delta2=28
```

This is the cleanest explanation of why the additive and product analogies both
felt right but should not be collapsed.  The additive lane is theorem-facing:
it follows the Stern-Brocot/Farey child chain linearly.  The product lane is
area-facing: it records the complete-bipartite incidence blow-up `K_{p,q}`.
The powers are too distorted for local proof order, but useful for exposing
fixed-rule tournament quotients that forgot magnitude.

## Frontier Confirmation

The S130 row bank through replacement `70` reproduces the same low frontier
seen in the larger S136 audit:

```text
M <= 3/41: AP, GW, 12->36
M <= 2/27: AP, GW, 12->36, 10->20, 13->26
```

So the exact router is still:

```text
2/27 branch: unit-visible C27 holes, petal/two-block rigidity.
3/41 branch: gcd3 -> gcd9 shell transfer, first K33 child.
```

The transfer quotient is necessary but not complete.  Its labels recur in loose
rows, so exact `M` and Farey branch must remain attached.

## Tournament Analysis

The conservative carrier order is a transitive tournament:

```text
q > C27 > p+q > K33 > p*q > octahedron > Clebsch/half-cube > PZ > powers.
```

But the majority-of-criteria gauge is more honest about the conflict between
sequence simplicity, graph-packet fit, state-lift fit, and magnitude stress.  It
has two directed 3-cycles and a nontrivial SCC:

```text
{p+q, p*q, octahedron, Clebsch/half-cube}.
```

That SCC is the new signal.  It says the middle packet layer is not linearly
ordered by nature.  Use a conservative order for proof routing, but do not
scalarize inside the SCC.  Keep the side-channel labels until either the C27
petal branch discharges, or the K33/octahedral/Clebsch packet is strong enough
to feed HYP-2908.

## Current Proof Shape

The sharpened lemma target is:

```text
After standard LRC14 reductions, every low-gap non-AP/GW atom either has a
unit-visible C27 hole, or has a sign-visible K33/octahedral/Clebsch packet whose
connected OCF state-lift would evaluate to 7.
```

This remains a proof-interface statement, not a proof.  The contribution is the
warning that the additive and packet carriers form a real tournament SCC under
reasonable gauges.  A finished proof should not pick one scalar from that SCC;
it should construct the labelled packet and only then pass to the forbidden-H
endpoint.
