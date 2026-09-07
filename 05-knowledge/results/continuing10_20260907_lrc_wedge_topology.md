# Actual graph topology deletes eleven of the fifteen projected 7200 words

**Status: PROVED relative to the two native-wedge theorems; FINITE-EXACT
consumer, independently audited.** This is a separate checkpoint from
`continuing10_20260907_lrc_composite_wedges.md`; its frozen producer is not
changed. The connected-complement hypothesis below is essential.

Let a primitive row of thirteen distinct positive integer speeds have a
selected six-label subset A with gcd(A)=7200 and a **connected strict-atlas
graph on its seven-label complement B**. Under hypothetical weak LRC(14)
failure, the inherited complete profile and minimum-tree consumer leaves the
fifteen words in the continuing8 certificate. Write d_i=gcd(7200,u_i).

The two new native-wedge theorems forbid, in the **actual** complement graph,
both margin paths `(4,24,18)` and `(16,24,18)`. At each positional vertex of
margin 24, either all edges to margin 18 are absent, or all edges to margins
4 and 16 are absent. This is the exact disjunction expressing that no two
such incident edges coexist. If both classes of edges are absent, either
branch may be used. Distinct equal-margin labels are retained as distinct
vertices.

Enumerate these binary choices for every 24 vertex. For a choice sigma let
H_sigma contain every locally compatible atlas edge not forbidden by that
choice, weighted by the inherited **minimum** uniform pair credit over all
ratios compatible with its two margins. Every actual graph without a forbidden
wedge is a subgraph of at least one H_sigma. Its connectivity forces that
H_sigma to be connected. Every actual spanning tree then has actual credit
at least the minimum weighted spanning tree of H_sigma. The inherited forest
inequality gives a safe lift if this minimum exceeds E_7200(d).

Consequently a word is impossible for a connected unsafe complement when
every branch is disconnected or has minimum-tree cost greater than E. This
uses MINIMUM trees on generous possible-edge graphs, not maximum trees of
independently invented edges. Branching preserves the needed quantifier over
all actual graph topologies; it does not assume that a previously selected
Kruskal tree is the actual graph.

There are exactly 29 branches on the fifteen inherited words. Their ordered
cost pairs below mean respectively “remove 24--18” and “remove 24--(4,16)”.
The sole word without 24 has one branch.

| Word | E | Branch minimum costs | Deleted |
|---|---:|---|---|
| 1,9,16,18,24,32,60 | 116 | 216,114 | no |
| 4,4,4,9,9,18,24 | 42 | 114,266 | yes |
| 4,4,8,9,9,18,24 | 42 | 114,168 | yes |
| 4,4,9,9,16,18,24 | 50 | 114,302 | yes |
| 4,4,9,9,18,24,32 | 66 | 114,190 | yes |
| 4,5,9,18,24,30,32 | 79 | 162,180 | yes |
| 4,8,8,9,9,18,24 | 42 | 114,92 | yes |
| 4,8,9,9,16,18,24 | 50 | 114,148 | yes |
| 4,8,9,9,18,24,32 | 66 | 114,92 | yes |
| 4,9,9,16,16,18,24 | 58 | 114,226 | yes |
| 4,9,9,16,18,24,32 | 74 | 114,112 | yes |
| 5,8,9,18,24,30,32 | 79 | 162,66 | no |
| 5,8,9,30,32,36,48 | 103 | 102 | no |
| 5,9,16,16,18,24,30 | 71 | 162,208 | yes |
| 5,9,16,18,24,30,32 | 87 | 162,66 | no |

Thus eleven whole projected words, rather than eleven chosen ratio assignments,
are deleted in this actual connected-complement scope. The four displayed
survivors are only survivors of this necessary-condition consumer. No actual
unsafe realization, joint low-credit ratio tree, or simultaneous minimizing
phase is asserted for them. No clock is deleted; 7200 remains among the
necessary clocks. The topology argument does not apply to disconnected B,
although the native-wedge safety theorems themselves do.

The standalone verifier pins the complete inherited minimum-tree certificate
and the two complete native-wedge certificates. It retains every positional
edge and every branch, reconstructs the unrestricted prior minima, and checks
each chosen edge by the independent minimum-cut property. Its JSON contains
all allowed weighted graphs and chosen trees, not just the surviving word
list. The source imports no mathematical producer.

```
python 04-computation/continuing10_20260907_lrc_wedge_topology.py
python -O 04-computation/continuing10_20260907_lrc_wedge_topology.py
```

Both modes pass 549 active exact gates, with identical actual LF output and
certificate bytes. Frozen pins:

```
source 6fc998e7596bfa2c4182bfe41423f3d4cbdd9d9a4f7ce77724afebdfdb653e6b
output 8ca50996eac493eaa9c9c98e543f73fff29758be008bafe7e6847b0af140082c
certificate 4008a6ae441672637becabb0f0ba17f8e0fe9e5975657fe9dddefe2bb5766e79
```

The [independent audit](continuing10_20260907_lrc_wedge_topology_audit.md) accepts the scoped proof.
Repository filing changes only routing and status prose; frozen source,
output and certificate bytes match the independently reviewed originals.
