---
id: HYP-3124
title: LRC14 tournament edge witness recursion
status: EVIDENCE / S268 labelled scout, recursive unlabeled packet scout, and S271 class-deck stress supplement; not a proof
source: codex-2026-06-27-S268; continued by codex-2026-06-27 recursive extension; extended by codex-2026-06-27-S271
tangent: T1198
technique: LTI-259
tournament_technique: LTT-157
script: 04-computation/lrc14_tournament_edge_witness_recursion_codex_20260627.py
result: 05-knowledge/results/lrc14_tournament_edge_witness_recursion_codex_20260627.out
related:
  - HYP-3131
  - HYP-3130
  - HYP-3129
  - HYP-3128
  - HYP-3127
  - HYP-3126
  - HYP-3125
  - HYP-3123
  - HYP-3122
  - HYP-3121
  - HYP-3120
  - HYP-3119
  - HYP-3118
  - HYP-3117
  - HYP-3116
  - HYP-3115
  - HYP-3112
  - HYP-3108
  - HYP-3106
  - HYP-3105
  - HYP-3054
  - HYP-3050
  - HYP-3049
  - HYP-2963
  - OPEN-Q-108
---

# HYP-3124: LRC14 Tournament Edge Witness Recursion

## Claim and Renumbering

This lane continues the S268 response to the edge-as-witness prompt.  The
provisional S268 namespace collided with live mainline work: S269 owns
HYP-3119/T1194/LTI-255/LTT-153 for the niche archive bridge, HYP-3121 is the
lift-and-decorrelate three-engine synthesis, S67 owns HYP-3122 for the
cap/phi4 quartic stabilizer, and S270 owns HYP-3123 for the chiral-stalk /
Cech finite-ruler proof angles.  The clean edge-witness namespace is therefore
HYP-3124/T1198/LTI-259/LTT-157.

The working object is not a tournament vertex and not a scalar edge count.  It
is an oriented edge `tail -> tip` whose proof value is the pair of recursive
witness obligations seen from both ends:

```text
edge_witness(tail -> tip) =
  (endpoint_role_word,
   outside_four_sector_deck,
   tail_deletion_child_packet,
   tip_deletion_child_packet,
   recursive_tail_child_edge_deck,
   recursive_tip_child_edge_deck,
   repair_sidecar_or_named_debt,
   terminal_exit)
```

In LRC14 language, a local edge quotient is legal only if the packet route is
constant on the quotient fiber, or the lost tail/tip coordinate is retained,
reconstructed, dual-annihilated, recursed into a smaller child, or named as
residual debt.

## Source Threads

The existing sources already point at this shape:

- HYP-3050 says an edge is naturally dualistic: it has a tail and a tip, and
  outside vertices split into four sector words around that directed edge.
- HYP-3054 names `edge_tail_tip_sector_word` as an observer-extension cut
  coordinate.
- HYP-3106 turns directed-edge sectors into a controlled-forgetting functor.
- HYP-3112 and HYP-3115 expose proof-relevant edge boundaries through
  ear-payload edges and Lee-Yang/Ising domain-wall edges.
- HYP-3116/HYP-3118 say any such edge witness is legal only if the destroyed
  coordinate is either retained, resurrected by a sidecar, or routed to named
  residual debt.
- HYP-3122 says the phi4 signal is a quartic wall stress, useful only after it
  can be localized to a retained packet or named cumulant debt.
- HYP-3123 says chiral/converse/rootless quotients need orientation sidecars;
  tail/tip recursion is the edge-local test for that guard.
- HYP-3128 says the Asano/Lee-Yang picture is a diagnostic, not a standalone
  joint-floor proof: the apex/tip block is zero-free, but the overcrowded
  R/tail block has interior Lee-Yang zeros, so any edge contraction that
  forgets the tail packet is illegal.
- HYP-3129 supplies the compatible repair: the multi-far `Rprime` floor is an
  elementary fixed-lattice SPEC bound on the retained `R-safe -> Q-lonely`
  edge, with a certified floor `Rprime >= 0.64178` on the tested family.
- HYP-3130 and HYP-3131 (the S69 Lee-Yang polydisk scout) say the far/apex tip side is
  not the binding obstruction: Gaussian/Beurling-Selberg minorants close the
  apex block with cap-matching constants, single/multi-far additions to a good
  base push PGF zeros outward (`rho >= 1.5589` in the HYP-3131 probes), and the
  remaining work is the bounded-core/tail coupling plus signed SPEC.

Challenged assumption: an edge is not just a relation between two vertices.
For LRC-style proof search it is a bidirectional proof obligation: the tail
asks what survives after pushing forward, while the tip asks what survives
after pulling back.  A legal witness must make both recursions compatible.

## S268 Labelled Scout

The incoming S268 scout
`04-computation/lrc14_tournament_edge_witness_recursion_codex_s268.py`
enumerates labelled tournaments through `n=5` and records each directed edge's
four-sector word, one-sided endpoint-deletion children, paired tail/tip child
signature, and one recursive layer inside both children.  Stored output:
`05-knowledge/results/lrc14_tournament_edge_witness_recursion_codex_s268.out`.

Exact census:

| n | labelled tournaments | directed-edge instances | sector words | sector + tail child | sector + tip child | sector + child pair | depth-1 signatures | sector groups split by child pair |
|---|----------------------|-------------------------|--------------|---------------------|--------------------|---------------------|--------------------|------------------------------------|
| 2 | 2 | 2 | 1 | 1 | 1 | 1 | 1 | 0 |
| 3 | 8 | 24 | 4 | 4 | 4 | 4 | 4 | 0 |
| 4 | 64 | 384 | 10 | 14 | 14 | 16 | 16 | 6 |
| 5 | 1024 | 10240 | 20 | 52 | 52 | 80 | 80 | 20 |

For `n=5`, every sector group is split by the paired endpoint-deletion child
object.  The four-sector word is therefore a strong local observable but not a
witness by itself; the first natural witness carrier is the sector word plus
both endpoint-deletion children.

## Extended Unlabeled Scout

The continuation scout
`04-computation/lrc14_tournament_edge_witness_recursion_codex_20260627.py`
implements the recursive edge object exactly through all unlabeled tournaments
on `n <= 6` vertices.  It compares score sequence, depth-0 edge-sector deck,
and depth-1/depth-2 recursive edge-witness deck:

```text
n  classes  score  edge_d0  edge_d1  edge_d2  d2_collisions
3        2      2        2        2        2              0
4        4      4        4        4        4              0
5       12      9       12       12       12              0
6       56     22       55       56       56              0
```

## S271 Class-Deck Stress Supplement

The S271 follow-up script
`04-computation/lrc14_tournament_edge_witness_recursion_codex_20260627.py`
keeps the same directed-edge carrier but shifts the finite test from labelled
edge instances through `n=5` to unlabelled class decks through `n=6`.  Stored output:
`05-knowledge/results/lrc14_tournament_edge_witness_recursion_codex_20260627.out`.

Class-deck readout:

| n | classes | sector counts | sector internal | roleless children | recursive children | full edge witness |
|---|---------|---------------|-----------------|-------------------|--------------------|-------------------|
| 3 | 2 | 2/2 | 2/2 | 2/2 | 2/2 | 2/2 |
| 4 | 4 | 4/4 | 4/4 | 4/4 | 4/4 | 4/4 |
| 5 | 12 | 12/12 | 12/12 | 12/12 | 12/12 | 12/12 |
| 6 | 56 | 55/56 | 55/56 | 56/56 | 56/56 | 56/56 |

The only `n=6` sector-count/internal collision is the converse pair
`344/345` with score sequence `(2,2,2,3,3,3)`, `c3=8`, and `H=43`.  Retaining
the paired endpoint children repairs that collision.  At the directed-edge
instance level, all `43` nontrivial `n=6` sector-internal fibers split by
tail/tip deletion children, and `16` recursive fibers split further by the
cross-sector orientation word.  This says the cross-sector word is not
cosmetic; it is the remaining covariance/phase payload inside recursive edge
fibers.

The covering-floor translation is now carried by HYP-3125/HYP-3127 rather than
by this HYP-3124 witness card.  Read this supplement as evidence that the
tail/tip deletion children and cross-sector word are genuine information
coordinates before an `R-safe packet -> Q-safe packet` is compressed into the
multi-far `Rprime` floor.

## Edge Information Rule After HYP-3131/HYP-3132

The useful abstraction is now an information split, not a pictorial edge
analogy.  In the lifted covering packet

```text
R-safe packet -> Q-lonely packet
```

the tip recursion carries stabilizing information: Q/apex loneliness,
Gaussian/Beurling-Selberg minorants, and HYP-3131's far-pushes-out Lee-Yang
radius signal.  The tail recursion carries binding information: the bounded
core, the retained SPEC lattice, exact low-frequency mass, and Parseval-tail
debt.  The cross-sector orientation word is the signed coupling between those
two child algebras; it is the phase information that is lost by a collapsed
Asano contraction.

This gives a sharper recursive proof obligation.  A legal edge quotient must
show that every tip child either has a zero-free/minorant/far-push certificate
or descends to a smaller Q child, and every tail child either has a
bounded-core SPEC certificate or descends to a smaller bounded-core packet.
If neither child closes, the first lost coordinate must be emitted as
observer-gluing, coordinate-resurrection, finite-ruler, phi4 cumulant, or
H7/F7 state-lift debt.  A future proof should treat this as an
edge-witness descent lemma rather than as another scalar invariant.

## Integration With Existing Threads

The first useful signal is at `n=6`: raw score data sees only `22` classes and
the one-step sector deck sees `55/56`, but retaining the recursive tail/tip
deletion children repairs the last collision and separates all `56` classes.
The same scout records high endpoint asymmetry at `n=6`: `648`
child-code-asymmetric directed edges, `744` root-perspective-asymmetric edges,
and `346` two-sided middle-sector edges across the unlabeled classes.  These
are exactly the places where an LRC quotient would be illegal unless it keeps
the tail child, tip child, observer-sector deck, or a named repair debt.

## Proof-Lens Tournament

The S268 meta-tournament over edge-witness reframes has one nontrivial SCC and
keeps `cross_sector_orientation_word` and `paired_tail_tip_deletion_recursion`
near the top of the Hamiltonian paths.  The continuation proof-lens tournament
uses proof lenses rather than runners or raw arcs and is transitive:

```text
recursive_edge_witness_packet
-> edge_coordinate_resurrection_guard
-> tail_tip_child_pair
-> decorrelation_edge_floor
-> h7_state_lift_edge_boundary
-> four_sector_observer_deck
-> phi4_ising_edge_wall
-> raw_H_value_shadow
-> raw_edge_count_shadow
```

Next executable test: attach this packet to HYP-3115 one-swap/domain-wall
edges and ask which walls become legal observer-gluing discharges, which
recurse to smaller tail/tip children, and which remain named HYP-2963/HYP-3098
debt.  When the same witness is used for covering-floor rows, use HYP-3125's
`edge_floor_packet` fields.

Fingerprint: `score_hist={0:1,...,8:1}`, `directed_3cycles=0`, singleton
SCCs, `hamiltonian_path_count=1`, and `4` edge flips against the fixed tie
path.  The combined reading is that sector orientation and child recursion are
not competing analogies: orientation is the local guard, and recursion is the
legal sidecar that decides whether the guard preserves proof information.

## LRC14 Use

The edge-witness packet should be added to HYP-2963/HYP-3098/HYP-3107 rows
with the fields:

```text
edge_witness_recursion_id
tail_child_packet
tip_child_packet
four_sector_observer_deck
child_deck_asymmetry
coordinate_resurrection_status
decorrelation_floor_status
asano_obstruction_status
spec_resonance_floor_status
minorant_apex_floor_status
bounded_core_binding_status
state_lift_boundary_status
phi4_edge_wall_status
terminal_exit_or_named_debt
```

The connection to HYP-3120 is that edge witnesses become one of the missing
packet-closure carriers: they can feed observer gluing, finite-address
normalization, coordinate resurrection, proof-circuit missing-input ledgers,
and Lee-Yang ear payloads without collapsing to a scalar.  The connection to
HYP-3121 is sharper: an edge cut can be asked to preserve both child event
algebras for the lift-and-decorrelate floor, while zero-mass edge cuts must
name the H=7/K33/F7 state-lift boundary.  HYP-3125/HYP-3126 sharpen the
positive-mass side of that field: `decorrelation_floor_status` should record
whether the edge cut is in the wide-V elementary decoupling regime, the
bounded-core `3/pi^2` floor regime, the finite `w0` check, or the remaining
minorant/constant-chase debt.  HYP-3127 adds the Asano contraction reading:
tail/tip recursion is a candidate contraction order for multi-far Lee-Yang
factors, but HYP-3128 narrows that reading: Asano cleanly packages the
apex/tip block and exposes the overcrowded tail's Lee-Yang zero obstruction;
it cannot certify the joint floor after the tail packet is collapsed to a
single zero-free factor.  HYP-3129 then supplies the proof-facing floor field:
`spec_resonance_floor_status` should record the resonance lattice, exact low
SPEC sum, Parseval-tail certificate, and constant-chase debt for the retained
two-ended edge packet.  HYP-3130 and HYP-3131 add the complementary orientation:
the tip/far side carries `minorant_apex_floor_status` and Lee-Yang outward-root
monotonicity, while `bounded_core_binding_status` names the tail-side core
that still needs the signed coupling proof.  HYP-3122/S67's `phi4` quartic
stabilizer is therefore only an edge-wall stress signal unless it points back
to the recursive tail/tip packets or to named quartic-cumulant debt.  HYP-3123
/S270 supplies the companion orientation test: a chiral, converse, or Cech
quotient using an edge boundary must either keep the cross-sector orientation
word and both endpoint children, or emit the first destroyed coordinate as
observer-gluing, coordinate-resurrection, or state-lift debt.
