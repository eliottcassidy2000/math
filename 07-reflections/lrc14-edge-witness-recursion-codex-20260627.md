# LRC14 Edge-Witness Recursion Reflection

**Session:** codex-2026-06-27
**Hypothesis:** HYP-3124
**Artifacts:** `04-computation/lrc14_tournament_edge_witness_recursion_codex_20260627.py`, `05-knowledge/results/lrc14_tournament_edge_witness_recursion_codex_20260627.out`

## Core Readout

The useful abstraction is not "an edge is a relation."  It is:

```text
edge_witness(tail -> tip) =
  (four_sector_observer_deck,
   tail_child_after_tip_deletion,
   tip_child_after_tail_deletion,
   recursive_tail_child_edge_deck,
   recursive_tip_child_edge_deck)
```

That object is information-rich in a way raw edge counts and raw H shadows are
not.  Exact enumeration through all unlabeled tournaments on `n <= 6` gives the
clean finite signal:

```text
n=6: score sequence = 22/56
n=6: depth-0 edge sector deck = 55/56
n=6: depth-1 recursive edge witness deck = 56/56
n=6: depth-2 recursive edge witness deck = 56/56
```

The first edge-sector collision is repaired by keeping both endpoint deletion
children.  That is exactly the kind of controlled-forgetting repair the LRC14
proof frontier keeps rediscovering: a quotient can be useful, but only after it
names the coordinate it would otherwise destroy.

## Connection To The Current Proof Frontier

HYP-3120 said the next useful work is a packet-closure router, not a scalar
invariant.  HYP-3124 supplies one concrete packet column family for that router:

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

HYP-3121 collapses the covering branch into three engines: q-witness/census,
lift-and-decorrelate, and H=7 state-lift.  The edge-witness language gives a
local way to ask which side of that split an edge cut belongs to.  Positive
mass cuts should retain both event algebras and feed `decorrelation_floor_status`.
Zero-mass cuts should name the edge boundary needed for K33/F7/H7 state-lift
construction.  Mixed cuts should go to coordinate resurrection or observer
gluing, not to raw H.

Incoming HYP-3125/HYP-3126 sharpen the positive-mass branch: a cut feeding
`decorrelation_floor_status` should say whether it is wide-V elementary
decoupling, bounded-core `3/pi^2` floor, finite `w0` check, or
minorant/constant-chase debt.

Incoming HYP-3127 adds the Asano contraction interpretation: the tail/tip
recursion can be treated as a contraction order for multi-far Lee-Yang
factors only when both endpoint children and the zero-free polydisk sidecar
survive the edge quotient.

Incoming HYP-3128/HYP-3129 revise that interpretation into a sharper
two-part rule.  The apex/tip block is Lee-Yang safe, but the overcrowded
R/tail block has interior Lee-Yang zeros, so a naive Asano contraction that
collapses the tail is blocked for the right reason.  The surviving proof route
is to keep the two-ended packet and attach the elementary SPEC
resonance-lattice certificate: exact low frequencies plus a Parseval-tail
bound, with the current certified floor `Rprime >= 0.64178` on the tested
multi-far family.  In packet terms, `asano_obstruction_status` is the alarm and
`spec_resonance_floor_status` is the positive repair field.

Incoming HYP-3130 and the S69 Lee-Yang polydisk scout add the orientation rule
for that repair.  The far/apex tips are stabilizing: the minorant closes the
apex block with cap-matching constants, and S69 finds multi-far nearest-zero
radius still bounded away from the unit disk (`rho >= 1.5589` in the probes).
The binding side is the bounded core/tail.  Edge recursion should therefore
ask whether the tail child preserves the bounded-core SPEC/minorant coupling,
not whether the far tips themselves create the obstruction.

HYP-3122/S67's `phi4` quartic-stabilizer signal fits as an edge-wall stress
test.  It is not yet a proof carrier by itself.  It becomes a proof-facing
sidecar only when the one-swap/Lee-Yang/Ising wall still points back to the
recursive tail/tip child packets or to named quartic-cumulant debt.  HYP-3123
/S270 adds the adjacent orientation guard: chiral, converse, rootless, and
Cech component quotients may use an edge boundary only after the
cross-sector orientation word and both endpoint recursions survive, or after
the first lost coordinate is named as repair debt.

## Challenged Assumption

The challenged assumption is that edge perspectives are already the right
unit.  The scout suggests a sharper rule: an edge perspective is almost enough,
but the legal unit is the edge plus its two deletion recursions.  The missing
`55 -> 56` bit at `n=6` is the warning label.  In LRC terms, an edge quotient
that forgets either child can silently identify different proof states.

## Next Move

Instantiate these fields on a small set of HYP-3098 observer-gluing rows and
HYP-2963/HYP-3107 residual packet representatives.  The decision test should
be finite and brutal:

1. Does the edge cut preserve the LRC predicate and endpoint ownership?
2. Do the tail and tip child packets land in the same terminal route, or does
   their asymmetry name a repair?
3. Is the cut positive-mass decorrelation, zero-mass state-lift, observer
   gluing, coordinate resurrection, or named debt?
4. If an Asano/Lee-Yang wall is used, does it keep the apex/tip zero-free
   side, expose the R/tail obstruction, and route the floor through the SPEC
   resonance-lattice certificate rather than through a collapsed tail factor?
5. If a minorant/Gaussian wall is used, does it certify the apex-tip floor or
   only an absolute-envelope failure, and does the signed coupling remain in
   the bounded-core tail packet?
6. If a phi4/Ising/Lee-Yang wall is used, does it point back to the child
   packets?
7. If a chiral/Cech quotient is used, does the edge boundary keep the
   cross-sector orientation word and both endpoint children?

Only after those answers are present should an edge, wall, or H value be
allowed to enter the proof route.
