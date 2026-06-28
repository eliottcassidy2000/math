---
id: HYP-3199
title: LRC14 n=4 Einheit/Erdos-870 tournament chart minimality
status: SYNTHESIS / executable n=4 chart scout and proof-packet proposal; not a proof
source: codex-2026-06-27
tangent: T1299
technique: LTI-299
tournament_technique: LTT-199
script: 04-computation/lrc14_n4_einheit_minimality_tournament_scout_codex_20260627.py
result: 05-knowledge/results/lrc14_n4_einheit_minimality_tournament_scout_codex_20260627.out
reflection: 07-reflections/lrc14-n4-einheit-minimality-tournament-models-codex-20260627.md
external_source:
  - https://github.com/davidturturean/erdos-870
  - https://github.com/davidturturean/erdos-870/blob/main/problem_statement.md
  - https://github.com/davidturturean/erdos-870/blob/main/Erdos870/MainTheorem.lean
related:
  - HYP-3160
  - HYP-3161
  - HYP-3153
  - HYP-3152
  - HYP-3151
  - HYP-3150
  - HYP-3149
  - HYP-3148
  - HYP-3147
  - HYP-3146
  - HYP-3145
  - HYP-3144
  - HYP-3143
  - HYP-3142
  - HYP-3141
  - HYP-3133
  - HYP-3124
  - HYP-3053
  - HYP-3049
  - HYP-3054
  - HYP-2129
  - HYP-2078
  - T1299
  - LTI-299
  - LTT-199
  - OPEN-Q-108
---

# HYP-3199: LRC14 n=4 Einheit/Erdos-870 tournament chart minimality

## Claim

The fixed-Hamiltonian-path n=4 tiling model is an abundant cover of score
classes, not a group quotient. The partial-score `(0,1,1,2)` model is an
exact two-coordinate Einheit/Klein-four chart. As a prompt-exact refinement
of HYP-3149's canary/filler quotient interface, HYP-3148's live-core
deletability audit, HYP-3145's filler-core interface, HYP-3146's
shift-package/canary packet, and HYP-3147's n=3 edge-flip kernel sidecar, the
proof-use lesson from Erdos #870 is local but sharp:
representation abundance cannot be used as a substitute for a minimality or
deletion sidecar.

For LRC14 packet work, every fixed-path tiling, score-sequence, or
tournament-class signal should carry these fields before it is used as proof
evidence:

```text
n4_model_scheme_id
n4_fixed_hamiltonian_path_word
n4_partial_score_seed
n4_free_arc_chart
n4_score_class_family
n4_fiber_multiplicity_by_class
n4_quotient_congruence_defect
n4_einheit_torsor_status
n4_deletable_arc_coordinate
n4_minimal_chart_status
n4_S_class_representation_multiplicity
erdos870_representation_abundance_status
erdos870_minimality_sidecar_status
lrc_edge_chart_section_status
n3_edge_flip_kernel
n3_minority_edge_gate
function_quartet_order_status
ordered_function_payload
symmetric_function_shadow
worpitzky_descent_word
compression_circuit_profile
scheme_A_cover_to_section_circuit
transformation_monoid_status
lee_yang_binomial_cap_status
phi4_off_circle_dip_status
large_radius_balance_q0_q6_R6
k8_bimodality_functional
bounded_core_resolvent_degree_ceiling
abel_ruffini_wall_status
ear_omega_sidecar_status
newton_maclaurin_quartic_AP_status
terminal_exit_or_named_debt
```

## Two n=4 Models

Use vertices `0,1,2,3`, fixed path `0->1,1->2,2->3`, and chord flips
`a=0->2`, `b=1->3`, `c=0->3`. Score classes are `T=(0,1,2,3)`,
`S=(1,1,2,2)`, `+=(0,2,2,2)`, and `-=(1,1,1,3)`.

Scheme A fixes the Hamiltonian path and exposes the three chord flips:

```text
*    E    a    b    c
E    T    +    -    S
a    +    T    S    S
b    -    S    T    S
c    S    S    S    T
```

Scheme B fixes four arcs `p01,p12,p23,c`, whose partial score sequence is
`(0,1,1,2)`, and leaves `x=a`, `y=b` free:

```text
*    E    x    y
E    T    +    -
x    +    T    S
y    -    S    T
```

The completed Scheme B table is the Klein-four torsor in score labels:

```text
T = Einheit = E
+ = x
- = y
S = x+y
```

## Computed Signals

The full fixed-path cube has these fibers:

```text
T: size=1 states=['E']
+: size=1 states=['a']
-: size=1 states=['b']
S: size=5 states=['c', 'ab', 'ac', 'bc', 'abc']
```

Thus `S` is the singular high-multiplicity class. The cover forgets enough
information that score-class multiplication is not well-defined:

```text
+*S can land in ('-', 'S')
-*S can land in ('+', 'S')
S*+ can land in ('-', 'S')
S*- can land in ('+', 'S')
S*S can land in ('+', '-', 'S', 'T')
```

The minimality/deletion audit separates cover abundance from essential chart
coordinates:

```text
tiling generators a,b,c reach=('T', '+', '-', 'S')
delete a: reach=('T', '-', 'S') deletable=False
delete b: reach=('T', '+', 'S') deletable=False
delete c: reach=('T', '+', '-', 'S') deletable=True
exact chart x,y reach=('T', '+', '-', 'S')
delete x: reach=('T', '-') deletable=False
delete y: reach=('T', '+') deletable=False
```

`c` behaves like filler in the fixed-path cover. `x,y` are the essential
order-two chart core.

## Function/Resolvent Continuation

The n=3 edge-flip graph is the smallest warning that the edge perspective is
not symmetric after quotienting.  With classes `C3=(1,1,1)` and
`T3=(0,1,2)`, the per-state edge multiplicity matrix is:

```text
C3 -> C3: 0, C3 -> T3: 3
T3 -> C3: 1, T3 -> T3: 2
```

The random-edge transition matrix is therefore
`[[0,1],[1/3,2/3]]`, with stationary distribution `(1/4,3/4)` and
nontrivial eigenvalue `-1/3`.  This is the user's two-node/three-edge graph:
from a cycle every flipped edge escapes to transitive, while from a transitive
tournament exactly one minority edge returns to cyclic and two edges remain
transitive.

Keep the Worpitzky degree-three word as a separate descent sidecar:

```text
x^3 = binom(x,3) + 4*binom(x+1,3) + binom(x+2,3)
Eulerian/descent weights = 1,4,1
```

The sidecar is suggestive because it has two extremes and a bulk, but it is
not the same object as the edge kernel.  The proof packet should record both.

The function quartet now becomes a test for what the quotient destroyed:

```text
symmetric shadows: a+b, a*b
ordered payloads:  a^b, b^a
```

Do not collapse a directed witness to a scalar until the ordered payloads have
been tested.  This is the same discipline as the n=4 chart: `S` can be
represented many ways, but those representatives are not equivalent proof
witnesses.

The exact compression from Scheme A to Scheme B is a tiny monotone circuit:

```text
x = a OR c
y = b OR c
T = not a and not b and not c
+ = a and not b and not c
- = not a and b and not c
S = c OR (a and b)
```

This explains the transformation-monoid behavior.  The projected canonical
`c` action is not a `V4` action on score classes: it sends `T,+,-` to `S` and
sends the canonical `S` representative back to `T`.  On the full cube, XOR is
invertible; after class projection, hidden representatives make the operation
multi-valued.  The compression is useful precisely because it is small, but it
must carry its lost preimage as witness debt.

For the LRC proof attempt, add a Lee-Yang/quartic/resolvent ledger rather than
using any single scalar:

```text
lee_yang_binomial_cap_status:
  circular-zero cap model = binomial/Pascal mass with de Moivre-Laplace bulk
phi4_off_circle_dip_status:
  dip = quartic off-circle correction packet, analogous to exp(-lambda*S^4-b*S^2) dS
large_radius_balance_q0_q6_R6:
  for a degree-6 reciprocal core, test q0=q6*R^6 before trusting a radius signal
k8_bimodality_functional:
  L_y = p0 + p6 + (1/10)*p3; consec is the k=8,9,10 benchmark with L_y(consec)<=cap
bounded_core_resolvent_degree_ceiling:
  expected hard core degree <= 4
abel_ruffini_wall_status:
  stay below the generic quintic wall; degree <= 4 keeps the Galois group in an S4-solvable window
ear_omega_sidecar_status:
  use SC ear / factor-critical odd ear / nested ear decompositions as recursive odd-cycle witnesses
newton_maclaurin_quartic_AP_status:
  test the quartic moment inequality against the arithmetic-progression section
```

This is a hypothesis map, not a theorem.  The recurring structure to measure
is: binomial cap supplies the circular-zero baseline, the `phi4` dip supplies
the off-circle correction, the large-radius balance supplies the reciprocal
root check, and the degree-four bounded core supplies the solvability
discipline.  The ear/Omega sidecar is the recursive tip/tail witness that
prevents a scalar cap from forgetting where the odd cycle or strongly
connected core actually lives.

## Incoming Mainline Integration

After rebasing over HYP-3160, HYP-3150/HYP-3151, HYP-3152, HYP-3153, and the
S72 Lee-Yang circle packet, read HYP-3199 as the prompt-exact n=4
deletion/minimality refinement, not as a competing function-compression or
Lee-Yang packet lane.

HYP-3150/HYP-3151 give the broader rule:

```text
compression q: X -> Y
observable  f: X -> Z
legal iff   f is constant on q-fibers, or the missing data is sidecarred
```

They sharpen this packet by proving that `x=a OR c`, `y=b OR c` has no affine
`GF(2)^3 -> GF(2)^2` substitute for the n=4 class target, and by recording
the k=8 centered quartic `u^4-5u^2+4` as the bounded-core degree-four wall.
Thus HYP-3199 should supply the exact Scheme A/Scheme B tables, `S`-fiber
multiplicity, quotient defect, and deletion audit to the HYP-3151
factor-through ledger.

The S72 Lee-Yang circle packet verifies the same scalar motifs in the current
hard-core language: `q0=q6*R^6` by Vieta, zero radii near a common circle,
`L_y=10*(q0+q6+0.1*q3)` as the k=8 bimodality functional, and the correction
that the flip action is a transformation monoid while solvability comes from
degree `<=4` / Galois group inside the solvable `S4` window.  This upgrades
the HYP-3199 Lee-Yang ledger from a wild signal list into explicit sidecars
that should be attached to HYP-3141/HYP-3142 rows.  HYP-3153 is the reserved
follow-on that should test these sidecars as one Lee-Yang/Worpitzky/quartic
compression packet; HYP-3199 supplies its exact n=4 cover/section and deletion
audit input.  HYP-3160 further sharpens the live k=8 target from a generic
quartic ledger to total covariance / variance extremality, with entropy ruled
out; the HYP-3199 contribution is to make sure any such covariance packet does
not silently quotient away the n=4 canary/fiber witnesses.

## Erdos-870 Transfer

The Erdos-870 repository gives a Lean-formalized negative answer to Erdos
Problem #870: for each `k >= 3` and each `C > 0`, there is an order-`k` basis
with at least `C log n` representations of large `n`, while no order-`k`
subbasis is minimal. The proof architecture uses an order-two source,
deterministic filler for `k >= 4`, and a clustered/parity route for `k=3`.

This does not transfer as a theorem to LRC14. The useful transfer is a proof
discipline:

```text
many representations / many tilings / many score-fiber states
  !=
minimal proof basis / essential chart coordinate / legal quotient
```

In the n=4 chart, the analogy is exact enough to be operational:

- order-two core -> exact `x,y` Einheit chart;
- filler -> deletable `c`;
- representation abundance -> the five-state `S` fiber;
- minimality sidecar -> deletion audit plus quotient-congruence defect.

## LRC14 Use

Before adding a fixed-path tiling or score-sequence signal to HYP-3141
edge-witness rows, require the packet to declare whether it is a cover or a
section. If it is only a cover, carry `n4_fiber_multiplicity_by_class` and
`n4_quotient_congruence_defect`. If it claims to be a chart, require
`n4_einheit_torsor_status=True` and a deletion audit showing which coordinates
are essential.

This is a small n=4 model, but it sharpens a recurring LRC14 guardrail:
abundance may be the easiest signal to measure and the wrong signal to trust.
The proof packet should prefer exact sections with deletion/minimality data
over large quotient fibers.

## Tournament Analysis

Vertices are model/proof carriers, not runners, raw arcs, or score labels.
The pairwise observable is retained exact quotient status plus
minimality/deletion information. The gauge is weighted retained proof payload;
ties use a lexical stable Hamiltonian path.

```text
score_hist={1: 1, 2: 1, 7: 1, 9: 1, 10: 1, 17: 2, 22: 2, 23: 2, 24: 1}
directed_3cycles=0
hamiltonian_path_count=1
priority_path=
  minimality_deletion_sidecar
  -> compression_circuit_canary_map
  -> lee_yang_quartic_resolvent_packet
  -> einheit_xy_exact_chart
  -> n3_edge_flip_worpitzky_kernel
  -> ear_omega_recursive_witness
  -> erdos870_abundance_nonminimality_gate
  -> S_class_multiplicity_alarm
  -> tiling_path_cube_cover
  -> partial_score_0112_seed
  -> score_sequence_family_shadow
  -> raw_flip_name_table
```

## Assumption Challenge

Candidate tournament vertex sets considered: raw arcs, runners, score classes,
flip coordinates, cube states, quotient fibers, deletion tests, exact charts,
fixed circle sections, section boundaries, wall-crossing events, residues,
cover arcs, Fourier modes, matroid circuits, and proof obligations.

Chosen vertex set: model/proof carriers. This preserves the LRC predicate only
when exact quotient status and minimality/deletion data are retained. It
destroys hidden `S`-fiber multiplicity, the fact that `c` is filler, and the
ambiguity of the fixed-path class product unless those fields are explicitly
attached.

Challenged assumption: the visible score-class table is not the algebra. The
algebra appears only after choosing a section; the fixed-path tiling cover is
not itself a quotient group.

## Next Work

1. Complete the Scheme B multiplication table with all four completed labels
   and test other partial-score seeds at n=4.
2. Join this with HYP-3053 fixed-path half-tiling fibers and HYP-3133/HYP-3134
   A000568 quotient discipline.
3. Attach `n4_quotient_congruence_defect`, `n4_deletable_arc_coordinate`, and
   `n4_minimal_chart_status` to HYP-3141 edge-witness rows.
4. Use the Erdos-870-style minimality audit as a sidecar for HYP-3142-style
   abundance and moment-majorization packets.
5. Add the K3 edge kernel, Worpitzky descent word, compression-circuit
   profile, Lee-Yang binomial cap, `phi4` dip, `q0=q6*R^6` balance,
   degree-four resolvent ceiling, and ear/Omega recursive witness to the next
   LRC14 hard-core runner.
