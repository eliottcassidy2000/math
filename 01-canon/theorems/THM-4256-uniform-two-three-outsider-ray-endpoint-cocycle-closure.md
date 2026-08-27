---
id: THM-4256
title: "Uniform two-to-three outsider-ray endpoint-cocycle closure"
status: >
  PROVED RELATIVE TO THM-4231/4252/4254 + FINITE-EXACT + INDEPENDENTLY
  AUDITED. For every integer g>=146 and every nine-body in the fixed
  thirty-label pool, adjoining (2g,3g) leaves Haar-safe mass at least 4/63.
  The exact new post-THM-4254 residual contribution is 73 edges; removing
  them leaves 180,991 residual edges. A 4,096-repair bridge proves all 870
  scales through g=1015, and a repair-specific oscillation bound proves every
  g>=1016. The latter threshold is sharp only for that sufficient predicate.
  No other ratio, below-boundary scale, physical entry, or LRC(14) follows.
source: root/cross-frontier-bridge/2026-08-26
depends_on:
  - THM-4231-arbitrary-pair-cofinal-depth-six-haar-repair-and-exact-outsider-lift
  - THM-4252-exact-signed-endpoint-cocycle-residual-edge-closure
  - THM-4254-fixed-ceiling-band-signed-endpoint-cocycle-cascade
related:
  - THM-4233-pair-specific-primitive-observable-oscillation-haar-charts
  - THM-4245-primitive-observable-component-floor-and-cofinal-gate-redundancy
common_header: 04-computation/lrc14_two_three_outsider_ray_thm4256/independent_common.hpp
common_header_sha256: 03420c21edbd95647503396bcb81d4e34645acbc9fb829d364495cf213e6637d
finite_bridge_script: 04-computation/lrc14_two_three_outsider_ray_thm4256/independent_finite_bridge.cpp
finite_bridge_O0_output: 05-knowledge/results/lrc14_two_three_outsider_ray_thm4256/independent_finite_bridge_O0.out
finite_bridge_O3_output: 05-knowledge/results/lrc14_two_three_outsider_ray_thm4256/independent_finite_bridge_O3.out
finite_bridge_script_sha256: 592896b7cb35f8103a3fcbbd99ab0d4ad7f5c4a320617a0141cb07a4512ef7db
finite_bridge_output_sha256: e73647e36093c619d9d4928a90840241be6fc6c009a1898533277f22717e7375
robust_tail_script: 04-computation/lrc14_two_three_outsider_ray_thm4256/independent_robust_tail.cpp
robust_tail_O0_output: 05-knowledge/results/lrc14_two_three_outsider_ray_thm4256/independent_robust_tail_O0.out
robust_tail_O3_output: 05-knowledge/results/lrc14_two_three_outsider_ray_thm4256/independent_robust_tail_O3.out
robust_tail_script_sha256: 38dc7dbf3e91df2ca74c0b49f2a465da175a71a9905746049b996dd9eb4eb7f1
robust_tail_output_sha256: 8a0395107ff76de54af84c6d40cf2493f2f8cd7ecac4ed9ee89aa6784190b9f6
complete_control_script: 04-computation/lrc14_two_three_outsider_ray_thm4256/independent_complete_deck_controls.cpp
complete_control_output: 05-knowledge/results/lrc14_two_three_outsider_ray_thm4256/independent_complete_deck_controls_O3.out
complete_control_script_sha256: 3de387ae18de96278dcf16f1dc667353694d83954d0b739c51d00765e0daf2d4
complete_control_output_sha256: 8d596773301a158c7f976eea78fa63d7b8ec702b4536399d906601761fc844b2
literal_hostile_script: 04-computation/lrc14_two_three_outsider_ray_thm4256/independent_joint_wall_hostile.py
literal_hostile_output: 05-knowledge/results/lrc14_two_three_outsider_ray_thm4256/independent_joint_wall_hostile.out
literal_hostile_script_sha256: d4cecbda6bf6845c5930fe5eb1895c37bd44da54a43591d9debba4e708a170d8
literal_hostile_output_sha256: 603cc470b11452e880f0d06b39964a227ea99cfecb93f1f1d4efe9eb92d9aaea
residual_postprocess_script: 04-computation/lrc14_two_three_outsider_ray_residual_postprocess_thm4256.py
residual_postprocess_output: 05-knowledge/results/lrc14_two_three_outsider_ray_residual_postprocess_thm4256.out
residual_postprocess_script_sha256: 45b58af8c567d68b3cde30e0947d7f5123f2bb3cd8c58fcd98a0f6afc7c7a7a2
residual_postprocess_output_sha256: 6992288d330b7e2d35dbbd4bde9a21f21b7aa3f9ad0e0f48c8d8bb38ecda9d04
packet_manifest_sha256: 200e40d2e66251e45864a39a22b3dee528d61608d42f457eba1973ac104330aa
hash_basis: raw LF bytes
audit: >
  PASS / ACCEPT. Independent O0/O3 bridge and robust-tail transcripts are
  byte-identical. The bridge checks 12,447,220,500 body-scale cases; the tail
  checks all 5,852,925 repairs and 14,307,150 bodies. A separate full-deck
  transform, literal joint-wall calculation, standalone primitive derivation,
  and ASan/UBSan smokes all pass. The semantic residual postprocessor freezes
  the exact 73-edge contribution and matches under normal and optimized Python.
---


# THM-4256 -- uniform two-to-three outsider-ray endpoint-cocycle closure

**PROVED RELATIVE TO THM-4231/4252/4254 + FINITE-EXACT +
INDEPENDENTLY AUDITED; LRC(14) REMAINS OPEN.**

## Statement and inheritance

For a finite positive set `A`, put

```text
G_A={x in R/Z:min_(a in A)||ax||>=1/14},
alpha=4/63,
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,264,286,290}.
```

The theorem is

```text
for every integer g>=146 and every B in binom(P,9),
mu(G_(B union {2g,3g}))>=4/63.                         (1)
```

The lower endpoint is the first integer scale with both newcomers strictly
above the pool maximum: `2g>max(P)=290` iff `g>=146`. This is not a claim of
theorem failure or optimality below 146; no claim is made for a smaller scale.

There is an important proof-graph overlap. In the maintained fixed-pool pair
census, THM-4254 closes the entire residual ceiling band and leaves maximum
unresolved endpoint `754`. Hence its inherited census supplies
`252<=g<=256`, while THM-4231's broad cofinal mechanism supplies every
`g>=257`. The block not covered en bloc by those cutoffs is `146<=g<=251`;
some individual points in that block were already closed.

The semantic post-THM-4254 residual contains exactly 73 ray edges, at scales

```text
146..206, 208..212, 214..217, 220, 221, 230.            (1a)
```

Their edge ledger has FNV `6be23222a3a20764` and SHA-256
`9d43ad3311533711e21c496141b05052c45ed68b8b03bcbce59737df3c7391ea`.
These 73 edges are the exact new proof-graph contribution of THM-4256.
Removing them leaves `180,991` edges, FNV `021bf0ed1581657f`, SHA-256
`9192c5d73aa5f123ddd10f0115dcaf7231fa518980610042e4cd3f8e73afd44f`,
and maximum endpoint `754`.

The proof-graph-minimal route uses the finite bridge on those 73 scales plus
the inherited cutoffs. Independently, the bridge continues through `g=1015`
and the robust tail from `g=1016` is a self-contained reusable mechanism.

## Exact primitive and component identity

For the primitive pair `(u,v)=(2,3)`, the THM-4252 notation gives

```text
N=14uv=84,       S=64,       beta=S/N=16/21,
DeltaC=max C_j-min C_j=536,
omega=DeltaC/N^2=67/882,
D=18,241,159,416,480.                                  (2)
```

For a repair `R in binom(P,8)`, let `U_R=G_(P\R)`. Its positive-length
components have total numerator length `m_R` on denominator `D` and component
count `c_R`. With THM-4252's signed endpoint integer `K_R(g)`, exact mass is

```text
mu(U_R intersect G_(2g) intersect G_(3g))
  =(g N S m_R+K_R(g))/(g D N^2),                       (3)
```

so the exact active predicate is

```text
g A_R+63K_R(g)>=0,
A_R=63 N S m_R-4D N^2.                                (4)
```

The centered primitive has range `DeltaC/N^2`, hence every oriented component
contributes at least `-D DeltaC` to `K_R(g)`. Therefore

```text
K_R(g)>=-c_R D DeltaC.                                 (5)
```

For every repair with `A_R>0`, define its phase-free robust threshold

```text
t_R=ceil(63 c_R D DeltaC/A_R).                         (6)
```

Equations `(4)--(6)` prove that the repair is active for every integer
`g>=t_R`. This is repair-specific and retains `m_R,c_R`; it is stronger than
substituting one worst pool-wide component count into the scalar gate.

## Uniform tail certificate

Exact fixed-pool geometry over all `binom(30,8)=5,852,925` repairs gives

```text
repairs with A_R>0                         5,654,814
repairs with t_R<=1016                       968,835
FNV of the latter deck              d809676561e7da8e.  (7)
```

The deck in `(7)`, ordered by

```text
(SplitMix64(mask xor 0x4245422842334245),mask),        (8)
```

has transversal number greater than nine: all `14,307,150` labelled
nine-bodies are disjoint from some repair in the deck. The scan records

```text
checks                  662,366,204
maximum prefix checks       761,091
prefix FNV              810f1ca450105a52.             (9)
```

Thus `(1)` holds for every `g>=1016`. As a boundary control, the complete set
of repairs satisfying this same phase-free bound at `g=1015` has `965,228`
members and misses exactly one body,

```text
{88,95,170,193,240,252,264,286,290}.                  (10)
```

This makes `1016` sharp only for the particular robust predicate `(6)`; it is
not an actual failure of `(1)` at `g=1015`.

## Finite exact bridge

Take the first `4,096` repairs in the global order `(8)`, before applying any
mass predicate. Their mask ledger is

```text
candidate FNV=225a0c689159c41e.                       (11)
```

For each of the `870` integers `146<=g<=1015`, evaluate `(3)` exactly for
these fixed candidates, retain the active subset, and enumerate every labelled
nine-body. Every scale closes. Only the 73 residual scales in `(1a)` are proof-graph-new; the remaining
bridge scales are an independently verified method control. The
aggregate transcript is

```text
scales                                  870
body-scale cases             12,447,220,500
repair/body checks           356,813,681,317
minimum active candidates              2,763 at g=256
maximum prefix checks                     873 at g=256
worst body at g=256
  {80,88,95,132,145,168,193,252,286}
transcript FNV                 266a67dd758f4e75.       (12)
```

The transcript ledger serializes, for each increasing `g`, the seven unsigned
little-endian words

```text
g, active count, active-mask FNV, failures,
total checks, maximum checks, worst-body mask.          (13)
```

At the hostile minimum-activity scale, the bundled complete-deck control finds
`3,879,889` active repairs and reproduces the same order-prefix length `873`
and the same worst body. The same exact control closes all bodies at
`g=146,1015,1016`; its four-scale transcript FNV is `93df134ded5e3ed4`.
These controls do not replace the all-scale bridge.

The unique body missed by the phase-free bound at `g=1015` is closed by the
finite bridge after 66 active repairs, with

```text
R={10,40,60,80,85,143,145,176},
84 D g (63 mu-4)=301,703,072,020,224,840>0.
```

A bundled full-joint-wall calculation, which does not use the endpoint
primitive, independently gives

```text
mu=607168742429/9120579708240,
mu-4/63=28084316509/9120579708240>0.
```

For the proof-graph-minimal conclusion, the bridge closes the 73 ray edges
listed in `(1a)`; existing canon supplies every other `g>=146`. Independently,
the full bridge proves through `g=1015`, and the robust tail proves every
`g>=1016`.

## Consequence and scope firewall

For each active repair disjoint from `B`,

```text
B subset P\R,
G_((P\R) union {2g,3g}) subset G_(B union {2g,3g}),
```

so the repair mass implies the body mass in the correct direction.

The theorem does not prove a primitive ratio other than `2:3`, a scale below
the stated strict-above-pool boundary, a smaller bridge candidate universe, a
globally minimal repair deck, physical entry, or LRC(14). The `g=1015`
control `(10)` is not a counterexample to the theorem. The independent
implementation reproduced the exact candidate order, all 870 body scans, the
robust thresholds, and the tail transversal audit. **QED.**
