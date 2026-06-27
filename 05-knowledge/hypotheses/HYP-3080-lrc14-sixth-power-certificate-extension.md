---
id: HYP-3080
status: SYNTHESIS / bounded Diophantine proof-carrier scout; not a proof of LRC14 or of the sixth-power equations
source: codex-2026-06-27-S248
tags: [lrc14, sixth-powers, equal-sums, fermat-catalan, roth-minkowski, residues, route-triples, controlled-forgetting, tournament-analysis]
related:
  - HYP-3076
  - HYP-3075
  - HYP-3074
  - HYP-3073
  - HYP-3072
  - HYP-3071
  - HYP-3070
  - HYP-3069
  - HYP-3066
  - HYP-3063
  - HYP-3062
  - HYP-3060
  - HYP-3058
  - HYP-3054
  - HYP-3048
  - HYP-2997
  - HYP-2963
  - THM-572
  - OPEN-Q-108
---

# HYP-3080: Sixth-power certificate extension for LRC14

This hypothesis extends the S244/HYP-3076 sixth-power collision sidecar with
an explicit LRC proof-interface certificate ledger.  S244 already says how the
two prompt equations enter the relation-lattice stack; S248 says what exact
tuple, residue, rank, and median-closure payload must survive before those
equations are usable as proof data:

```text
a^6 + b^6 + c^6 = d^6 + e^6 + f^6
a^6 + b^6       = d^6 + e^6
```

The inherited S244 split is:

```text
3-vs-3 sixth-power equality -> native support-six relation data
2-vs-2 sixth-power equality -> rank-lowered square-cube shadow, padded if used
```

S248 adds the certificate coordinates needed before that split is used inside a
route-state median, controlled-forgetting quotient, Roth-Minkowski height
fence, or Fermat-Catalan pressure argument.  It also keeps the S243/HYP-3075
Hurwitz-Markov-Pell lesson nearby: scalar coincidences are proof data only
after their arithmetic addresses are retained.

The S248 scout checks positive unordered pairs through `250` and positive
unordered triples through `80`.

```text
pair inputs checked:   31375
pair collisions:       0
triple inputs checked: 88560
triple collision sums: 5
triple certificates:   5
primitive sample certs: 3/5
shared-term certs:     0/5
```

The first primitive triple certificate is:

```text
3^6 + 19^6 + 22^6 = 10^6 + 15^6 + 23^6
```

The bounded absence of a positive two-lane collision is not a proof that none
exist.  It is a proof-carrier warning: if a quotient sees a two-lane
sixth-power equality and has forgotten the unordered pair, primitive gcd, and
residue payload, it is not theorem-safe.  The existence of small three-lane
collisions gives the complementary warning: route triples can balance
sixth-power mass without identifying the lanes, so raw equal-sum mass cannot
be used as a median center.

## Packet fields

Add the following sidecar vocabulary when a sixth-power equality is used in an
LRC packet, route, certificate, or discharge state:

```text
binary_sixth_gaussian_owner_gate
ternary_sixth_diagonal_current
sixth_power_residue_signature
sixth_power_collision_rank        # 2-lane or 3-lane
sixth_power_collision_sum
sixth_power_left_lane_tuple
sixth_power_right_lane_tuple
primitive_all_terms_gcd
left_lane_gcd
right_lane_gcd
shared_term_filter
mod14_sixth_power_word
mod27_sixth_power_word
mod41_sixth_power_word
two_lane_rigidity_gate
three_lane_resonance_graph_id
raw_equal_sum_scalar_shadow
sixth_power_collision_exit
```

The discharge rule is the same controlled-forgetting rule as HYP-2990 and the
S240 route-state interface.  Forgetting from a collision certificate to the
raw scalar sum destroys lane rank, primitive gcd, shared-term cancellation,
and CRT residue words.  That loss is legal only when the missing coordinate is
reconstructed, dual-annihilated, descended to a smaller packet, stopped at
AP/GW boundary equality, or routed to named THM-572/F7 debt.

## Relation to prior carriers

HYP-3076/S244 supplies the upstream support-arity rule.  The ternary equation
is native because it already has the three-positive/three-negative support-six
shape; the binary equation is rank-lowered because `x^6=(x^3)^2` and needs an
explicit canceling pair before it can be lifted into a six-slot ledger.
HYP-3060/S242 remains the older Desargues/Beal finalizer context: binary
sixth-power equality behaves like a Gaussian-owner / rigidity gate, while
ternary equality behaves like a diagonal-current / route-triple carrier.
HYP-3080 is the packet-language extension: it says what certificate data must
survive before either side is medianized, quotient-compressed, or treated as a
route-state discharge.

HYP-3058 imported Fermat-Catalan through the reciprocal-sum / hyperbolic
pressure sidecar.  Sixth powers sit in that hyperbolic zone, but these
equal-sum equations are additive-energy identities, not direct
Fermat-Catalan equations.  The right import is a certificate sidecar:
hyperbolic pressure can warn that the lane has large curvature debt, while the
collision certificate names the actual equality, residues, and lane rank.

HYP-3062 imported Roth and Minkowski through relation-lattice and height
sidecars.  The sixth-power carrier supplies a concrete low-height wall family:
equal sums of sixth powers define short integer relations among monomials
`x^6`.  Before a height or lattice estimate is used, the packet should retain
the relation support, primitive gcd, covolume or height class, CRT residue
word, and exceptional-approximant exit.

HYP-3063, HYP-3069, and HYP-3074 give the medianization interface.  A
three-lane collision should be treated like a route triple: after legal
sidecars are attached, its center is usable only if it remains a legal
packet/route/certificate/sidecar/discharge state.  The rank-2 gate is the
counterpart of the no-hidden-pair check; the rank-3 graph is the counterpart
of a serious route triple whose center may need extra sidecars.

## Tournament Analysis

Vertices are proof obligations and sidecar carriers, not runners and not the
integers `a,b,c,d,e,f`.

```text
labelled_lrc_packet_sheaf
sixth_power_collision_certificate
two_lane_rigidity_gate
three_lane_resonance_graph
route_state_median_sidecar
ramanujan_crt_residue_word
roth_minkowski_height_fence
fermat_catalan_hyperbolic_pressure
partial_cube_simplex_lane
raw_equal_sum_scalar_shadow
```

The S248 pairwise observable is weighted retained payload axes.  Under the
rank-2 rigidity gauge and the rank-3 resonance gauge the tournament is
transitive with one Hamiltonian path, but the two gauges have `5` edge flips.
The important flip is local:

```text
rank-2 gauge: two_lane_rigidity_gate > three_lane_resonance_graph
rank-3 gauge: three_lane_resonance_graph > two_lane_rigidity_gate
```

That is the intended proof use.  The carrier must know whether the active
packet is using a two-lane equality as a rigidity gate or a three-lane equality
as a resonance graph.  The raw equal-sum scalar is last in both gauges.

## Next pull

Run this sidecar over HYP-2963 route triples and any packet family that uses a
power-lift, Fermat-Catalan, Roth-Minkowski, or partial-cube/simplex analogy.
For each use, emit the rank-2 or rank-3 collision certificate and test whether
route-state closure keeps the median center legal.  A failed center should be
classified as missing collision certificate, missing CRT residue word, missing
height/lattice fence, missing gated route sidecar, or explicit THM-572/F7 debt.

Script:

```text
04-computation/lrc14_sixth_power_certificate_extension_s248.py
```

Result:

```text
05-knowledge/results/lrc14_sixth_power_certificate_extension_s248.out
```
