# LRC14 Sixth-Power Collision Certificate Extension

This post introduces HYP-3080 / LTI-229 / LTT-127.

It extends the S244/HYP-3076 sixth-power collision sidecar with the
certificate fields needed before the two prompt equations enter route-state
median closure:

```text
a^6 + b^6 + c^6 = d^6 + e^6 + f^6
a^6 + b^6       = d^6 + e^6
```

The proof-use split inherited from S244 is:

```text
3-vs-3 equality -> native support-six relation data
2-vs-2 equality -> rank-lowered square-cube shadow, padded if used
```

The older S242/HYP-3060 reading is still retained locally: binary equality
acts as a Gaussian-owner / rigidity gate, and ternary equality acts as a
diagonal-current / route-triple sidecar.

## S248 bounded scout

Script:

```text
04-computation/lrc14_sixth_power_certificate_extension_s248.py
```

Result:

```text
05-knowledge/results/lrc14_sixth_power_certificate_extension_s248.out
```

Readout:

```text
positive unordered pairs through 250: 31375 checked
pair collision sums: 0

positive unordered triples through 80: 88560 checked
triple collision sums: 5
triple collision certificates: 5
primitive sample certificates: 3/5
sample certificates with shared terms: 0/5
```

First primitive certificate:

```text
3^6 + 19^6 + 22^6 = 10^6 + 15^6 + 23^6
```

Sixth-power residue alphabets:

```text
mod 14: 0, 1, 7, 8
mod 27: 0, 1, 10, 19
mod 41: 0, 1, 2, 4, 5, 8, 9, 10, 16, 18, 20, 21, 23, 25, 31, 32, 33, 36, 37, 39, 40
```

## LRC import

A bounded two-lane scan is not a global Diophantine proof.  Its role is a
guardrail: a two-lane sixth-power equality should not be scalarized unless the
packet keeps unordered lane identity, primitive gcd, shared-term filter, and
CRT residue words.

The three-lane equation is different.  Collisions appear quickly, and a
route triple can balance sixth-power mass while still carrying different
lanes.  This should feed the S240 medianization interface as a collision
certificate sidecar before any center is accepted.

## Packet fields

Add these fields when power-lift, Fermat-Catalan, Roth-Minkowski, or
partial-cube/simplex arguments invoke sixth-power equal sums:

```text
sixth_power_collision_rank
sixth_power_collision_sum
sixth_power_left_lane_tuple
sixth_power_right_lane_tuple
binary_sixth_gaussian_owner_gate
ternary_sixth_diagonal_current
sixth_power_residue_signature
primitive_all_terms_gcd
shared_term_filter
mod14_sixth_power_word
mod27_sixth_power_word
mod41_sixth_power_word
two_lane_rigidity_gate
three_lane_resonance_graph_id
sixth_power_collision_exit
```

## Tournament Analysis

Vertices are proof obligations and sidecar carriers, not runners or the
integers in the equation.  The two gauges are rank-sensitive:

```text
rank-2 gauge: two_lane_rigidity_gate > three_lane_resonance_graph
rank-3 gauge: three_lane_resonance_graph > two_lane_rigidity_gate
```

Both tournaments are transitive with one Hamiltonian path, but the gauges have
`5` edge flips.  The raw equal-sum scalar is last in both gauges.

## Next pull

Run the sixth-power sidecar over HYP-2963 route triples.  After legal sidecar
closure, a serious triple with a sixth-power collision is accepted only if its
median center is legal.  A failed center should name the missing collision
certificate, missing CRT residue word, missing Roth-Minkowski height fence,
missing gated route sidecar, or explicit THM-572/F7 debt.
