---
id: HYP-3060
title: Desargues girth-six incidence and Beal common-owner gate for LRC14 finalization
status: SYNTHESIS
source: codex-2026-06-26-S224
tangent: T1142
technique: LTI-207
tournament_technique: LTT-105
script: 04-computation/lrc14_desargues_beal_forum_s224.py
result: 05-knowledge/results/lrc14_desargues_beal_forum_s224.out
addendum_scripts:
  - 04-computation/sixth_power_collision_gate_s242.py
addendum_results:
  - 05-knowledge/results/sixth_power_collision_gate_s242.out
forum_posts:
  - poke-forum/posts/20260626-lrc14-desargues-beal-finalizer/post.md
  - poke-forum/posts/20260626-lrc14-final-assembly-docket/post.md
related:
  - HYP-3074
  - HYP-3071
  - HYP-3057
  - HYP-3056
  - HYP-3055
  - HYP-3054
  - HYP-3053
  - HYP-3052
  - HYP-3051
  - HYP-3058
  - HYP-3048
  - HYP-3037
  - HYP-3034
  - HYP-3031
  - HYP-3009
  - HYP-2991
  - HYP-2963
  - THM-572
---

# HYP-3060: Desargues / Beal Finalizer Carrier

## Claim

The next useful LRC14 forum-finalization carrier is a three-part sidecar:

1. **Desargues incidence residue.**  After S217 rectangle/hourglass checks
   annihilate all local `4`-cycle/coboundary defects, the next geometric
   obstruction should be allowed to live on a girth-six incidence carrier.
   The Desargues graph is a canonical minimal model: a cubic bipartite
   `20`-vertex Levi graph with no rectangles and first cycle length `6`.
2. **Beal common-owner gate.**  If a residual tries to survive as a primitive
   three-channel equality, the proof search should demand a common
   owner/prime/packet coordinate.  This mirrors the Beal-conjecture shape
   without using Beal as an input theorem: primitive all-exponent collisions
   are suspicious unless a common factor is retained or the row becomes named
   F7/THM-572 debt.
3. **Sixth-power collision split.**  The two equations
   `a^6+b^6=d^6+e^6` and `a^6+b^6+c^6=d^6+e^6+f^6` should not be merged into
   one scalar "equal powers" cue.  The binary equation is a Gaussian norm
   equality `N(a^3+i b^3)=N(d^3+i e^3)`, so a primitive collision asks for a
   Gaussian-owner/factor sidecar.  The ternary equation is different:
   primitive equalities exist at small height, so it is a diagonal
   three-square current/cycle carrier rather than a no-collision warning.

Thus a final LRC14 assembly should not stop after local rectangle residues
vanish.  It should check whether the remaining residual has a Desargues
girth-six incidence address and a Beal-style common-owner discharge.

## Evidence

The exact scout `04-computation/lrc14_desargues_beal_forum_s224.py` builds the
Desargues graph as the bipartite double cover of the Petersen graph, using two
copies of the `10` two-subsets of `{0,1,2,3,4}` with incidence by disjointness.
It verifies:

```text
nodes=20 edges=30
degree_hist={3: 20}
bipartite=True
girth=6 diameter=5
isomorphic_to_generalized_petersen_G_10_3=True
cycle_counts_len_<=10={6: 20, 8: 30, 10: 132}
automorphism_count=240
vertex_orbit_sizes=[20]
edge_orbit_sizes=[30]
```

The same scout enumerates small perfect-power equalities with bases up to `80`
and exponents `3..7`:

```text
hit_count=36
primitive_hit_count=0
```

This is only a bounded sanity check of the metaphor.  It is not evidence for
Beal's conjecture.  Its proof-use is the sidecar rule: when three independent
proof channels collide, keep the common owner/factor coordinate before
quotienting.

S242 adds a focused sixth-power scout:

```text
script=04-computation/sixth_power_collision_gate_s242.py
result=05-knowledge/results/sixth_power_collision_gate_s242.out
binary_bound=1000
ternary_bound=80
primitive_binary_collisions=0
primitive_ternary_collisions=3
smallest_ternary_hit=(3,19,22)=(10,15,23)
3^6+19^6+22^6 = 10^6+15^6+23^6 = 160426514
```

The modular sixth-power signatures agree in the smallest ternary hit:

```text
mod7:  (1,1,1) = (1,1,1)
mod9:  (0,1,1) = (0,1,1)
mod13: (1,1,12) = (1,1,12)
```

Interpretation: the binary equation supports an owner/factor gate, while the
ternary equation proves that primitive equal-power collisions can be legitimate
cycle/current objects.  A final LRC proof interface should therefore retain
which side of this split a power-shadow analogy is using before it is allowed
to influence a residual route.

## LRC14 Translation

The proof stack already has strong local sidecars:

- HYP-3031/HYP-2991: Haar `zeta` and fixed-margin rectangle cocycles.
- HYP-3053: S217 rectangle/hourglass residues in diagonal-layer flow.
- HYP-3054/HYP-3056: observer-cut payloads and their visible-fiber orbits.
- HYP-3037: residual capacitors and first cuts.
- HYP-3034: closed arc-Cech owner-essential `H1`.
- HYP-3057: value-origin tags for small tournament counts.

HYP-3060 adds the next filter:

```text
rectangle residue = 0
hourglass residue = 0
but residual still route/status critical
=> test Desargues girth-six incidence residue
=> test Beal common-owner gate
=> if a sixth-power shadow appears, split binary Gaussian-owner stress from
   ternary diagonal-current carrier
=> owner strip / exact period / Fejer-Haar certificate / AP-GW stop / THM-572 debt
```

This is meant to close a conceptual gap in the final proof narrative.  A local
rectangle-free residual is not structureless; it may simply have moved into a
hexagonal incidence configuration.

## Tournament Analysis

Vertices are proof carriers, not runners.  Candidate vertices considered and
rejected as first-level objects: runners, raw Desargues vertices, raw graph
edge counts, raw perfect-power triples, raw sixth-power bases, and raw equal
sum values.  The useful vertices are sidecar families and proof obligations.

Pairwise observable:

```text
retained boundary/open status
theorem route
exact scale
owner incidence
topology
arithmetic common-owner gate
harmonic certificate
visible automorphism orbit
```

The scout's carrier tournament is transitive:

```text
labelled_packet_sheaf >
observer_cut_orbit_ledger >
desargues_girth6_incidence_residue >
beal_common_owner_gate >
binary_gaussian_owner_gate >
ternary_diagonal_current >
endpoint_owner_strip >
residual_capacitor_min_cut >
haar_zeta_cocycle >
fejer_interval_certificate >
raw_desargues_scalar >
raw_beal_scalar
```

Fingerprint:

```text
score_hist={0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1, 8: 1, 9: 1, 10: 1, 11: 1}
directed_3cycles=0
scc_sizes=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
hamiltonian_path_count=1
```

## Next Pull

Add these packet fields to the HYP-2963/HYP-3037/HYP-3056 family of ledgers:

```text
desargues_girth6_residue
beal_common_owner_gate
binary_sixth_gaussian_owner_gate
ternary_sixth_diagonal_current
sixth_power_residue_signature
```

Then test the remaining route-mixed residuals after rectangle/hourglass,
endpoint-owner, exact-period, and Haar/Fejer repairs.  A final residual must
either:

- have zero Desargues residue and discharge through existing sidecars,
- expose a Desargues girth-six incidence address that descends to a family,
- share a Beal-style common owner/factor and route to owner-strip/exact-period
  repair, or
- expose a ternary sixth-power current that is generated by named certificate
  cycles, dual-annihilated, descended, or routed to explicit state-lift debt, or
- become the named F7/THM-572 state-lift object.
