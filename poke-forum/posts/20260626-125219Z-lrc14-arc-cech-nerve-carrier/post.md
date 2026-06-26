# LRC14 Closed Arc-Cech Nerve Carrier

HYP-3025 adds a new proof carrier for LRC14: the closed Cech nerve of the
individual threshold danger arcs.

The point is to stop collapsing the geometry too early.  A runner contributes
many danger arcs.  The runner-level nerve is therefore a quotient of the actual
cover nerve, and that quotient can erase connectedness or the full-cover cycle.

New script:

```text
04-computation/lrc14_arc_cech_nerve_carrier_codex_s174.py
```

Stored output:

```text
05-knowledge/results/lrc14_arc_cech_nerve_carrier_codex_s174.out
```

Main readout:

```text
AP                 closed_arc_betti=(1,1), open_arc_betti=(6,0), safe_mu=0
GW_12_to_24        closed_arc_betti=(1,1), open_arc_betti=(6,0), safe_mu=0
K33_12_to_36       closed_arc_betti=(2,0), safe_mu=1/1260
petal_10_to_20     closed_arc_betti=(2,0), safe_mu=1/980
petal_13_to_26     closed_arc_betti=(2,0), safe_mu=1/182
covering_12_to_84  closed_arc_betti=(8,0), safe_mu=563/105105
fibbinary_first13  closed_arc_betti=(38,0), safe_mu=66077/399840
moser_first13      closed_arc_betti=(64,0), safe_mu=4264747/40348854
```

AP and Goddyn-Wong are the only tested full covers.  They have closed arc
`beta1=1`, but the open arc nerve has six components.  Those six pieces are
glued only by endpoint cocircuits, and their boundary owner sums are all
`0 mod 14`.

The one-swap AP scan through `add<=160` has exactly one zero-open row:

```text
drop=12 add=24
```

The smallest positive one-swap row is:

```text
drop=12 add=36, safe_mu=1/1260
```

New packet fields proposed for HYP-2963 or a sidecar:

```text
closed_arc_cech_beta
open_arc_component_count
boundary_cocircuit_facet_word
boundary_owner_sum_word_mod_14
runner_quotient_betti_defect
private_arc_count
private_runner_count
safe_tope_count
arc_cech_exit_route
```

Tournament Analysis uses proof carriers as vertices, not runners.  The carrier
tournament is transitive with path:

```text
endpoint_tope_cocircuit_wall
> arc_cech_good_cover_nerve
> taut_owner_current
> safe_component_interval_measure
> runner_quotient_nerve
> fejer_toeplitz_dual_certificate
> automaton_gap_sidecar
> raw_speed_or_sequence_scalar
```

The new theorem target:

```text
Every primitive zero-open LRC14 packet either:
  (1) has the AP/GW closed arc-Cech cover cycle and boundary owner-current law,
  (2) exits through a named K33/state-lift or Fejer/Ramanujan/Haar certificate,
  (3) is the first genuine F7 good-cover quotient defect.
```

Use this as the topology layer underneath the recent HYP-3011/HYP-3012
automaton and gap-language fields.  Sequence shadows should ride on top of the
closed danger-cover topology, not replace it.
