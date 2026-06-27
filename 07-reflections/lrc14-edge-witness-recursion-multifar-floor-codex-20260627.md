# LRC14 Edge Witness Class-Deck Stress and the Multi-Far Handoff

*codex-2026-06-27.  S271 class-deck supplement to the
HYP-3124/T1198/LTI-259/LTT-157 edge-witness lane, now feeding the
HYP-3125/HYP-3127 multi-far floor program.*

## What changed

The edge-as-witness idea now has a finite scout:

```text
04-computation/lrc14_tournament_edge_witness_recursion_codex_20260627.py
05-knowledge/results/lrc14_tournament_edge_witness_recursion_codex_20260627.out
```

The computation treats a directed tournament edge `tail -> tip` as a
two-ended observer.  S268 first established the labelled `n<=5` packet under
HYP-3124; this S271 addendum tests the same packet against unlabelled class
decks through `n=6` and the multi-far floor.  The outside vertices are split
into four sectors, then the edge keeps both recursive deletion children:

```text
edge_witness =
  sector_counts
  + sector_internal_deck
  + tail_deletion_child
  + tip_deletion_child
  + cross_sector_orientation_word
```

The class-level result is tight and useful:

```text
n=3,4,5: sector counts/internal decks already separate all unrooted classes.
n=6: sector counts/internal decks separate 55/56 classes.
     The lone collision is the known converse pair 344/345.
     Adding tail/tip deletion children separates 56/56.
```

At the edge-instance level, the payload story is stricter.  For `n=6`,
all `43` nontrivial sector-internal edge fibers are split by tail/tip deletion
children, and `16` recursive fibers are split further by cross-sector
orientation.  So the cross-sector word is not redundant; it is just not needed
for the particular `n=6` class-deck collision once the two deletion children
are retained.

## HYP-3125 handoff

HYP-3121 reduces the open covering branch to the few-apex/multi-far floor:

```text
Rprime = mu(R-safe and Q-lonely)/(mu(R-safe) mu(Q-lonely)) >= c
```

The edge-witness reading is:

```text
tail = R-safe marginal
tip  = Q-lonely marginal
cross-sector word = covariance payload
tail/tip children = recursive marginal tests
```

The remote HYP-3125 scout is now the canonical measured `Rprime` packet.  This
class-deck supplement contributes a narrower point: a quotient may forget the
exact witness time only after it keeps enough two-ended payload, especially
tail/tip deletion children and cross-sector orientation, or it must name the
lost coordinate before HYP-3125 compresses it into an `edge_floor_packet`.

## Four proof routes worth testing

1. **Edge-recursion packet.**  Add `edge_witness_id`, `tail_child`,
   `tip_child`, and `cross_sector_orientation_word` to HYP-2963/HYP-3098
   few-apex packets.  This is the finite combinatorial carrier.

2. **L2 floor certificate.**  Promote
   `lrc14_spectrum_L2tail_synthesis_kpswf12.py` from a spectral scout into a
   row certificate for HYP-3121: exact low modes plus Parseval tail must output
   a floor and the destroyed-coordinate ledger.

3. **Gaussian smoothing.**  The incoming S254 scale-separation work makes this
   more than a metaphor: prove a heat-smoothed correlation floor, then
   desmooth with finite-ruler margins and an explicit channel-discrepancy
   budget.  This makes the Fourier tail less brittle, but it must retain
   boundary atoms and denominator thresholds.

4. **Asano endpoint contraction.**  The incoming S68 work supplies the right
   shape: encode tail/tip interactions in a multiaffine polynomial and
   contract endpoint variables only if the zero-free statement still reads
   back to endpoint `Phi/P` activation.  This is a Lee-Yang-style route, not a
   substitute for metric arc data.

The Elliott-Halberstam cue fits as an averaged-distribution route: define
finite channel errors over `(q, residue, edge_witness_id, Fourier shell)` and
prove an average bound.  Incoming work now makes the guardrail sharper: EH is
not load-bearing for the proof unless that finite inequality is actually
stated and checked.

## Working conjecture

The uniform multi-far floor should be attacked as a coupled finite/analytic
packet:

```text
edge_witness_packet
+ L2/Gaussian Fourier floor
+ endpoint Phi/P activation
+ finite-ruler desmoothing
+ named residual exit
```

The scout's carrier tournament is transitive, with this priority path:

```text
edge_tail_tip_recursive_witness
> l2_cauchy_schwarz_tail_floor
> cross_sector_orientation_word
> gaussian_heat_kernel_smoothing
> endpoint_phi_activation_gate
> normal_fan_cech_finite_ruler
> asano_endpoint_contraction
> elliott_halberstam_level_packet
> raw_bonferroni_scalar
> raw_edge_count_shadow
```

The negative controls matter: Bonferroni already fails on the few-apex examples
from HYP-3121, and raw edge counts do not preserve the LRC predicate.  Any
route that does not keep a two-ended payload, an analytic floor, or an endpoint
activation gate is only telemetry.

## Assumption challenge

The vertices in this session are not runners.  The scout explicitly considered
runners, gaps, fixed circle sections, section boundaries, wall crossings,
residues, cover arcs, Fourier modes, matroid circuits/topes, directed edges,
and proof obligations.  For the multi-far floor, the useful vertices are
edge witnesses, residue/Fourier channels, and proof-carrier packets because
they preserve the predicate `R-safe intersect Q-lonely nonempty`.

The destroyed information is exact time, endpoint owner, Fourier phase, and
sometimes tail/tip role.  Those are precisely the coordinates the next packet
schema has to retain, average with a finite error budget, resurrect by
Gaussian/Asano/endpoint sidecars, or route to HYP-2968/HYP-3121 debt.
