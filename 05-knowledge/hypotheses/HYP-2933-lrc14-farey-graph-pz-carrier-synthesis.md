---
id: HYP-2933
title: LRC14 Farey graph/PZ carrier synthesis
status: PROOF-INTERFACE / carrier hierarchy; not a proof of LRC14
source: codex-2026-06-23-S132
related:
  - HYP-2932
  - HYP-2931
  - HYP-2930
  - HYP-2928
  - HYP-2908
  - HYP-2891
  - HYP-2887
  - HYP-2823
  - HYP-2183
  - HYP-2899
  - HYP-2900
  - THM-572
---

# HYP-2933: LRC14 Farey Graph/PZ Carrier Synthesis

HYP-2931's denominator verdict and HYP-2932's `K_{p,q}` product split survive
the wider synthesis: for
`M(S)=p/q`,

```text
M(S)-1/14 = (14p-q)/(14q).
```

Thus `q` is still the theorem coordinate.  The prompt's four mutated Farey
payloads are best treated as labelled side channels:

- `p+q` is the additive Stern-Brocot / `n+2` ledger.
- `p*q` is the HYP-2932 complete-bipartite area/coimage ledger
  `|E(K_{p,q})|`.
- `q^p` and `p^q` are magnitude-leak stress tests.

The unit-excess child chain `p/(14p-1)` makes the role split exact.  Along this
chain, `q` advances by `+14`, `p+q` advances by `+15`, and `p*q` is quadratic
with constant second difference `28`.  This explains how the earlier `n+2`
and `n*2` recursion pictures can coexist without replacing the binding scale:
the additive lane tracks the local Farey theorem coordinate, while the product
lane tracks two-dimensional packet/coimage growth.

## Graph Carriers

The octahedron, Clebsch graph, and halved cube requested in the prompt have
different proof roles and should not be collapsed to one runner tournament.

Exact fingerprints from `lrc14_farey_graph_pz_carriers_codex_s132.py`:

- Octahedron `L(K4)`: `6` vertices, `12` edges, degree `4`, `8` triangles,
  cycle rank `7`.  Read as the support-six repeated-packet current/curl
  carrier; its cycle rank is the apex-7 module.
- Clebsch folded 5-cube: `16` vertices, `40` edges, degree `5`, triangle-free,
  cycle rank `25`.  Adjacent pairs have `0` common neighbors and nonadjacent
  pairs have `2`; closed neighborhoods form a `2-(16,6,2)` design.  Read as
  the folded residual-mask covariance carrier.
- Halved 5-cube: `16` vertices, `80` edges, degree `10`, `160` triangles,
  cycle rank `65`; under even representatives it is exactly
  `complement(Clebsch)`.  Read as the cut/covariance complement of the
  Clebsch carrier.

This matches HYP-2887 and HYP-2891: octahedral structure belongs to the
support-six current/face-curl side, while Clebsch and half-cube structure
belongs to signed mask covariance and folded cut data.

## Paley-Zygmund Role

Paley-Zygmund is useful as an existence gateway, but it is too lossy for the
tight cap comparison by itself.  In the toy independent six-sector empty model
with `a=(6/7)^k`, the second-moment lower bound loses about `0.16` at `k=8`
and about `0.11` at `k=12` against the exact union probability.

This reinforces HYP-2823: the live cap route needs the labelled degree-4
factorial moment / Delsarte region, not a raw second-moment scalar.

## Tournament Analysis

Tournament vertices are proof carriers, not runners.  The pairwise observable
is

```text
(theorem scale, label retention, LRC relevance, scalarization risk).
```

The switch is lexicographic role score, with lower scalarization risk better.
The tie Hamiltonian path is the listed carrier order.  The resulting role
tournament is transitive (`c3=0`) with Hamiltonian role order

```text
q_binding_scale
> p_plus_q_additive
> octahedron_LK4
> p_times_q_product
> Clebsch_halfcube
> PZ_second_moment
> power_payloads.
```

The small inversion relative to the script's printed carrier list is
intentional: octahedral support-six current data is more LRC-relevant than the
raw product payload, even though both preserve labelled side-channel data.

## Proof Target

The reduced `|M14|<=6` atom should carry a labelled address

```text
(e, q, p+q, p*q), where e=14p-q,
```

plus octahedral/Clebsch packet data.  The next theorem should force such an
atom into either:

1. HYP-2930's non-tight Farey class, proving `M(S)>1/14`; or
2. HYP-2908 / THM-572's forbidden `H=7` tournament-state lift.

This is not an LRC14 proof.  It is a sharper interface: keep the Farey fraction
address, keep packet graph labels, and postpone scalarization until the exact
non-tightness or forbidden-state endpoint is available.

## Assumption Challenge

Vertices considered: runners, Farey fractions, mutated fraction payloads,
residual masks, `K4` edges, half-cube cuts, sector-empty events, and proof
obligations.  The selected quotient preserves which carrier retains the data
needed for the LRC14 binding-scale / state-lift target.  It destroys raw
row-specific witnesses, which must remain in HYP-2930/HYP-2931 and the
exact-period packet ledgers.
