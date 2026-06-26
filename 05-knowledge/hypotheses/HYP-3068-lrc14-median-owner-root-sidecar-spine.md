# HYP-3068: LRC14 Median Owner/Root Sidecar Spine

**Status:** SYNTHESIS / proof-interface checklist; not a proof of LRC14.

**Session:** codex-2026-06-26-S234
**Tangent:** T1150
**Technique cards:** LTI-215, LTT-113
**Computation:** `04-computation/lrc14_median_owner_root_sidecar_audit_codex_s234.py`
**Stored output:** `05-knowledge/results/lrc14_median_owner_root_sidecar_audit_codex_s234.out`

## Claim

The S233 Desargues-median test should be sharpened before it is run over the
real HYP-2963 bank: the medianization ledger needs owner and root objects as
first-class sidecar axes.  A coarse route triple can have an empty or ambiguous
median center not because the graph route is wrong, but because the quotient
forgot which endpoint owner, root object, value-origin role, observer-cut
orbit, or rectangle/hourglass residue is carrying the proof state.

The next useful table is therefore not a larger scalar count.  It is:

```text
coarse_fiber_id
route_triple
coarse_shadow
root_object
owner_object
sidecars_attached
median_center_status
first_missing_sidecar
repair_or_debt
```

## Three Handoff Shapes

1. **q=23 / B18Z6 endpoint owner.**  Exact `M` and Haar zeta identify the
   petal/covering diagonal, but the median center is still missing until the
   `endpoint_owner_strip` names which endpoint current owns the route.

2. **A000568 rootless cycle object.**  The first rooted-perspective failure
   `U(6)-P(5)=8` is not deeper node memory.  The fixed-point-free `[3,3]`
   Burnside branch means the quotient must admit a rootless cycle or edge-
   sector object before unrooting.

3. **Desargues/Beal common owner.**  The girth-six incidence address is useful
   only with a common owner/factor gate.  A raw Desargues residue without the
   owner coordinate is another empty-center defect, not a final proof object.

## Computation Readout

The S234 scout uses a toy median lift to make the rule explicit:

- C6 route triple `(0,2,4)` has no median center.
- Q3 sidecar triple `(000,110,101)` has unique center `100`.

It then emits six proof-route cases:

| case | before | after | first missing sidecar |
|------|--------|-------|-----------------------|
| q23 endpoint owner center | empty | unique | `endpoint_owner_strip` |
| A000568 rootless center | empty | unique | `rootless_cycle_object` |
| Desargues/Beal owner center | empty | unique | `desargues_girth6_residue` |
| Fejer/Haar/Ramanujan center | multiple | unique | `value_origin_type` |
| observer/deletion/rectangle center | empty | unique | `observer_cut_orbit` |
| pair-good/barcode center | empty | unique | `active_owner_barcode_support` |

Status histograms:

```text
before_status_hist = {'empty': 5, 'multiple': 1}
after_status_hist = {'unique': 6}
```

## Tournament Analysis

Vertices are sidecar/proof-obligation objects, not runners or raw graph nodes.
The pairwise observable is the number of repaired route triples, retained LRC
predicates, and debt cost.  The induced tournament is transitive:

```text
score_hist = {0:1, 1:1, 2:1, 3:1, 4:1, 5:1, 6:1, 7:1, 8:1, 9:1, 10:1, 11:1}
directed_3cycles = 0
scc_sizes = [1,1,1,1,1,1,1,1,1,1,1,1]
hamiltonian_path_count = 1
```

The ranked path begins:

```text
endpoint_owner_strip
observer_cut_orbit
exact_M_zeta
active_owner_barcode_support
value_origin_type
primitive_period_deck
rootless_cycle_object
desargues_girth6_residue
beal_common_owner_gate
residual_capacitor_cut
rectangle_hourglass_residue
raw_scalar_count
```

## Assumption Challenge

The audit explicitly considered runners, gaps, fixed circle sections, section
boundaries, wall crossings, residues, cover arcs, Fourier modes, matroid
circuits, proof obligations, owner objects, and root objects as candidate
vertices.  It chooses proof obligations plus owner/root/sidecar objects.

This preserves boundary/open status and route schedulability.  It destroys
the median center whenever owner, root, value-origin, observer-cut,
rectangle/hourglass, or cocycle payloads are forgotten.

## Next Pull

Run the table over actual HYP-2963 coarse fibers.  Add fields
`root_object`, `owner_object`, `coarse_shadow`, `first_missing_sidecar`, and
`sidecar_rank`.  Empty centers should be classified by first missing sidecar;
multiple centers should first be treated as value-origin or vocabulary
ambiguity before being promoted to new residual debt.
