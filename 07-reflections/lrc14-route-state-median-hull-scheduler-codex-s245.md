# LRC14 Route-State Median-Hull Scheduler - S245

The useful angle from this pass is that the final proof interface should not
look like one more ranking of routes.  It should look like a median graph
check on completed proof states.

The completed state has five layers:

```text
packet / route / certificate / sidecar / discharge
```

The sidecar laws in this scout are unary implications.  If a route label is used, exact `M`,
endpoint owner, safe topology, and magnitude cocycle are required.  If an
observer-cut payload is used, value origin, deletion fiber, and cross-sector
orientation are required.  If a partial-cube word is used, bit phase and LRC
packet data are required.  If a Toeplitz or square-peg witness is used, the
positive-scale/noncollapse coordinate is required.  If Desargues/Beal is used,
the girth-six residue and common-owner gate are required.

That unary implication form is the structural point: legal states are closed
under coordinatewise majority.  So a serious triple can have a unique center
without pretending one route has already won.  Any future conjunctive sidecar
guard has to be compiled into its own named coordinate or separately checked;
arbitrary Horn theories do not get this median property for free.

## Computed Readout

The S245 scout defines `41` coordinates and `34` unary Horn sidecar rules.  Starting
from ten route/sidecar states, the median hull has `31` states.  It checks
`29,791` triples and finds `0` raw illegal majorities,
`closure_added_features_hist={0: 29791}`, `0` interval-intersection failures,
and `0` illegal centers after closure.

That is not a proof of LRC14.  It is a consistency check for the final proof
interface: the recent sidecar vocabulary can be organized as a finite
median-closed state space.

## Angle 1: Low-Frontier Scheduler

The triple

```text
AP/GW boundary route
C27 petal route
K33 state-lift route
```

has a legal median center, but it drops the specific AP/GW, C27, and K33
terminal atoms.  The center keeps the common packet/route/certificate/discharge
type.

This is exactly the right failure mode.  It says the low frontier is not a
place to choose by scalar pressure.  It is a scheduler center.  To terminate,
we need a separating sidecar: C27 shell transfer, K33 owner/common-factor
incidence, AP/GW endpoint boundary equality, or named residual debt.

## Angle 2: Noncollapse Sidecar Centers

The triples involving Toeplitz scale, Fejer/Toeplitz certificates,
Desargues/Beal girth-six residue, Moser partial cubes, and hyperbolic
reciprocal pressure also produce scheduler centers.  The common lesson is
that a certificate-like or geometry-like route should not be trusted unless
the median center still contains:

```text
exact_M
endpoint_owner
safe_topology
certificate_payload
sidecar_payload
```

When a specific terminal atom drops out of the median, the proof has not
failed.  It has located the next sidecar that must split the fiber.

## Proof-Interface Proposal

Add these fields to packet ledgers when doing final assembly:

```text
median_state_id
median_hull_id
median_center_kind
median_dropped_atoms
median_required_refinement
specific_discharge_atom
```

The terminal proof rule becomes:

```text
Every serious route triple has a unique legal median center.
If the center is terminal, discharge it.
If the center is a scheduler, add the named separating sidecar.
If no separating sidecar exists, name residual debt.
```

This reframes the remaining work.  We are not looking for another universal
scalar.  We are checking whether the current sidecar sheaf is median-complete
over packet fibers.

## Tournament Analysis

Vertices are proof-interface coordinates, not runners:

```text
labelled_packet_sheaf >
route_state_median_center >
horn_sidecar_closure >
discharge_atom_type >
observer_cut_payload >
partial_cube_cut_payload >
toeplitz_noncollapse_gate >
hyperbolic_triple_pressure >
raw_route_label
```

Pairwise observable: retained LRC predicate, route specificity, certificate
legality, sidecar closure, and discharge refinement.

The script's tournament is transitive, with one Hamiltonian path.  That is not
evidence the problem is easy; it says the proof-interface order is coherent
once the vertices are chosen at the sidecar level.

## Next Move

Run this on actual HYP-2963 packet fibers.  The first targets should be:

- AP/GW, C27 petal, and K33 low-frontier packets.
- q=23 residual capacitor packets.
- Moser/fibbinary automatic fibers after S231 bridge-rank sidecars.
- Fejer/Toeplitz certificate packets against Desargues/Beal finalizers.

For each triple, record the median center, dropped terminal atoms, and the
first sidecar that can split the scheduler center.
