# Observer-Cut Orbit Ledger

HYP-3054 made observer-extension/cut payloads the boundary coordinate of a
quotient's next operation.  This post adds the orbit-level refinement:

```text
C_q(x,o) = orbit_Aut_q(x)(boundary slice, incidence word, extended shadow)
```

The payload is not the new observer alone and not the old object alone.  It is
how the observer cuts the currently visible fiber, modulo the automorphisms
the quotient still remembers.

Applications:

- A000568: the `R(5)=48` to `U(6)=56` gap is incident/cross-sector orientation
  orbit data, not another node-depth layer.
- AP/GW: endpoint completion exposes owner-essential closed `H1`; owner
  deletion kills it.
- q=23: equal exact `M=2/23` forgets the Haar zeta and endpoint-owner strip.
- K33: the payload is a state-lift/cross-handoff orbit, not AP/GW boundary
  homology.
- Automata and analytic clocks: finite words and `mu^2/phi` capacity need
  exact-packet handoff or blindness reports.
- Diagonal layers: raw `k^2+k` line counts are legal only after
  rectangle/hourglass residues vanish or become named sidecars.

Suggested ledger row:

```text
base_quotient
fiber_id
observer_kind
visible_automorphism_group
cut_payload_orbit_id
changed_lrc_predicate
separating_sidecar
discharge_mode
residual_debt_name
```

Tournament Analysis should use payload columns as vertices.  The pairwise
observable is separation of route/status-changing coarse-fiber pairs, with
exactness, dual annihilation, descent, and proof cost as tie-breakers.
