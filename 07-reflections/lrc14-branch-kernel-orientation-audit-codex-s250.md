# LRC14 Branch-Kernel Orientation Audit

S250 turns the S249/HYP-3081 Robbins-Fermat-Catalan branch-orientation idea,
downstream of the S247/HYP-3079 Lean q-cusp finite-tail ledger, S246/HYP-3078
q-cusp principal-part gate, S245/HYP-3077 median-hull scheduling, and
S240/HYP-3074 route-state closure, into an executable
packet-bank audit.  The test is deliberately bounded: it runs on the existing
HYP-2963/S168 bank and asks whether known sidecars protect every bridge in the
branch graph.

The useful result is that the raw scalar view is exactly as bad as expected:
it forms a star with `5` bridges and all `5` are naked.  Once the section map,
Haar grid exits, endpoint-owner strips, C27 petal discharge, K33 state-lift
debt, power-lift guards, and Desargues/Beal finalizer gate are attached, the
protected graph has `80` nodes, `83` edges, `69` bridges, and `0` naked
bridges.  Contracting the protected bridge edges leaves a one-node core, so
Robbins' strong-orientation obstruction vanishes inside this bounded ledger.

This should change the shape of the LRC14 proof target.  The next missing
lemma is not a new scalar invariant on packets.  It is a global-reduction and
debt-discharge statement:

```text
primitive residual
  -> S250 branch-kernel ledger
  -> section/grid exit or named state-lift debt
  -> no new naked bridge after contraction
```

The audit also gives a useful negative lesson.  Coarse quotients still mix
routes: automatic words, branch labels, q-factorizations, and power-lift guards
each create mixed fibers.  The branch graph becomes safe only after the exits
are named and protected.  That is the controlled-forgetting rule in a more
graph-theoretic form.

Tournament Analysis used branch carriers and proof debts as vertices, not
runners.  The retained-payload tournament is transitive with one Hamiltonian
path:

```text
observer_cut_payload_orbit
> Robbins_no_bridge_assembly
> labelled_packet_bank
> residual_section_exit
> Haar_grid_exit
> endpoint_owner_strip
> C27_petal_branch
> K33_state_lift_branch
> covering_moment_branch
> Fermat_Catalan_no_lift_guard
> Desargues_Beal_finalizer_gate
> named_residual_debt
> raw_scalar_shadow
```

The challenged assumption is still load-bearing: the vertices are not runners,
arcs, or scalar speeds.  They are branch kernels, observer-cut payloads,
section exits, grid exits, no-lift guards, finalizer gates, and residual debts.
This quotient preserves boundary/open status, exact scale, q-threshold,
endpoint owner, route handoff, section exit, and residual naming.  It destroys
raw runner identity and raw automatic/power numerology unless those coordinates
are retained as sidecars.

The honest caveat remains sharp.  HYP-3082 is not LRC14.  It does not prove the
global reduction to the packet bank, THM-572, K33 discharge, or the covering
family theorem.  It makes those the next clean proof obligations.
