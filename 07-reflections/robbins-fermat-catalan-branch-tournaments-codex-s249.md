# Robbins-Fermat-Catalan Branch Tournaments

The useful synthesis is that Robbins gives the shape of a safe proof network,
while Fermat-Catalan gives a discipline for branch endpoints that try to hide
inside perfect powers.

After the rebased S233-S245 median, center-control, observability, carrier,
energy, route-state closure, arithmetic-sidecar, and median-hull scheduler
work, this S248 lens sits one layer later:
HYP-3067 asks for unique median centers in route triples, HYP-3068 requires
owner/root sidecars for empty centers, HYP-3069 closes named route states under
Boolean median centers, HYP-3070 checks raw-vs-legal route-triple center
control, HYP-3071 names the cycle-class observability/F7 coordinate, HYP-3074
applies legal sidecar closure before medians, HYP-3075 and HYP-3076 keep
arithmetic coincidences as sidecars, HYP-3078 keeps q-Pochhammer/modular-cusp
tails legal only through finite named polar debt, HYP-3079 gives the
Lean-facing finite-tail and padding ledger, HYP-3077 shows legal centers may
still be nonterminal schedulers, and HYP-3081 asks whether the resulting proof
graph has any naked bridge left before contraction.

Robbins' theorem says a connected graph can be strongly oriented exactly when
it has no bridges.  In the LRC proof stack, a bridge is a dependency edge that
the quotient has made load-bearing: delete it and one side no longer knows how
to verify the other.  That is the graph version of controlled forgetting.  A
quotient may forget a coordinate only if the proof graph remains strongly
oriented after reconstruction, exactness, dual annihilation, descent, AP/GW
boundary stop, or named residual debt.

Once oriented, the proof graph wants to decompose into branch corridors between
small tournaments.  The small tournaments are not runner tournaments by
default.  They are local switchboards: AP/GW boundary atoms, K33 state-lift
exits, C27 petal transfers, q=23 Haar squares, A000568 edge-sector decks,
residual capacitor plates, diagonal-layer rectangle/hourglass cells, Fejer
certificate routes, and automaton or power-lift guards.  The branches are
typed corridors between them.

This is the missing shape behind several older notes:

```text
observer-cut orbit ledger = what a branch interface sees
value-origin ledger = what the small kernel count actually means
hyperbolic reciprocal sidecar = when a triple lane has curvature debt
Robbins no-bridge rule = whether that interface is load-bearing
small tournament kernel = finite local exit switchboard
Fermat-Catalan no-lift guard = finite exception rule for power branches
```

The Fermat-Catalan side should stay negative and local.  The unit-excess lane
has stress points like `p=2, q=27=3^3`; those points are not magic.  They are
places where an infinite-looking branch could be a hidden power lift unless
the exponent vector, p-adic valuation, cyclotomic factor, and finite exception
id are retained.  In a branch graph, that means no power branch may be
contracted unless it has a no-lift certificate or lands in named residual
debt.

The proof angle this suggests is not a new scalar invariant.  Build the proof
graph over HYP-2963 packet fibers, attach small tournament kernels at every
local switch, tag their value origin using HYP-3057, and test bridgelessness
after the known sidecars are added.  If all non-AP/GW branches orient into
strict-open or certificate kernels, AP/GW remain the only closed boundary seed
cycle.  If some branch remains naked, it has to become a named residual rather
than a vague obstruction.

Tournament Analysis for this pass uses branch carriers:

```text
observer_cut_payload_orbit
> Robbins_no_bridge_assembly
> small_tournament_kernel
> endpoint_owner_closed_H1
> residual_capacitor_cut
> K33_state_lift_branch
> C27_petal_transfer_branch
> q23_Haar_square_branch
> diagonal_rectangle_hourglass_branch
> Fermat_Catalan_no_lift_guard
> automaton_gap_branch
> Fejer_Ramanujan_certificate_branch
> raw_scalar_shadow
```

The assumption challenged here is the old reflex that tournament vertices are
runners or even arcs.  For this lens the vertices are proof branches, kernels,
observer cuts, residual plates, power-lift guards, and certificate obligations.
The quotient preserves LRC boundary/open status, route handoff, owner current,
observer-cut payload, and named residual exits.  It destroys raw runner names
and scalar magnitudes unless they are carried as sidecars.

Next implementation target: extend the S220 observer-cut ledger with
`branch_id`, endpoint kernel iso classes, `bridge_status`,
`reverse_verification_mode`, and `power_lift_guard`.  A route/status fiber that
still has a naked bridge after that extension is a real theorem obligation.
