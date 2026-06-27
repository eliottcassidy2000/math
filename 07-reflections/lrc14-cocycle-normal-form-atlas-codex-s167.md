# LRC14 Cocycle Normal-Form Atlas

The useful reframe is that the repo has been doing cohomology even when it did
not call it that.

The early scalar failures now look predictable.  Residues mod `14`, raw
tournament classes, component counts, product values, and divisor summaries all
forget some cochain.  Once the forgotten cochain is not a coboundary, the
quotient starts mixing AP/GW boundary packets with loose rows.  The familiar
warning "residues lie" becomes a cohomological warning: a scalar shadow is safe
only after every vertical cocycle on its fiber has been killed, reconstructed,
or named.

HYP-2990 made this a zipper rule.  HYP-2991 gave the first explicit local
coordinate, the mixed Haar fixed-margin curl

```text
zeta(T)=T00-T01-T10+T11.
```

HYP-2997 says to stop treating that as a special case.  Endpoint owner
currents, C27 carries, Farey excess, exact-period Ramanujan projectors,
Toeplitz PSD duals, boundary-moment chart transitions, tournament `H^1`, and
the state-lift residual are all cocycle channels in one packet complex.

That makes AP/GW cleaner.  They are not just two tight examples with shared
residue DNA.  They are the zero-open boundary packets where the current
channels close:

```text
Farey excess zero
endpoint debt as boundary cocircuit
Haar zeta stopped at boundary
C27/Jacobsthal branch closed at AP or unique D3 acceleration
no Fejer/Toeplitz PSD failure because no open safe interval exists
no F7/THM-572 residual invoked
```

The counterexample search also becomes more precise.  A primitive bad packet
cannot be "mysterious"; it has to exhibit a first nonzero class: Haar curl,
wall current, Farey excess, C27 carry, PSD/dual obstruction, exact-period
survivor, tournament pressure, chart transition, curried section derivative,
or the named state lift.  If it exhibits none of those, the packet has already
collapsed to AP/GW boundary accounting or a certificate.

The Tournament Analysis vertex choice matters.  For this pass, runners and
arcs are the wrong vertices.  The vertices are obstruction classes.  Under the
conservative retention gauge, the carrier tournament is transitive:

```text
labelled_packet_total_cocycle
> haar_zipper_2cocycle
> farey_excess_mediant_1cocycle
> endpoint_owner_boundary_cocycle
> tope_cocircuit_wall_cocycle
> c27_carry_lift_1cocycle
> state_lift_obstruction_class
> fejer_toeplitz_dual_coboundary
> ramanujan_exact_period_character_cocycle
> boundary_moment_multichart_cocycle
> curried_section_derivative_cocycle
> tournament_path_h1_cocycle
> raw_scalar_shadow
```

That order should not be overread as mathematical importance.  It is a
proof-engineering priority order: start with the full packet cochain, retain
the local Haar curl and exact scale, then let endpoint/topology/C27/state-lift
and harmonic exits do their work.  Raw scalar shadows are allowed only as
diagnostics after the cocycle ledger is closed.

The next computation should not add another global ranker.  It should build a
packet-level cocycle ledger for the HYP-2963 bank and tag every row by its
first nonzero class.  The hoped-for result is that AP/GW are exactly the
all-closed zero-open packets, loose rows expose an early certificate or open
witness, and any real survivor is forced into the named `F_lift` sector where
THM-572 can be applied.
