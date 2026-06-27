# Unit-Distance n=21 Endpoint Universality (Codex S129)

The endpoint-ear route has now done something useful, but not what I first
hoped.  It does not close `u(22)=60`.  It tells us why the proof cannot be
only a Hamiltonian-spine proof.

The five exact `n=21`, `57`-edge graph6 cores are endpoint-universal:

```text
every core vertex can be a Hamiltonian-path endpoint
every vertex deletion remains traceable
every incident edge can serve as a reattachment edge after deletion
```

The endpoint-options histogram equals the degree histogram in all five cores.
That is almost comically strong.  The unit spine is not fragile; it is
overdetermined.

The immediate consequence is a redirect:

```text
any positive-degree graph-only one-vertex extension of a 57-edge core
preserves a unit Hamiltonian spine
```

In particular, every degree-4 neighbor set over a 57-edge core is compatible
with some spine endpoint at the graph level.  Since a `61`-edge `n=22` graph
would be exactly a degree-4 extension of a `57`-edge `n=21` core in the
min-degree-4 branch, traceability cannot rule it out.

So the proof target has to be geometric:

```text
forbid a realizable unit-cocyclic 4-set over every exact 57-edge core
```

or use the totally-unfaithful obstruction library that already kills graph
quotients too loose to realize in the plane.

This is the same scalar-plus-side-channel lesson again.  The graph quotient
keeps adjacency, edge count, deletion decks, and Hamiltonian endpoints, but it
forgets the new point's circle center.  For `u(22)`, the missing datum is not
"can the new vertex attach to a path endpoint?"  It is "can one point in the
plane be at distance one from exactly those four core vertices without forcing
a fifth neighbor or an unfaithful subgraph?"

The next computation should therefore not count more Hamiltonian paths.  It
should enumerate, for each of the five exact cores and each degree-4 subset,
whether that subset is unit-cocyclic in some faithful embedding.  If the answer
is no for all subsets, the `57+4=61` branch dies.  If some subset is
unit-cocyclic, the proof must show its extension forces either a fifth edge
(`62`, impossible by the known upper bound), an unfaithful graph obstruction,
or a collision with the exact `n=21` embedding constraints.

Tournament Analysis after S129 has a new order:

```text
unit-cocyclic 4-set geometry
> M_L exact extension census
> endpoint universality audit
> degree-deletion core ledger
> totally-unfaithful obstruction library
> raw graph traceability
> raw edge count
```

The spine proof is now a solved side-channel for the five exact cores.  The
remaining proof lives on the circle-center/cocyclicity side channel.
