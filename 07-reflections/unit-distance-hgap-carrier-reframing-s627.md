# Unit-Distance H-Gap Carrier Reframing S627

The best way to make the unit-distance tournament analysis as meaningful as
the LRC tournament analysis is to stop making the points do all the work.
Points are too raw.  They know the metric, but they do not remember which proof
obligation a unit edge is serving.

S627 uses carrier packets instead:

- a unit Hamiltonian spine;
- the tile/bulk unit edges beyond that spine;
- the boundary/frontier recursion;
- the exact non-lattice/Moser carrier;
- the OCF `H=7`/`H=21` guardrail;
- the round-LRC worry-set channel;
- the literal unit-flip tournament `H`;
- the raw edge count quotient;
- the equidecomposability ledger.

The incoming S626/HYP-2203 Moser unit-spine packet fits this perfectly.  It
says the unit spine and tie-order witness are real side channels in the
non-lattice Moser lane.  S627 asks what happens when that kind of carrier data
is crushed into a scalar like `7` or `21`.

That reframing makes the H=7 question click.  Exact unit distances have
`u(5)=7`, but the literal unit-tournament count in the S599z spine model is
`H=15`, not `7`.  The visible `7` is legal because it is an additive edge
count already split as

```text
7 = 4 spine edges + 3 tile/bulk edges.
```

The tournament impossibility is about a different carrier: a single
Hamiltonian-path/OCF evaluation cannot be `7`.  So `u(5)=7` is not a
counterexample or even a near miss.  It is a warning label: the quotient that
sees only the scalar `7` has forgotten the scissors decomposition.

The `21` echo is subtler.  The triangular/Harborth lattice row has `21` unit
edges at `n=11`, split as

```text
21 = 10 spine edges + 11 tile/bulk edges.
```

But the exact small planar value is already `u(11)=23`.  So `21` is a
triangular-lattice carrier echo, not the global unit-distance optimum and not a
literal tournament with `H=21`.  That fits the current repo guardrail:
`H=7` and `H=21` are permanent tournament gaps, while unit-distance edge counts
live in an additive geometric carrier where those same integers can appear
after decomposition.

This is also why the equinumerosity/equidecomposability language is useful.
The raw equality `7=7` is equinumerosity.  The proof content is
equidecomposability: which pieces, side channels, and allowed operations
realize that number?  Tournament `H`, LRC worry sets, and unit-distance edge
counts can share integers while refusing to share decompositions.

The Tournament Analysis bears this out.  Under pure geometry, the unit spine
and frontier shell lead.  Under the H-gap gauge, the OCF guardrail,
equidecomposability ledger, and round-LRC channel lead.  Under scaling, the
Moser/non-lattice exact carrier leads.  Each single gauge is transitive, but
the majority tournament has a large SCC, five directed 3-cycles, and 25
Hamiltonian paths.  That is not noise; it is the proof-route tension we want.

So the useful theorem-shape is not "the unit-distance tournament has forbidden
H values."  It is:

> A unit-distance quotient that collapses spine, tile, frontier, and exact
> carrier data into a raw scalar can falsely ask for a forbidden tournament
> H-value.  The repair is to retain the carrier decomposition.

That is close to the LRC lesson: raw residues are not enough; the owner,
carry, pinch, shell, and witness side channels decide whether the quotient is
conservative.
