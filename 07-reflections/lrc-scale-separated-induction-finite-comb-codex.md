# LRC scale-separated induction is finite-comb topology, not pure size

The useful inductive reduction is real but narrower than the tempting slogan.
Given a seed safe set `A` at threshold `1/n`, a new speed `v` has an unsafe
comb with exact density `2/n`.  If `A` has `c` interval components and measure
`mu`, the boundary loss is at most `2c/v`, so

`meas(A cap Safe_v) >= (1-2/n)mu - 2c/v`.

This makes the large-speed branch rigorous: for a fixed smaller seed, LRC at
the previous threshold gives an open seed-safe interval at the new threshold,
and sufficiently large `v` preserves a safe point.

The important warning is equally rigorous.  Scaling a seed by `q` preserves its
safe measure but multiplies its number of components.  Therefore no induction
depending only on runner count can be uniform.  The induction invariant has to
remember safe-set topology: interval count, arc count, exact-period packet
structure, or a bounded-core atlas.

For LRC14 this makes the current proof split cleaner.  Node 2 should compact
or classify bounded/scale-normalized seeds.  Node 3 can then add committed
large speeds by finite-comb equidistribution.  The AP-core seed
`{1,...,11,13}` is the model case: at threshold `1/14` it has measure
`426/35035`, only four components, and the crude comb estimate already
certifies every added speed `v>=768`, including the radical/lcm committed
speeds `30030`, `60060`, and `510510`.

Incoming KPS-S31v is the same argument one Bonferroni layer higher.  For a
bounded core `G_C`, each large danger comb removes at most
`(1/7)meas(G_C)+arcCount/(7v)` at threshold `1/14`; the union bound closes
`r<=6` large speeds because the leading coefficient `1-r/7` stays positive.
So the one-speed induction atom and the multi-large `r<=6` lemma are not
separate mechanisms.  They are the same finite-comb topology estimate with a
positive seed floor and an arc budget.  The next genuinely different target is
`r>=7`, where the union bound turns vacuous and the second-moment /
resonant-pair defect ledger has to take over.

Incoming KPS-S31w supplies the global map around that local atom.  Remove a
scale-separated large speed, use the finite-comb estimate to keep positive
measure, and descend to a smaller LRC seed.  If a small prime is omitted, use
`t=1/p` directly.  Normalize by dilation.  What survives all three reductions
is exactly the bounded covering core.  That is the sharp limit of size
induction: it handles the unbounded and non-covering cases, but the bounded
core still needs the finite three-gap / AP-hull / Legendre-Venn extremality.

Assumption challenged: an induction quotient should have vertices equal to
runner sets of a given size.  The proof-relevant vertices are safer as
topological/analytic carriers: finite comb budget, scale-normalized seed,
Node-3 large-speed branch, Node-2 bounded core, pure size induction, and raw
runner vertices.  In that tournament, pure size induction sits below the
component-aware carriers because it destroys the only term controlling the
Weyl error.
