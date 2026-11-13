# Minimal incidence closure and signed converse transport

**PROVED IN NOTES / FINITE-EXACT after the independent audits listed in
[AUDIT.md](AUDIT.md).** No unrestricted Hamiltonian spectrum, compact new
A000568 kernel, all-order black-core connectivity, or Bott-periodic tournament
sequence is claimed. Chan–Pak and LRC(14) remain OPEN.

The full slot polynomial Q loses exactly the information relevant to where
native-complement paths land. We now have both a canonical closure operation
and concrete minimum measurement results. The native move complements every
nonpath arc of a retained Hamilton path. All counts preserve source automorphism
weights, converse fibre sizes, parallel lines and doubled loops.

## The smallest extra incidence, with the cost model declared

[The closure theorem and complete census](closure/PROOF.md) start with equality
of full Q and repeatedly split cells by their actual row sums into current
cells. The limit is the unique coarsest stable refinement. For the stochastic
matrix it is strong lumpability; for colored kernels it preserves every word
of color-and-target observations, for every source state.

| Order | Merged classes | Q cells | Raw incidence closure | Uncolored probability closure | Minimum named anchor columns |
|---:|---:|---:|---|---|---:|
| 3 | 2 | 2 | 2 | 2 | 0 |
| 4 | 3 | 3 | 3 | 3 | 0 |
| 5 | 10 | 8 | 8→10 | 8→10 | 1 |
| 6 | 34 | 24 | 24→34 | 24→34 | 2 |
| 7 | 272 | 156 | 156→272 | 156→271→272 | 9 |

The last column restricts measurements to literal columns A(u,v), or P(u,v),
for a fixed set of named target classes v. Both colored versions have the same
minimum, counting a two-output colored column as one anchor. These are exact
finite optima, not an all-order pattern or an unrestricted information lower
bound. The n7 minimum has a separate exhaustive budget-eight rejection.

At order five the intrinsic bit **blue degree equals one** already closes Q.
It splits both collision pairs: blue degrees (0,1) for codes2,5 and (1,3)
for codes8,10. Alternatively the single bit A(u,8)>0 works. At order seven,
codes5 and11 are the delayed pair: raw masses2 and1 enter distinct targets
301 and291 having the same initial Q. Normalization makes both transitions
probability one. The targets split first, then their predecessors. Retaining
source degree or edge color avoids the delay.

This answers the minimum-state question through order seven: there is no
smaller strongly lumpable partition refining Q at these orders. Linear
realizations, history-dependent models and changing-order recurrences remain
separate possibilities; [ROOT_OBSERVERS.md](ROOT_OBSERVERS.md) states the exact
algebraic distinctions.

## The signed sidecar closes the ordinary wiring

Let C be the ordinary-class dart matrix, R converse, A the merged matrix,
and choose one representative x_u from each chiral pair. The missing matrix is

    D(u,v)=C(x_u,x_v)-C(x_u,Rx_v).

For two chiral fibres the ordinary 2-by-2 block is [[a,b],[b,a]], with

    a=A(u,v)/4+D(u,v)/2,   b=A(u,v)/4-D(u,v)/2.

Thus A, fibre sizes and D reconstruct C exactly. Changing representatives
conjugates D by diagonal signs. Nonnegativity, parity and actual tournament
realization remain explicit guards. See [SIGNED_LIFT.md](pairings/SIGNED_LIFT.md).

Let Sigma join distinct chiral fibres wherever D is nonzero. A class permutation
fixing the merged vertices preserves ordinary wiring iff its flip choice is
constant on each component of Sigma. Its kernel is C2^c, where c is the number
of these components. If the merged weighted graph is rigid, the exact fixed-point
polynomial is

    I_C(x)=x^s product_i (1+x^(2m_i)),

where s counts self-converse classes and m_i are signed-component sizes.
The broader base-symmetry lift problem is a guarded sign-coboundary test:
zero support, magnitudes and diagonal entries must agree, and all ratio-sign
products around support cycles must be positive.

| Order | Signed component sizes | Ordinary-wiring-compatible involutions |
|---:|---|---|
| 5 | 1,1 | x^8(1+x^2)^2 |
| 6 | 22 | x^12+x^56 |
| 7 | 184 | x^88+x^456 |

The merged graph is independently certified rigid at these three orders.
At six and seven, ordinary incidence forces all chiral flips to move together:
only identity and global converse survive. The merged-only freedom C2^22 or
C2^184 was real freedom of the weaker data, not freedom of the full matrix.
The cycle-index formula also counts arbitrary class colorings under this action;
it does not count new tournament orders.

The [black-core boundary](closure/BLACK_BOUNDARY.md) identifies a separate
connectivity obligation. Collapse each signed component but retain positive-black-
degree self-converse nodes and every black edge. The quotient is connected iff
the full black core is connected. At order five the core is connected while
signed support is disconnected: the black path 2--8--4 joins the two chiral
fibres through a self-converse class. At six and seven the positive self-converse
sector has five and six components, each attached to the single signed component.
These are finite attachment certificates; all-order black-core connectivity
remains open.

## What actually unifies the other lanes

[Connected completion fibres and periods](bridges/PROOF.md) recover the intended
Neményi lane: isothermal principal stress nets, with separate stress and
displacement period conditions. A connected exchange graph likewise provides
paths but does not force transport around loops to be trivial. For bijective
edge transports, global sections are exactly root states fixed by every loop.
The small fixed-upper-tail graceful fibre is K4; its 64 sign transports have
8 balanced assignments and 56 obstructed ones. The same coboundary equation
controls a tested tournament symmetry lift, on its own distinct support graph.

A sheaf of sections of a twisted double cover can already be a sheaf and have
no global section. A genuine two-pass plus example needs the first pass to
repair empty-overlap data before the second glues arbitrary local labels.
Quillen plus preserves sign characters because they kill every perfect normal
subgroup. Trivial tensor amplification preserves specified flat holonomy,
although forgetting the connection may leave a trivial underlying bundle.
These are exact maps and obstructions, with primary references and hostile
examples in the bridge note.

[The Clifford boundary](CLIFFORD_BOUNDARY.md) supplies an additional exact map:
converse gives an ungraded Cl_(1,0) representation with sector dimensions
s+t and t. Its dimension defect is precisely s, the unavoidable fixed-class
count above. An invertible anticommuting negative-square generator exists
only after balancing the sectors; preserving incidence additionally requires
an intertwiner between their operators. No nonnegative class permutation can
supply that anticommuting generator. This identifies both a real connection
to Clifford representation theory and the obstruction to reading it as Bott
periodicity of tournament counts.

The [incoming reconciliation](INCOMING_RECONCILIATION.md) merges Session317's
physical consecutive-odd triples, exact one-return factorization, target inverse
and fixed-bank capacity bound. All four replays pass. Those constructions
advance the spectrum lane with coupled terminal coordinates; the universal
constructor and target-selection step remain open.

## Reproduction and next decisive questions

[The audit](AUDIT.md) distinguishes complete universes, independent algorithms,
proof reviews, candidate corrections and inherited citations. Run
`python3 05-knowledge/results/incidence_closure_20261113/verify_packet.py`
for the closed-manifest and receipt checks. It lists the heavier reproduction
commands without silently rerunning them. All checker gates survive `python -O`.

The [late incoming merge](LATE_INCOMING.md) credits Session318's concurrent
signed-lift theorem and order-six collapse. It also supplies a richer character
operator palette, a literal-pullback formulation of reusable gluing, and a
fixed-flat product theorem for graceful completions. Equal-H reservoir switches
can change Q(2), whereas our full-Q twins change incidence: the two obstructions
locate successive information losses rather than competing descriptions.

The next sharp tests are whether signed support stays connected at higher orders,
whether full-Q stable closure ever ceases to be discrete, and whether an inter-order
constructor transports the signed matrix and joint terminal response without
requiring the full class bank as input. None follows from the present finite
census. The existing earlier-value A000568 recurrence remains intact; the
new rational recurrence here closes walks on an already constructed fixed-order
incidence matrix.
