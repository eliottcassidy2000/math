# LRC14 forbidden-H7 state lift

The prompt "digraphs, forbidden `H=7`, binary arcs" becomes sharp only after
one correction: arbitrary digraphs are too loose.  A 4-vertex present/absent
digraph can already have exactly seven Hamiltonian paths.  The obstruction is
not binary data by itself.  It is complete binary orientation data, or the
equivalent odd-cycle conflict packet of a tournament.

The exact atom is very small.  If a connected graph `G` has `I(G,2)=7`, then
`G=K_3`: four vertices already give `1+2n>=9`, while any connected non-clique
on three vertices has an independent pair and value at least `11`.  So a
tournament with `H=7` would have to have conflict graph `Omega=K_3`.  THM-343
blocks `H=7` for all tournaments, and THM-201 blocks `K_3` as a connected
component of `Omega`.  This is not a subgraph ban: THM-344 realizes `H=63` at
`n=8` with `Omega=K31`, which contains many `K_3` subgraphs.

That is the usable bridge to LRC14:

```text
counterexample -> binary packet conflict graph with I(.,2)=7 -> K3 -> impossible.
```

The whole issue is the first arrow.  The vertices cannot be raw runners.  They
should be the finite atom's primitive obligations: sector-pair covers,
wall-crossing events, exact-period phi packets, support-six relation packets,
or cover-arc packets.  The quotient must keep exactly the information that
OCF keeps: conflict incidence plus two side states per primitive packet.

This ties into HYP-2906 and HYP-2905 cleanly.  HYP-2906 peels any one speed
larger than `13` times the second speed, so a counterexample entering this
finite atom is top-balanced.  HYP-2905 says the remaining object is a boundary
state, not a runner set.  HYP-2908 now says what a successful finite atom
theorem should look like: convert the apex-7 over-cover boundary state into
the forbidden connected `I=7` conflict atom.

The route is therefore proof-oriented, not disproof-oriented.  If the state
lift lands in arbitrary digraphs or even graphs, it fails because those
categories realize `7`.  If it lands in tournament conflict packets, it proves
that an LRC14 counterexample cannot exist.
