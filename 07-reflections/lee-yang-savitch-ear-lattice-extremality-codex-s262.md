# Lee-Yang, Savitch, Bravais, And Ear-Lattice Extremality

HYP-3108 is useful because it refuses to turn the new PGF-root signal into a
new scalar.  The S262 scout produced two maps.

The first map is a packet map.  Each named row carries its full miss-count law
`q_t`, the roots of `G_N(z)=sum q_t z^t`, a quartic `phi4` fit, `q0`
component geometry, co-emptiness Perron gap, Bravais q-lattice address, and a
cell-state transition graph.  The main readout is sharp: `consec_8` remains in
the all-complex PGF stratum with nearest root `|z|=1.489`, while `break_8`
has a real-root collision with nearest `|z|=0.121`.  Root confinement is not a
renaming of `p0`; it preserves a coordinate that scalar coverage forgets.

The Bravais readout is equally important.  The named packets have full affine
rank `6` in the `q_t` simplex.  That means a proof that uses only one moment
or one cap value is discarding real coordinates.  The discard may be legal, but
it has to name the sidecar that restores the lost information.  This matches
the lesson of HYP-3106: controlled forgetting is legal only when the next
operation's required sidecar is retained.

The bounded-bank appendix corrects the first naive lattice guess.  On all
`{0}+7` rows from `1..13`, high `p0` correlates negatively with the largest
mod-7 Bravais peak and positively with residue entropy.  The lattice signal is
therefore reciprocal flatness plus nonlinear sector correlation, not a
crystalline Bragg peak.

The second map is a proof-frontier map.  HYP-3107 exposes open Lean fields;
Savitch's theorem suggests looking at reachability in that proof-state graph
through recursive midpoint certificates.  In the scout, paths from raw packets
or PGF roots to the terminal LRC14 interface repeatedly have midpoints at
`observer_gluing_certificate`, `bravais_lattice_address`, or
`finite_address_packet`.  That is a useful compression of what remains: the
root signal can suggest where to go, but the proof still needs a certified
observer/coverage packet at the midpoint.

The ear-decomposition lens gives a new tournament use that is less brittle
than raw H=7/H=21 analogy.  Instead of asking whether a constructed LRC
subproblem "is" a tournament with forbidden H, record the cell-state transition
digraph and ask what ear-growth data is needed to glue it into an
observer-certificate.  The state graph is naturally directed and cyclic; its
ear excess is a sidecar for gluing complexity, not a terminal contradiction by
itself.

The bolder guesses are still guesses.  The quartic `lambda` from the
`exp(-lambda S^4 - b S^2)` fit is sometimes positive and sometimes negative,
so it is not yet a monotone extremality parameter.  The close apex-7 root-angle
gaps in the wide k=11/12 rows are suggestive, but they are not a zero-free
region theorem.  The Savitch route is a map of proof-state reachability, not a
replacement for the missing coverage extremality and observer-gluing theorem.

Tournament Analysis used sidecars as vertices.  It ranked
`savitch_reachability`, `ear_decomposition_state_graph`, and `pgf_zero_curve`
ahead of `raw_scalar_p0`, with no directed 3-cycles and one Hamiltonian path.
That transitivity should not be overread as proof.  It says only that, for the
current frontier, these sidecars retain more of the payload that HYP-3107 says
future formal proof nodes will need.
