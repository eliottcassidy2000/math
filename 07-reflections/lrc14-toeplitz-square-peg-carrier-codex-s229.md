# S229 Reflection: Toeplitz Square-Peg Scale Gate

The important merge was realizing that the prompt's "Toeplitz conjecture" is
not the same as the repo's existing Toeplitz PSD matrix thread, but they meet
at the same warning: a quotient can prove the wrong thing if it forgets the
coordinate that keeps the witness nondegenerate.

For square pegs, the lost coordinate is side length.  For LRC14, it is strict
safe mass / positive scale / non-boundary status.  The square constraints give
four clean packet channels: midpoint balance, equal diagonal radius,
quarter-turn orthogonality, and positive scale.  D4 symmetry says the object is
not four unordered points but two antipodal pairs with a cyclic-order sidecar.

This fits neatly after S224, S226, and S227.  Desargues/Beal says a
rectangle-free or three-channel residual still needs incidence/common-owner
payload; Roth-Minkowski says a volume/height estimate still needs lattice and
exceptional-approximant sidecars; Moser/fibbinary says sequence/cube/simplex
motifs need native-transition, Theta, and gated-subcube sidecars; Toeplitz
square-peg says a four-witness residual still needs a noncollapse gate.  The
next real implementation should add `ordered_quad_collapse_mode` and
`toeplitz_square_scale_gate` to HYP-2963-style packet ledgers and audit every
final residual before it can be called proved.
