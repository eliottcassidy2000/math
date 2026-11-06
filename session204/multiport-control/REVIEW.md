# Independent audit of the multiport-control release

Status: PASS. No substantive mathematical or implementation issue found.

Reviewed `session204-multiport-control/release/RESULT.md`, `verify.py`, and `literal_forest.py`. Both normal and optimized verifier runs pass, reporting formal dimension 210, full reachable dimension 8, observable quotient dimension 4, and seven literal controls. All three manifest hashes pass.

## Physical semantics and exact certificate

The seed rows define an actual tournament. Every digit introduces one vertex with exactly the declared incidences; forgetting/relabeling changes no graph edges. Both closing introductions beat ports zero and one and lose to all other vertices, including the previously introduced closing vertex. Root port two persists.

The unordered endpoint-pattern transfer correctly enumerates predecessor/successor component choices, excludes joining a component to itself, and counts multiplicities. Discarding more than two bulk starts is safe because future vertices cannot precede a bulk start. The final readout counts exactly h, root-start s, and one/two-component root-ending L; there is a unique root-ending component, so no ordering multiplicity is lost.

The independent forest initializer uses canonical subset partitions with the least remaining vertex to avoid overcounting unordered components. Its misleading legacy docstring mentions two active ports, but its implementation is general and the actual seed call passes three. This is cosmetic only.

The rational span algorithm includes each independent encountered state and schedules both children. Dependent children need not be expanded, by linearity. Hence closure at termination is an all-word invariant-subspace certificate, not a finite-word extrapolation. Exact rational identities `(R Mi - Ai R)V=0` prove the quotient dynamics; rank four does not imply the full eight-dimensional reachable forest state is recoverable from four readouts. The report states this distinction correctly.

The fourth coordinate is the closing start count s, initially one and preserved by both quotient matrices. Therefore the constants in the displayed affine formulas are valid. Integrality follows because D,u,v,s are integer forest readouts, notwithstanding fractional displayed coefficients.

## Mixed words and correction obstruction

The verifier proves annihilation by M_i^2-3M_i+2I on the entire repeated-digit invariant span, for each mixed closing row D M_j. Together with the two checked initial terms, this rigorously establishes both mixed-run formulas.

For nominal N0=90 and ordinary updates N->2N+d, word 0^n1 has N=180*2^n+1. Subtracting 6N+23 from 1125*2^n-13 gives 45*2^n-42. Word 1^n0 has N=182*2^n-2; subtracting its nominal readout from 1050*2^n+53 gives 42*(1-2^n). These agree with the report.

A finite-state additive correction with fixed real value per state is bounded. The first discrepancy is unbounded, so such correction cannot make these unchanged physical operations and unchanged closing readout implement the declared nominal encoding on all words. This does not exclude changed closing contexts, changed operations, unbounded corrections, or alternative/nonlinear encodings. The release states these restrictions sufficiently clearly.

## Additional independent controls

`independent.py` uses a separate ordinary subset endpoint DP, with L obtained by summing `end_root(S)*H(complement(S))`, rather than unordered forest partition enumeration. It validates all orientations and compares B,s with the transferred closing readouts on every word of length at most three: 15 cases, all PASS.

This packet proves an exact infinite physical family and a candidate-specific failed-encoding obstruction. It proves no tail occupancy, F(8,137), or eventual spectrum threshold.

## Release packaging changes

After this mathematical audit, the legacy initializer docstring was corrected to describe arbitrary supplied retained ports. The independent checker now locates verify.py relative to its own file, so the package is portable. Neither change alters mathematical logic. SHA256SUMS refers to the packaged files, not the original pre-packaging manifest.
