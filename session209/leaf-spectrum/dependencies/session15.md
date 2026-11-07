# The deficit-two polynomial is D(z)=z^4

Research date: October 20, 2026. Status: PROVED with complete finite certificates and two independent audits.

The unique strongly connected, self-converse tournament with exactly two black Hamiltonian-path orbits is the strong four-vertex tournament S4. Thus the coefficient of z^4 is 1, every other coefficient is 0, and the polynomial's degree is 4.

The proof does not extrapolate from a small-order census. It first proves an effective all-order reduction to seven possible new orders, then reconstructs and excludes every one of 2,121,984 necessary symmetry-normalized presentations by an explicit prefix obstruction. The finite spaces and all individual certificates were independently reconstructed and replayed.

Start here:

- [Complete theorem and exact census](session15-deficit-centers/DEFICIT_TWO_COMPLETE_CENSUS.md)
- [All-order cutoff and structural proof](session15-deficit-centers/PREFIX_CONNECTIVITY_CUTOFF.md)
- [Short hand proof of the local configuration](session15-deficit-centers/HAND_LOCAL_CONFIGURATION_SUPPLEMENT.md)
- [Exact leaf series and minimal recurrences](session15-deficit-centers/EXACT_LEAF_SERIES.md)
- [Independent cutoff audit](session15-deficit-cutoff-audit/AUDIT.md)
- [Independent full census audit](session15-deficit-census-audit/AUDIT.md)
- [Sources and hostile boundaries](session15-deficit-centers/SOURCE_AND_SCOPE_AUDIT.md)

The preserved theorem manuscripts were frozen while their audits were being written, so their opening historical status lines refer to pending review. The approval reports and this navigation record give the current final status; the mathematical text is exactly the audited version.

The complete paired-black leaf generating function is now

    L(z)=(2z^4-z^5-z^7)/[(1-z-z^3)(1-z^2-z^6)].

Its reduced denominator has degree 9, so the minimal constant-coefficient eventual recurrence has order 9. Mixed leaves have generating function z^4/(1-z^2-z^6), with minimal recurrence order 6. This closes the mixed-leaf center classification, not the separate classification of nonleaf bridges.

Reproduction: use the producer directory's README. Before independently replaying the compressed certificates, either run the producer scanner to create the raw five-byte streams or decompress them with gzip. The independent affine reconstruction expects the three retained directories to be siblings, as packaged here. No compiled executable is required; C++ sources and Python scripts are included.
