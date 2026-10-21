# Independent review of actual-H multiplicities

APPROVED mathematically and computationally. This review concerns actual H, not the quotient H/|Aut|. No full odd-spectrum completeness theorem is assumed.

The SCC-word coefficient formula, square-divisor center formula and merged Burnside formula have been checked independently. Highest poles cannot cancel because their core coefficients are positive. The exact degree is maximum strong factorization length, not arithmetic Omega in general. Every self-converse half contributes twice to a full core, establishing the lower-order symmetry bound. Parity and absent-center qualifications are necessary and correctly retained.

Busch's primary publisher abstract independently confirms H(S)>=5^((n-1)/3):
https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3
The DOI endpoint was inaccessible to the lookup tool, but the publisher search result and primary PDF mirror agreed:
https://www.kurims.kyoto-u.ac.jp/EMIS/journals/EJC/Volume_13/PDF/v13i1n3.pdf
Only the stated lower bound is used, not an equality classification.

The logarithmic onsets were independently derived: a shifted binomial agrees with its zero-truncated version already at n=M-r. Summing Busch bounds over factors gives M-r<=3log_5H. The fixed-center parity version adds at most one. These are safe upper bounds, not minimal-onset assertions. The resulting coarse H^{O(log H)} algorithm is for a fixed actual-H multiplicity series, not unrestricted A000568.

The exact two-degree loss along7*3^a, the leading coefficient2/(a-2)!, the7-power rays, and the odd-exponent self-converse refinement were checked. The even ray is only an upper bound, as required. The top-core metagraph-degree consequence is correctly restricted to the inherited fundamental paired-black convention and does not claim adjacency. Final text distinguishes automorphism-free atoms from path-rigid h=1 atoms.

Computational review:
- Read the complete C++ canonicalization, augmentation, automorphism counting, reachability and Hamilton DP implementation.
- Recompiled from source and reran all58,368 augmentations; resulting6880-row file is byte-identical.
- Independently checked every output is a valid tournament, all canonical codes distinct, all group orders odd and dividing8!, and all Hamilton counts divisible by their group orders.
- The sum8!/Aut=2^28 certifies coverage of the entire labelled universe once canonical uniqueness and valid automorphism multiplicities are established. The fixed-path mass sumH/Aut=2^21 is a separate control.
- Producer checker replayed normally and under-O:2063 checks, byte-identical receipts. Added mass/validity verifier also matches in both modes.
- Exact cutoff5^8>73^3 makes the strong atom inventory globally complete for oddH<=73. No unknown large-order atom is silently omitted from those rational functions.

The H39 row correctly includes the3*13 decomposition and is quadratic. H63 has two automorphism-free non-self-converse strong atoms at order8; therefore ordinary multiplicity2(n-7), self-converse multiplicity0, merged multiplicityn-7. H49 has eight strong atoms and the stated linear law. All fixed-H laws are order-indexed: every attained H occurs infinitely often across unbounded order by singleton padding.

Reproduce the independent full augmentation with a standard C++17 compiler, then run review.py; run the packet's check.py normally and under-O. The regenerated n8 CSV is deliberately not republished as a duplicate; its digest is recorded in mass-receipt.json.
