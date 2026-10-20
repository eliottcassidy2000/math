# Strong tournaments of bounded black defect form a finite family

Publication update, 2026-10-20: This is a preserved precursor. The later [complete deficit-two center classification](../deficit_two_centers_20261020/README.md) proves D(z)=z^4, and a separate [effective order-bound proof](../effective_black_defect_order_bound_20261020/README.txt) supplies a cubic bound. The unknown-center and non-effective language below records this earlier packet, not the current theorem boundary.

PROVED from cited published theorems: for each fixed d>=0, only finitely many strong tournaments have Delta=H/|Aut|-b<=d, where b counts literal-blue Hamiltonian-path orbits and is zero for non-self-converse tournaments.

The proof restores the missing degree coordinate from cycle-rotation orbits, then compares Thomassen's factorial cycle lower bound with Dixon's exponential automorphism upper bound. It is not a classification or a census extrapolation.

Consequences:

- Strong self-converse black defect tends to infinity uniformly with order.
- The strong deficit-two center series D(z) is a finite but still unknown polynomial.
- MixedLeaf(z)=D(z)/(1-z^2-z^6), and the exact total leaf series is rational.
- Total leaf counts are asymptotic to half the number of words in atoms of sizes 1 and 3, independently of the unknown finite center alphabet.
- Bounded paired-black quotient degree forces bounded order within the strong core, including non-self-converse nodes.
- None of this classifies nonleaf bridges or proves D(z)=z^4.

Read `BOUNDED_BLACK_DEFECT_FINITENESS.md` for the full proof and `SOURCE_AUDIT.md` for precise source scope and retrieval limits. Independent exhaustive checks and their reproduction command are in `independent/INDEPENDENT_REVIEW.md`.

The computations use no external package and are validations only. Publication destinations and hashes are in `PUBLICATION_MANIFEST.json`.

Formal-series controls: `python3 check_series_consequences.py` verifies the recurrence, concatenation inequality, and asymptotic error bound through order 1000 for the known S4 control and for a hypothetical finite alphabet containing extra even and odd terms. The hypothetical terms are tests, not newly claimed centers.

The retained census payload is compressed as `independent/verification_results.json.gz`; decompression yields the exact JSON produced by `python3 independent/verify_cycle_bounds.py --max-order 6`. The console output and full independent review are uncompressed.
