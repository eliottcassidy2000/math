# The deficit-two center polynomial

Main conclusion, independently APPROVED: D(z)=z^4. The unique strong self-converse tournament whose black Hamiltonian-path-orbit deficit is exactly two is S4. The degree is four; coefficient d4 is one and every other coefficient is zero.

Read in this order:

1. DEFICIT_TWO_COMPLETE_CENSUS.md: complete argument, exact finite spaces, coverage proof, certificate specification and result.
2. PREFIX_CONNECTIVITY_CUTOFF.md: independently approved all-order reduction, half-supported permutation obstruction, regular orbit connectivity and degree28 cutoff.
3. EXACT_LEAF_SERIES.md: exact rational leaf series, coefficient formulas and recurrence orders/onsets.
4. Independent cutoff and residual audit reports in the two sibling audit directories.

The decisive new idea is not a larger unrestricted tournament census. A triangle of half-supported even-order reroutings already requires at least four black path orbits. This forces any hypothetical deficit-two center to have connectivity at most four. Hamilton-cycle symmetry supplies large regular rotation blocks and the opposite connectivity bound, leaving just seven orders. Symmetry-normalized finite presentations at those orders all have explicit four-cycle prefix obstructions.

Exact finite coverage: 2,121,984 presentations, no survivors. They are not claimed to be pairwise nonisomorphic classes. No random sample replaces an exhaustive exclusion. All presentation witnesses fit into about188KB compressed.

Reproduce the finite computation from this directory:

    python3 verify_local_configuration.py
    python3 generate_residual_frames.py
    g++ -O3 -std=c++20 scan_residual_frames.cpp -o scan_residual_frames
    for f in frames_n*.txt; do ./scan_residual_frames "$f"; done
    python3 check_residual_samples.py
    python3 check_leaf_series.py

The raw witness files are emitted by the C++ scan. Their compressed checked versions and SHA256 index are retained. residual_all_scan.jsonl records the complete census; residual_sample_checks.json records separate invariant/predicate controls, not the main exhaustion.

Inherited inputs are the already audited cycle-cut lemma and the exact orders through nine. The new external source is Moon's vertex-pancyclic theorem, verified in its primary1966 paper. No factorial bound, unpublished theorem, assumed connectivity of the metagraph, or whole-fiber isomorphism enumeration is used.

Scope: this determines the strong deficit-two center alphabet and hence all mixed paired-black leaves. It does not classify all nonleaf bridges or prove global metagraph connectivity. Staging and publication are add-only; the historical session14 proof packets remain unchanged.
