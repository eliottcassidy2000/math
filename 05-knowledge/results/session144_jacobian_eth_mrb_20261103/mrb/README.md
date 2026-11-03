# MRB: what survives a representation, and what determines a constant

**Session 144, 2026-11-03.** Source snapshot `origin/main=5e1cedf` in the read-only source checkout. This packet is a research result, not a new canon theorem number or a publication action.

**Main result:** even a complete branch-jump profile, any fixed finite prefix, and every algebraic asymptotic coefficient uniformly on the closed right half-plane fail to determine the alternating constant. The new smooth-Laplace perturbation strengthens Session31's obstruction beyond the global decay sidecar that earlier note left available.

**Other results:** a complete periodic-filter convergence criterion, an exact conversion between cutoff and Abel finite parts, and an explicit precision-controlled evaluator of the fixed MRB constant. The latter does not decide SAT or imply any ETH consequence. The first two classify a specified finite-dimensional family of summation operators, not all formulas for MRB and not polynomial maps.

Files:

- `RESEARCH.md`: proofs, inheritance, source links, status and connection ledger.
- `check.py`, `results.json`, `check.out`: independent exact/interval/numerical tests, each labelled.
- `test_rejections.py`, `rejections.out`: mutation controls that remain live under Python `-O`.
- `audit/`: independent analysis and precision-evaluator audit with its own reproducible script and receipts.
- `inherited/`: the two original Session31 scripts fetched by targeted `git show`; replay receipt for the first. Their historical assertions are Python assertions and are not claimed to remain live under `-O`.
- `compare_certificates.py`, `certificate_comparison.json`: exact rational containment checks against the freshly replayed inherited interval.
- `SOURCES.json`, `MANIFEST.json`: source/status ledger and file hashes.

Run from this directory:

    python3 check.py > check.out
    cp results.json results.normal.json
    python3 -O check.py > check.optimized.out
    cmp results.normal.json results.json
    cmp check.out check.optimized.out
    python3 test_rejections.py > rejections.out
    python3 -O test_rejections.py > rejections.optimized.out
    cmp rejections.out rejections.optimized.out
    python3 compare_certificates.py
    python3 -O compare_certificates.py
    (cd audit && python3 exact_precision_audit.py --precisions 1 8 16 32 --output exact_precision_results.json)
    (cd audit && python3 -O exact_precision_audit.py --precisions 1 8 16 32 --output exact_precision_results_optimized.json)
    python3 audit/compare_receipts.py

Proof novelty: the results are independently derived, with prior art explicitly separated. The new obstruction and universal mask/finite-part synthesis are new relative to the recovered Session31/52 lineage; no literature-wide originality claim is made.

Final executed scope: the new checks pass normally and with `-O`; the precision evaluator certifies p=1,8,16,32 with exact rational elementary enclosures. All-p precision complexity and the smooth-bump theorem are analytic proofs, not finite-run extrapolations. No p=64 execution certificate is claimed.
