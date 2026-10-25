# Tangent-Legendre CM uniqueness certificate

Read RESEARCH.md for the complete proof and scope. The claimed all-integer result is: for integer m>4, the Legendre curve with parameter tan^2(pi/m) has CM iff m=8. This uses standard cited CM results plus bounded exact verification; it does not concern a modular formula for the MRB constant or defect.

Run:

    python check.py > results.normal.json
    python -O check.py > results.optimized.json
    cmp results.normal.json results.optimized.json
    python test_rejections.py > rejections.normal.json
    python -O test_rejections.py > rejections.optimized.json
    cmp rejections.normal.json rejections.optimized.json
    python verify_manifest.py

Requires Python and SymPy. Both normal and optimized checks use exact arithmetic. Main certificate includes 30 proven-complete candidate conductors, all minimal polynomials, a separate requested m5..24 box, rational Taylor bounds, and every primitive reduced form for discriminant absolute value below32. Mutation controls independently damage the conductor universe, rational map, analytic threshold, and CM positive control and require rejection.

SOURCES.json identifies primary literature, precise claims used, non-used claims, and SHA256 receipts for privately inspected source PDFs. This text-only publication edition does not include third-party articles; bibliography entries are provenance, not runtime dependencies. MANIFEST.json freezes the exact recursive payload file set, including inherited source bytes and their hashes. The verifier rejects missing, altered, and unlisted nested files. Independent audit is separate; producer cannot self-certify the independent verdict.
