# Independent audit: PASS, conditional on inherited census

Reviewed RESULT.md and final check.py. All algebraic claims and scope qualifications are sound. This audit does not replay tournament-class enumeration; it independently checks the cited receipts, transformations, and coefficient identities. No packages were installed, and neither deps nor recover_inventory.py was executed by this audit.

## Source checks

- Independently extracted A17=2z^6, A31=2z^6+2z^7, A89=14z^7+7z^8+2z^9 and self-converse C89=z^8 from the inherited inventory375.
- Confirmed index76 of the inherited odd375..623 totals is184.
- Read the targeted earlier b3_ear_census_625.out through git show because that tracked file is absent from the sparse worktree. It explicitly confirms all125 totals entry by entry and the strong-order cutoff.
- Independently counted22 rows with H=527 in n8_classes.csv; all22 are strong.
- Read the inherited SCC proof and checked the exact ordinary-isomorphism convention. Neither automorphism quotient nor converse merging is silently substituted.

## Mathematical checks

Unique ordered SCC decomposition yields exactly one527 atom or an ordered17/31 pair, giving A527/(1-z)^2+8(z^12+z^13)/(1-z)^3. The complete census assumption gives S=184 but leaves the first moment W unknown. No exact constant is claimed without W.

Coefficient extraction proves N527=8n²+8n+1152-W for n>=11. The zero boundary contributions at n=11 are correct; positive decomposable contribution starts at12. Bounds1634<=W<=2120, forward difference16(n+1), and second difference16 are correct. The onset is guaranteed, not claimed minimal.

The prime89 specialization gives N89=23n-149 for n>=8, with an actual mismatch at7 proving minimality. The self-converse parity and converse-merged formulas follow from the unique prime center.

For the quadratic sequence modulo72, any eventual period d forces72|16d by subtracting adjacent increment conditions, hence9|d; d=9 works. The linear89 sequence has least period72 because23 is a unit modulo72. The unknown W only changes a constant offset, so neither period depends on recovering it.

## Tests

Final author script:2863 checks, normal and -O outputs byte-identical.

Additional independent standard-library script builds generating-series coefficients by repeated cumulative summation, rather than reusing the author's binomial formula. It checks1000 seeded admissible coefficient distributions,89000 coefficient equalities, forward differences, the89 boundary, and every candidate period1..72. Both execution modes pass identically. See independent_audit.py and independent_audit.json.

The randomized coefficient distributions are algebraic test data, not claims that such census distributions are realized by tournaments. The proof uses linearity, so testing four extremal placements in the author's script is sufficient as a control and is not described as an enumeration of all distributions.

No material correction is required. The exact W527 remains an explicitly open data-extraction task.
