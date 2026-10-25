# Session55 final integration audit

25 October 2026. Verdict: PASS within the stated mathematical and publication scope. No unresolved mathematical blocker. Frozen producer files were not edited. This audit performs no remote publication and makes no priority claim.

## New dyadic theorem: independent proof and controls

The literal array identities for M4, M8 and M16 are correct. Modulo 4 the binary-valued line sums are integers between zero and three, so modular equality is integer equality and the opposite-pair identity forces constancy. Reduction modulo 4 and division of line-sum differences then give exactly 4T or constant 1 modulo 8, and 1+8T or constants 0,4 modulo 16. This proves equality of the FULL realizable collision-relation sets at 2,8,16, not merely the minimal antichains.

Consequently CRT proves Parker(m) iff Parker(4m), and Parker(2m) iff Parker(8m) iff Parker(16m), for every odd m, including those outside the finite odd Parker list. Neither assertion requires arbitrary lifting across noncoprime factors.

The AP constructions are valid: odd m>=5 cannot divide 24 or 48; 1,9,17 are distinct square units modulo 2^a for a>=5. Complementary fixed-label three-block partitions have singleton intersections. The root-bit argument establishes existence of square roots at every dyadic level; it need not preserve an arbitrary chosen old root. The Lo Shu construction handles a>=7. The finite 32/64 negative tables and the odd classification complete the seven-layer proof.

The displayed modulus-50 witness and all nine roots check directly. Any square-magic lift to 100 would have constant mod-4 reduction, hence constant parity, contradicting this witness. Nonexistence of any other distinct witness modulo 100 additionally uses Parker(25); this is properly distinguished from failure to lift one array. The six divisor counterpairs follow from the complete list, not from a bounded heuristic search.

Our audit.py independently enumerates ALL ambient triples (c,u,v), rather than the producer's square top-row/center reconstruction. Counts at 2,4,8,16,32,64 are 8,2,9,10,73,522, with mask counts 4,1,4,4,10,20. Literal M4/M8/M16 identities, full relation equality, 32/64 obstruction and the 50-to-100 parity failure all pass. Normal and optimized outputs agree. Producer replay also passes in both modes, checking its 203 candidates and six counterpairs.

## New ordinary generating function

The gcd-to-indicator expansion is exact. Substitution into the standard power-sum series gives the nine rational summands with the correct substitution factors. Independently reducing them with SymPy gives the stated numerator 8z^9 Q(z) of degree 26 and denominator of degree 27, with coprime numerator and denominator. We independently recomputed all 247 rank-two Smith-subset inclusion-exclusion contributions, recovering the source count coefficients.

The poles have exact orders 4,3,2,2,2,2 at primitive root orders 1,2,3,4,5,6. In particular contributions at primitive third roots add, rather than cancel. Primitive fourth-, fifth-, and sixth-root poles force every quasipolynomial period to be divisible by 60; the indicator expression supplies period 60. No primitive 60th-root pole is needed or present. The source theorem counts unrestricted residue arrays, not square-entry arrays, roots or D4 classes. The new proof and root summaries preserve that distinction.

The producer's standard-library rational-arithmetic replay passes normally and with -O. Our independent symbolic reduction uses SymPy 1.14.0; this is an additional dependency of this audit only, not of the producer replay.

## Summary and dependency review

The final README and SYNTHESIS were checked against producer theorems and the separate finite-field, integer-search, lifting/score, odd-ring, p-adic, modular-count, weighted-density, and collision-sieve audits. Their quantitative statements and logical scopes agree. In particular:

- The 136 total and 28 odd Parker-ring counts restrict to n>=2. Including modulus 1 gives 137 and 29. Square residues include zero; labels and distinct residue values are retained.
- The all-modulus proof uses prime-field classification only at exponent one, elementary prime-power constructions/obstructions, the odd AP CRT construction, and the 203-candidate exhaustion. It does not infer global classification from a scan to 3216.
- The finite-field classification relies on a cited Weil bound and separately audited bounded census. This audit does not rerun the 1730-field census or reprove Weil.
- Cain's Section 9 was freshly retrieved: Conjecture 9.3 states odd maximum 201, and Conjecture 9.4 states unrestricted maximum 3216. The root correctly identifies these dependencies and makes no claim of first proof. Source: https://arxiv.org/html/1908.03236v2#S9 . The audit does not certify novelty from an unsuccessful literature search.
- Prime 59 concerns the weak canonical residue/no-wrap predicate. The integer problem remains open. The root's final minimum-prime wording correctly says the smallest integer center determines the smallest faithful prime, without claiming injectivity or an unproved short-prime-gap theorem.
- Square-count series and densities, ordinary counts, weighted CRT masses and position orbits remain separate. The field asymptotic is uniform over odd fields and follows from character sums, not probabilistic independence. The ratios 64/63 and 1728/91 are arithmetic comparisons, not MRB/CM/Monster identifications.
- Residual score feasibility preserves named endpoints and has the sharp n=4 cycle exception. No all-size native metagraph wiring, A000568 prediction, integer-square solution, search record, or external priority is asserted.
- The collision filter's bounded rejection count is accurate, and the root explicitly declines a speedup claim.

## Diagrams

Both SVGs were rendered with Inkscape and independently viewed. Labels, legends, contrasts and warnings are readable with no observed overlap or clipping. In the displayed CRT grid the mod-5 partition A is the fibre partition of i-j, and mod-7 partition B is that of i+j+1, up to colour relabelling. Our checker verifies equivalence relations on every pair of fixed cells, so swapping the attributions would fail. The nine ordered coordinate pairs are distinct. The integer row sums 10,45,80 are correctly contrasted with their common residue 10 modulo 35.

The field/ring diagram correctly distinguishes orders 27,81,243, the characteristic-two field obstruction, the dyadic ring threshold 128, and 50 non-Parker / 100 Parker / 200 non-Parker. Its independent 103-modulus controls, extension-field controls, deterministic replay and self-tests pass normally and with -O. Universal statements depend on pinned proofs, not the plotted finite examples.

## Integration and reproducibility

The definitive preaudit map has SHA256 1c3503e3c85f45b22cf078207e505d167faad2f8a02a17ec7ddb907207cc29cc and 82 files. Every mapped source byte length and SHA256 passes. Every relative Markdown link in both root documents resolves against that map plus the first/second publication plans, with referenced prior frozen files independently hash-checked. The parent's suggested clean-repository snapshot did not contain the Session55 directory; the check therefore uses the actual frozen-first/frozen-second source files named by the verified publication plans. This was a local path mismatch, not a broken publication link.

All seven dependency pins of the two new proof packets match their exact source files. Source hashes, including the definitive root files, appear in SOURCE_PINS.json. Existing published targets are not changed by this audit. Publisher path-only relocations, final tree preservation, commit creation and remote readback remain the publisher's responsibility and are not claimed complete here.

Commands from the workspace root:

    python3 -B session55-audit-final-integration/audit.py
    python3 -O -B session55-audit-final-integration/audit.py
    python3 -B session55-audit-final-integration/check_integration.py
    python3 -O -B session55-audit-final-integration/check_integration.py
    python3 -B session55-dyadic-reentry/verify.py --verify
    python3 -O -B session55-dyadic-reentry/verify.py --verify
    python3 -B session55-ordinary-count-gf/verify.py
    python3 -O -B session55-ordinary-count-gf/verify.py
    python3 -B session55-diagrams/replay.py --self-test
    python3 -O -B session55-diagrams/replay.py

Normal/optimized arithmetic results are byte-identical. audit.py also supports the published sibling name diagrams. check_integration.py is explicitly a prepublication workspace check using the retained publication maps, not a standalone public-tree checker. Publication should replay audit.py after assembly and verify relocated links/pins with its own staged-tree checks.
