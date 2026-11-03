# Independent synthesis audit

**Verdict: PASS. Both precision edits and the direct source link were applied and verified.**

Compared NOTE.md against the exact five source packets:

- `/tmp/s173_clocks/PROOF.md` and its independent audit;
- `/tmp/s173_squares/NOTE.md` and `/tmp/s173_squares_audit/AUDIT.md`;
- `/tmp/s173_runner/result.md`, including its independent audit section;
- `/tmp/s173_spectrum/NOTE.md` and its independent audit;
- `/tmp/s173_527/RESULT.md` and its independent audit.

This is a cross-packet statement/provenance audit, not another execution of all five finite searches. The clock and spectrum packets were separately developed or independently tested by this auditor earlier in the session.

## Findings

1. **Marked-cycle provenance:** the exact inherited note section3 agrees: 13→2→5→3→7→4 requires a declared closing identification, and the full valuation output fails to descend modulo9 at11 and2. The chosen five-cycle extension and affine clock are correctly identified as constructions rather than recovered canonical dynamics.
2. **Group statements:** A9, C2×A9, five sector sizes, affine cycle type and positive-integer nonperiodicity agree with the clock proof. No unsupported Monster or tournament conclusion is imported.
3. **Square-sum statements:** the 216 displacement, height44 injectivity, nonsquare28 control, sharp placement bounds, 36-vertex witness, connector paths, minimum number of additions and conditional minimum maximum addition317 all match the proved packet. No all-n search termination is inferred from the finite path certificates.
4. **Runner statements:** masks, tree/cycle defects, phase translation condition and composite moduli agree with the proved proposition. The eight-sheet cycle bound is only inconclusive, so the text should say the alternative tree succeeds rather than that changing to a tree is mathematically required.
5. **External LRC source:** independently reopened [Allikvere, arXiv:2609.02604v2](https://arxiv.org/html/2609.02604v2) during this audit. The page identifies Jaan Allikvere, version2 dated24 September2026, states the k-nonzero-speed/k+1-total-runner convention, and gives Theorem1.1 for fourteen and fifteen total runners. The synthesis correctly labels it a computer-assisted preprint whose global computation was not rerun. This audit confirms the source statement, not its full proof or computation.
6. **Storage:** uniform substitution, full-profile blindness, exact r-connectivity, intact-bit recovery, nonuniform read head and bit-erasure limitation are preserved correctly. Exact v3=M excludes7 and21 within this family without implying a general classification.
7. **Composition counts:** ordinary-isomorphism convention, inherited inventories, A527 generating series, eventual quadratic, unknown W, normalized44→89 interpretation and exact least residue periods9/72 all match the conditional algebraic packet. No missing census constant was silently supplied.
8. **Open targets:** recovering W, changing constructor arithmetic, building an H-preserving clock map, and optimizing occupied seams are appropriately posed as further work. No implication between independently declared models is asserted.

## Applied wording improvements

- Replace the claim that the eight-sheet example "requires changing" trees with: the prescribed cycle bound is inconclusive, while a better tree succeeds.
- Replace "extra vertices12 and317 are necessary in number and sufficient" with: two additions are necessary, and {12,317} suffice. This avoids suggesting that each displayed added vertex is individually compulsory in every minimum-cardinality solution.
- Link the primary arXiv source directly for ease of inspection.

The final NOTE.md SHA256 is `cc3d138d0939f659b3c543776d7c2c0d5f8aa3381e9f8d6d53fd55f4695d306c`.

The sibling relative links are designed for the publication layout `synthesis/` alongside `clocks/`, `squares/`, `runner/`, `spectrum/` and `counts527/`; they are not expected to resolve from the separate temporary working directories. Verify that layout when publishing.
