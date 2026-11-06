# Final exact promotion and packaging verification

**Verdict: PASS. The exact final packet inherits all three independently audited conclusions in their stated scopes.**

The final versions are byte-preserved in `promoted-inputs/`. Original pre-audit files remain in `frozen-inputs/`; the intermediate dyadic terminology revision remains in `terminology-inputs/`. The three substantive audit reports remain unchanged by this final packaging check.

## Exact approved changes

The promotion verifier reconstructs each approved final source directly from its frozen original and requires byte-for-byte equality:

- PROOF.md: audit status promotion and the correct clarification that the U<=63 inventory has 55 profiles, 34 of which also have D<=33.
- DIGIT-ANSATZ.md: audit status promotion only.
- DYADIC-SPAN.md: audit status promotion and the reviewed heading correction to “Three linearly independent homogenized vectors in every tail.”
- check_dyadic_span.py: the receipt scope text changes from “actual affine span3” to “actual homogenized linear span3.” No algorithm changes.

The main and digit implementations and their four receipts are identical to their original audited bytes. Both physical dependencies are also identical. The regenerated dyadic receipts differ only in scope text and verifier hash; every numerical field is unchanged.

README.md was independently reviewed. It accurately separates the three scopes, preserves the C1 membership condition and open universal-occupancy problem, and explicitly distinguishes the general compiler's optimal prefix with a nine-letter target-33 suffix from the separately proved optimal seven-letter family suffix.

## Exact final text hashes

- PROOF.md: `fac1169ad772839fd80be856c67d05d01d01f703404456f24a139b8d738c54e0`.
- DIGIT-ANSATZ.md: `cd7b3fcd2a11e9542e2312d02faade85c508c35f503452fd9c85c780d67d33c7`.
- DYADIC-SPAN.md: `9f37c511c6dcd1cda8fb292826470d9218e55fd9b02bb0e93232c872e1c4d859`.
- README.md: `06c316a94534acd652e5e28255f915b91822659bde05ed417e2d8bf969b562fd`.

All final implementation and receipt hashes are pinned in `PROMOTION-VERIFICATION.json` and enforced by `check_promotion.py`.

## Portable release-tree test

The verifier materialized a fresh temporary packet with the intended published layout, including `audit/` and its preserved inputs. It did not rely on the original scratch directory name or neighboring research files.

All twelve replay runs matched their exact expected receipts:

- Main author verifier, normal and optimized.
- Digit author verifier, normal and optimized.
- Dyadic author verifier, normal and optimized.
- Main independent verifier under `audit/`, normal and optimized.
- Digit independent verifier under `audit/`, normal and optimized.
- Dyadic independent verifier under `audit/`, normal and optimized.

Both README repair examples also succeeded from the clean release tree:

- `--repair-word ZOOZ` retains Z, deletes three letters, and appends the valid nine-letter suffix ZOZJOJJJJ at target 33.
- `--repair-word OZJJ --target 5` retains OZJ and appends ZJJJ at target 5.

Reproduce from the published packet root with:

```bash
python3 audit/check_promotion.py
```

The resulting JSON should match `audit/PROMOTION-VERIFICATION.json` byte-for-byte. This check includes its own exact final source checks and the twelve portable replays.

## Final scope

1. Constructive maximal-prefix replacement is exact conditional on actual C1 target membership. It is not a universal occupancy theorem and does not globally minimize suffix length.
2. The digit result excludes the stated coefficient-wise separate-index even-scaling identities.
3. The dyadic result extends that exclusion to actual identities on full products of independent dyadic residue-class tails. Coupled or thinner subsets remain outside the result.

No mathematical or implementation correction remains outstanding. No GitHub writes or Git-index changes were performed by this audit.
