# Independent audit: persistent physical carrier pairs

**Verdict: PASS.** Signed by the independent Session215 mathematical auditor. Exact checked release files are recorded in CARRIERS.sha256.

## Mathematical review

The new pairs use actual rooted seeds, followed by JO or ZJO. Independent recursive start-path counts, converse counts, subset convolutions and literal context counts reproduce all nine coordinates for all four resulting graphs. At base 4225 the starts are 25 and 11, hence j=8 and 3; at base 4357 the starts are 13 and 14, hence j=4 and 4. Both remain in the required [0,8] window.

The shared states (X,e,V,d) are (4225,0,8910,195) and (4357,0,9836,251). The independent compiled graph counts are respectively 12677,12679 and 13073,13075. The common O/Z projection is closed, so these are genuine all-length binary-language families, not numerical collisions alone. O resets the previously mismatched e without changing X or d and doubles equal V; the recorded pre-reset values correctly explain the mechanism.

The necessity proof for universal O/Z equality uses empty word, Z, OZ and ZZ, extracting successively X, V/e and d. It is valid without an abstract realizability assumption. Solving pure-Z gives exactly 8715·2^k−4295 and 9585·2^k−4977 for k>=1. Common binary words retain the fixed starts and labels. Earlier finite collisions with mismatched d correctly fail to imply persistent synchronization.

The old-scalar-orbit exclusions at 4225 and 4357 use the inherited target-sensitive order bound, not an arbitrary search depth. Since B=X−2j<=X<=4357, an old eight-vertex seed needs at most six moves after removing terminal O. The finite depth-six check is sufficient. The stronger full-state separation holds forever since d is invariant and differs from all four old d values. The proof correctly does not claim that every new ray base avoids the old scalar image.

## Reproduction

- The new selected carrier receipt was replayed normally and under Python -O with byte-identical output.
- check_new_carriers.py independently verifies the four new complete profiles and actual compiled graph counts; its normal and optimized outputs agree. audit_support.py supplies a memoized start-recursion engine separate from the authored iterative endpoint engine.
- check_selected_carriers.py additionally independently checks all four earlier selected finite collision records and their physical contexts; SELECTED_CARRIERS.json records the result.
- The deterministic C++ seed-bank search was compiled and replayed. Its stdout matches NEW_BANK.txt; it reports 12893 roots, 341727 stored keys and 20 selected paired keys. These are stored-representative counts, not isomorphism classes or all recipes.
- Dependency Git blobs and SHA256 were checked at 0f8750677dd20b594a0e032a3b56aa23f58a950c. The bank is used only for discovery, not as an assumed complete class census.

The auxiliary length-nine collision and length-ten ray searches are explicitly scoped to their retained representatives. No all-word completeness, new scalar-spectrum record, F(8,137) completion or all-odd threshold follows. The approved positive result is the two physical persistent pairs and their exact infinite continuation languages.

Reproduce from this audit folder after publication:

    python3 check_new_carriers.py
    python3 -O check_new_carriers.py
    python3 check_selected_carriers.py
    python3 ../carriers/verify_new_carriers.py

The proof's universal statements rest on the exact transfer and stated counting arguments; bounded replay checks are independent controls.
