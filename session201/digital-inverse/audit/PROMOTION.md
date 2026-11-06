# Independent repair, promotion, and relocation readback

**PASS.** The numerical inverse, finite-template complexity theorem, and exact targeted certificates retain their audited scope. The release now enforces the original-index precondition vmin>=2. Internal scalar h=1 remains legal. Universal occupancy remains **OPEN**.

## Exact release changes

The independently regenerated `PROMOTION.diff` matches exactly. It contains only the proof/README status updates, the targeted-certificate status update retaining author-only long-search statistics, and the minimal API repair: one docstring sentence and a `core.need(vmin>=2, ...)` guard.

- Final `PROOF.md` SHA256: `970114f1d0b3fbfa02c4420dd69d5efcc597b858a8637f11ae3127033080d915`.
- Final `verify.py` SHA256: `712870cd1b069db1206b4ccc82ebd57fddbacf55116a22ec477813a0fd92b21d`.
- Exact `PROMOTION.diff` SHA256: `efb4c90bea2dbb897e3ea01b10ccba9b1b1d1d63e0c9ca6bbdebbca0b6577f10`.

No interval, gcd, recursion, leaf, index-map, completeness, complexity, witness, or occupancy claim changed. All original changed files are preserved under `originals/`. The twenty frozen artifacts match their original hashes either in those archives or as unchanged release files. The source freeze list is preserved exactly. Every one of the sixteen previous independent audit artifacts and its original manifest was copied byte-for-byte.

## API repair verification

`check_api_guard.py` was read and independently run normally and under optimization. It requires the exact expected error for original vmin=1 and verifies that target three for row (1,1,0,0) still returns legal original v=2 through the scalar singleton boundary. Its two included receipts match both replays.

The existing independent API comparison was also replayed normally and under optimization against the promoted verifier. All 58,725 complete/window controls and 13,500 first-result controls still pass. The only receipt change is the intended boundary status:

    unguarded vmin=1 returns scalar index one
    -> vmin>=2 API precondition enforced

The original comparison receipts remain unchanged. The promoted receipts are saved separately as `PROMOTED_COMPARISON_NORMAL.json` and `PROMOTED_COMPARISON_OPTIMIZED.json`; their common SHA256 is `c1988e50ab28d35837d0f98694d10eea1799be5c55a0fa0a2a2803fa9677393b`. This explicitly closes the original report's API-domain note rather than rewriting its history.

## Relocated execution

The following were independently replayed normally and with `-O` from the release layout:

- `verify.py`: both outputs exactly match the original validation receipt.
- `check_targets.py`: both outputs exactly match the original prime/profile certificate receipt.
- `check_api_guard.py`: both outputs match the included guard receipts.
- `verify.py --target 2134117 --all`: both outputs exactly match the original exhaustive receipt, SHA256 `4aa27c7cc49672a329245c125cc884278d16c6a6efa6221dc61164fb7ae215bd`.
- From `audit/`, `independent.py`: both outputs exactly match the original independent receipt.
- From `audit/`, `compare_author.py`: both outputs match each other with only the documented repaired boundary status.

The independent complete finite-menu exhaustion and the original exhaustive node statistics are therefore preserved. The three long positive first-witness searches were not repeated and their statistics remain author-only. Exact positive certificates and Lucas chains were independently reverified; no literal large-graph path enumeration is claimed.

The audit manifest is refreshed to include this report and the two promoted comparison receipts. Refresh any aggregate release manifest after these additions. No GitHub write or remote publication is certified by this local approval.
