# Fixed symmetry-slack generating functions

Research date: 2026-10-21. Status: PROVED in the stated scopes, with independent proof audits and FINITE-EXACT controls. Publication is established by the enclosing verified repository checkpoint, not by preparation metadata.

The [canonical compiler proof](session23-connections/SYMMETRY_SLACK_RATIONALITY.md) proves rationality at every fixed symmetry slack. Read its [producer overview](session23-connections/README.md), [independent audit](session23-slack-review/AUDIT.md), and [final exact-byte binding](session23-slack-review/final-review-receipt.json). The inherited bounded exceptional-symmetry theorem remains a dependency; no efficient general inventory is claimed.

The separately audited [renewal transfer](session23-root/SLACK_RENEWAL_TRANSFER.md) is bundled here. Its [review](session23-slack-renewal-review/REVIEW.md) and [final binding](session23-slack-renewal-review/final-review-receipt.json) distinguish synthetic compiler controls from actual low-slack tournaments and do not claim a complete general numerator inventory.

## Reproduction

From this packet directory:

    python3 -B session23-connections/verify_symmetry_slack.py 7
    python3 -B -O session23-connections/verify_symmetry_slack.py 7
    python3 -B session23-slack-review/check_actual_tournaments.py
    python3 -B session23-slack-renewal-review/check_slack_transfer.py
    python3 -B -O session23-slack-renewal-review/check_slack_transfer.py

The actual-tournament checker rewrites its receipt with elapsed runtime. Use a disposable copy to preserve archived bytes. Frozen source and audit files retain their original path descriptions and historical receipts. The producer order-seven controls, independent actual-tournament audit, and synthetic transfer checks are separate universes, not interchangeable censuses or formal verification.

No metagraph adjacency classification, unrestricted tournament-enumeration efficiency, or resolution of the repository's major open conjectures follows.
