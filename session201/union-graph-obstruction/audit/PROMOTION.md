# Union witness graph promotion audit

Verdict: **PASS**. The promoted packet preserves the independently audited mathematical statement and proof. Status promotion is supported within the exact scope below.

## Compared release

Packet: `session201/union-graph-obstruction/`.

Promoted `PROOF.md` SHA-256:
`91886431c15894428cadd2b50c6f1ffe11309b4997c3f9d70a19a11a999e8cfe`.

- Read `PROMOTION.diff` and `RELEASE.md`.
- Archived `originals/PROOF.md`, `originals/README.md`, and `originals/MANIFEST.sha256` match the frozen independently audited input byte for byte.
- Current `verify.py`, `INPUT-MENU.json`, and both `VALIDATION-*.json` files match that frozen input byte for byte.
- Current proof differs only by changing the status to PROVED / INDEPENDENTLY AUDITED and replacing the two reproduction commands with packet-local commands.
- Current README differs only by changing “unreviewed” to “independently audited”.
- Independently regenerated unified diffs match `PROMOTION.diff` exactly.
- All 15 files listed in the historical independent audit manifest match their original audit files and listed hashes. The copied audit manifest itself also matched byte for byte before refresh.

## Relocated reproduction

Ran the release's current `verify.py` in its own packet directory under normal Python and `python3 -O`. Both exit successfully, with no stderr, and output matches the corresponding frozen author receipt exactly.

Ran the copied `audit/independent_verify.py` in its own audit directory under normal Python and `python3 -O`. Both exit successfully, with no stderr, and output matches the corresponding frozen independent receipt exactly.

New replay outputs are retained as `PROMOTION-replay-*.json`; exact checks and receipt hashes are in `PROMOTION-receipts.json`. These are relocation checks, not replacements for the historical evidence.

## Preservation and release handoff

Historical audit files and frozen input snapshots remain unchanged. The historical audit manifest is additionally retained byte for byte in `PROMOTION-original-audit-manifest.sha256.txt`. The current audit manifest has been refreshed to include the new promotion artifacts. No mathematical source or verifier code was edited, and no GitHub write was made.

The publisher must regenerate the aggregate release `MANIFEST.sha256` after incorporating this promotion audit, then verify it from the release packet root. At review time the aggregate release manifest had not yet been created; this does not affect the validated input comparisons or replays.

## Exact approved scope

The proof establishes non-binary-automaticity of the index-preserving ternary union witness relation (u,v,n), using strict physical-word increase and selection of the empty-template fiber minimum. It does not establish nonregularity of projected target support, disprove occupancy, or settle cofiniteness. The twelve-template target occupancy assertion remains **OPEN**. The separate digital-graph companion and other counting/approximation packets are outside this promotion audit.
