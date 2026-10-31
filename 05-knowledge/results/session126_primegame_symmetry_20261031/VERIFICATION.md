# Session126 publication verification

**Result:** owned routes, byte-integrity checks, and whitespace checks pass.
The repository documentation guard retains one unchanged baseline failure.

## Packet and route checks

All **five** `SHA256SUMS` manifests were rechecked successfully:

| Component | Verified manifest entries |
|---|---:|
| Balance dynamics | 27 |
| Five families | 18 |
| Inheritance | 6 |
| Metagraph CA | 13 |
| PRIMEGAME, including both arithmetic and CA subdirectories | 23 |

These are **87 manifest entries**, not six separate manifests.

The reviewed route set, before adding this receipt, contains **21 Markdown
documents and 38 local links**. Every destination resolves; sparse baseline
targets were checked against Git objects. The recent-session index reaches
the Session126 route, synthesis, and component packets. The owned
`git diff --check` passes. No external-link availability check is claimed.

Exact TSV bytes are retained through the packet's `*.tsv -text` attribute.
Git reports `text: unset` for all three TSV files, and their checksums pass.

## Documentation-guard boundary

The sparse validation checkout received **54 missing baseline inputs**, without
replacing existing files. Its bounded theorem fixture contains **34 files**:
33 required frontier theorems and one unproved router-smoke candidate. This is
not a full theorem-graph audit.

`agents/check_docs.py` reports only the inherited
`05-knowledge/hypotheses/INDEX.md` line-budget failure: **121 maintained lines
versus a limit of 120**. That file was verified unchanged from the baseline;
no repair or mathematical-text edit was applied. All other checks exercised
by this bounded invocation pass.

These are local verification results. **No CI run or CI status is asserted.**
