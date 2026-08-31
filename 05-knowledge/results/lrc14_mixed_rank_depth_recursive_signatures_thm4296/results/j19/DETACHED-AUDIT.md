# Detached J19 common-deck audit

Status: `FINITE-EXACT` detached computation; no maintained files edited.

## Typed transfer

- Source: exact inactive-signature ideal `J={19}` in the post-THM4287 live
  residual, plus the 421-mask common deck.
- Target: delete deck atom 19 and repair its eight private rank-9 body
  obligations with masks that stay active on every row of `J`.
- Preserved predicate: literal fixed-pool activity on all 36 ideal rows and
  rank-9 body hitting by the rebuilt deck.
- Destroyed data: this quotient forgets physical-entry information and all
  live rows outside `J`; it is not an LRC(14) theorem.
- Required sidecars: the live residual, full signature atlas, deck order, and
  the private-body ledger of deleted atom 19.
- Hostile test: enumerate every rank-8 mask and every response to the eight
  private bodies, reject any one-mask full response, then replay every rank-9
  body against the rebuilt deck.

## Result

The ideal has 36 rows, FNV `5c8af37cf2f002e7`.  Deleted mask `1804aa01`
is strictly inactive on all 36.  There are 2,212,775 masks active on every
row (FNV `8133b80b94da565c`), forming 41 response classes with no `ff`
class.  Thus one-mask common-deck surgery is impossible.  Masks `0000e649`
and `0184a205` have responses `0c` and `f3`, so together cover all eight
obligations and prove the exact minimum is two.

The rebuilt 422-mask deck has FNV `dc50478119bc6c12`.  All 15,192
mask-row activity cells are strict (zero equalities), and exhaustive replay
of 14,307,150 rank-9 bodies performs 406,270,061 deck checks with zero
failures; body ledger FNV is `10c4c3ed46d44bf1`.

Literal weakest-gap diagnostics are:

- `0000e649`: row `(6,522)`, `4215823/7059272220`;
- `0184a205`: row `(179,347)`, `141485665543/104830971024780`.

The earlier support/cocycle weakest-row annotation for the second witness is
representation-scale diagnostic data, not a canonical theorem claim.  It
must not replace the literal result above.  This discrepancy does not affect
activity signs, the no-one-mask result, the two-mask witness, or body cover.

## Reproduction

From the repository root in PowerShell:

```powershell
g++ -std=c++20 -O3 -DNDEBUG -pthread .scratch/incoming_lrc_signal_agent/detached_j19_audit.cpp -o .scratch/incoming_lrc_signal_agent/detached_j19_audit.exe
g++ -std=c++20 -O0 -pthread .scratch/incoming_lrc_signal_agent/detached_j19_audit.cpp -o .scratch/incoming_lrc_signal_agent/detached_j19_audit_O0.exe

$j19Inputs = @(
  '05-knowledge/results/lrc14_repaired_carrier_endpoint637_descent_thm4287/inputs/joint421_masks.txt',
  '05-knowledge/results/lrc14_repaired_carrier_endpoint637_descent_thm4287/inputs/full_signatures_primary.csv',
  '05-knowledge/results/lrc14_repaired_carrier_endpoint637_descent_thm4287/results/proof_graph/final_residual.csv'
)
& .\.scratch\incoming_lrc_signal_agent\detached_j19_audit.exe @j19Inputs 4
& .\.scratch\incoming_lrc_signal_agent\detached_j19_audit_O0.exe @j19Inputs 8
```

Both optimization levels reproduce `detached_j19_audit.out` after CRLF/LF
normalization.  Normalized stdout SHA-256:
`4ea0f1d5ef6f23aabdaf43b744fa39e7dae17942e459b25198391ded32afbcd5`.

Artifact SHA-256 values:

- source: `b8f83eb24cc187aff3aaf11ed63d218168dbdc28d700ae6eb30790021393824f`;
- frozen output: `4ea0f1d5ef6f23aabdaf43b744fa39e7dae17942e459b25198391ded32afbcd5`.
