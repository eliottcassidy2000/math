# Publication and reproducibility gates

This first checkpoint contains only the independently accepted local enumeration
lane. Other Session67 lanes remain ongoing until separately audited.

The root enumerator was rebuilt with `-O0` and `-O3 -DNDEBUG`, C++17,
`-Wall -Wextra -pedantic`. Both independently regenerated all 17 census/table
files through order 8; every byte matched (order-8 line bytes after decompression).
Explicit exceptions, not removable assertions, enforce its invariants.

The independent auditor checks all production LUT destinations and legality
against complete arc relations through order 8; an independent all-permutation
class oracle through 6; complete node and coloured-line tables through 7; and
all order-8 line/node mass and degree accounting. It separately checks the
connected-plus-closed criterion on all 4,140 partitions of the order-4 cube,
and the mod-two Laplacian hostile. The audit report states exact scope.

Root exported-ledger verification passes in normal and optimized Python modes:
444,711 line rows and 7,412 ordinary nodes over orders 1–8. These totals count
stored rows/nodes across orders, not ordinary classes of one order.

`FIRST_REPLAY_RECEIPT.json` records the local replays. `FIRST_DOCS_CHECK.txt`
records repository documentation guards. `FILE_MANIFEST.json` seals the exact
curated file set, excluding only itself. Temporary binaries, caches and the
exploratory path-move probe are excluded. `verify_packet.py` rejects changed,
extra, missing and symbolic-link files and checks frozen inherited-source pins.

Remote verification is performed after publication: exact target commit,
ancestry, owned-path Git blob hashes/modes and checks for that exact SHA.
Local receipts do not imply publication. Empty GitHub status/check lists are
absence of CI evidence, not a CI pass.
