# S203 audit of the inherited 2^32 coverage endpoint

Status: **PASS AFTER CORRECTION.** The source defect is confirmed, repaired in the isolated replay, and the full odd prefix has been recertified using F1+F2 alone.

## Precise inherited claim

THM-4875 (D), in `01-canon/theorems/THM-4875-old-mark-digit-lower-bound-lorentzian-gadget-adjoint-and-one-apex-joins-realise-every-odd-number-to-2-32-except-7-and-21.md`, claims FINITE-EXACT positive existence for every odd count in `[1, 2^32]` except 7 and 21. The last odd target is 4,294,967,295. This is a full prefix, not only its last dyadic block, and not an eventual-spectrum theorem. The relevant adopted audit is commit `baba886cbe54a4eb47c7d8a8ceafdee8e4fd5b98` (2026-11-04 S78), not the earlier unaudited Session193 description.

The independent audit's `aud_cover_32.out` reports exactly the misses `{7,21}` after enumerating 47,964,383,242 F1 pairs and 134,010,742,236 F2 evaluations. Its source and a bank of 2,732 actual rooted triples are committed. No final bitset or complete per-target witness certificate is saved in that directory. Consequently a source/output check cannot substitute for reconstructing the coverage map if the source has undefined behavior.

## Boundary defect and its exact mathematical effect

In inherited `aud_cover.c`, the table `GK` has `2^16` entries and `gfast(v)` reads `GK[v >> 16]`. This is in bounds when `v < 2^32`. The F0 pass nevertheless loops over `0 <= v <= 2^32`. Its last iteration reads `GK[65536]`, one beyond the table. ASan and UBSan independently diagnose a global buffer overrun in the isolated deterministic probe `boundary_probe.c`. The original source is preserved unchanged as `aud_cover_inherited.c`.

The legacy one-apex function in this program is **not** the session201 unit-source function with the same letter. Its definition is

    g(0)=1; g(2v+1)=g(v); g(2v)=g(v)+v+(v mod 2).

For every k>=1 it satisfies g(2^k)=2^k+1. Indeed g(2)=3; for k>=2, induction gives g(2^k)=(2^(k-1)+1)+2^(k-1). Thus g(2^32)=4,294,967,297 exceeds the target cap and this endpoint is mathematically irrelevant to F0 coverage. Deleting that iteration or evaluating it correctly leaves the intended F0 set unchanged.

That mathematical observation alone does not prove the old compiled run was unaffected: an out-of-bounds value could be spuriously marked. The present local optimized build returned the even value 4,294,901,760, which is rejected by the odd-count filter. This is a local observation, not a claim about the archived executable's layout. The separate author implementation `chanpak_one_apex_join_cover_20261104.c` uses recursive `gval` and has no corresponding table-bound error.

## Corrected replay and completed result

`aud_cover_corrected.c` makes the minimal general repair: if `v >> 16 >= 2^16`, fall back to the direct closed formula. It adds deterministic comparisons among fast, direct, and recursive evaluations at `2^32-2`, `2^32-1`, and `2^32`.

The full replay is compiled with `-DNOF0`. It reconstructs the F1 and F2 maps and checks every odd target, not just a previously reported miss list. If its misses are exactly 7 and 21, those two families alone certify the prefix, so no F0 pass or redundancy assumption is needed. If extra misses occur, they will be reported before a separate corrected F0 check.

Commands from the workspace root:

    gcc -O3 -march=native -fopenmp -DNOF0 session203-independent-audit/aud_cover_corrected.c -o /tmp/s203-cover-corrected
    OMP_NUM_THREADS=4 /tmp/s203-cover-corrected 32 session203-independent-audit/aud_bank7_inherited.txt > session203-independent-audit/CORRECTED_NOF0_32.out
    gcc -O3 -march=native session203-independent-audit/aud_bank_inherited.c -o /tmp/s203-bank-replay
    /tmp/s203-bank-replay 7 > session203-independent-audit/BANK_REPLAY.txt 2> session203-independent-audit/BANK_REPLAY.log
    cmp session203-independent-audit/BANK_REPLAY.txt session203-independent-audit/aud_bank7_inherited.txt

The bank replay enumerates all labeled tournaments of orders 1 through 7 rooted at vertex 0. Relabeling makes this complete for rooted isomorphism types. Its exact subset path counts construct the triples used by the coverage pass.

The corrected pass is a replay of the independent S78 implementation, not a newly authored third independent sieve. No shared dirty checkout or original author's proof has been modified. The endpoint correction does not claim new numerical coverage or any all-order result.

## Completed replay result

The corrected full replay exited0 and found exactly7 and21 missing from all odd targets in `[1,2^32]`. It used **F1+F2 alone**; F0 was skipped. Therefore the earlier F0 endpoint defect is no longer part of this finite coverage certificate. The elapsed time was465.1seconds with4 threads. Counts:47,964,383,242 F1 pairs and134,010,742,236 F2 evaluations. All22 per-block F1 miss counts match the archived independent run. `COVERAGE_REPLAY_RESULT.json` records the exact scope and totals.

The independent S78 bank generator was rebuilt and rerun on every labeled tournament through order7. Its output matched all2,732 archived triples byte for byte. `BANK_REPLAY.log` records the per-order class of inputs and distinct triples.

The endpoint controls return exactly:

    g(2^32-2)=2,147,483,649; g(2^32-1)=1; g(2^32)=4,294,967,297.

No claim beyond the finite prefix follows. The largest odd certified target is4,294,967,295.
