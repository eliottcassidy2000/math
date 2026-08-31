# Detached audit: r629 mixed optimum and post-9017 descent

Status: `FINITE-EXACT`.  This is the maintained detached-audit packet.

## Typed statement

- Source: the 28 failed rank-9 bodies at `(100,629)` after the r630 repair.
- Candidate universe: active rank-8 or rank-9 fixed-pool masks at that pair.
- Response: the labelled subset of the 28 bodies disjoint from a candidate.
- Preserved predicate: literal fixed-pair activity and body disjointness.
- Destroyed data: the response quotient forgets behavior at other pairs and
  physical-entry information.
- Sidecars: the complete mixed response atlas, ordered carrier inputs, and the
  THM4287 residual.
- Hostile tests: exhaustive atlas reconstruction under `-O0`; exact dual load
  on every response class; direct literal replay of every selected witness;
  list-level carrier reconstruction; and full raw descent replay.

## Exact optimum

The exhaustive atlas contains 419 distinct response classes. After canonical
LF normalization, its independent `-O0` reconstruction is byte-identical to
the frozen TSV, raw-LF SHA-256
`9d60fcfec24de98843212d2dff15821dbad2e641d3bde43ecf811c9323066182`.

The exact dual uses quarter units with nonzero weights
`0:1, 1:1, 3:1, 4:1, 5:2, 7:1, 8:1, 9:2, 23:2, 24:1, 26:1`.
The total is `14/4 = 7/2`; every one of the 419 class loads is at most
`4/4`, with 42 tight classes.  Therefore any integer cover has size at least
four.

The matching four-mask upper certificate is:

| mask | rank | colex ordinal | response | cover | mass | margin ticks |
|---|---:|---:|---:|---:|---:|---:|
| `002ac4c0` | 8 | 270405 | `0c00000` | 2 | 214760471161330 | 31451714968590 |
| `3882a082` | 9 | 14119948 | `0000c27` | 6 | 215227333274140 | 60864028075620 |
| `0041c325` | 9 | 519935 | `0687198` | 10 | 214424024475422 | 10255573756386 |
| `08c28e40` | 9 | 5363554 | `f178240` | 11 | 215167289205628 | 57081251759364 |

Their union is all 28 obligations, with multiplicities 27 once and one twice.
Mask FNV is `5601ef0ed7c3ecaa`; response/activity FNV is
`e83f7f763d141b63`.  Hence the exact mixed response-cover optimum is four.

## Carrier and raw descent

Independent ordered-list reconstruction gives a 9,017-mask carrier with FNV
`07689a1534ce7327`, comprising 9,011 rank-8 and six rank-9 masks.  It is
distinct and contains the full 421-mask common deck.

An independent `-O0` literal scan exactly matches the frozen endpoint
transcript and failure CSV.  The normalized transcript SHA-256 is
`342f7fae712fe566db0b9cdde57de93edc98351976c84df3d889311961fded28`;
the normalized failure CSV SHA-256 is
`52ab54aca80e5308538597b0685a6134e69ae75108f2dc1a80df47800ea8885b`.

There are 48 raw pair rows: 47 successes and one failure row.  Every layer
from 636 through 629 closes.  The first failed pair is `(100,628)`, with four
bodies:

`05346408 15306408 17581400 27d01008`

Their FNV is `2e7f9ddccc403b9d`; the 48-row scan ledger is
`be2cbab7fb58d91e`.  This is a finite descent boundary, not an LRC(14) proof.

## Reproduction

Compile and run the independent optimum consumer:

```powershell
g++ -std=c++20 -O3 -DNDEBUG 04-computation/lrc14_mixed_rank_depth_recursive_signatures_thm4296/r629_mixed_optimum_detached_audit.cpp -o .scratch/r629_mixed_optimum_detached_audit.exe
g++ -std=c++20 -O0 04-computation/lrc14_mixed_rank_depth_recursive_signatures_thm4296/r629_mixed_optimum_detached_audit.cpp -o .scratch/r629_mixed_optimum_detached_audit_O0.exe
& .\.scratch\r629_mixed_optimum_detached_audit.exe 05-knowledge/results/lrc14_mixed_rank_depth_recursive_signatures_thm4296/inputs/r629_failures28.csv 05-knowledge/results/lrc14_mixed_rank_depth_recursive_signatures_thm4296/results/endpoint/r629_mixed_atlas.tsv
```

Reconstruct the carrier and consume raw rows:

```powershell
python 04-computation/lrc14_mixed_rank_depth_recursive_signatures_thm4296/r629_mixed_carrier_identity_audit.py 05-knowledge/results/lrc14_mixed_rank_depth_recursive_signatures_thm4296/inputs/reconstructed_final8951.txt 05-knowledge/results/lrc14_mixed_rank_depth_recursive_signatures_thm4296/inputs/additions45.txt 05-knowledge/results/lrc14_mixed_rank_depth_recursive_signatures_thm4296/inputs/endpoint638_response_witness9.txt 05-knowledge/results/lrc14_mixed_rank_depth_recursive_signatures_thm4296/inputs/joint421_masks.txt
python 04-computation/lrc14_mixed_rank_depth_recursive_signatures_thm4296/r629_post9017_descent_consumer.py DESCENT_OUT FAILURE_CSV
```

`DESCENT_OUT` and `FAILURE_CSV` are explicit arguments: they are the raw
stdout transcript and failure CSV produced by
`mixed_carrier_descending_detached.cpp` for the 9,017-mask carrier with lower
endpoint 620 and appended masks
`0010e125 002ac4c0 3882a082 0041c325 08c28e40`.  The maintained packet keeps
the concise consumer report; its final raw scan has subsequently advanced to
the 9,019-mask carrier and is stored as `final_mixed9019_descent.out`.

Fresh promotion replay on 2026-08-31: the optimum report matched under both
`-O3` and `-O0`; the carrier and raw-row consumers each matched their
maintained concise report exactly.  No distinct O0 report is stored because
its bytes equal the ordinary optimum report.

Independent artifact SHA-256 values:

- optimum source/output: `5eaa57512bdf47e27a7be5bf567600ec7973bff683e3eca94bdaafa8b046d27d` / `c0339471d1e915a0ff9ec3aad4731feef923b19f812b363f923ef383b1a3b879`;
- carrier source/output: `d6ea9ff1bb6adbc8638e0b0b90155c18ba9c40ce509b3e8a3eda800376e6a0f3` / `c44f7f0382fdff7743f51b78fe659ad7e6a71b68fe28742f8ed37a105064026a`;
- row consumer source/output: `0fcb93bdc6a5792f4a51befb24daadf56203a24f41bcb27909f1adc22a247554` / `48d38cafc662596efb5274ab6d4c5795930a8aec7155e3b002d151a18c68bd07`.
