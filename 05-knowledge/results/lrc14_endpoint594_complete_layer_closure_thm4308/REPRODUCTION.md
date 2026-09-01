# THM-4308 reproduction

Run from the repository root. These commands reconstruct the THM-4307 carrier
from its canonical component ledgers and consume the complete 25-row
endpoint-594 residual layer of THM-4307. The typed consumer separately marks
the 22 rows not already supplied by THM-4306. These operations do not merge
proof objects or construct a physical lonely-runner entry.

```powershell
$code = '04-computation/lrc14_endpoint594_complete_layer_closure_thm4308'
$packet = '05-knowledge/results/lrc14_endpoint594_complete_layer_closure_thm4308'
$old = '05-knowledge/results/lrc14_mixed_rank_depth_recursive_signatures_thm4296'
$t4300 = '05-knowledge/results/lrc14_size_preserving_response_staircase_thm4300'
$t4302 = '05-knowledge/results/lrc14_endpoint596_response_exchange_thm4302'
$t4306 = '05-knowledge/results/lrc14_index265_recursive_ideal_thm4306'
$t4307 = '05-knowledge/results/lrc14_endpoint595_support_threshold_residual_hypergraph_compression_thm4307'
$scratch = '.scratch/lrc14_thm4308_replay'
New-Item -ItemType Directory -Force $scratch, "$scratch/typed" | Out-Null

$common = @(
  "$old/inputs/joint421_masks.txt",
  "$old/inputs/reconstructed_final8951.txt",
  "$old/inputs/additions45.txt",
  "$old/inputs/endpoint638_response_witness9.txt",
  "$old/inputs/current_residual22647.csv",
  "$t4300/results/proof_graph/typed_union1624.csv",
  "$t4300/inputs/repairs76.txt",
  "$t4302/inputs/additions4.txt",
  "$t4302/inputs/delete73.txt",
  "$t4307/inputs/additions10.txt",
  "$t4307/inputs/final_delete5104.txt",
  "$t4307/results/proof_graph/residual_top594.csv"
)
```

## Direct layer replay

```powershell
g++ -std=c++20 -O3 -DNDEBUG -pthread -I. `
  "$code/endpoint594_carrier_audit.cpp" `
  -o "$scratch/endpoint594_carrier_audit_O3.exe"
g++ -std=c++20 -O2 -pthread -I. `
  "$code/endpoint594_carrier_audit.cpp" `
  -o "$scratch/endpoint594_carrier_audit_O2.exe"

& "$scratch/endpoint594_carrier_audit_O3.exe" @common `
  "$scratch/endpoint594_pair_audit_O3.csv" `
  "$scratch/endpoint594_failures_O3.csv" 32 `
  > "$scratch/endpoint594_raw_O3.out"
if ($LASTEXITCODE -ne 0) { throw 'O3 raw replay failed' }

& "$scratch/endpoint594_carrier_audit_O2.exe" @common `
  "$scratch/endpoint594_pair_audit_O2.csv" `
  "$scratch/endpoint594_failures_O2.csv" 32 `
  > "$scratch/endpoint594_raw_O2.out"
if ($LASTEXITCODE -ne 0) { throw 'O2 raw replay failed' }
```

The transcripts, pair CSVs, and empty failure CSVs must agree byte-for-byte.
Each transcript reports 3,925 masks, 25 rows, 357,678,750 body tests, zero
failures, and pair-ledger FNV `996ebeec37e58e98`.

## Exact all-mask singleton boundary

```powershell
g++ -std=c++20 -O3 -DNDEBUG -pthread -I. `
  "$code/endpoint594_all_singleton_quotient.cpp" `
  -o "$scratch/endpoint594_all_singleton_quotient_O3.exe"
g++ -std=c++20 -O2 -pthread -I. `
  "$code/endpoint594_all_singleton_quotient.cpp" `
  -o "$scratch/endpoint594_all_singleton_quotient_O2.exe"

& "$scratch/endpoint594_all_singleton_quotient_O3.exe" @common `
  "$scratch/endpoint594_all_singletons_O3.csv" `
  "$scratch/endpoint594_all_protected_masks_O3.csv" 25 `
  > "$scratch/endpoint594_all_singleton_quotient_O3.out"
if ($LASTEXITCODE -ne 0) { throw 'O3 singleton quotient failed' }

& "$scratch/endpoint594_all_singleton_quotient_O2.exe" @common `
  "$scratch/endpoint594_all_singletons_O2.csv" `
  "$scratch/endpoint594_all_protected_masks_O2.csv" 25 `
  > "$scratch/endpoint594_all_singleton_quotient_O2.out"
if ($LASTEXITCODE -ne 0) { throw 'O2 singleton quotient failed' }
```

The O3/O2 transcripts and both pairs of ledgers must be byte-identical. They
must report 17 singleton obligations, 14 protected masks, and the exact
one-deletion split `3,911/14`.

The narrower zero-joint implementation is a separate control for the
nonjoint branch:

```powershell
g++ -std=c++20 -O3 -DNDEBUG -I. `
  "$code/endpoint594_singleton_quotient.cpp" `
  -o "$scratch/endpoint594_singleton_quotient.exe"
& "$scratch/endpoint594_singleton_quotient.exe" @common `
  "$scratch/endpoint594_singletons.csv" `
  "$scratch/endpoint594_protected_masks.csv" `
  > "$scratch/endpoint594_singleton_quotient.out"
```

It must report 14 nonjoint singleton obligations on 11 masks, with 3,493 safe
and 11 unsafe nonjoint single deletions. It deliberately makes no claim about
joint masks.

Run the hostile control:

```powershell
& "$scratch/endpoint594_carrier_audit_O3.exe" @common `
  "$scratch/hostile_00b0832c_pair.csv" `
  "$scratch/hostile_00b0832c_failures.csv" 32 00b0832c `
  > "$scratch/hostile_00b0832c_raw.out"
if ($LASTEXITCODE -ne 2) { throw 'hostile deletion did not fail exactly' }
```

It must produce exactly two failures, `(96,594,054e5001)` and
`(96,594,0d4c5001)`.

## Typed consequence

```powershell
python -B "$code/typed_endpoint594_consumer.py" `
  --universe "$old/inputs/current_residual22647.csv" `
  --prior-union "$t4306/results/proof_graph/typed_union2014.csv" `
  --prior-residual "$t4306/results/proof_graph/final_residual20633.csv" `
  --top594 "$t4306/results/proof_graph/residual_top594.csv" `
  --carrier-target "$t4307/results/proof_graph/residual_top594.csv" `
  --pair-audit "$scratch/endpoint594_pair_audit_O3.csv" `
  --failures "$scratch/endpoint594_failures_O3.csv" `
  --output-dir "$scratch/typed" `
  > "$scratch/typed_endpoint594_consumer.out"
if ($LASTEXITCODE -ne 0) { throw 'typed consumer failed' }
```

The consumer ends in `VERDICT PASS`, writes the 2,036-row union and
20,611-row residual, and freezes all 16 rows of the endpoint-593 frontier.

## Frozen-byte audit

`SHA256SUMS` uses lower-case SHA-256, two spaces, and paths relative to this
packet. Verify every listed file before relying on the theorem.
