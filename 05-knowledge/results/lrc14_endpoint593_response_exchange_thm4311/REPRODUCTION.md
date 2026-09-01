# THM-4311 reproduction

Run from the repository root after THM-4309 and THM-4310 are present.  These
commands reconstruct carriers from canonical component ledgers rather than
trusting a serialized final carrier.

```powershell
$code = '04-computation/lrc14_endpoint593_response_exchange_thm4311'
$packet = '05-knowledge/results/lrc14_endpoint593_response_exchange_thm4311'
$scratch = '.scratch/lrc14_thm4311_replay'
$old = '05-knowledge/results/lrc14_mixed_rank_depth_recursive_signatures_thm4296'
$t4300 = '05-knowledge/results/lrc14_size_preserving_response_staircase_thm4300'
$t4302 = '05-knowledge/results/lrc14_endpoint596_response_exchange_thm4302'
$t4309 = '05-knowledge/results/lrc14_endpoint595_support_threshold_residual_hypergraph_compression_thm4309'
$t4310 = '05-knowledge/results/lrc14_endpoint594_complete_layer_closure_thm4310'
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
  "$t4309/inputs/additions10.txt",
  "$t4309/inputs/final_delete5104.txt"
)
```

## Baseline hostile replay

```powershell
g++ -std=c++20 -O3 -DNDEBUG -pthread -I. `
  "$code/endpoint593_carrier_audit.cpp" `
  -o "$scratch/endpoint593_carrier_audit_O3.exe"
g++ -std=c++20 -O2 -pthread -I. `
  "$code/endpoint593_carrier_audit.cpp" `
  -o "$scratch/endpoint593_carrier_audit_O2.exe"

foreach ($tag in @('O3','O2')) {
  & "$scratch/endpoint593_carrier_audit_$tag.exe" @common `
    "$packet/residual_top593.csv" `
    "$scratch/endpoint593_pair_audit_$tag.csv" `
    "$scratch/endpoint593_failures_$tag.csv" 32 `
    > "$scratch/endpoint593_raw_$tag.out"
  if ($LASTEXITCODE -ne 2) { throw "baseline $tag did not hostile-fail" }
}

if ((Compare-Object (Get-Content "$scratch/endpoint593_raw_O3.out") `
                    (Get-Content "$scratch/endpoint593_raw_O2.out")) -or
    (Compare-Object (Get-Content "$scratch/endpoint593_pair_audit_O3.csv") `
                    (Get-Content "$scratch/endpoint593_pair_audit_O2.csv")) -or
    (Compare-Object (Get-Content "$scratch/endpoint593_failures_O3.csv") `
                    (Get-Content "$scratch/endpoint593_failures_O2.csv"))) {
  throw 'baseline O2/O3 mismatch'
}
```

Expected: exactly `96,593,34087401`, with 228,914,400 total cases.

## Complete response census and strict inactivity

```powershell
g++ -std=c++20 -O3 -DNDEBUG -I. `
  "$code/endpoint593_response_capacity.cpp" `
  -o "$scratch/endpoint593_response_capacity.exe"
g++ -std=c++20 -O3 -DNDEBUG -I. `
  "$code/independent_complement_response_audit.cpp" `
  -o "$scratch/independent_complement_response_audit.exe"

& "$scratch/endpoint593_response_capacity.exe" @common `
  "$t4309/results/final_full_prefix_pair_audit.csv" `
  "$t4310/results/direct/endpoint594_pair_audit_O3.csv" `
  "$packet/residual_top593.csv" `
  "$scratch/addition1.txt" "$scratch/common_inactive432.txt" `
  > "$scratch/endpoint593_response_capacity.out"
if ($LASTEXITCODE -ne 0) { throw 'response/capacity audit failed' }

& "$scratch/independent_complement_response_audit.exe" `
  > "$scratch/independent_complement_response_audit.out"
if ($LASTEXITCODE -ne 0) { throw 'independent response audit failed' }
```

Expected: response minimum one; selected addition `0036092c`; zero common
inactive masks on both the prior 416 rows and the full 432 rows.

## Exact one-deletion quotient

```powershell
foreach ($spec in @(@('O3','-O3','-DNDEBUG'), @('O2','-O2',''))) {
  $tag,$opt,$ndebug = $spec
  $flags = @('-std=c++20',$opt,'-pthread','-I.')
  if ($ndebug) { $flags += $ndebug }
  & g++ @flags "$code/endpoint593_singleton_deletion_quotient.cpp" `
    -o "$scratch/endpoint593_singleton_deletion_quotient_$tag.exe"
  if ($LASTEXITCODE -ne 0) { throw "quotient $tag build failed" }

  & "$scratch/endpoint593_singleton_deletion_quotient_$tag.exe" @common `
    "$t4309/results/final_full_prefix_pair_audit.csv" `
    "$t4310/results/direct/endpoint594_pair_audit_O3.csv" `
    "$packet/residual_top593.csv" "$scratch/addition1.txt" `
    "$scratch/singleton432_$tag.csv" `
    "$scratch/protected_nonjoint432_$tag.txt" `
    "$scratch/safe_delete1_$tag.txt" 32 `
    > "$scratch/endpoint593_singleton_deletion_quotient_$tag.out"
  if ($LASTEXITCODE -ne 0) { throw "quotient $tag failed" }
}
```

The O2/O3 transcript and all three ledgers must byte-match.  Expected:
1,520 private obligations, 412 protected and 3,093 safe nonjoint masks, with
selected safe deletion `0006e281`.

## Full 432-row exchange replay

```powershell
foreach ($spec in @(@('O3','-O3','-DNDEBUG'), @('O2','-O2',''))) {
  $tag,$opt,$ndebug = $spec
  $flags = @('-std=c++20',$opt,'-pthread','-I.')
  if ($ndebug) { $flags += $ndebug }
  & g++ @flags "$code/endpoint593_exchange_full432_audit.cpp" `
    -o "$scratch/endpoint593_exchange_full432_audit_$tag.exe"
  if ($LASTEXITCODE -ne 0) { throw "exchange $tag build failed" }

  & "$scratch/endpoint593_exchange_full432_audit_$tag.exe" @common `
    "$t4309/results/final_full_prefix_pair_audit.csv" `
    "$t4310/results/direct/endpoint594_pair_audit_O3.csv" `
    "$packet/residual_top593.csv" "$scratch/addition1.txt" `
    "$scratch/safe_delete1_$tag.txt" `
    "$scratch/exchange_full432_pair_audit_$tag.csv" `
    "$scratch/exchange_full432_failures_$tag.csv" 32 `
    > "$scratch/exchange_full432_raw_$tag.out"
  if ($LASTEXITCODE -ne 0) { throw "exchange $tag failed" }
}
```

The O2/O3 transcript, 432-row pair ledger, and empty failure ledger must
byte-match.  Expected: 6,180,688,800 tests and zero failures.

## Typed consequence and packet verification

```powershell
python -B "$code/typed_endpoint593_consumer.py" `
  --universe "$old/inputs/current_residual22647.csv" `
  --prior-union "$t4310/results/proof_graph/typed_union2036.csv" `
  --prior-residual "$t4310/results/proof_graph/final_residual20611.csv" `
  --top593 "$t4310/results/proof_graph/residual_top593.csv" `
  --prefix391-audit "$t4309/results/final_full_prefix_pair_audit.csv" `
  --top594-audit "$t4310/results/direct/endpoint594_pair_audit_O3.csv" `
  --pair-audit "$scratch/exchange_full432_pair_audit_O3.csv" `
  --failures "$scratch/exchange_full432_failures_O3.csv" `
  --output-dir "$scratch/typed" `
  > "$scratch/typed_endpoint593_consumer.out"
if ($LASTEXITCODE -ne 0) { throw 'typed consumer failed' }

python -B "$code/verify_endpoint593_packet.py" `
  --repo . --scratch "$scratch" `
  > "$scratch/verify_endpoint593_packet.out"
if ($LASTEXITCODE -ne 0) { throw 'packet verifier failed' }
```

The canonical packet records the frozen outputs. Reproductions are written
only below `$scratch`; compare them byte-for-byte with the corresponding
packet files before relying on the result.
