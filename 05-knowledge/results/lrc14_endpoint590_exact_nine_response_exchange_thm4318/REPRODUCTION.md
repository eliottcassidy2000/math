# Reproducing the endpoint-590 nine-for-nine packet

Run from the repository root. Canonical sources are in the computation
directory and frozen outputs are in the result packet. The commands below
reproduce the expensive audits in a separate output directory so the frozen
packet is not modified.

```powershell
$c='04-computation/lrc14_endpoint590_exact_nine_response_exchange_thm4318'
$p='05-knowledge/results/lrc14_endpoint590_exact_nine_response_exchange_thm4318'
$r="$p/response"
$x="$p/exchange"
$o='.scratch/lrc14_endpoint590_thm4318_reproduced'
New-Item -ItemType Directory -Force -Path $o | Out-Null

$old='05-knowledge/results/lrc14_mixed_rank_depth_recursive_signatures_thm4296'
$t4300='05-knowledge/results/lrc14_size_preserving_response_staircase_thm4300'
$t4302='05-knowledge/results/lrc14_endpoint596_response_exchange_thm4302'
$t4309='05-knowledge/results/lrc14_endpoint595_support_threshold_residual_hypergraph_compression_thm4309'
$t4311='05-knowledge/results/lrc14_endpoint593_response_exchange_thm4311'
$t4313='05-knowledge/results/lrc14_endpoint592_fortythree_response_exchange_thm4313'
$t4314='05-knowledge/results/lrc14_endpoint591_complete_layer_closure_deletion_boundary_thm4314'

$base=@(
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
  "$t4309/inputs/final_delete5104.txt",
  "$t4311/addition1.txt",
  "$t4311/safe_delete1_O3.txt"
)

foreach($spec in @(@('O3','-O3','-DNDEBUG'),@('O2','-O2',''))){
  $tag,$opt,$ndebug=$spec
  $flags=@('-std=c++20',$opt,'-pthread','-I.')
  if($ndebug){$flags+=$ndebug}

  & g++ @flags "$c/endpoint590_final_carrier_audit.cpp" -o "$o/baseline_$tag.exe"
  & "$o/baseline_$tag.exe" @base "$t4314/typed/residual_top590.csv" `
    "$t4313/cover43.csv" "$t4313/delete43_low_activity.txt" `
    "$o/endpoint590_pair_$tag.csv" "$o/endpoint590_failures_$tag.csv" 16 `
    > "$o/endpoint590_raw_$tag.out"
  if($LASTEXITCODE -ne 0){throw "baseline $tag failed"}

  & g++ @flags "$c/endpoint590_response_structure.cpp" -o "$o/response_$tag.exe"
  & "$o/response_$tag.exe" @base "$t4313/cover43.csv" `
    "$t4313/delete43_low_activity.txt" "$o/endpoint590_failures_$tag.csv" `
    "$o/endpoint590_response_signatures_$tag.csv" "$o/greedy_$tag.csv" `
    > "$o/endpoint590_response_structure_$tag.out"
  if($LASTEXITCODE -ne 0){throw "response census $tag failed"}

  & g++ @flags "$c/endpoint590_cover_direct_audit.cpp" -o "$o/direct_$tag.exe"
  & "$o/direct_$tag.exe" @base "$t4313/cover43.csv" `
    "$t4313/delete43_low_activity.txt" "$o/endpoint590_failures_$tag.csv" `
    "$r/endpoint590_solver_cover.csv" "$o/endpoint590_cover_incidences_$tag.csv" `
    > "$o/endpoint590_cover_direct_$tag.out"
  if($LASTEXITCODE -ne 0){throw "direct cover audit $tag failed"}

  & g++ @flags "$c/endpoint590_low_witness_census.cpp" -o "$o/low_$tag.exe"
  & "$o/low_$tag.exe" @base "$t4314/typed/residual_top590.csv" `
    "$t4313/cover43.csv" "$t4313/delete43_low_activity.txt" `
    "$o/endpoint590_all_low_$tag.csv" "$o/endpoint590_protected_low_$tag.csv" 16 `
    > "$o/endpoint590_low_witness_$tag.out"
  if($LASTEXITCODE -ne 0){throw "low-witness census $tag failed"}

  & g++ @flags "$c/endpoint590_protected_deletion_quotient.cpp" -o "$o/quotient_$tag.exe"
  & "$o/quotient_$tag.exe" @base "$t4313/cover43.csv" `
    "$t4313/delete43_low_activity.txt" "$t4313/exchange43_final_pair_$tag.csv" `
    "$t4313/typed/residual_top591.csv" `
    "$t4314/endpoint591_all_witness_low_$tag.csv" `
    "$t4314/typed/residual_top590.csv" "$o/endpoint590_all_low_$tag.csv" `
    "$r/endpoint590_solver_cover.csv" `
    "$o/endpoint590_singleton493_$tag.csv" "$o/endpoint590_protected493_$tag.csv" `
    "$o/endpoint590_safe_old493_$tag.csv" "$o/endpoint590_delete9_quotient_$tag.txt" 16 `
    > "$o/endpoint590_deletion_quotient_$tag.out"
  if($LASTEXITCODE -ne 0){throw "deletion quotient $tag failed"}

  & g++ @flags "$c/endpoint590_exchange493_audit.cpp" -o "$o/exchange_$tag.exe"
  & "$o/exchange_$tag.exe" @base "$t4313/cover43.csv" `
    "$t4313/delete43_low_activity.txt" "$t4313/exchange43_final_pair_$tag.csv" `
    "$t4313/typed/residual_top591.csv" "$t4314/typed/residual_top590.csv" `
    "$r/endpoint590_solver_cover.csv" `
    "$x/endpoint590_delete9_capacity_cover_candidate_O3.txt" `
    "$o/endpoint590_exchange493_pair_$tag.csv" `
    "$o/endpoint590_exchange493_failures_$tag.csv" 16 `
    > "$o/endpoint590_exchange493_$tag.out"
  if($LASTEXITCODE -ne 0){throw "simultaneous exchange $tag failed"}
}

foreach($spec in @(@('O3','-O3','-DNDEBUG'),@('O2','-O2',''))){
  $tag,$opt,$ndebug=$spec
  $flags=@('-std=c++20',$opt)
  if($ndebug){$flags+=$ndebug}
  & g++ @flags "$c/endpoint590_exact_no8_search.cpp" -o "$o/no8_$tag.exe"
  & "$o/no8_$tag.exe" "$o/endpoint590_response_signatures_$tag.csv" `
    "$r/endpoint590_cover_dual_weights.csv" 8 > "$o/endpoint590_exact_no8_$tag.out"
  if($LASTEXITCODE -ne 0){throw "exact no-eight search $tag failed"}
}

python -B "$c/typed_endpoint590_consumer.py" `
  --universe "$old/inputs/current_residual22647.csv" `
  --prior-union "$t4314/typed/typed_union2100.csv" `
  --prior-residual "$t4314/typed/final_residual20547.csv" `
  --top590 "$t4314/typed/residual_top590.csv" `
  --pair-audit "$o/endpoint590_exchange493_pair_O3.csv" `
  --failures "$o/endpoint590_exchange493_failures_O3.csv" `
  --output-dir "$o/typed" `
  > "$o/typed_endpoint590_consumer.out"
python -B -O "$c/typed_endpoint590_consumer.py" `
  --universe "$old/inputs/current_residual22647.csv" `
  --prior-union "$t4314/typed/typed_union2100.csv" `
  --prior-residual "$t4314/typed/final_residual20547.csv" `
  --top590 "$t4314/typed/residual_top590.csv" `
  --pair-audit "$o/endpoint590_exchange493_pair_O3.csv" `
  --failures "$o/endpoint590_exchange493_failures_O3.csv" `
  --output-dir "$o/typed_opt" `
  > "$o/typed_endpoint590_consumer_opt.out"
if($LASTEXITCODE -ne 0){throw 'typed consumer failed'}
```

The full simultaneous replay evaluates 7,053,424,950 row-body cases and can
take substantial CPU time. To verify the frozen packet without rerunning those
scans:

```powershell
$repo=(Get-Location).Path
$packet=(Resolve-Path '05-knowledge/results/lrc14_endpoint590_exact_nine_response_exchange_thm4318').Path
$verifier='04-computation/lrc14_endpoint590_exact_nine_response_exchange_thm4318/verify_endpoint590_packet.py'
python -B $verifier --repo $repo --packet $packet
python -B -O $verifier --repo $repo --packet $packet
```

Both verifier modes contain no Python `assert` dependency and must print the
same seven-line `PASS` transcript. `SHA256SUMS` is generated only after all text
files are normalized to LF with no UTF-8 BOM.
