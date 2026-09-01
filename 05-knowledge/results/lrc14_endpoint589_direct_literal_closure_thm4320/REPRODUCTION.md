# Reproduction

Run from the repository root in PowerShell. The commands rebuild the two typed
successors, the exact endpoint-589 carrier audit, two independent direct-wall
audits, the fixed-fifty bridge, and the hardened packet verifier.

```powershell
$p='05-knowledge/results/lrc14_endpoint589_direct_literal_closure_thm4320'
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

python -B "$p/derive_post590_typed_successor.py" `
  --universe "$old/inputs/current_residual22647.csv" `
  --prior-union "$t4314/typed/typed_union2100.csv" `
  --prior-residual "$t4314/typed/final_residual20547.csv" `
  --top590 "$t4314/typed/residual_top590.csv" `
  --output-dir "$p/post590_typed" > "$p/post590_typed.out"
python -B -O "$p/derive_post590_typed_successor.py" `
  --universe "$old/inputs/current_residual22647.csv" `
  --prior-union "$t4314/typed/typed_union2100.csv" `
  --prior-residual "$t4314/typed/final_residual20547.csv" `
  --top590 "$t4314/typed/residual_top590.csv" `
  --output-dir "$p/post590_typed_opt" > "$p/post590_typed_opt.out"

foreach($spec in @(@('O2','-O2',''),@('O3','-O3','-DNDEBUG'))){
  $tag,$opt,$ndebug=$spec
  $flags=@('-std=c++20',$opt,'-pthread','-I.')
  if($ndebug){$flags+=$ndebug}
  & g++ @flags "$p/endpoint589_exchanged_carrier_audit.cpp" `
    -o "$p/endpoint589_exchanged_carrier_audit_$tag.exe"
  & "$p/endpoint589_exchanged_carrier_audit_$tag.exe" @base `
    "$t4313/cover43.csv" "$t4313/delete43_low_activity.txt" `
    "$p/post590_typed/residual_top589.csv" "$p/endpoint590_add9.csv" `
    "$p/endpoint590_delete9.txt" "$p/endpoint589_pair_$tag.csv" `
    "$p/endpoint589_failures_$tag.csv" 16 `
    > "$p/endpoint589_raw_$tag.out"
  if($LASTEXITCODE -ne 0){throw "endpoint589 carrier audit $tag failed"}

  & g++ -std=c++20 $opt $ndebug -Wall -Wextra -Wpedantic -Werror `
    "$p/endpoint589_direct_literal_primary.cpp" `
    -o "$p/direct_primary_$tag.exe"
  & "$p/direct_primary_$tag.exe" "$p/endpoint589_failures_$tag.csv" `
    "$p/direct_primary_$tag.csv" > "$p/direct_primary_$tag.out"
  if($LASTEXITCODE -ne 0){throw "endpoint589 direct audit $tag failed"}

  & g++ -std=c++20 $opt $ndebug -Wall -Wextra -Wpedantic -Werror `
    "$p/endpoint589_literal_lower_bound_independent.cpp" `
    -o "$p/literal_lower_bound_independent_$tag.exe"
  & "$p/literal_lower_bound_independent_$tag.exe" `
    "$p/endpoint589_failures_$tag.csv" `
    "$p/literal_body_masses_independent_$tag.csv" `
    > "$p/literal_lower_bound_independent_$tag.out"
  if($LASTEXITCODE -ne 0){throw "independent literal audit $tag failed"}

  & g++ @flags "$p/q50_active_carrier.cpp" `
    -o "$p/q50_active_carrier_$tag.exe"
  & "$p/q50_active_carrier_$tag.exe" @base `
    "$t4313/cover43.csv" "$t4313/delete43_low_activity.txt" `
    "$p/endpoint590_add9.csv" "$p/endpoint590_delete9.txt" `
    "$p/q50_active_carrier_$tag.csv" > "$p/q50_active_carrier_$tag.out"
  if($LASTEXITCODE -ne 0){throw "q50 active-carrier dump $tag failed"}
}

python -B "$p/endpoint589_direct_literal_independent.py" `
  "$p/endpoint589_failures_O2.csv" "$p/direct_independent.csv" `
  > "$p/direct_independent.out"
python -B -O "$p/endpoint589_direct_literal_independent.py" `
  "$p/endpoint589_failures_O2.csv" "$p/direct_independent_opt.csv" `
  > "$p/direct_independent_opt.out"

python -B "$p/analyze_fixed50_petal_bridge.py" `
  "$p/direct_primary_O2.csv" > "$p/petal_bridge.out"
python -B -O "$p/analyze_fixed50_petal_bridge.py" `
  "$p/direct_primary_O2.csv" > "$p/petal_bridge_opt.out"

python -B "$p/q50_failure_structure.py" `
  "$p/endpoint589_failures_O2.csv" "$p/q50_active_carrier_O2.csv" `
  "$p/q50_vertex_degrees.csv" "$p/q50_hub_fibres.csv" `
  "$p/q50_neither_bodies.csv" > "$p/q50_structure.out"
python -B -O "$p/q50_failure_structure.py" `
  "$p/endpoint589_failures_O2.csv" "$p/q50_active_carrier_O2.csv" `
  "$p/q50_vertex_degrees_opt.csv" "$p/q50_hub_fibres_opt.csv" `
  "$p/q50_neither_bodies_opt.csv" > "$p/q50_structure_opt.out"

python -B "$p/derive_post589_typed_successor.py" `
  --universe "$old/inputs/current_residual22647.csv" `
  --prior-union "$p/post590_typed/typed_union2113.csv" `
  --prior-residual "$p/post590_typed/final_residual20534.csv" `
  --top589 "$p/post590_typed/residual_top589.csv" `
  --output-dir "$p/post589_typed" > "$p/post589_typed.out"
python -B -O "$p/derive_post589_typed_successor.py" `
  --universe "$old/inputs/current_residual22647.csv" `
  --prior-union "$p/post590_typed/typed_union2113.csv" `
  --prior-residual "$p/post590_typed/final_residual20534.csv" `
  --top589 "$p/post590_typed/residual_top589.csv" `
  --output-dir "$p/post589_typed_opt" > "$p/post589_typed_opt.out"

$repo=(Get-Location).Path
$packet=(Resolve-Path -LiteralPath $p).Path
python -B "$p/verify_endpoint589_direct_literal_closure_packet.py" `
  --repo "$repo" --packet "$packet" `
  > "$p/verify_endpoint589_direct_literal_closure_packet.out"
python -B -O "$p/verify_endpoint589_direct_literal_closure_packet.py" `
  --repo "$repo" --packet "$packet" `
  > "$p/verify_endpoint589_direct_literal_closure_packet_opt.out"
if($LASTEXITCODE -ne 0){throw 'packet verifier failed'}
```

The verifier uses explicit runtime gates rather than `assert`, checks manifest
closure, and must report
`ENDPOINT589_DIRECT_LITERAL_CLOSURE_PACKET_VERIFY PASS` identically in normal
and optimized Python modes.
