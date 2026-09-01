# Reproduction

Run from the repository root in PowerShell. These commands rebuild the
inherited carrier serialization, two independent compiler-mode carrier scans,
two wall-mass scans, and the independent typed successor.

```powershell
$d='.scratch/lrc14_endpoint588_independent_agent_20260901/packet'
$old='05-knowledge/results/lrc14_mixed_rank_depth_recursive_signatures_thm4296'
$t4300='05-knowledge/results/lrc14_size_preserving_response_staircase_thm4300'
$t4302='05-knowledge/results/lrc14_endpoint596_response_exchange_thm4302'
$t4309='05-knowledge/results/lrc14_endpoint595_support_threshold_residual_hypergraph_compression_thm4309'
$t4311='05-knowledge/results/lrc14_endpoint593_response_exchange_thm4311'
$t4313='05-knowledge/results/lrc14_endpoint592_fortythree_response_exchange_thm4313'
$p='05-knowledge/results/lrc14_endpoint589_direct_literal_closure_thm4320'

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
  "$t4311/safe_delete1_O3.txt",
  "$t4313/cover43.csv",
  "$t4313/delete43_low_activity.txt",
  "$p/endpoint590_add9.csv",
  "$p/endpoint590_delete9.txt"
)

g++ -std=c++20 -O2 -I. "$d/export_inherited_carrier.cpp" `
  -o "$d/export_inherited_carrier.exe"
& "$d/export_inherited_carrier.exe" @base "$d/inherited_carrier3925.txt" `
  > "$d/export_inherited_carrier.out"

foreach($spec in @(@('O2','-O2',''),@('O3','-O3','-DNDEBUG'))){
  $tag,$opt,$ndebug=$spec
  $flags=@('-std=c++20',$opt,'-pthread','-Wall','-Wextra','-Wpedantic','-Werror')
  if($ndebug){$flags+=$ndebug}
  & g++ @flags "$d/endpoint588_cleanroom_carrier_audit.cpp" `
    -o "$d/endpoint588_cleanroom_carrier_audit_$tag.exe"
  & "$d/endpoint588_cleanroom_carrier_audit_$tag.exe" `
    "$old/inputs/joint421_masks.txt" "$d/inherited_carrier3925.txt" `
    "$p/post589_typed/residual_top588.csv" `
    "$d/endpoint588_cleanroom_pair_$tag.csv" `
    "$d/endpoint588_cleanroom_failures_$tag.csv" 16 `
    > "$d/endpoint588_cleanroom_carrier_$tag.out"

  & g++ -std=c++20 $opt $ndebug -Wall -Wextra -Wpedantic -Werror `
    "$d/endpoint588_cleanroom_literal_audit.cpp" `
    -o "$d/endpoint588_cleanroom_literal_audit_$tag.exe"
  & "$d/endpoint588_cleanroom_literal_audit_$tag.exe" `
    "$d/endpoint588_cleanroom_failures_$tag.csv" `
    "$d/endpoint588_cleanroom_literal_$tag.csv" `
    > "$d/endpoint588_cleanroom_literal_$tag.out"
}

foreach($spec in @(@('normal',''),@('opt','-O'))){
  $tag,$opt=$spec
  python -B $opt "$d/endpoint588_cleanroom_typed_successor.py" `
    --universe "$old/inputs/current_residual22647.csv" `
    --prior-union "$p/post589_typed/typed_union2141.csv" `
    --prior-residual "$p/post589_typed/final_residual20506.csv" `
    --top588 "$p/post589_typed/residual_top588.csv" `
    --output-dir "$d/typed_$tag" `
    > "$d/endpoint588_cleanroom_typed_$tag.out"
}

python -B "$d/verify_independent_endpoint588_packet.py" --packet "$d"
python -B -O "$d/verify_independent_endpoint588_packet.py" --packet "$d"
```

The carrier scan tests 944,271,900 labelled row-body cases. The packet verifier
must print the same `PASS` transcript in normal and optimized Python modes.
