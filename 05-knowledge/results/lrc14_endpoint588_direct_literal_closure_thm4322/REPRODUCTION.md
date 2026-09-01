# Reproduction

Run from the repository root in PowerShell.  The packet freezes all data
inputs; only the maintained THM-4318 source include is read outside it.

```powershell
$s='.scratch/lrc14_endpoint588_scout_agent_20260901/packet'
$i="$s/inputs/carrier"
$t="$s/inputs/typed"
$base=@(
  "$i/joint421_masks.txt", "$i/reconstructed_final8951.txt",
  "$i/additions45.txt", "$i/endpoint638_response_witness9.txt",
  "$t/universe22647.csv", "$i/typed_union1624.csv", "$i/repairs76.txt",
  "$i/additions4.txt", "$i/delete73.txt", "$i/additions10.txt",
  "$i/final_delete5104.txt", "$i/addition1.txt", "$i/safe_delete1.txt"
)

foreach($spec in @(@('O2','-O2',''),@('O3','-O3','-DNDEBUG'))){
  $tag,$opt,$ndebug=$spec
  $flags=@('-std=c++20',$opt,'-pthread','-I.')
  if($ndebug){$flags+=$ndebug}
  & g++ @flags "$s/endpoint588_exchanged_carrier_audit.cpp" `
    -o "$s/endpoint588_exchanged_carrier_audit_$tag.exe"
  & "$s/endpoint588_exchanged_carrier_audit_$tag.exe" @base `
    "$i/cover43.csv" "$i/delete43_low_activity.txt" `
    "$t/residual_top588.csv" "$i/endpoint590_add9.csv" `
    "$i/endpoint590_delete9.txt" "$s/endpoint588_pair_$tag.csv" `
    "$s/endpoint588_failures_$tag.csv" 12 > "$s/endpoint588_raw_$tag.out"
  if($LASTEXITCODE -ne 0){throw "carrier audit $tag failed"}

  & g++ -std=c++20 $opt $ndebug -Wall -Wextra -Wpedantic -Werror `
    "$s/endpoint588_direct_literal_primary.cpp" `
    -o "$s/endpoint588_direct_literal_primary_$tag.exe"
  & "$s/endpoint588_direct_literal_primary_$tag.exe" `
    "$s/endpoint588_failures_$tag.csv" `
    "$s/endpoint588_direct_primary_$tag.csv" `
    > "$s/endpoint588_direct_primary_$tag.out"
  if($LASTEXITCODE -ne 0){throw "literal audit $tag failed"}
}

& g++ -std=c++20 -O3 -DNDEBUG -pthread -Wall -Wextra -Wpedantic -Werror `
  "$s/endpoint588_direct_literal_rawcell_independent.cpp" `
  -o "$s/endpoint588_direct_literal_rawcell_independent_O3.exe"
& "$s/endpoint588_direct_literal_rawcell_independent_O3.exe" `
  "$s/endpoint588_failures_O3.csv" `
  "$s/endpoint588_direct_rawcell_independent_O3.csv" 12 `
  > "$s/endpoint588_direct_rawcell_independent_O3.out"

python -B "$s/derive_post588_typed_successor.py" `
  --universe "$t/universe22647.csv" `
  --prior-union "$t/typed_union2141.csv" `
  --prior-residual "$t/final_residual20506.csv" `
  --top588 "$t/residual_top588.csv" --output-dir "$s/post588_typed" `
  > "$s/post588_typed.out"
python -B -O "$s/derive_post588_typed_successor.py" `
  --universe "$t/universe22647.csv" `
  --prior-union "$t/typed_union2141.csv" `
  --prior-residual "$t/final_residual20506.csv" `
  --top588 "$t/residual_top588.csv" --output-dir "$s/post588_typed_opt" `
  > "$s/post588_typed_opt.out"

python -B "$s/verify_endpoint588_literal_closure_packet.py" `
  --packet (Resolve-Path $s).Path
python -B -O "$s/verify_endpoint588_literal_closure_packet.py" `
  --packet (Resolve-Path $s).Path
```
