# Reproduction

Run from the repository root in PowerShell.  All data inputs are frozen in
the packet; the carrier source imports the maintained THM-4318 code.

```powershell
$s='.scratch/lrc14_endpoint586_scout_agent_20260901/packet'
$i="$s/inputs/carrier"; $t="$s/inputs/typed"
$base=@(
  "$i/joint421_masks.txt", "$i/reconstructed_final8951.txt",
  "$i/additions45.txt", "$i/endpoint638_response_witness9.txt",
  "$t/universe22647.csv", "$i/typed_union1624.csv", "$i/repairs76.txt",
  "$i/additions4.txt", "$i/delete73.txt", "$i/additions10.txt",
  "$i/final_delete5104.txt", "$i/addition1.txt", "$i/safe_delete1.txt"
)
foreach($spec in @(@('O2','-O2',''),@('O3','-O3','-DNDEBUG'))){
  $tag,$opt,$ndebug=$spec; $flags=@('-std=c++20',$opt,'-pthread','-I.')
  if($ndebug){$flags+=$ndebug}
  & g++ @flags "$s/generic_endpoint_exchanged_carrier_audit.cpp" `
    -o "$s/carrier_$tag.exe"
  & "$s/carrier_$tag.exe" @base "$i/cover43.csv" `
    "$i/delete43_low_activity.txt" "$t/residual_top586.csv" `
    "$i/endpoint590_add9.csv" "$i/endpoint590_delete9.txt" `
    586 12 a1b617faa2e7f63f "$s/endpoint586_pair_$tag.csv" `
    "$s/endpoint586_failures_$tag.csv" 12 > "$s/endpoint586_raw_$tag.out"
  if($LASTEXITCODE -ne 0){throw "carrier $tag failed"}

  & g++ -std=c++20 $opt $ndebug -Wall -Wextra -Wpedantic -Werror `
    "$s/generic_endpoint_direct_literal_primary.cpp" `
    -o "$s/direct_$tag.exe"
  & "$s/direct_$tag.exe" 586 "$s/endpoint586_failures_$tag.csv" `
    "$s/endpoint586_direct_primary_$tag.csv" `
    > "$s/endpoint586_direct_primary_$tag.out"
}

& g++ -std=c++20 -O3 -DNDEBUG -pthread -Wall -Wextra -Wpedantic -Werror `
  "$s/generic_endpoint_direct_literal_rawcell_independent.cpp" `
  -o "$s/rawcell_O3.exe"
& "$s/rawcell_O3.exe" 586 4090 ffb884b2b17e6ef4 `
  "$s/endpoint586_failures_O3.csv" `
  "$s/endpoint586_direct_rawcell_independent_O3.csv" 12 `
  > "$s/endpoint586_direct_rawcell_independent_O3.out"

python -B "$s/generic_derive_typed_successor.py" `
  --universe "$t/universe22647.csv" --prior-union "$t/typed_union2217.csv" `
  --prior-residual "$t/final_residual20430.csv" `
  --top "$t/residual_top586.csv" --output-dir "$s/post586_typed" `
  > "$s/post586_typed.out"
python -B -O "$s/generic_derive_typed_successor.py" `
  --universe "$t/universe22647.csv" --prior-union "$t/typed_union2217.csv" `
  --prior-residual "$t/final_residual20430.csv" `
  --top "$t/residual_top586.csv" --output-dir "$s/post586_typed_opt" `
  > "$s/post586_typed_opt.out"

python -B "$s/verify_endpoint586_literal_closure_packet.py" `
  --packet (Resolve-Path $s).Path
python -B -O "$s/verify_endpoint586_literal_closure_packet.py" `
  --packet (Resolve-Path $s).Path
```
