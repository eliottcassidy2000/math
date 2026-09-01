# Reproduction

Run from this packet directory in PowerShell.  All mathematical inputs are
frozen under `inputs/`; no maintained endpoint-587 source is imported.

```powershell
$d=(Get-Location).Path

foreach($spec in @(@('O2','-O2',''),@('O3','-O3','-DNDEBUG'))){
  $tag,$opt,$ndebug=$spec
  $flags=@('-std=c++20',$opt,'-pthread','-Wall','-Wextra','-Wpedantic','-Werror')
  if($ndebug){$flags += $ndebug}
  & g++ @flags "$d/endpoint587_cleanroom_carrier_audit.cpp" `
    -o "$d/endpoint587_cleanroom_carrier_$tag.exe"
  & "$d/endpoint587_cleanroom_carrier_$tag.exe" `
    "$d/inputs/joint421_masks.txt" `
    "$d/inputs/inherited_carrier3925.txt" `
    "$d/typed_normal/residual_top587.csv" `
    "$d/carrier_pair_$tag.csv" "$d/carrier_failures_$tag.csv" 10 `
    > "$d/carrier_$tag.out"
  if($LASTEXITCODE -ne 0){throw "carrier audit $tag failed"}
}

python -B "$d/endpoint587_cleanroom_wall_audit.py" `
  --failures "$d/carrier_failures_O2.csv" `
  --detail "$d/wall_detail_normal.csv" > "$d/wall_normal.out"
python -B -O "$d/endpoint587_cleanroom_wall_audit.py" `
  --failures "$d/carrier_failures_O3.csv" `
  --detail "$d/wall_detail_opt.csv" > "$d/wall_opt.out"

python -B "$d/endpoint587_cleanroom_typed_successor.py" `
  --universe "$d/inputs/universe22647.csv" `
  --prior-union "$d/inputs/typed_union2207.csv" `
  --prior-residual "$d/inputs/final_residual20440.csv" `
  --output-dir "$d/typed_normal" > "$d/typed_normal.out"
python -B -O "$d/endpoint587_cleanroom_typed_successor.py" `
  --universe "$d/inputs/universe22647.csv" `
  --prior-union "$d/inputs/typed_union2207.csv" `
  --prior-residual "$d/inputs/final_residual20440.csv" `
  --output-dir "$d/typed_opt" > "$d/typed_opt.out"

python -B "$d/verify_endpoint587_cleanroom_packet.py" --packet "$d"
python -B -O "$d/verify_endpoint587_cleanroom_packet.py" --packet "$d"
```

The verifier must print `PASS` in both Python modes.  Compiled executables are
reproducible build products and are intentionally excluded from `SHA256SUMS`.

