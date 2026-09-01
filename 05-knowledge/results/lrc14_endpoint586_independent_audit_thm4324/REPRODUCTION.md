# Reproduction

Run from the repository root in PowerShell. The packet freezes the inherited
carrier, protected deck, and prior typed partition, so the endpoint audit has
no source imports.

```powershell
$p='.scratch/lrc14_endpoint586_independent_audit_20260901/packet'

python -B "$p/cleanroom_typed_transition.py" `
  --universe "$p/inputs/universe22647.csv" `
  --prior-union "$p/inputs/typed_union2217.csv" `
  --prior-residual "$p/inputs/final_residual20430.csv" `
  --output "$p/typed_normal" > "$p/typed_normal.out"
python -B -O "$p/cleanroom_typed_transition.py" `
  --universe "$p/inputs/universe22647.csv" `
  --prior-union "$p/inputs/typed_union2217.csv" `
  --prior-residual "$p/inputs/final_residual20430.csv" `
  --output "$p/typed_opt" > "$p/typed_opt.out"

foreach($spec in @(@('O2','-O2',''), @('O3','-O3','-DNDEBUG'))) {
  $tag,$opt,$ndebug=$spec
  $flags=@('-std=c++20',$opt,'-pthread','-Wall','-Wextra','-Wpedantic','-Werror')
  if($ndebug){$flags+=$ndebug}
  & g++ @flags "$p/cleanroom_carrier_audit.cpp" -o "$p/carrier_$tag.exe"
  & "$p/carrier_$tag.exe" "$p/inputs/joint421_masks.txt" `
    "$p/inputs/inherited_carrier3925.txt" `
    "$p/typed_normal/reconstructed_top586.csv" `
    "$p/carrier_pair_$tag.csv" "$p/carrier_failures_$tag.csv" 12 `
    > "$p/carrier_$tag.out"
}

python -B "$p/cleanroom_literal_mass.py" `
  --failures "$p/carrier_failures_O2.csv" `
  --detail "$p/literal_normal.csv" > "$p/literal_normal.out"
python -B -O "$p/cleanroom_literal_mass.py" `
  --failures "$p/carrier_failures_O3.csv" `
  --detail "$p/literal_opt.csv" > "$p/literal_opt.out"

python -B "$p/verify_independent_endpoint586_packet.py" `
  --packet (Resolve-Path $p).Path
python -B -O "$p/verify_independent_endpoint586_packet.py" `
  --packet (Resolve-Path $p).Path
```

The frozen carrier was serialized once from the maintained THM-4318
construction by the already-audited adapter
`04-computation/lrc14_endpoint588_direct_literal_closure_thm4322/independent/export_inherited_carrier.cpp`.
Its output is frozen as `inputs/inherited_carrier3925.txt`; the endpoint-586
program itself reads only that ledger and checks its count, ranks, FNV64, and
protected-deck containment.
