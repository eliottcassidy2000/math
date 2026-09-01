# THM-4314 reproduction

Run from the repository root with THM-4313 present. Generated files remain in
scratch. The explicit LF normalization mirrors the repository's exact-artifact
checkout policy on Windows.

```powershell
$code='04-computation/lrc14_endpoint591_complete_layer_closure_deletion_boundary_thm4314'
$packet='05-knowledge/results/lrc14_endpoint591_complete_layer_closure_deletion_boundary_thm4314'
$scratch='.scratch/lrc14_thm4314_replay'
$old='05-knowledge/results/lrc14_mixed_rank_depth_recursive_signatures_thm4296'
$t4300='05-knowledge/results/lrc14_size_preserving_response_staircase_thm4300'
$t4302='05-knowledge/results/lrc14_endpoint596_response_exchange_thm4302'
$t4309='05-knowledge/results/lrc14_endpoint595_support_threshold_residual_hypergraph_compression_thm4309'
$t4311='05-knowledge/results/lrc14_endpoint593_response_exchange_thm4311'
$t4313='05-knowledge/results/lrc14_endpoint592_fortythree_response_exchange_thm4313'
New-Item -ItemType Directory -Force $scratch,"$scratch/typed","$scratch/typed_opt","$scratch/derived" | Out-Null

$common=@(
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
$tail=@(
  "$t4311/addition1.txt",
  "$t4311/safe_delete1_O3.txt",
  "$t4313/typed/residual_top591.csv",
  "$t4313/cover43.csv",
  "$t4313/delete43_low_activity.txt"
)

foreach($spec in @(@('O2','-O2',''),@('O3','-O3','-DNDEBUG'))){
  $tag,$opt,$ndebug=$spec
  $flags=@('-std=c++20',$opt,'-pthread','-I.')
  if($ndebug){$flags+=$ndebug}

  & g++ @flags "$code/endpoint591_final_carrier_audit.cpp" -o "$scratch/endpoint591_final_carrier_audit_$tag.exe"
  & "$scratch/endpoint591_final_carrier_audit_$tag.exe" @common @tail "$scratch/endpoint591_pair_$tag.csv" "$scratch/endpoint591_failures_$tag.csv" 16 > "$scratch/endpoint591_raw_$tag.out"
  if($LASTEXITCODE -ne 0){throw "capacity audit $tag failed"}

  & g++ @flags "$code/endpoint591_all_witness_census.cpp" -o "$scratch/endpoint591_all_witness_census_$tag.exe"
  & "$scratch/endpoint591_all_witness_census_$tag.exe" @common @tail "$scratch/endpoint591_all_witness_low_$tag.csv" 16 > "$scratch/endpoint591_all_witness_$tag.out"
  if($LASTEXITCODE -ne 0){throw "all-witness audit $tag failed"}

  & g++ @flags "$code/endpoint591_twohit_census.cpp" -o "$scratch/endpoint591_twohit_census_$tag.exe"
  & "$scratch/endpoint591_twohit_census_$tag.exe" @common @tail "$scratch/endpoint591_twohit_$tag.csv" "$scratch/endpoint591_twohit_pairs_$tag.csv" 16 > "$scratch/endpoint591_twohit_$tag.out"
  if($LASTEXITCODE -ne 0){throw "protected-joint audit $tag failed"}
}

# Canonical packet bytes are LF-normalized.
Get-ChildItem $scratch -File | Where-Object {$_.Extension -in '.out','.csv'} | ForEach-Object {
  $text=[System.IO.File]::ReadAllText($_.FullName).Replace("`r`n","`n").Replace("`r","`n")
  [System.IO.File]::WriteAllText($_.FullName,$text,[System.Text.UTF8Encoding]::new($false))
}

python -B "$code/derive_endpoint591_deletion_boundary.py" --input-dir $scratch --output-dir "$scratch/derived" --joint "$old/inputs/joint421_masks.txt"
if($LASTEXITCODE -ne 0){throw 'deletion-boundary derivation failed'}

$typedCommon=@(
  '--universe',"$old/inputs/current_residual22647.csv",
  '--prior-union',"$t4313/typed/typed_union2087.csv",
  '--prior-residual',"$t4313/typed/final_residual20560.csv",
  '--top591',"$t4313/typed/residual_top591.csv",
  '--pair-audit',"$scratch/endpoint591_pair_O3.csv",
  '--failures',"$scratch/endpoint591_failures_O3.csv"
)
python "$code/typed_endpoint591_consumer.py" @typedCommon --output-dir "$scratch/typed" > "$scratch/typed_endpoint591_consumer.out"
python -O "$code/typed_endpoint591_consumer.py" @typedCommon --output-dir "$scratch/typed_opt" > "$scratch/typed_endpoint591_consumer_opt.out"
if($LASTEXITCODE -ne 0){throw 'typed successor failed'}

python "$code/verify_endpoint591_packet.py" --repo . --packet $packet > "$scratch/verify_endpoint591_packet.out"
python -O "$code/verify_endpoint591_packet.py" --repo . --packet $packet > "$scratch/verify_endpoint591_packet_opt.out"
if($LASTEXITCODE -ne 0){throw 'frozen packet verification failed'}

Get-ChildItem $scratch -File -Recurse | Where-Object {$_.Extension -in '.out','.csv'} | ForEach-Object {
  $text=[System.IO.File]::ReadAllText($_.FullName).Replace("`r`n","`n").Replace("`r","`n")
  [System.IO.File]::WriteAllText($_.FullName,$text,[System.Text.UTF8Encoding]::new($false))
}

$comparisons=@()
foreach($name in @(
  'endpoint591_raw_O2.out','endpoint591_raw_O3.out',
  'endpoint591_pair_O2.csv','endpoint591_pair_O3.csv',
  'endpoint591_failures_O2.csv','endpoint591_failures_O3.csv',
  'endpoint591_all_witness_O2.out','endpoint591_all_witness_O3.out',
  'endpoint591_all_witness_low_O2.csv','endpoint591_all_witness_low_O3.csv',
  'endpoint591_twohit_O2.out','endpoint591_twohit_O3.out',
  'endpoint591_twohit_O2.csv','endpoint591_twohit_O3.csv',
  'endpoint591_twohit_pairs_O2.csv','endpoint591_twohit_pairs_O3.csv',
  'typed_endpoint591_consumer.out','typed_endpoint591_consumer_opt.out',
  'verify_endpoint591_packet.out','verify_endpoint591_packet_opt.out'
)){$comparisons+=,@("$scratch/$name","$packet/$name")}
foreach($name in @(
  'endpoint591_critical_singletons.csv',
  'endpoint591_distinct_twohit_witness_pairs.csv',
  'endpoint591_deletion_boundary.out','endpoint591_O2_O3_identity.out'
)){$comparisons+=,@("$scratch/derived/$name","$packet/$name")}
foreach($dir in @('typed','typed_opt')){
  foreach($name in @('typed_union2100.csv','final_residual20547.csv','residual_top590.csv')){
    $comparisons+=,@("$scratch/$dir/$name","$packet/$dir/$name")
  }
}
foreach($pair in $comparisons){
  if((Get-FileHash -Algorithm SHA256 $pair[0]).Hash -cne
     (Get-FileHash -Algorithm SHA256 $pair[1]).Hash){
    throw "frozen artifact mismatch: $($pair[0])"
  }
}
```

The packet verifier contains no Python `assert` statements. Both verifier
paths must report:

```text
ENDPOINT591_PACKET_VERIFY PASS
target_rows=13 body_tests=185992950 failures=0
singleton_unsafe=2 pair_unsafe=7869 pair_safe=7692981
typed_union=2100 residual=20547 next_endpoint=590 next_rows=13
```

The final loop compares every regenerated O2/O3, derived, typed, and verifier
artifact byte-for-byte with its namesake in the frozen packet.
