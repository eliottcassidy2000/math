# THM-4309 reproduction

Run from the repository root.  The commands below reconstruct the inherited
carrier rather than trusting a serialized carrier file.  They operate only on
the fixed thirty-label pool and do not construct a physical lonely-runner
entry.

```powershell
$code = '04-computation/lrc14_endpoint595_support_threshold_residual_hypergraph_compression_thm4309'
$packet = '05-knowledge/results/lrc14_endpoint595_support_threshold_residual_hypergraph_compression_thm4309'
$old = '05-knowledge/results/lrc14_mixed_rank_depth_recursive_signatures_thm4296'
$t4300 = '05-knowledge/results/lrc14_size_preserving_response_staircase_thm4300'
$t4302 = '05-knowledge/results/lrc14_endpoint596_response_exchange_thm4302'
$t4303 = '05-knowledge/results/lrc14_endpoint595_twentyfive_closure_thm4303'
$scratch = '.scratch/lrc14_thm4309_replay'
New-Item -ItemType Directory -Force $scratch | Out-Null
```

## Complete response quotient and exact minimum ten

Build the primary full-universe scanner and the independent
complement-generated scanner:

```powershell
g++ -std=c++20 -O3 -DNDEBUG -I. "$code/response_atlas.cpp" `
  -o "$scratch/response_atlas.exe"
g++ -std=c++20 -O3 -DNDEBUG -I. `
  "$code/independent_complement_response_audit.cpp" `
  -o "$scratch/independent_complement_response_audit.exe"

& "$scratch/response_atlas.exe" `
  "$old/inputs/joint421_masks.txt" `
  "$old/inputs/reconstructed_final8951.txt" `
  "$old/inputs/additions45.txt" `
  "$old/inputs/endpoint638_response_witness9.txt" `
  "$old/inputs/current_residual22647.csv" `
  "$t4300/results/proof_graph/typed_union1624.csv" `
  "$t4300/inputs/repairs76.txt" `
  "$t4302/inputs/additions4.txt" `
  "$t4302/inputs/delete73.txt" `
  "$t4303/results/endpoint595_failures.csv" `
  "$scratch/response_atlas.csv" "$scratch/inactive.csv" `
  > "$scratch/response_atlas.out"

& "$scratch/independent_complement_response_audit.exe" `
  "$t4303/results/endpoint595_failures.csv" `
  > "$scratch/independent_complement_response_audit.out"

python -B "$code/solve_response_cover.py" `
  "$scratch/response_atlas.csv" > "$scratch/solve_response_cover.out"
python -B "$code/joint_certificate_audit.py" `
  "$scratch/response_atlas.csv" "$packet/inputs/additions10.txt" `
  > "$scratch/joint_certificate_audit.out"
```

The generated atlas and inactive CSVs must byte-match `results/response_atlas.csv`
and `results/inactive.csv`.  The path-bearing `OUTPUT_ATLAS` line is the only
expected textual difference in the primary transcript.  The independent
transcript and both exact-solver transcripts must otherwise match their frozen
counterparts line-for-line.  The joint certificate checks integer dual total
`15030`, response load at most `1667`, and the ten-mask cover.

## Activity support, hostile threshold, and residual repair

Build the support census, direct raw verifier, residual response quotient, and
the structurally independent repaired-threshold verifier:

```powershell
g++ -std=c++20 -O3 -DNDEBUG -I. "$code/activity_support_census.cpp" `
  -o "$scratch/activity_support_census.exe"
g++ -std=c++20 -O3 -DNDEBUG -pthread -I. `
  "$code/final_full_prefix_raw_audit.cpp" `
  -o "$scratch/final_full_prefix_raw_audit_O3.exe"
g++ -std=c++20 -O2 -pthread -I. `
  "$code/final_full_prefix_raw_audit.cpp" `
  -o "$scratch/final_full_prefix_raw_audit_O2.exe"
g++ -std=c++20 -O3 -DNDEBUG -I. "$code/residual_retention_atlas.cpp" `
  -o "$scratch/residual_retention_atlas.exe"
g++ -std=c++20 -O3 -DNDEBUG -pthread -I. `
  "$code/k350_retention_certificate.cpp" `
  -o "$scratch/k350_retention_certificate_O3.exe"
g++ -std=c++20 -O2 -pthread -I. `
  "$code/k350_retention_certificate.cpp" `
  -o "$scratch/k350_retention_certificate_O2.exe"
g++ -std=c++20 -O3 -DNDEBUG -pthread -I. `
  "$code/independent_repaired_threshold_replay.cpp" `
  -o "$scratch/independent_repaired_threshold_replay.exe"
```

Recompute all `9,019*391` activity cells:

```powershell
& "$scratch/activity_support_census.exe" `
  "$old/inputs/joint421_masks.txt" `
  "$old/inputs/reconstructed_final8951.txt" `
  "$old/inputs/additions45.txt" `
  "$old/inputs/endpoint638_response_witness9.txt" `
  "$old/inputs/current_residual22647.csv" `
  "$t4300/results/proof_graph/typed_union1624.csv" `
  "$t4300/inputs/repairs76.txt" `
  "$t4302/inputs/additions4.txt" `
  "$t4302/inputs/delete73.txt" `
  "$scratch/activity_support_census.csv" `
  > "$scratch/activity_support_census.out"
```

The CSV and transcript must match `results/activity_support_census.*`.

The following common argument list reconstructs `C_596`:

```powershell
$common = @(
  "$old/inputs/joint421_masks.txt",
  "$old/inputs/reconstructed_final8951.txt",
  "$old/inputs/additions45.txt",
  "$old/inputs/endpoint638_response_witness9.txt",
  "$old/inputs/current_residual22647.csv",
  "$t4300/results/proof_graph/typed_union1624.csv",
  "$t4300/inputs/repairs76.txt",
  "$t4302/inputs/additions4.txt",
  "$t4302/inputs/delete73.txt"
)
```

Delete the entire support-`<=350` pool.  Exit code two is expected, and the
failure file must regenerate the frozen 84-obligation ledger exactly:

```powershell
& "$scratch/final_full_prefix_raw_audit_O3.exe" @common `
  "$packet/inputs/additions10.txt" `
  "$packet/inputs/threshold350_delete5141.txt" `
  "$scratch/threshold350_pair.csv" `
  "$scratch/threshold350_failures.csv" 32 `
  > "$scratch/threshold350_raw.out"
if ($LASTEXITCODE -ne 2) { throw 'D350 hostile replay did not fail exactly' }
```

Build and solve the residual quotient:

```powershell
& "$scratch/residual_retention_atlas.exe" `
  "$packet/inputs/threshold350_delete5141.txt" `
  "$scratch/threshold350_failures.csv" `
  "$scratch/residual_retention_atlas.csv" `
  > "$scratch/residual_retention_atlas.out"
python -B "$code/solve_residual_retention.py" `
  "$scratch/residual_retention_atlas.csv" `
  > "$scratch/solve_residual_retention.out"
```

The quotient has 53 types on 84 obligations.  The matching cover and packing
values must both be 37.

## Final 3,925-mask raw replay

Replay `D_final=D_350\H_37` at two optimization levels.  Each invocation scans
all `391*binom(30,9)=5,594,095,650` row-body cases.

```powershell
& "$scratch/final_full_prefix_raw_audit_O3.exe" @common `
  "$packet/inputs/additions10.txt" "$packet/inputs/final_delete5104.txt" `
  "$scratch/final_pair_O3.csv" "$scratch/final_failures_O3.csv" 32 `
  > "$scratch/final_raw_O3.out"
if ($LASTEXITCODE -ne 0) { throw 'O3 final replay failed' }

& "$scratch/final_full_prefix_raw_audit_O2.exe" @common `
  "$packet/inputs/additions10.txt" "$packet/inputs/final_delete5104.txt" `
  "$scratch/final_pair_O2.csv" "$scratch/final_failures_O2.csv" 32 `
  > "$scratch/final_raw_O2.out"
if ($LASTEXITCODE -ne 0) { throw 'O2 final replay failed' }
```

Both pair ledgers and both empty failure ledgers must be byte-identical.  The
transcripts report carrier size 3,925, ranks `3858/67`, carrier FNV
`6fbd0bffcf0ed78b`, zero failures, and pair-ledger FNV `bb28f7e567c4a4b0`.

Run the independent threshold/hypergraph certificate at O3 and O2:

```powershell
$certargs = @common + @(
  "$packet/inputs/additions10.txt",
  "$packet/inputs/threshold350_delete5141.txt",
  "$packet/results/controls/threshold350_failures84.csv",
  "$packet/inputs/retain37.txt",
  "$packet/inputs/final_delete5104.txt",
  "$packet/inputs/retention_packing37_indices.txt"
)

& "$scratch/k350_retention_certificate_O3.exe" @certargs `
  "$scratch/k350_atlas_O3.csv" "$scratch/k350_edges_O3.csv" `
  > "$scratch/k350_certificate_O3.out"
& "$scratch/k350_retention_certificate_O2.exe" @certargs `
  "$scratch/k350_atlas_O2.csv" "$scratch/k350_edges_O2.csv" `
  > "$scratch/k350_certificate_O2.out"
```

The two certificate transcripts, atlases, and edge ledgers must be
byte-identical.  Finally, verify the exact set difference, packing, cover,
edge alignment, hostile replay, final replay, and all frozen identities:

```powershell
python -B "$code/verify_k350_retention.py" `
  "$scratch/k350_certificate_O3.out" "$scratch/k350_certificate_O2.out" `
  "$scratch/k350_atlas_O3.csv" "$scratch/k350_atlas_O2.csv" `
  "$scratch/k350_edges_O3.csv" "$scratch/k350_edges_O2.csv" `
  "$scratch/final_raw_O3.out" "$scratch/final_raw_O2.out" `
  "$scratch/final_pair_O3.csv" "$scratch/final_pair_O2.csv" `
  "$scratch/final_failures_O3.csv" "$scratch/final_failures_O2.csv" `
  "$packet/inputs/threshold350_delete5141.txt" `
  "$packet/inputs/retain37.txt" "$packet/inputs/final_delete5104.txt" `
  "$packet/inputs/retention_packing37_indices.txt" `
  "$packet/results/controls/threshold350_failures84.csv" `
  "$scratch/threshold350_raw.out" "$scratch/threshold350_pair.csv" `
  "$scratch/threshold350_failures.csv" `
  > "$scratch/verify_k350_retention.out"
```

Successful output ends with the `FULL391 ... FAILURES=0` identity line and an
exit code of zero.

The second independent implementation also reconstructs `D_350`, the
residual quotient, both optimality certificates, and `D_final`, then enumerates
every body whose witness set can change under that deletion:

```powershell
& "$scratch/independent_repaired_threshold_replay.exe" @common `
  "$packet/inputs/additions10.txt" 350 `
  "$packet/inputs/threshold350_delete5141.txt" `
  "$packet/inputs/retain37.txt" "$packet/inputs/final_delete5104.txt" `
  "$packet/results/controls/threshold350_failures84.csv" `
  "$packet/inputs/retention_packing37_indices.txt" 32 `
  > "$scratch/independent_repaired_threshold_replay.out"
```

It must test 5,588,707,348 relevant body-row cases with zero failures and
match `results/independent/independent_repaired_threshold_replay.out`.

## Typed consequence and stopping control

```powershell
python -B "$code/typed_union_consumer.py" `
  "$old/inputs/current_residual22647.csv" `
  "$t4300/results/proof_graph/typed_union1624.csv" `
  "$t4300/results/proof_graph/final_residual21023.csv" `
  "$packet/results/final_full_prefix_pair_audit.csv" `
  "$scratch/typed_union1661.csv" "$scratch/final_residual20986.csv" `
  "$scratch/residual_top594.csv" > "$scratch/typed_union_consumer.out"
```

The three CSVs and transcript must match `results/proof_graph/`.  As a hostile
stopping control, substitute `inputs/threshold360_delete5548.txt` in the direct
raw command.  It must exit two with 391 failures on 36 rows and pair-ledger FNV
`6b4fa1bf8aab4347`, matching `results/controls/threshold360_*`.

## Frozen-byte audit

The packet manifest uses lower-case SHA-256 followed by two spaces and the
path relative to this directory.  Verify every listed artifact before relying
on the theorem.
