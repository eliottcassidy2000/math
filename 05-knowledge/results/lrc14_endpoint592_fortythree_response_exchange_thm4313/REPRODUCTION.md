# THM-4313 reproduction

Run from the repository root after THM-4311 is present. All generated files
remain below the scratch directory.

~~~powershell
$code = '04-computation/lrc14_endpoint592_fortythree_response_exchange_thm4313'
$packet = '05-knowledge/results/lrc14_endpoint592_fortythree_response_exchange_thm4313'
$scratch = '.scratch/lrc14_thm4313_replay'
$old = '05-knowledge/results/lrc14_mixed_rank_depth_recursive_signatures_thm4296'
$t4300 = '05-knowledge/results/lrc14_size_preserving_response_staircase_thm4300'
$t4302 = '05-knowledge/results/lrc14_endpoint596_response_exchange_thm4302'
$t4309 = '05-knowledge/results/lrc14_endpoint595_support_threshold_residual_hypergraph_compression_thm4309'
$t4310 = '05-knowledge/results/lrc14_endpoint594_complete_layer_closure_thm4310'
$t4311 = '05-knowledge/results/lrc14_endpoint593_response_exchange_thm4311'
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
~~~

## Hostile baseline

Compile endpoint592_exchanged_carrier_scout.cpp at both O3 and O2, then run:

~~~powershell
& "$scratch/endpoint592_exchanged_carrier_scout_TAG.exe" @common "$t4311/addition1.txt" "$t4311/safe_delete1_O3.txt" "$packet/residual_top592.csv" "$scratch/endpoint592_pair_TAG.csv" "$scratch/endpoint592_failures_TAG.csv" 32 > "$scratch/endpoint592_raw_TAG.out"
~~~

Replace TAG by O3 and O2. Each executable must exit 2 because this is a
hostile control. Expected: 500,750,250 cases, 2,468 failures on seven rows,
and obligation FNV 2209b8d6760280cc. The transcripts and ledgers must
byte-match.

## Complete responses, packing, and direct cover

Compile endpoint592_response_multiplicity.cpp, verify_global_packing8.cpp, and
endpoint592_cover_direct_audit.cpp at both O3 and O2. For each TAG run:

~~~powershell
& "$scratch/endpoint592_response_multiplicity_TAG.exe" @common "$t4311/addition1.txt" "$t4311/safe_delete1_O3.txt" "$scratch/endpoint592_failures_TAG.csv" "$scratch/greedy51_TAG.csv" "$scratch/greedy51_responses_TAG.csv" "$scratch/endpoint592_obligation_multiplicity_TAG.csv" > "$scratch/endpoint592_response_multiplicity_TAG.out"
& "$scratch/verify_global_packing8_TAG.exe" "$scratch/endpoint592_failures_TAG.csv" > "$scratch/verify_global_packing8_TAG.out"
& "$scratch/endpoint592_cover_direct_audit_TAG.exe" "$scratch/endpoint592_failures_TAG.csv" "$packet/cover43.csv" > "$scratch/endpoint592_cover_direct_audit_TAG.out"
~~~

Expected complete responder counts are 188,462 at rank eight and 2,205,776
at rank nine. The packing verifier proves lower bound eight. The explicit
43-mask cover has FNV ca3cb80f471f2e7e and misses zero obligations.

## Retained-pool dual and hostile complete pricing

Compile endpoint592_retained_pool_export.cpp at O3, then run:

~~~powershell
& "$scratch/endpoint592_retained_pool_export.exe" @common "$t4311/addition1.txt" "$t4311/safe_delete1_O3.txt" "$scratch/endpoint592_failures_O3.csv" "$scratch/pool_greedy.csv" "$scratch/pool_greedy_responses.csv" "$scratch/candidate_pool_bits.csv" > "$scratch/pool_export.out"
python -B "$code/certify_pool_dual_integer.py" --pool "$scratch/candidate_pool_bits.csv" --dual "$packet/pool_lp_dual.txt" --output "$scratch/pool_dual_integer_certificate.json" --weights-output "$scratch/pool_dual_integer_weights.txt" --masks-output "$scratch/retained_pool_masks.csv"
~~~

Compile full_universe_integer_dual_pricing.cpp at both O3 and O2 and run:

~~~powershell
& "$scratch/full_universe_integer_dual_pricing_TAG.exe" "$scratch/endpoint592_failures_O3.csv" "$scratch/pool_dual_integer_weights.txt" "$scratch/candidate_pool_bits.csv" "$scratch/full_universe_integer_dual_violations_TAG.csv" > "$scratch/full_universe_integer_dual_pricing_TAG.out"
~~~

The retained-pool certificate proves lower bound 36 only inside that
37,497-mask pool. Complete pricing must find 20,986 omitted-mask violations,
FNV 9f5c2aabbea12034. This is the hostile control preventing globalization.
The certificate pass also exports the compact ordered pool-mask ledger; compare
it byte-for-byte with `retained_pool_masks.csv`. Its mask FNV is
`557039eeb8db4ed4`, and the packet verifier checks that all 43 cover masks occur
in this ledger.

## Singleton quotient and complete exchange

Compile endpoint592_full467_singleton_capacity.cpp and
endpoint592_full467_exchange_audit.cpp at O3 and O2 with -pthread. For each
TAG run:

~~~powershell
& "$scratch/endpoint592_full467_singleton_capacity_TAG.exe" @common "$t4309/results/final_full_prefix_pair_audit.csv" "$t4310/results/direct/endpoint594_pair_audit_O3.csv" "$t4311/residual_top593.csv" "$t4311/addition1.txt" "$t4311/safe_delete1_O3.txt" "$packet/residual_top592.csv" "$packet/cover43.csv" "$t4311/singleton432_O3.csv" "$scratch/singleton467_cover43_final_TAG.csv" "$scratch/protected467_cover43_final_TAG.csv" "$scratch/safe_old_stats467_cover43_final_TAG.csv" 32 > "$scratch/singleton_capacity_cover43_final_TAG.out"
& "$scratch/endpoint592_full467_exchange_audit_TAG.exe" @common "$t4309/results/final_full_prefix_pair_audit.csv" "$t4310/results/direct/endpoint594_pair_audit_O3.csv" "$t4311/residual_top593.csv" "$t4311/addition1.txt" "$t4311/safe_delete1_O3.txt" "$packet/residual_top592.csv" "$packet/cover43.csv" "$packet/delete43_low_activity.txt" "$scratch/exchange43_final_pair_TAG.csv" "$scratch/exchange43_final_failures_TAG.csv" 32 > "$scratch/exchange43_final_raw_TAG.out"
~~~

Expected: 1,510 private obligations and 3,135 individually safe old masks.
The separate raw exchange must replay 6,681,439,050 cases with zero failures
and final carrier FNV a0d08a38c10bdab7. All O2/O3 products must byte-match.

## Typed consequence and frozen-packet check

~~~powershell
python -B "$code/typed_endpoint592_consumer.py" --universe "$old/inputs/current_residual22647.csv" --prior-union "$t4311/typed/typed_union2052.csv" --prior-residual "$t4311/typed/final_residual20595.csv" --top592 "$packet/residual_top592.csv" --prefix391-audit "$t4309/results/final_full_prefix_pair_audit.csv" --top594-audit "$t4310/results/direct/endpoint594_pair_audit_O3.csv" --top593 "$t4311/residual_top593.csv" --pair-audit "$scratch/exchange43_final_pair_O3.csv" --failures "$scratch/exchange43_final_failures_O3.csv" --output-dir "$scratch/typed" > "$scratch/typed_endpoint592_consumer.out"
python -B "$code/verify_endpoint592_packet.py" --packet "$packet"
~~~

The packet verifier must report a 2,087-row typed union, a 20,560-row
residual, and next endpoint 591 on 13 rows. Reproductions stay in scratch;
compare them byte-for-byte with the frozen packet before relying on them.
