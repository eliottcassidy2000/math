# Reproduction

Run from this packet directory on an arm64 macOS host with Apple clang 17 and
Python 3.  A full run enumerates 90 copies of the 14,307,150-body universe and
seven complete rank-eight response quotients.

```bash
sha256sum -c SHA256SUMS

packet_root="$PWD"
audit_tmp="$(mktemp -d /private/tmp/thm4282-final8996-repro.XXXXXX)"

clang++ -std=c++20 -O3 -DNDEBUG -pthread \
  -I "$packet_root" \
  "$packet_root/source/detached_literal_continuation_audit.cpp" \
  -o "$audit_tmp/detached_literal_audit"

"$audit_tmp/detached_literal_audit" \
  "$packet_root/inputs/reconstructed_final8951.txt" \
  "$packet_root/inputs/joint421_masks.txt" \
  "$packet_root/inputs/combined_residual_520_367.csv" \
  "$packet_root/inputs/additions45.txt" \
  "$audit_tmp/pair_rows.tsv" \
  > "$audit_tmp/audit.stdout" \
  2> "$audit_tmp/audit.progress"

cmp "$audit_tmp/audit.stdout" \
    "$packet_root/results/detached_literal_audit.stdout"
cmp "$audit_tmp/pair_rows.tsv" \
    "$packet_root/results/pair_rows_663_645.tsv"
test "$(wc -l < "$audit_tmp/audit.progress")" -eq 90
test "$(rg -c '^PROGRESS [0-9]+/90 PAIR [0-9]+,[0-9]+$' \
    "$audit_tmp/audit.progress")" -eq 90

python3 "$packet_root/source/derive_final8996_consequence.py" \
  "$packet_root/inputs/combined_residual_520_367.csv" \
  "$audit_tmp/pair_rows.tsv" \
  "$audit_tmp/post_final8996_residual.csv" \
  > "$audit_tmp/consequence.stdout"

cmp "$audit_tmp/consequence.stdout" \
    "$packet_root/results/set_consequence.stdout"
cmp "$audit_tmp/post_final8996_residual.csv" \
    "$packet_root/results/post_final8996_residual.csv"
```

The base carrier reconstruction transcript is frozen in
`controls/reconstruct_final8951.stdout`.  Its Python source is included, but
the 59 historical replay-band transcripts and five historical prefix outputs
are intentionally not duplicated in this compact packet; they remain in the
canonical repository.  The detached carrier/body audit itself is fully
self-contained here from the frozen reconstructed input onward.

`results/detached_literal_audit.progress` is evidence of the frozen run but is
not byte-compared because worker completion order is nondeterministic.
