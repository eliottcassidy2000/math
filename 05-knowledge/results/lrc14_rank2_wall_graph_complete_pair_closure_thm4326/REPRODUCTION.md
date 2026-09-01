# Reproduction

Run from the repository root.  A full degree pass over 181,194 pairs takes
several minutes on the reference machine; the exact 107-row branch replay is
much faster than the flat `C(30,9)` control.

```powershell
$p = '.scratch/lrc14_literal_exit_theory_20260901/packet'

python -B "$p/export_thm4231_remainder.py" "$p/rebuilt_remainder.csv"

g++ -std=c++20 -O2 "$p/rank2_degree_screen.cpp" -o "$p/degree_O2.exe"
g++ -std=c++20 -O3 "$p/rank2_degree_screen.cpp" -o "$p/degree_O3.exe"
& "$p/degree_O2.exe" "$p/thm4231_remainder181194.csv" "$p/rebuilt_degree_O2.csv" > "$p/rebuilt_degree_O2.out"
& "$p/degree_O3.exe" "$p/thm4231_remainder181194.csv" "$p/rebuilt_degree_O3.csv" > "$p/rebuilt_degree_O3.out"

g++ -std=c++20 -O2 "$p/rank2_allbody_audit.cpp" -o "$p/allbody_O2.exe"
g++ -std=c++20 -O3 "$p/rank2_allbody_audit.cpp" -o "$p/allbody_O3.exe"
& "$p/allbody_O2.exe" --degree-csv "$p/thm4231_remainder181194_degree_O3.csv" > "$p/rebuilt_bad107_O2.out"
& "$p/allbody_O3.exe" --degree-csv "$p/thm4231_remainder181194_degree_O3.csv" > "$p/rebuilt_bad107_O3.out"

g++ -std=c++20 -O2 "$p/rank2_branch_bound_independent.cpp" -o "$p/branch_O2.exe"
g++ -std=c++20 -O3 "$p/rank2_branch_bound_independent.cpp" -o "$p/branch_O3.exe"
& "$p/branch_O2.exe" "$p/thm4231_remainder181194_degree_O3.csv" > "$p/rebuilt_branch_O2.out"
& "$p/branch_O3.exe" "$p/thm4231_remainder181194_degree_O3.csv" > "$p/rebuilt_branch_O3.out"

& "$p/allbody_O2.exe" '50,212' '50,274' '100,110' > "$p/rebuilt_threshold3_O2.out"
& "$p/allbody_O3.exe" '50,212' '50,274' '100,110' > "$p/rebuilt_threshold3_O3.out"

& "$p/degree_O3.exe" "$p/universe22647.csv" "$p/rebuilt_universe_degree.csv" > "$p/rebuilt_universe_degree.out"

python -B "$p/verify_rank2_arbitrary_pair_packet.py" --packet "$p"
python -B -O "$p/verify_rank2_arbitrary_pair_packet.py" --packet "$p"

Push-Location "$p/independent_cleanroom"
python -B verify_packet.py
python -B -O verify_packet.py
Pop-Location
```

The rebuilt full CSVs must match each other and the frozen full CSV.  The
flat and branch exact ledgers must agree on all 107 `(minimum,body)` pairs.
The verifier additionally reconstructs the proof-graph set arithmetic,
checks the exact source-ledger order and the unique signed-64-bit overflow
row, validates the finite-null-wall convention and normalized three-row
threshold argument, rejects the obsolete signed-64-bit scanner/FNV path, and
closes the manifest.  The nested clean-room verifier independently checks its
event-sweep construction, 107-row optimizer, 110 raw-cell replays, overflow
control, and primary fieldwise crosscheck.
