# Reproduction

Run from the repository root in PowerShell.  Set `$p` to this packet.

```powershell
$p='.scratch/lrc14_endpoint589_q50_cover_theory_agent_20260901/capacity_19_95_independent_audit_packet'

g++ -std=c++20 -O2 -march=native -DNDEBUG -I"$p" `
  "$p/audit_capacity_residual_19_95.cpp" -o "$p/audit_O2.exe"
g++ -std=c++20 -O3 -march=native -DNDEBUG -I"$p" `
  "$p/audit_capacity_residual_19_95.cpp" -o "$p/audit_O3.exe"

& "$p/audit_O2.exe" "$p/input_endpoint589_failures.csv" `
  "$p/input_packing19.csv" "$p/input_cover95.csv" > "$p/replay_O2.out"
if($LASTEXITCODE -ne 0){throw 'O2 audit failed'}

& "$p/audit_O3.exe" "$p/input_endpoint589_failures.csv" `
  "$p/input_packing19.csv" "$p/input_cover95.csv" > "$p/replay_O3.out"
if($LASTEXITCODE -ne 0){throw 'O3 audit failed'}

if(-not ((Get-Content -Raw "$p/replay_O2.out") -ceq `
         (Get-Content -Raw "$p/replay_O3.out"))){
  throw 'O2/O3 replay mismatch'
}

python -B "$p/q50_cover_two_for_one_direct.py" $p `
  > "$p/q50_cover_two_for_one_direct.out"
python -B -O "$p/q50_cover_two_for_one_direct.py" $p `
  > "$p/q50_cover_two_for_one_direct_opt.out"
if($LASTEXITCODE -ne 0){throw 'local-exchange audit failed'}

python -B "$p/verify_packet.py"
python -B -O "$p/verify_packet.py"
if($LASTEXITCODE -ne 0){throw 'packet verifier failed'}
```

The C++ audit takes roughly 25--40 seconds per build on the reference machine.
The Python verifier independently rebuilds the wall geometry and directly
replays the 95-mask upper certificate; the complete 20,160,075-mask lower
stream is the C++ audit's proof-critical exhaustive step.
The local-exchange audit enumerates every possible single replacement after
removing each of the `4,465` cover-mask pairs.

`freeze_manifest.py` regenerates `SHA256SUMS` only when intentionally freezing
a new packet version; do not run it as part of ordinary verification.
