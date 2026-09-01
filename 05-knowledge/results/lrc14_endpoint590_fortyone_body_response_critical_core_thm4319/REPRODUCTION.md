# Reproducing the endpoint-590 core41 packet

Run from the packet directory. The packet is standalone and contains the
frozen endpoint-590 failure ledger, complete response-signature atlas, and
integer dual used by the exact search.

## Exact search, two compiler variants

```powershell
$packet = (Get-Location).Path
$tmpBase = Join-Path $env:TEMP "endpoint590-core41-repro"
New-Item -ItemType Directory -Force -Path $tmpBase | Out-Null

g++ -std=c++20 -O2 -DNDEBUG endpoint590_obstruction41_no8.cpp `
  -o (Join-Path $tmpBase "core41-O2.exe")
g++ -std=c++20 -O3 -DNDEBUG endpoint590_obstruction41_no8.cpp `
  -o (Join-Path $tmpBase "core41-O3.exe")

$o2 = (& (Join-Path $tmpBase "core41-O2.exe") `
  endpoint590_response_signatures.csv endpoint590_cover_dual_weights.csv) -join "`n"
$o3 = (& (Join-Path $tmpBase "core41-O3.exe") `
  endpoint590_response_signatures.csv endpoint590_cover_dual_weights.csv) -join "`n"
if ($o2 -cne $o3) { throw "O2/O3 transcript mismatch" }
if (($o2 + "`n") -cne ((Get-Content -Raw endpoint590_obstruction41_no8_O2.out) `
  -replace "`r`n", "`n")) { throw "frozen exact-search transcript mismatch" }
$o2
```

The transcript must report no depth-eight cover of the 41-point core and no
depth-seven cover after any of its 41 single deletions.

## Hardened packet verifier

```powershell
python verify_endpoint590_core41_packet.py .
python -O verify_endpoint590_core41_packet.py .
```

The two outputs must be byte-identical and end in `VERDICT PASS`. The verifier
recomputes the response restriction, inclusion quotient, all 41 direct
eight-cover witnesses, the direct nine-cover witness, the dual capacity, both
pairwise-incompatibility clique censuses, semantic input identities, transcript
identity, the explicit proper four-coloring of the core incompatibility graph,
and manifest closure. It contains no Python `assert`, so `-O` does not weaken
its checks.

## Regenerating deterministic ledgers

The core and restricted-trace ledgers can be rebuilt into a temporary
directory without changing the packet:

```powershell
$regen = Join-Path $env:TEMP "endpoint590-core41-ledgers"
New-Item -ItemType Directory -Force -Path $regen | Out-Null
python build_core41_ledgers.py `
  endpoint590_failures.csv endpoint590_response_signatures.csv `
  (Join-Path $regen "endpoint590_core41.csv") `
  (Join-Path $regen "endpoint590_core41_restricted_traces.csv") `
  (Join-Path $regen "endpoint590_core41_structure.out")

foreach ($name in @("endpoint590_core41.csv", `
                    "endpoint590_core41_restricted_traces.csv", `
                    "endpoint590_core41_structure.out")) {
  $left = (Get-FileHash -Algorithm SHA256 (Join-Path $regen $name)).Hash
  $right = (Get-FileHash -Algorithm SHA256 $name).Hash
  if ($left -cne $right) { throw "regeneration mismatch: $name" }
}
```

The exploratory MILP used to discover the criticality witnesses is deliberately
not part of this proof packet. The hardened verifier checks every frozen
witness directly against the complete response atlas.
