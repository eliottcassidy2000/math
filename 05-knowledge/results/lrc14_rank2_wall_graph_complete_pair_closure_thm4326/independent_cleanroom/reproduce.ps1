$ErrorActionPreference = "Stop"
$Packet = Split-Path -Parent $MyInvocation.MyCommand.Path
$BuildRoot = Join-Path ([IO.Path]::GetTempPath()) ("lrc14-rank2-independent-" + [guid]::NewGuid().ToString("N"))

function Assert-Exit([string]$Label) {
    if ($LASTEXITCODE -ne 0) { throw "$Label failed with exit code $LASTEXITCODE" }
}

function File-Hash([string]$Path) {
    return (Get-FileHash -LiteralPath $Path -Algorithm SHA256).Hash.ToLowerInvariant()
}

function Assert-Same([string]$Left, [string]$Right, [string]$Label) {
    $leftHash = File-Hash $Left
    $rightHash = File-Hash $Right
    if ($leftHash -ne $rightHash) { throw "$Label mismatch: $leftHash != $rightHash" }
}

function Invoke-Audit(
    [string]$Exe, [string]$Input, [string]$Screen, [string]$Exact, [string]$Summary
) {
    $lines = & $Exe $Input $Screen $Exact
    Assert-Exit "audit $Exe"
    $text = (($lines -join "`r`n") + "`r`n")
    [IO.File]::WriteAllText($Summary, $text, [Text.UTF8Encoding]::new($false))
    if ($text -notmatch "VERDICT PASS") { throw "audit did not pass: $Exe" }
}

New-Item -ItemType Directory -LiteralPath $BuildRoot | Out-Null
try {
    $Source = Join-Path $Packet "src/rank2_event_sweep_cleanroom.cpp"
    $WideSource = Join-Path $Packet "src/rank2_event_sweep_wide_cleanroom.cpp"
    $SelectedSource = Join-Path $Packet "src/rank2_selected_ratio_exact.cpp"

    $NarrowO2 = Join-Path $BuildRoot "narrow_O2.exe"
    $NarrowO3 = Join-Path $BuildRoot "narrow_O3.exe"
    $WideO2 = Join-Path $BuildRoot "wide_O2.exe"
    $WideO3 = Join-Path $BuildRoot "wide_O3.exe"
    $SelectedO2 = Join-Path $BuildRoot "selected_O2.exe"
    $SelectedO3 = Join-Path $BuildRoot "selected_O3.exe"

    & g++ -std=c++20 -O2 -Wall -Wextra -pedantic $Source -o $NarrowO2
    Assert-Exit "compile narrow O2"
    & g++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -pedantic $Source -o $NarrowO3
    Assert-Exit "compile narrow O3"
    & g++ -std=c++20 -O2 -Wall -Wextra -pedantic $WideSource -o $WideO2
    Assert-Exit "compile wide O2"
    & g++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -pedantic $WideSource -o $WideO3
    Assert-Exit "compile wide O3"
    & g++ -std=c++20 -O2 -Wall -Wextra -pedantic $SelectedSource -o $SelectedO2
    Assert-Exit "compile selected O2"
    & g++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -pedantic $SelectedSource -o $SelectedO3
    Assert-Exit "compile selected O3"

    $NarrowInput = Join-Path $Packet "inputs/universe22647.csv"
    $BroadInput = Join-Path $Packet "inputs/thm4231_remainder181194.csv"
    Invoke-Audit $NarrowO2 $NarrowInput (Join-Path $BuildRoot "narrow_screen_O2.csv") (Join-Path $BuildRoot "narrow_exact_O2.csv") (Join-Path $BuildRoot "narrow_O2.out")
    Invoke-Audit $NarrowO3 $NarrowInput (Join-Path $BuildRoot "narrow_screen_O3.csv") (Join-Path $BuildRoot "narrow_exact_O3.csv") (Join-Path $BuildRoot "narrow_O3.out")
    Invoke-Audit $WideO2 $BroadInput (Join-Path $BuildRoot "broad_screen_O2.csv") (Join-Path $BuildRoot "broad_exact_O2.csv") (Join-Path $BuildRoot "broad_O2.out")
    Invoke-Audit $WideO3 $BroadInput (Join-Path $BuildRoot "broad_screen_O3.csv") (Join-Path $BuildRoot "broad_exact_O3.csv") (Join-Path $BuildRoot "broad_O3.out")

    Assert-Same (Join-Path $BuildRoot "narrow_screen_O2.csv") (Join-Path $BuildRoot "narrow_screen_O3.csv") "narrow O2/O3 screen"
    Assert-Same (Join-Path $BuildRoot "narrow_screen_O3.csv") (Join-Path $Packet "results/coarse_screen_O3.csv") "narrow frozen screen"
    Assert-Same (Join-Path $BuildRoot "narrow_exact_O2.csv") (Join-Path $BuildRoot "narrow_exact_O3.csv") "narrow O2/O3 exact"
    Assert-Same (Join-Path $BuildRoot "narrow_exact_O3.csv") (Join-Path $Packet "results/hostile107_exact_O3.csv") "narrow frozen exact"
    Assert-Same (Join-Path $BuildRoot "broad_screen_O2.csv") (Join-Path $BuildRoot "broad_screen_O3.csv") "broad O2/O3 screen"
    Assert-Same (Join-Path $BuildRoot "broad_screen_O3.csv") (Join-Path $Packet "results/thm4231_coarse_screen_O3.csv") "broad frozen screen"
    Assert-Same (Join-Path $BuildRoot "broad_exact_O2.csv") (Join-Path $BuildRoot "broad_exact_O3.csv") "broad O2/O3 exact"
    Assert-Same (Join-Path $BuildRoot "broad_exact_O3.csv") (Join-Path $Packet "results/thm4231_hostile_exact_O3.csv") "broad frozen exact"

    foreach ($Opt in @("O2", "O3")) {
        $Exe = if ($Opt -eq "O2") { $SelectedO2 } else { $SelectedO3 }
        $SelectedLines = & $Exe
        Assert-Exit "selected ratio $Opt"
        $SelectedText = (($SelectedLines -join "`r`n") + "`r`n")
        $SelectedPath = Join-Path $BuildRoot "selected_$Opt.out"
        [IO.File]::WriteAllText($SelectedPath, $SelectedText, [Text.UTF8Encoding]::new($false))
        $FrozenText = (Get-Content -LiteralPath (Join-Path $Packet "results/selected_ratio_exact_$Opt.out") -Raw).Replace("`r`n", "`n")
        if ($SelectedText.Replace("`r`n", "`n") -ne $FrozenText) { throw "selected $Opt output mismatch" }
    }

    $RawLines = & python (Join-Path $Packet "src/rawcell_winner_replay.py") (Join-Path $BuildRoot "broad_exact_O3.csv") (Join-Path $BuildRoot "selected_O3.out") (Join-Path $BuildRoot "rawcell.csv")
    Assert-Exit "raw-cell replay"
    if (($RawLines -join "`n") -notmatch "VERDICT PASS") { throw "raw-cell replay did not pass" }
    Assert-Same (Join-Path $BuildRoot "rawcell.csv") (Join-Path $Packet "results/rawcell_winner_replay.csv") "raw-cell frozen ledger"

    $FilterLines = & python (Join-Path $Packet "src/global_minimum_filter.py") (Join-Path $BuildRoot "broad_screen_O3.csv") (Join-Path $BuildRoot "broad_exact_O3.csv") (Join-Path $BuildRoot "selected_O3.out")
    Assert-Exit "global minimum filter"
    if (($FilterLines -join "`n") -notmatch "VERDICT PASS") { throw "global filter did not pass" }

    & python (Join-Path $Packet "verify_packet.py") --no-rerun-rawcell
    Assert-Exit "packet verifier"
    Write-Output "REPRODUCTION VERDICT PASS"
}
finally {
    $ResolvedBuild = [IO.Path]::GetFullPath($BuildRoot)
    $ResolvedTemp = [IO.Path]::GetFullPath([IO.Path]::GetTempPath())
    $Leaf = Split-Path -Leaf $ResolvedBuild
    if (-not $ResolvedBuild.StartsWith($ResolvedTemp, [StringComparison]::OrdinalIgnoreCase) -or
        -not $Leaf.StartsWith("lrc14-rank2-independent-", [StringComparison]::Ordinal)) {
        throw "refusing to remove unverified build path: $ResolvedBuild"
    }
    if (Test-Path -LiteralPath $ResolvedBuild) {
        Remove-Item -LiteralPath $ResolvedBuild -Recurse -Force
    }
}
