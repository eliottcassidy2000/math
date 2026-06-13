param(
    [string]$ScratchDir = ".scratch\lean",
    [switch]$SmokeTest
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$RepoRoot = Resolve-Path (Join-Path $ScriptDir "..")
$ScratchPath = Join-Path $RepoRoot $ScratchDir
$ElanHome = if ($env:ELAN_HOME) { $env:ELAN_HOME } else { Join-Path $env:USERPROFILE ".elan" }
$ElanBin = Join-Path $ElanHome "bin"

function Invoke-VersionCheck {
    param(
        [string]$Tool,
        [int]$TimeoutSeconds = 60
    )

    $OutputPath = Join-Path ([IO.Path]::GetTempPath()) "$Tool-version.out"
    $ErrorPath = Join-Path ([IO.Path]::GetTempPath()) "$Tool-version.err"
    Remove-Item -LiteralPath $OutputPath, $ErrorPath -Force -ErrorAction SilentlyContinue

    $Process = Start-Process `
        -FilePath $Tool `
        -ArgumentList "--version" `
        -NoNewWindow `
        -PassThru `
        -RedirectStandardOutput $OutputPath `
        -RedirectStandardError $ErrorPath

    if (-not $Process.WaitForExit($TimeoutSeconds * 1000)) {
        Stop-Process -Id $Process.Id -Force -ErrorAction SilentlyContinue
        throw "$Tool --version timed out after ${TimeoutSeconds}s. The Lean toolchain may still be downloading."
    }
    $Process.Refresh()

    $Stdout = ""
    $Stderr = ""
    if (Test-Path $OutputPath) {
        $Stdout = [string](Get-Content $OutputPath -Raw)
    }
    if (Test-Path $ErrorPath) {
        $Stderr = [string](Get-Content $ErrorPath -Raw)
    }
    $StdoutText = "$Stdout"
    $StderrText = "$Stderr"
    if ($StdoutText.Trim()) { Write-Host $StdoutText.Trim() }
    if ($StderrText.Trim()) { Write-Host $StderrText.Trim() }
    if ($null -ne $Process.ExitCode -and $Process.ExitCode -ne 0) {
        throw "$Tool --version failed with exit code $($Process.ExitCode)."
    }
}

if (($env:Path -split [IO.Path]::PathSeparator) -notcontains $ElanBin) {
    $env:Path = "$ElanBin$([IO.Path]::PathSeparator)$env:Path"
}

Write-Host "Repo root: $RepoRoot"
Write-Host "ELAN_HOME: $ElanHome"
Write-Host "Scratch: $ScratchPath"
Write-Host "PATH contains elan bin: $(($env:Path -split [IO.Path]::PathSeparator) -contains $ElanBin)"
Write-Host ""

$Missing = @()
foreach ($Tool in @("elan", "lean", "lake")) {
    $Command = Get-Command $Tool -ErrorAction SilentlyContinue
    if (-not $Command) {
        $Missing += $Tool
        Write-Host "${Tool}: MISSING"
    } else {
        Write-Host "${Tool}: $($Command.Source)"
        Invoke-VersionCheck -Tool $Tool
    }
    Write-Host ""
}

if ($Missing.Count -gt 0) {
    throw "Missing Lean tools: $($Missing -join ', '). Run scripts/bootstrap_lean.ps1 first."
}

if ($SmokeTest) {
    New-Item -ItemType Directory -Force $ScratchPath | Out-Null
    $SmokeFile = Join-Path $ScratchPath "smoke.lean"
    Set-Content -Path $SmokeFile -Value "theorem lean_bootstrap_smoke : True := by`n  trivial`n" -Encoding ASCII
    Write-Host "Running smoke test: lean $SmokeFile"
    lean $SmokeFile
}
