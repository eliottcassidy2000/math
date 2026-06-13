param(
    [string]$ScratchDir = ".scratch\lean",
    [string]$ToolchainVersion = "4.30.0"
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

function Write-Step {
    param([string]$Message)
    Write-Host "==> $Message"
}

function Download-File {
    param(
        [string]$Uri,
        [string]$OutFile
    )

    [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12
    $LastError = $null
    for ($Attempt = 1; $Attempt -le 3; $Attempt++) {
        try {
            Invoke-WebRequest -Uri $Uri -OutFile $OutFile -UseBasicParsing -Headers @{"User-Agent" = "codex-lean-bootstrap"}
            return
        } catch {
            $LastError = $_
            Start-Sleep -Seconds $Attempt
        }
    }

    if (Get-Command curl.exe -ErrorAction SilentlyContinue) {
        & curl.exe -L --fail --retry 3 --output $OutFile $Uri
        if ($LASTEXITCODE -eq 0) {
            return
        }
    }

    throw $LastError
}

function Install-LeanToolchainArchive {
    param([string]$Version)

    $ArchiveName = "lean-$Version-windows.tar.zst"
    $ArchiveUri = "https://releases.lean-lang.org/lean4/v$Version/$ArchiveName"
    $ArchivePath = Join-Path ([IO.Path]::GetTempPath()) $ArchiveName
    $ExtractPath = Join-Path ([IO.Path]::GetTempPath()) "lean-$Version-extract"
    $ToolchainPath = Join-Path $ElanHome "toolchains\leanprover--lean4---v$Version"

    if (Test-Path $ToolchainPath) {
        Write-Step "Lean $Version toolchain already extracted"
        return
    }

    Write-Step "Downloading Lean $Version toolchain archive"
    $Aria = Get-Command aria2c -ErrorAction SilentlyContinue
    if (-not $Aria) {
        $Aria = Get-ChildItem "$env:LOCALAPPDATA\Microsoft\WinGet\Packages" -Recurse -Filter aria2c.exe -ErrorAction SilentlyContinue | Select-Object -First 1
    }

    if ($Aria) {
        $AriaPath = if ($Aria.Source) { $Aria.Source } else { $Aria.FullName }
        & $AriaPath `
            --continue=true `
            --max-connection-per-server=16 `
            --split=16 `
            --min-split-size=1M `
            --dir=([IO.Path]::GetTempPath()) `
            --out=$ArchiveName `
            $ArchiveUri
        if ($LASTEXITCODE -ne 0) {
            throw "aria2c failed to download $ArchiveUri"
        }
    } else {
        Download-File -Uri $ArchiveUri -OutFile $ArchivePath
    }

    Write-Step "Extracting Lean $Version toolchain"
    if (Test-Path $ExtractPath) {
        Remove-Item -Recurse -Force $ExtractPath
    }
    New-Item -ItemType Directory -Force $ExtractPath | Out-Null
    New-Item -ItemType Directory -Force (Join-Path $ElanHome "toolchains") | Out-Null
    tar -xf $ArchivePath -C $ExtractPath

    $ExtractedRoot = Join-Path $ExtractPath "lean-$Version-windows"
    if (-not (Test-Path $ExtractedRoot)) {
        throw "Lean archive did not contain expected directory: $ExtractedRoot"
    }

    if (Test-Path $ToolchainPath) {
        Remove-Item -Recurse -Force $ToolchainPath
    }
    Move-Item -LiteralPath $ExtractedRoot -Destination $ToolchainPath
    Remove-Item -Recurse -Force $ExtractPath
}

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$RepoRoot = Resolve-Path (Join-Path $ScriptDir "..")
$ScratchPath = Join-Path $RepoRoot $ScratchDir
$ElanHome = if ($env:ELAN_HOME) { $env:ELAN_HOME } else { Join-Path $env:USERPROFILE ".elan" }
$ElanBin = Join-Path $ElanHome "bin"

Write-Step "Preparing Lean scratch directory"
New-Item -ItemType Directory -Force $ScratchPath | Out-Null

if (($env:Path -split [IO.Path]::PathSeparator) -notcontains $ElanBin) {
    $env:Path = "$ElanBin$([IO.Path]::PathSeparator)$env:Path"
}

if (-not (Get-Command elan -ErrorAction SilentlyContinue)) {
    Write-Step "Installing elan"
    New-Item -ItemType Directory -Force $ElanHome | Out-Null
    $ZipPath = Join-Path ([IO.Path]::GetTempPath()) "elan-windows.zip"
    $ExtractPath = Join-Path ([IO.Path]::GetTempPath()) "elan-windows"
    Download-File `
        -Uri "https://github.com/leanprover/elan/releases/latest/download/elan-x86_64-pc-windows-msvc.zip" `
        -OutFile $ZipPath
    if (Test-Path $ExtractPath) {
        Remove-Item -Recurse -Force $ExtractPath
    }
    Expand-Archive -Path $ZipPath -DestinationPath $ExtractPath -Force
    $Installer = Get-ChildItem -Path $ExtractPath -Filter "elan-init.exe" -Recurse | Select-Object -First 1
    if (-not $Installer) {
        throw "Downloaded elan archive did not contain elan-init.exe"
    }
    & $Installer.FullName -y --default-toolchain stable
    if (($env:Path -split [IO.Path]::PathSeparator) -notcontains $ElanBin) {
        $env:Path = "$ElanBin$([IO.Path]::PathSeparator)$env:Path"
    }
} else {
    Write-Step "elan already installed"
}

Write-Step "Installing Lean $ToolchainVersion toolchain"
Install-LeanToolchainArchive -Version $ToolchainVersion
elan default $ToolchainVersion

Write-Step "Verifying Lean tools"
elan --version
lean --version
lake --version

Write-Host ""
Write-Host "Lean environment ready."
Write-Host "Scratch directory: $ScratchPath"
Write-Host "For this shell, elan bin is on PATH: $ElanBin"
