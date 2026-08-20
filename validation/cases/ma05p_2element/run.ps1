[CmdletBinding()]
param(
    [string]$FebioExe = "",
    [string]$PluginPath = "",
    [switch]$SkipSync
)

$ErrorActionPreference = "Stop"

$caseDir = $PSScriptRoot
$repoRoot = (Resolve-Path (Join-Path $caseDir "..\..\..")).Path
$validationRoot = Join-Path $repoRoot "validation"
$modelSource = Join-Path $caseDir "input\model\ma05p.feb"
$vfmScriptSource = Join-Path $caseDir "scripts\vfm_config.py"
$validator = Join-Path $caseDir "scripts\validate.py"
$expectedMetrics = Join-Path $caseDir "expected\metrics.json"

if (-not $FebioExe) {
    $FebioExe = if ($env:FEBIO_EXE) {
        $env:FEBIO_EXE
    } else {
        "C:\Program Files\FEBioStudio\bin\febio4.exe"
    }
}

if (-not $PluginPath) {
    if ($env:VFM_PLUGIN) {
        $PluginPath = $env:VFM_PLUGIN
    } else {
        $pluginCandidates = @(
            (Join-Path $repoRoot "out\build\local-vs-release\FEBio_VFM_Task\Release\FEBio_VFM_Task.dll"),
            (Join-Path $repoRoot "out\build\x64-release\FEBio_VFM_Task\Release\FEBio_VFM_Task.dll"),
            (Join-Path $repoRoot "out\build\x64-release\FEBio_VFM_Task\FEBio_VFM_Task.dll")
        )
        $PluginPath = $pluginCandidates | Where-Object { Test-Path -LiteralPath $_ } | Select-Object -First 1
    }
}

$requiredFiles = [ordered]@{
    "FEBio executable" = $FebioExe
    "VFM plugin" = $PluginPath
    "FEBio model" = $modelSource
    "VFM configuration" = $vfmScriptSource
    "metrics validator" = $validator
}
foreach ($item in $requiredFiles.GetEnumerator()) {
    if (-not $item.Value -or -not (Test-Path -LiteralPath $item.Value -PathType Leaf)) {
        throw "$($item.Key) was not found: $($item.Value)"
    }
}

$uv = (Get-Command uv -ErrorAction SilentlyContinue).Source
if (-not $uv) {
    throw "uv was not found on PATH"
}

if (-not $SkipSync) {
    & $uv sync --project $validationRoot --python 3.13 --frozen
    if ($LASTEXITCODE -ne 0) {
        throw "uv failed to restore the validation environment"
    }
}

$pythonExe = Join-Path $validationRoot ".venv\Scripts\python.exe"
if (-not (Test-Path -LiteralPath $pythonExe -PathType Leaf)) {
    throw "Validation Python environment not found: $pythonExe"
}

$pythonInfoJson = & $pythonExe -c "import json, site, sys; print(json.dumps({'base_prefix': sys.base_prefix, 'site_packages': site.getsitepackages()}))"
if ($LASTEXITCODE -ne 0) {
    throw "Could not inspect the validation Python environment"
}
$pythonInfo = $pythonInfoJson | ConvertFrom-Json

$runId = Get-Date -Format "yyyyMMdd-HHmmss"
$runDir = Join-Path $repoRoot "out\validation\ma05p_2element\$runId"
New-Item -ItemType Directory -Path $runDir -Force | Out-Null
foreach ($relativeDir in @("temp\debug\result", "temp\debug\tecplot")) {
    New-Item -ItemType Directory -Path (Join-Path $runDir $relativeDir) -Force | Out-Null
}
$runModel = Join-Path $runDir "model.feb"
$runVfmScript = Join-Path $runDir "vfm_config.py"
$metricsPath = Join-Path $runDir "validation_metrics.json"
$runLog = Join-Path $runDir "run.log"
Copy-Item -LiteralPath $modelSource -Destination $runModel
Copy-Item -LiteralPath $vfmScriptSource -Destination $runVfmScript

$pluginDir = Split-Path -Parent $PluginPath
$febioDir = Split-Path -Parent $FebioExe
$vcpkgBin = Join-Path $repoRoot "out\vcpkg_installed\x64-windows\bin"
$pathParts = @($pythonInfo.base_prefix, $pluginDir, $febioDir)
if (Test-Path -LiteralPath $vcpkgBin -PathType Container) {
    $pathParts += $vcpkgBin
}
$env:PATH = (($pathParts + $env:PATH) -join [IO.Path]::PathSeparator)
$env:PYTHONHOME = $pythonInfo.base_prefix
$env:PYTHONPATH = ($pythonInfo.site_packages -join [IO.Path]::PathSeparator)
$env:VFM_VALIDATION_METRICS = $metricsPath
$env:VFM_VALIDATION_SEED = "20260820"

# fun_for_optim_T currently accumulates the multi-virtual-field loss into one
# shared scalar inside an OpenMP loop.  Keep this validation deterministic until
# that core reduction is fixed without altering the legacy source encoding.
$previousOmpNumThreads = $env:OMP_NUM_THREADS
$env:OMP_NUM_THREADS = "1"

$febioArguments = @(
    "-i", $runModel,
    "-dump",
    "-import", $PluginPath,
    "-task=VFM", $runVfmScript
)

Write-Host "Case directory : $caseDir"
Write-Host "Run directory  : $runDir"
Write-Host "FEBio          : $FebioExe"
Write-Host "Plugin         : $PluginPath"
Write-Host "Python home    : $($pythonInfo.base_prefix)"

Push-Location $runDir
try {
    & $FebioExe @febioArguments *> $runLog
    $febioExitCode = $LASTEXITCODE
}
finally {
    Pop-Location
    if ($null -eq $previousOmpNumThreads) {
        Remove-Item Env:OMP_NUM_THREADS -ErrorAction SilentlyContinue
    } else {
        $env:OMP_NUM_THREADS = $previousOmpNumThreads
    }
}

$febioReportedErrorTermination = Select-String `
    -LiteralPath $runLog `
    -Pattern "E R R O R   T E R M I N A T I O N" `
    -SimpleMatch `
    -Quiet

if ($febioExitCode -ne 0 -or $febioReportedErrorTermination) {
    Get-Content -LiteralPath $runLog -Tail 100
    if ($febioReportedErrorTermination) {
        throw "FEBio reported error termination in $runLog"
    }
    throw "FEBio failed with exit code $febioExitCode. See $runLog"
}

$febioSolveSkipped = Select-String `
    -LiteralPath $runLog `
    -Pattern "Skipping FEBio solve by configuration." `
    -SimpleMatch `
    -Quiet
if ($febioSolveSkipped) {
    throw "The MA05P validation unexpectedly skipped the FEBio forward solve."
}

if (-not (Test-Path -LiteralPath $metricsPath -PathType Leaf)) {
    throw "The MA05P VFM configuration did not produce metrics: $metricsPath"
}

& $pythonExe $validator --actual $metricsPath --expected $expectedMetrics
if ($LASTEXITCODE -ne 0) {
    exit $LASTEXITCODE
}

Write-Host "Validation run completed: $runDir"
