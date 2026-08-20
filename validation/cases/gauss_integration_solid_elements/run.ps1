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
$modelDir = Join-Path $caseDir "input\models"
$vfmScriptSource = Join-Path $caseDir "scripts\vfm_config.py"
$validator = Join-Path $caseDir "scripts\validate.py"
$expected = Join-Path $caseDir "expected\cases.json"

if (-not $FebioExe) {
    $FebioExe = if ($env:FEBIO_EXE) { $env:FEBIO_EXE } else {
        "C:\Program Files\FEBioStudio\bin\febio4.exe"
    }
}

if (-not $PluginPath) {
    $pluginCandidates = @(
        (Join-Path $repoRoot "out\build\local-vs-release\FEBio_VFM_Task\Release\FEBio_VFM_Task.dll"),
        (Join-Path $repoRoot "out\build\x64-release\FEBio_VFM_Task\Release\FEBio_VFM_Task.dll"),
        (Join-Path $repoRoot "out\build\x64-release\FEBio_VFM_Task\FEBio_VFM_Task.dll")
    )
    $PluginPath = $pluginCandidates | Where-Object { Test-Path -LiteralPath $_ } | Select-Object -First 1
}

foreach ($path in @($FebioExe, $PluginPath, $vfmScriptSource, $validator, $expected)) {
    if (-not $path -or -not (Test-Path -LiteralPath $path -PathType Leaf)) {
        throw "Required file was not found: $path"
    }
}

$uv = (Get-Command uv -ErrorAction Stop).Source
if (-not $SkipSync) {
    & $uv sync --project $validationRoot --python 3.13 --frozen
    if ($LASTEXITCODE -ne 0) { throw "uv failed to restore the validation environment" }
}
$pythonExe = Join-Path $validationRoot ".venv\Scripts\python.exe"
$pythonInfo = (& $pythonExe -c "import json, site, sys; print(json.dumps({'base_prefix': sys.base_prefix, 'site_packages': site.getsitepackages()}))") | ConvertFrom-Json

$pluginDir = Split-Path -Parent $PluginPath
$febioDir = Split-Path -Parent $FebioExe
$vcpkgBin = Join-Path $repoRoot "out\vcpkg_installed\x64-windows\bin"
$pathParts = @($pythonInfo.base_prefix, $pluginDir, $febioDir)
if (Test-Path -LiteralPath $vcpkgBin -PathType Container) { $pathParts += $vcpkgBin }
$env:PATH = (($pathParts + $env:PATH) -join [IO.Path]::PathSeparator)
$env:PYTHONHOME = $pythonInfo.base_prefix
$env:PYTHONPATH = ($pythonInfo.site_packages -join [IO.Path]::PathSeparator)

$runRoot = Join-Path $repoRoot ("out\validation\gauss_integration_solid_elements\" + (Get-Date -Format "yyyyMMdd-HHmmss"))
$models = @("hex8", "tet4", "penta6", "pyra5")
$validationFailures = @()

foreach ($elementType in $models) {
    $modelSource = Join-Path $modelDir "$elementType.feb"
    if (-not (Test-Path -LiteralPath $modelSource -PathType Leaf)) {
        throw "Model was not found: $modelSource"
    }

    $runDir = Join-Path $runRoot $elementType
    New-Item -ItemType Directory -Path (Join-Path $runDir "temp\debug\result") -Force | Out-Null
    $runModel = Join-Path $runDir "model.feb"
    $runVfmScript = Join-Path $runDir "vfm_config.py"
    $metrics = Join-Path $runDir "metrics.json"
    $runLog = Join-Path $runDir "run.log"
    Copy-Item -LiteralPath $modelSource -Destination $runModel
    Copy-Item -LiteralPath $vfmScriptSource -Destination $runVfmScript

    $env:VFM_GAUSS_ELEMENT_TYPE = $elementType
    $env:VFM_GAUSS_METRICS = $metrics
    $arguments = @("-i", $runModel, "-dump", "-import", $PluginPath, "-task=VFM", $runVfmScript)

    Write-Host "Running Gaussian integration validation: $elementType"
    Push-Location $runDir
    try {
        & $FebioExe @arguments *> $runLog
        $febioExitCode = $LASTEXITCODE
    }
    finally {
        Pop-Location
    }

    if ($febioExitCode -ne 0 -or -not (Test-Path -LiteralPath $metrics -PathType Leaf)) {
        Get-Content -LiteralPath $runLog -Tail 100
        throw "FEBio/VFM Gaussian integration run failed for $elementType (exit code $febioExitCode)"
    }
    & $pythonExe $validator --actual $metrics --expected $expected
    if ($LASTEXITCODE -ne 0) {
        $validationFailures += $elementType
    }
}

if ($validationFailures.Count -gt 0) {
    throw "Gaussian integration validation failed for: $($validationFailures -join ', ')"
}
Write-Host "Gaussian integration validation completed: $runRoot"
