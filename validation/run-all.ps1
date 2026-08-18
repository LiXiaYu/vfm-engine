[CmdletBinding()]
param(
    [string]$FebioExe = "",
    [string]$PluginPath = "",
    [switch]$SkipSync
)

$ErrorActionPreference = "Stop"
$caseRoot = Join-Path $PSScriptRoot "cases"
$runners = Get-ChildItem -LiteralPath $caseRoot -Filter "run.ps1" -File -Recurse | Sort-Object FullName

if ($runners.Count -eq 0) {
    throw "No validation case runners were found below $caseRoot"
}

foreach ($runner in $runners) {
    Write-Host "Running validation case: $($runner.Directory.Name)"
    & $runner.FullName -FebioExe $FebioExe -PluginPath $PluginPath -SkipSync:$SkipSync
    if ($LASTEXITCODE -ne 0) {
        exit $LASTEXITCODE
    }
}
