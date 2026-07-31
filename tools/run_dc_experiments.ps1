param(
    [string]$Executable = ".\build\netlist_dc_solver.exe",
    [string]$NetlistDirectory = ".\netlist",
    [string]$OutputDirectory = ".\result\dc_experiments",
    [double]$Temperature = 26.0
)

# Run all four netlist DC algorithms and give each run its own CSV directory.
$modes = @(
    @{ Name = "direct_newton"; Arguments = @("--nr") },
    @{ Name = "source_stepping"; Arguments = @() },
    @{ Name = "sequential_source_stepping"; Arguments = @("--sequential-sources") },
    @{ Name = "gmin_stepping"; Arguments = @("--gmin-stepping") }
)

# Every Netlist.txt below the input directory is one independent experiment.
$netlists = @(Get-ChildItem -Path $NetlistDirectory -Recurse -Filter "Netlist.txt")
if ($netlists.Count -eq 0) {
    throw "No Netlist.txt files were found under $NetlistDirectory."
}

New-Item -ItemType Directory -Force -Path $OutputDirectory | Out-Null
$comparison = @()

foreach ($netlist in $netlists) {
    $circuitName = $netlist.Directory.Name
    foreach ($mode in $modes) {
        $reportDirectory = Join-Path $OutputDirectory "$circuitName\$($mode.Name)"
        New-Item -ItemType Directory -Force -Path $reportDirectory | Out-Null

        # The solver writes dc_steps.csv and dc_summary.csv for this one run.
        $arguments = @($netlist.FullName, "--temp", $Temperature) +
                     @($mode.Arguments) +
                     @("--report-dir", $reportDirectory)
        & $Executable @arguments
        $exitCode = $LASTEXITCODE
        $summaryPath = Join-Path $reportDirectory "dc_summary.csv"

        if (Test-Path $summaryPath) {
            $row = Import-Csv $summaryPath | Select-Object -First 1
            $row | Add-Member -NotePropertyName circuit -NotePropertyValue $circuitName
            $row | Add-Member -NotePropertyName process_exit_code -NotePropertyValue $exitCode
            $comparison += $row
        }
    }
}

# This file is suitable for direct comparison in VS Code, Excel, or Python.
$comparison |
    Select-Object circuit, algorithm, converged, process_exit_code, `
        total_newton_iterations, accepted_newton_iterations, `
        failed_probe_newton_iterations, accepted_steps, `
        accepted_line_search_reductions, total_runtime_ms, final_residual_inf |
    Export-Csv -NoTypeInformation -Encoding utf8 `
        (Join-Path $OutputDirectory "dc_algorithm_comparison.csv")

Write-Host "Experiment reports were written to: $OutputDirectory"
